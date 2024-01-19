from openmm.app import *
from openmm import *
from openmm.unit import *
from pdbfixer import PDBFixer
from sys import stdout

import time
import os
import logging
import argparse
import datetime
import sys

def run_simulation(pdb, params=None):
    """
    requires a PDBFile object
    """
    total_simulation_time = 1*nanoseconds # production run time
    total_equilibriation_time = 100*picoseconds # equilibriation run time
    temperature = 300 # kelvin
    step_size = 0.002 # of a picosecond

    report_every = 1000 # every X steps
    coordinates_output_file = 'output.pdb'
    state_output_file = 'state.csv'

    if params is not None:
        # total_simulation_time     = params['total_simulation_time']
        # total_equilibriation_time = params['total_equilibriation_time']
        # temperature  = params['temperature']
        # step_size    = params['step_size']
        # report_every = params['report_every']

        coordinates_output_file = params['coordinates_output_file']
        state_output_file       = params['state_output_file']

    using_pbc = False # use periodic boundary conditions
    restrain_backbone = True # restrain the backbone during NVT and NPT equilibration, but not during production

    total_n_equilibriation_steps = int(total_equilibriation_time / (step_size*picoseconds) ) # 100 ps
    total_n_steps = int(total_simulation_time / (step_size*picoseconds))+1
    logging.info(f"total equilibriation time (NVT+NPT)    = {total_equilibriation_time*2}")
    logging.info(f"total_n_equilibriation_steps (NVT+NPT) = {total_n_equilibriation_steps*2}")
    logging.info(f"total_simulation_time = {total_simulation_time}")
    logging.info(f"total_n_steps = {total_n_steps}")

    ##############
    # specify forcefield
    ##############
    forcefield = ForceField('charmm36.xml', 'charmm36/water.xml') # uses 'default' forcefields from Jorgensen et al. (1983)
    # forcefield = ForceField('amber14-all.xml', 'amber14/tip3p.xml')
    
    ##############
    # remove crystal water and add hydrogen atoms
    ##############
    modeller = Modeller(pdb.topology, pdb.positions)
    modeller.deleteWater()
    residues = modeller.addHydrogens(forcefield)

    ##############
    # add solvent
    ##############
    # t0 = time.time()
    modeller.addSolvent(forcefield, padding=1.0*nanometer)
    # print(f"took = {round(time.time() - t0, 4)}s")

    ##############
    # setup system and integrator
    ##############
    system = forcefield.createSystem(modeller.topology, 
                                    nonbondedMethod=PME, 
                                    nonbondedCutoff=1.0*nanometer, 
                                    constraints=HBonds
    )
    integrator = LangevinMiddleIntegrator(temperature*kelvin,
                                        1/picosecond, 
                                        step_size*picoseconds
    )

    simulation = Simulation(modeller.topology, system, integrator)
    simulation.context.setPositions(modeller.positions)

    ##############
    # local energy minimization
    ############## 
    logging.info("Minimizing energy...")
    t0 = time.time()
    simulation.minimizeEnergy()
    logging.info(f"Minimized energy in {round(time.time() - t0, 4)}s")

    ##############
    # Reporters
    ##############
    # position reporter
    simulation.reporters.append(
        PDBReporter(coordinates_output_file, report_every)
    )

    # state reporter 
    simulation.reporters.append(
        StateDataReporter(
            state_output_file, 
            report_every, 
            step=True, 
            potentialEnergy=True, 
            temperature=True,
            volume=True
        )
    )

    ##############
    # Restrain protein backbone
    ##############
    if restrain_backbone:
        logging.info("Restrain protein backbone...")
        if using_pbc:
            restraint = CustomExternalForce('k*periodicdistance(x, y, z, x0, y0, z0)^2')
        else:
            restraint = CustomExternalForce('k*((x-x0)^2+(y-y0)^2+(z-z0)^2)')
        restraint_force_system_idx = system.addForce(restraint)
        restraint.addGlobalParameter('k', 100.0*kilojoules_per_mole/nanometer)
        restraint.addPerParticleParameter('x0')
        restraint.addPerParticleParameter('y0')
        restraint.addPerParticleParameter('z0')

        # restrain all Calpha atoms
        for atom in pdb.topology.atoms():
            if atom.name == 'CA':
                restraint.addParticle(atom.index, pdb.positions[atom.index])


    ##############
    # run NVT equilibriation
    ##############
    logging.info(f"Running NVT equilibriation for {total_n_equilibriation_steps} steps...")
    simulation.step(total_n_equilibriation_steps)
    
    ##############
    # run NPT equilibriation
    ##############
    logging.info(f"Running NPT equilibriation for {total_n_equilibriation_steps} steps...")
    system.addForce(
        MonteCarloBarostat(1*bar, temperature*kelvin)
    )
    simulation.context.reinitialize(preserveState=True)

    simulation.step(total_n_equilibriation_steps)

    ##############
    # remove restraint if there is one
    ##############
    if restrain_backbone:
        logging.info("Removing backbone restraint...")
        system.removeForce(restraint_force_system_idx)
    
    ##############
    # run production
    ##############
    logging.info(f"Running production for {total_n_steps} steps...")
    simulation.step(total_n_steps)

if __name__=="__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--pdb_file', type=str, required=True, help='PDB file to simulate.')
    parser.add_argument('--pdb_dir',  type=str, required=True, help='directory containing pdb file(s).')
    parser.add_argument("--prefix",   default="", type=str, required=False, help="prefix for output files")
    parser.add_argument("--device",   default="cpu", type=str,choices=["cpu","cuda"], required=False, help="device to run on")
    parser.add_argument("--slurm_id", default="", type=str, required=False, help="slurm id (for logging)")

    args = parser.parse_args()
    pdb_file = args.pdb_file
    pdb_dir  = args.pdb_dir
    prefix   = args.prefix
    slurm_id = args.slurm_id
    device   = args.device
    device   = device.upper()
    if prefix != "":
        prefix += "_"
    if slurm_id != "":
        slurm_id += "_"

    today = datetime.datetime.now()
    logfilename = f"{slurm_id}{prefix}sim_protein_in_water.log"
    logging.basicConfig(
        filename="logs/"+logfilename,
        filemode='a',
        format='%(asctime)s,%(msecs)d %(name)s %(levelname)s %(message)s',
        datefmt='%H:%M:%S',
        level=logging.INFO
    )
    logging.info(f"Job ID was {slurm_id}")
    logging.info(f"n CPU cores = {os.cpu_count()}")
    
    if pdb_dir[-1] == "/":
        pdb_dir = pdb_dir[:-1]
    pdb_fpath = pdb_dir + "/" + pdb_file.split("/")[-1]
    
    # pdb = PDBFile("starPep_06810_test_pdbfixer.pdb")
    pdb = PDBFile(pdb_fpath)

    logging.info(f"using pdb file: {pdb_file}")
    logging.info(f"num residues in pdb = {pdb.topology.getNumResidues()}")
    
    params = {}
    params["coordinates_output_file"] = f"outputs/{slurm_id}{prefix}coordinates.pdb" 
    params["state_output_file"]       = f"outputs/{slurm_id}{prefix}state.csv"

    ##############################
    # setup and run simulation
    ##############################
    logging.info(f"Running simulation on {pdb_file}")
    t0 = time.time()
    run_simulation(pdb, params=params)
    logging.info(f"Simulation took {round(time.time() - t0, 4)}s")