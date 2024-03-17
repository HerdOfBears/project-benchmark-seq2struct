from openmm.app import *
from openmm import *
from openmm.unit import *
from pdbfixer import PDBFixer
from sys import stdout
from helpers import make_dir_for_starpepid

import time
import os
import logging
import argparse
import datetime
import sys

# Function to add backbone position restraints
def add_backbone_posres(system, pdb, restraint_force, periodic_boundaries=True):
    if periodic_boundaries:
        force = CustomExternalForce("k*periodicdistance(x, y, z, x0, y0, z0)^2")
    else:
        force = CustomExternalForce('k*((x-x0)^2+(y-y0)^2+(z-z0)^2)')        

    force_amount = restraint_force * kilocalories_per_mole/angstroms**2
    force.addGlobalParameter("k", force_amount)
    force.addPerParticleParameter("x0")
    force.addPerParticleParameter("y0")
    force.addPerParticleParameter("z0")

    # restrain all Calpha atoms
    for atom in pdb.topology.atoms():
        if atom.name in ('CA', 'C', 'N'): # all heavy atoms
            force.addParticle(atom.index, pdb.positions[atom.index])

    system.addForce(force)
    

def run_simulation(pdb, params=None):
    """
    requires a PDBFile object
    """
    total_simulation_time = 10*nanoseconds # production run time
    total_equilibriation_time = 100*picoseconds # equilibriation run time
    temperature = 300 # kelvin
    step_size = 0.002 # of a picosecond
    energy_threshold = None

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
        energy_threshold        = params["energy_threshold"]

    using_pbc = True # use periodic boundary conditions
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
    # add solvent (neutralizes the system)
    ##############
    modeller.addSolvent(forcefield, padding=1.0*nanometer)

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
    if energy_threshold is not None:
        simulation.minimizeEnergy(tolerance=energy_threshold*nanometer)
    else:
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
            volume=True,
            density=True
        )
    )

    ##############
    # Restrain protein backbone
    ##############
    if restrain_backbone:
        logging.info("Restrain protein backbone...")
        add_backbone_posres(system, pdb, 100.0, periodic_boundaries=using_pbc)                
        simulation.context.reinitialize(preserveState=True) # reinitialize context with additional force

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
    simulation.context.reinitialize(preserveState=True) # reinitialize context with additional force

    simulation.step(total_n_equilibriation_steps)

    ##############
    # remove restraint if there is one
    ##############
    if restrain_backbone:
        logging.info("Removing backbone restraint...")
        simulation.context.setParameter('k', 0.0) # remove restraint
    
    ##############
    # run production
    ##############
    logging.info(f"Running production for {total_n_steps} steps...")
    simulation.step(total_n_steps)

if __name__=="__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--pdb_file', type=str, required=True, help='PDB file to simulate.')
    parser.add_argument('--input_dir',  type=str, required=True, help='directory containing pdb file(s).')
    parser.add_argument("--output_dir", default="outputs", type=str, required=False, help="output directory")
    parser.add_argument("--prefix",     default="", type=str, required=False, help="prefix for output files")
    parser.add_argument("--device",     default="cpu", type=str,choices=["cpu","cuda"], required=False, help="device to run on")
    parser.add_argument("--slurm_id",   default="", type=str, required=False, help="slurm id (for logging)")
    parser.add_argument("--E_threshold",default=None,type=int,required=False, help="energy threshold in kJ/(mole * nm) for minimizeEnergy().")

    args = parser.parse_args()
    pdb_file    = args.pdb_file
    input_dir   = args.input_dir
    output_dir  = args.output_dir
    prefix      = args.prefix
    slurm_id    = args.slurm_id
    energy_threshold = args.E_threshold
    device      = args.device
    device      = device.upper()
    if prefix != "":
        prefix += "_"
    if slurm_id != "":
        slurm_id += "_"

    if energy_threshold is not None:
        if energy_threshold <0:
            raise ValueError(f"E_threshold must be positive integer. Got {energy_threshold}")
    
    if input_dir[-1] == "/":
        input_dir = input_dir[:-1]
    pdb_fpath = input_dir + "/" + pdb_file.split("/")[-1]

    starpep_id = "_".join(pdb_file.split("/")[-1].split("_")[:2]) + "_"
    output_dir = make_dir_for_starpepid(starpep_id.strip("_"), input_dir)

    today = datetime.datetime.now()
    logfilename = f"{starpep_id}{slurm_id}{prefix}sim.log"
    logging.basicConfig(
        filename=output_dir+logfilename,
        filemode='a',
        format='%(asctime)s,%(msecs)d %(name)s %(levelname)s %(message)s',
        datefmt='%H:%M:%S',
        level=logging.INFO
    )
    logging.info(f"=========Job ID was {slurm_id}============")
    logging.info(f"output directory = {output_dir}")
    logging.info(f"input arguments =\n{vars(args)}")
    
    # check if simulation already complete
    if os.path.exists(output_dir + "simulation_complete.txt"):
        logging.info("Simulation already complete. Exiting...")
        sys.exit(0)
    
    pdb = PDBFile(pdb_fpath)

    logging.info(f"using pdb file: {pdb_file}")
    logging.info(f"num residues in pdb = {pdb.topology.getNumResidues()}")
    
    # energy_threshold=50
    # logging.warning(f"HARD CODED energy_threshold=50")

    params = {}
    params["coordinates_output_file"] = f"{output_dir}{starpep_id}{slurm_id}{prefix}coordinates.pdb" 
    params["state_output_file"]       = f"{output_dir}{starpep_id}{slurm_id}{prefix}state.csv"
    params["energy_threshold"]        = energy_threshold
    ##############################
    # setup and run simulation
    ##############################
    t0 = time.time()
    run_simulation(pdb, params=params)
    logging.info(f"Simulation took {round(time.time() - t0, 4)}s")
    os.mknod(output_dir + "simulation_complete.txt")
