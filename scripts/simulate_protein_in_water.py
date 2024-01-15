from openmm.app import *
from openmm import *
from openmm.unit import *
from pdbfixer import PDBFixer
from sys import stdout

import time
import logging
import argparse

def run_simulation(pdb):
    """
    requires a PDBFile object
    """
    temperature = 300 # kelvin
    step_size = 0.004 # of a picosecond
    report_every = 1000 # every X steps

    coordinates_output_file = 'output.pdb'
    state_output_file = 'state.csv'

    total_simulation_time = 1*nanoseconds
    total_n_steps = int(total_simulation_time / (step_size*picoseconds))+1
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
    # run (NVT) simulation
    ##############
    simulation.step(total_n_steps)
        

if __name__=="__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--pdb_file', type=str, required=True, help='PDB file to simulate.')
    parser.add_argument('--pdb_dir',  type=str, required=True, help='directory containing pdb file(s).')

    args = parser.parse_args()
    pdb_file = args.pdb_file
    pdb_dir  = args.pdb_dir

    logging.basicConfig(
        out="logs/test_simulate_protein_in_water.log", 
        level=logging.INFO
    )
    if pdb_dir[-1] == "/":
        pdb_dir = pdb_dir[:-1]
    pdb_fpath = pdb_dir + "/" + pdb_file.split("/")[-1]
    
    # pdb = PDBFile("starPep_06810_test_pdbfixer.pdb")
    pdb = PDBFile(pdb_fpath)

    logging.info(f"num residues in pdb = {len(pdb.topology.residues())}")

    logging.info(f"Running simulation on {pdb_file}")
    t0 = time.time()
    run_simulation(pdb)
    logging.info(f"Simulation took {round(time.time() - t0, 4)}s")