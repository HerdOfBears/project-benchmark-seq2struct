# run this script from project-benchmark-seq2struct/ directory

from pdbfixer import PDBFixer
from openmm.app import *

import os
import sys
import argparse
import logging
import time

# function to fix a given pdb file with pdbfixer, and saves the fixed 
# pdb to a provided output file name. Logs information along the way.
def fix_pdb(pdb_file, output_file):
    logging.info('Fixing pdb file: {}'.format(pdb_file))
    fixer = PDBFixer(filename=pdb_file)
    fixer.findMissingResidues()
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.0)
    PDBFile.writeFile(fixer.topology, fixer.positions, open(output_file, 'w'))
    logging.info('Fixed pdb file saved as: {}'.format(output_file))

def main(pdb_dir, is_experimental=False):

    # check directory exists
    if not os.path.isdir(pdb_dir):
        logging.error(f'Directory {pdb_dir} does not exist. Current dir is {os.getcwd()}')
        sys.exit(1)

    # get all pdb files in the given directory
    pdb_files = []
    for file_ in os.listdir(pdb_dir):
        condn1 = file_.endswith('.pdb')
        condn2 = "fixed" not in file_ 
        if condn1 and condn2:
            if is_experimental:
                if "model" in file_:
                    pdb_files.append(file_)
            else:
                pdb_files.append(file_)
    
    # iterate through all pdb files and fix them
    # SAVES fixed pdbs to same directory
    output_suffix = "_pdbfixed.pdb"
    counter = 0
    t0 = time.time()
    t00 = time.time()
    for pdb_file in pdb_files:
        pdb_fpath = pdb_dir + "/" + pdb_file
        output_file = pdb_dir + "/" + pdb_file.replace('.pdb', output_suffix)
        fix_pdb(pdb_fpath, output_file)
        counter += 1
        if counter % 100 == 0:
            logging.info(f'Fixed {counter} pdb files in {time.time() - t00} seconds.')
            t00 = time.time()
    logging.info(f'Fixed {len(pdb_files)} pdb files in {time.time() - t0} seconds.')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--pdb_dir', type=str, required=True, help='Directory containing pdb files to fix.') 
    parser.add_argument('--experimental', type=str, required=False, choices=["y", "n"], help='If the pdb files are experimental models, use "y".')

    args = parser.parse_args()
    pdb_dir = args.pdb_dir
    experimental = args.experimental

    if experimental is None:
        experimental = "n"
        
    if experimental == "y":
        model_name = "experimental"
    else:
        model_name = pdb_dir.split('/')[-1]        
        if model_name == "":
            model_name = pdb_dir.split('/')[-2]
        elif model_name == " ": # if model_name is empty
            model_name = pdb_dir.split('/')[-2]
    log_name = f"logs/fix_pdb_{model_name}.log"
    logging.basicConfig(filename=log_name, level=logging.INFO)

    if experimental == "y":
        is_experimental = True
    else:
        is_experimental = False
    main(pdb_dir, is_experimental)