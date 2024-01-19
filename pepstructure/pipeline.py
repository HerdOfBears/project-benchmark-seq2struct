
from helpers_pymol import align_pdb_files, get_second_structure_from_pdb_file

import os
import sys
import argparse
import logging
import time
import pickle as pkl


def compute_rmsd_of_model_outputs(model1, model2):
    """
    takes two model names and computes RMSD between each pdb file 
    in the model output directories
    e.g. model1 = "esmfold", model2 = "omegafold"
    this script will then look in directories with same names
    """

    model1_dir = f"outputs/{model1}"
    model2_dir = f"outputs/{model2}"

    logging.warning("ignoring files with pdbfixed in the name")

    # grab all of the starPep IDs
    starpep_ids = []
    for file in os.listdir(model1_dir):
        if file.endswith(".pdb") and "pdbfixed" not in file:
            starpep_ids.append("starPep_" + file.split("_")[1])
    
    # sort the starpep ids
    starpep_ids.sort()
    data = {}
    data["starpep_id"] = []
    data[f"{model1}_compared_to_{model2}_rmsd"] = []
    data[f"{model1}_compared_to_{model2}_secondary_structure_mismatches"] = []
    t0 = time.time()
    for starpep_id in starpep_ids:
        file1 = f"{model1_dir}/{starpep_id}_{model1}_prediction.pdb"
        file2 = f"{model2_dir}/{starpep_id}_{model2}_prediction.pdb"

        # perform alignment
        rmsd, n_aligned_residues, raw_alignment_score = align_pdb_files(file1, file2)

        # get secondary structure of each
        second_struct1 = get_second_structure_from_pdb_file(file1)
        second_struct2 = get_second_structure_from_pdb_file(file2)

        # compute the difference b/w secondary structures
        n_mismatches = 0
        for i,c in enumerate(second_struct1):
            if c != second_struct2[i]:
                n_mismatches += 1
        
        # add to data
        data["starpep_id"].append(starpep_id)
        data[f"{model1}_compared_to_{model2}_rmsd"].append(rmsd)
        data[f"{model1}_compared_to_{model2}_n_aligned_residues"].append(n_aligned_residues)
        data[f"{model1}_compared_to_{model2}_raw_alignment_score"].append(raw_alignment_score)
        data[f"{model1}_compared_to_{model2}_secondary_structure_mismatches"].append(n_mismatches)
    
    logging.info(f"took {round(time.time() - t0,4)}s to process") 

    return data

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--model1', type=str, required=True, help='name of  first model')
    parser.add_argument('--model2', type=str, required=True, help='name of second model')

    args = parser.parse_args()
    model1 = args.model1
    model2 = args.model2

    if not os.path.exists("outputs/"):
        raise FileNotFoundError("outputs/ directory does not exist.")

    logging.basicConfig(
        filename=f"logs/data_process_{model1}_compared_to_{model2}.log",
        filemode='a',
        format='%(asctime)s,%(msecs)d %(name)s %(levelname)s %(message)s',
        datefmt='%H:%M:%S',
        level=logging.INFO
    )

    data = compute_rmsd_of_model_outputs(model1, model2)
    with open(f"outputs/{model1}_compared_to_{model2}.pkl", "wb") as f:
        pkl.dump(data, f)
    logging.info(f"saved data to outputs/{model1}_compared_to_{model2}.pkl")