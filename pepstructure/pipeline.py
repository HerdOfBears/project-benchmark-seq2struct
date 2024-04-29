import os
import argparse
import logging
import warnings
import time

import numpy as np
import pandas as pd
import MDAnalysis as mda


def compute_rmsd(structure_uni1:mda.Universe, structure_uni2:mda.Universe) -> float:
    """
    computes the RMSD between two MDAnalysis Universes
    """
    rmsd= mda.analysis.rms.rmsd(
        structure_uni1.select_atoms("backbone"), 
        structure_uni2.select_atoms("backbone")
    )
    return rmsd

def compare_2models(model1:str, model2:str) -> pd.DataFrame:
    """
    compares two models by computing the RMSD between them
    """
    choices = ["pepfold4", "alphafold", "esmfold", "omegafold", "rosettafold", "experimental"]
    if model1.lower() not in choices:
        raise ValueError(f"model1 {model1} not recognized. Choose from {choices}")
    if model2.lower() not in choices:
        raise ValueError(f"model2 {model2} not recognized. Choose from {choices}")


    warnings.warn("Using hardcoded path to intersection files.")
    intersection_all_path = "intersection_all.txt"
    intersection_all_except_pepfold = "intersection_of_models_except_pepfold.txt"
    if model1=="pepfold4" or model2=="pepfold4":
        intersection_path = intersection_all_except_pepfold
    else:
        intersection_path = intersection_all_path
    
    with open(intersection_path, "r") as f:
        overlapping_starpep_ids = f.readlines()
    
    ###########################################
    # construct file tree
    ###########################################
    logging.info("constructing file tree")
    rmsd_values = []
    model1_file_tree = {}
    model2_file_tree = {}
    for starpep_id in overlapping_starpep_ids:
        starpep_id = starpep_id.strip()
        
        model1_file_tree[starpep_id] = {}
        model2_file_tree[starpep_id] = {}

        # construct tree of files
        for model in [model1, model2]:
            for model_file in os.listdir(f"preliminary_analysis/{model}/{starpep_id}"):
                if model == model1:
                    if model_file.endswith("pbcFixed.pdb"):
                        model1_file_tree[starpep_id]["trajectory"] = model_file
                    else:
                        model1_file_tree[starpep_id][   "initial"] = model_file
                else:
                    if model_file.endswith("pbcFixed.pdb"):
                        model2_file_tree[starpep_id]["trajectory"] = model_file
                    else:
                        model2_file_tree[starpep_id][   "initial"] = model_file

    ###########################################
    # now for each starpep_id, compute the RMSD
    ###########################################
    options = [("b","b"), ("a","a"), ("b","a"), ("a","b")]
    # if compare_option not in options:
    #     raise ValueError(f"compare_option {compare_option} not recognized. Choose from {options}")
    options_dict = {}
    options_dict["before"] = "initial"
    options_dict["b"     ] = "initial"
    options_dict["after" ] = "trajectory"
    options_dict["a"     ] = "trajectory"
    
    outputs = {}
    outputs["starpep_id"] = []
    for option in options:
        new_name = f"{model1}_{option[0]}_with_{model2}_{option[1]}"
        outputs[new_name] = []

    logging.info(f"computing RMSD of {model1} with {model2}")
    for starpep_id in overlapping_starpep_ids:
        outputs["starpep_id"].append(starpep_id)
        t0 = time.time()
        for option in options:
            model1_choice = model1_file_tree[starpep_id][ options_dict[option[0]] ]
            model2_choice = model2_file_tree[starpep_id][ options_dict[option[1]] ]

            # load the structures
            model1_uni = mda.Universe(f"preliminary_analysis/{model1}/{starpep_id}/{model1_choice}")
            model2_uni = mda.Universe(f"preliminary_analysis/{model2}/{starpep_id}/{model2_choice}")

            # set to last frame (if trajectory)
            if model1_choice == "trajectory":
                model1_uni.trajectory[-1]
            if model2_choice == "trajectory":
                model2_uni.trajectory[-1]

            # compute the RMSD
            rmsd = compute_rmsd(model1_uni, model2_uni)
            rmsd_values.append(rmsd)

            outputs[f"{model1}_{option[0]}_with_{model2}_{option[1]}"].append(rmsd)
        t1 = time.time()
        logging.info(f"starpep_id {starpep_id} took {t1-t0} seconds to compute RMSD of each option")

    logging.info("finished computing RMSD")
    outputs = pd.DataFrame(outputs)
    return outputs


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--model1', type=str, required=True, help='name of  first model')
    parser.add_argument('--model2', type=str, required=True, help='name of second model')

    args = parser.parse_args()
    model1 = args.model1
    model2 = args.model2

    logging.basicConfig(
        filename=f"logs/rmsd_{model1}_compared_to_{model2}.log",
        filemode='a',
        format='%(asctime)s,%(msecs)d %(name)s %(levelname)s %(message)s',
        datefmt='%H:%M:%S',
        level=logging.INFO
    )

    data = compare_2models(model1, model2)
    logging.info(f"outputting results to csv")
    data.to_csv(f"preliminary_analysis/rmsd_{model1}_with_{model2}.csv", index=False)