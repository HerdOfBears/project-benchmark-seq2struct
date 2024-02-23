import MDAnalysis as mda
import numpy as np

from MDAnalysis.analysis import rms

import os
import argparse

def compute_rmsd(uni1, uni2, sel="backbone"):
    """
    Compute the RMSD between two MDAnalysis universes.

    Parameters
    ----------
    uni1 : MDAnalysis.Universe
        The first universe.
    uni2 : MDAnalysis.Universe
        The second universe.
    sel : str
        The selection string to use for the RMSD calculation.

    Returns
    -------
    rmsd : float
        The RMSD between the two universes.
    """
    # align the two structures
    ref = uni1.select_atoms(sel)
    mobile = uni2.select_atoms(sel)
    rmsd = rms.rmsd(mobile.positions, ref.positions)
    return rmsd

def get_model_most_similar_to_others(starPep_dir):
    if starPep_dir[-1] != "/":
        best_structure_dir = starPep_dir + "/Best_Structure/"
    else:
        best_structure_dir = starPep_dir + "Best_Structure/"
    
    models = []
    for file in os.listdir(best_structure_dir):
        if file.startswith("model"):
            models.append(file)
    
    rmsd_matrix = np.zeros((len(models), len(models)))
    for i, model1 in enumerate(models):
        for j, model2 in enumerate(models):
            if i == j:
                continue
            uni1 = mda.Universe(best_structure_dir + model1)
            uni2 = mda.Universe(best_structure_dir + model2)
            rmsd = compute_rmsd(uni1, uni2)
            rmsd_matrix[i, j] = rmsd
    
    most_similar_model_mean   = np.argmin(np.mean(  rmsd_matrix, axis=1))
    most_similar_model_median = np.argmin(np.median(rmsd_matrix, axis=1))

    return models[most_similar_model_mean], models[most_similar_model_median]


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--input_dir", type=str, required=True)

    args = parser.parse_args()
    input_dir = args.input_dir

    if input_dir[-1] != "/":
        input_dir += "/"

    starpep_directories = []
    for directory in os.listdir(input_dir):
        if os.path.isdir(input_dir + directory):
            starpep_directories.append(directory)

        
    for starpep_dir in starpep_directories:
        print(starpep_dir)
        best_structure_mean, best_structure_median = get_model_most_similar_to_others(input_dir + starpep_dir)
        print(best_structure_mean, best_structure_median)
        print("\n")
        break

if __name__ == "__main__":
    main()