import pymol
import Bio.PDB as bPDB

import os

def align_pdb_files(file1, file2):
    """
    takes two pdb files and aligns them using pymol's align function
    returns the rmsd and number of aligned atoms
    """
    for file in [file1, file2]:
        if not os.path.exists(file):
            raise FileNotFoundError(f"{file} does not exist")
        if not file.endswith(".pdb"):
            raise ValueError(f"{file} is not a pdb file")

    # load the two pdb files
    pdb1 = pymol.cmd.load(file1)
    pdb2 = pymol.cmd.load(file2)

    allobjects = pymol.cmd.get_object_list('all')
    # pymol's align returns a list of 7 items:
    # 0: RMSD
    # 1: number of aligned atoms
    # 2: number of refinement ccles
    # 3: RMSD before refinement
    # 4: Number of aligned atoms before refinement
    # 5: raw alignment score
    # 6: number of residues aligned
    alignment = pymol.cmd.align(allobjects[0], allobjects[1], cycles=5) 
    rmsd                = alignment[0]
    raw_alignment_score = alignment[5]
    n_aligned_residues  = alignment[6]
    
    return rmsd, n_aligned_residues, raw_alignment_score

def get_dssp_from_pdb_file(file):
    """
    takes a pdb file and returns a dssp object
    """
    if not os.path.exists(file):
        raise FileNotFoundError(f"{file} does not exist")
    if not file.endswith(".pdb"):
        raise ValueError(f"{file} is not a pdb file")

    # load the pdb file into a bio.PDB structure

    p = bPDB.PDBParser()
    structure = p.get_structure("", file)
    model = structure[0]
    dssp = bPDB.DSSP(model, file, dssp="mkdssp")

    keys = list(dssp.keys())
    second_structs = []
    for key in keys:
        second_structs.append( dssp[key][2] )
    return second_structs

if __name__ == "__main__":
    file1 = "outputs/esmfold/starPep_44878_esmfold_prediction.pdb"
    file2 = "outputs/omegafold/starPep_44878_omegafold_prediction.pdb"
    
    rmsd, n_aligned_residues, raw_alignment_score = align_pdb_files(file1, file2)
    print(f"rmsd = {rmsd}")
    print(f"n_aligned_residues = {n_aligned_residues}")
    print(f"raw_alignment_score = {raw_alignment_score}")

    dssp1 = get_dssp_from_pdb_file(file1)
    dssp2 = get_dssp_from_pdb_file(file2)
    print(f"dssp1 = {''.join(dssp1)}")
    print(f"dssp2 = {''.join(dssp2)}")

    print("file has no default behavior")