from __future__ import absolute_import, division, print_function
from .due import due, Doi

import numpy as np
import torch
import time
import logging
import os
import sys


from transformers import AutoTokenizer, EsmForProteinFolding

__all__ = [
    "get_sequences_from_fasta", 
    "get_model_esm", 
    "sequence_to_pdb_esm_batch", 
    "write_pdb_outputs_"
]


# Use duecredit (duecredit.org) to provide a citation to relevant work to
# be cited. This does nothing, unless the user has duecredit installed,
# And calls this with duecredit (as in `python -m duecredit script.py`):
due.cite(Doi("10.1167/13.9.30"),
         description="Template project for small scientific Python projects",
         tags=["reference-implementation"],
         path='')


def get_sequences_from_fasta(fasta_file):
    """
    Reads a fasta file and returns a dict of sequences
    Assumes that information after ">" is relevant.

    inputs:
        fasta_file: path to fasta file
    outputs:
        starpepid_to_sequence: dict of starpep ids and their sequences 
    """
    logging.info("Reading fasta file using readlines...")
    with open(fasta_file, "r") as f:
        lines = f.readlines()
    
    starpepid_to_sequence = {}
    for line in lines:
        if line.startswith(">"):
            starpep_id = line[1:].strip()
            starpepid_to_sequence[starpep_id] = ""
        else:
            starpepid_to_sequence[starpep_id] += line.strip()
    
    return starpepid_to_sequence

def get_model_esm(model_path):
    """
    Loads the ESM model from the specified path
    
    inputs:
        model_path: path to the model directory of parameters and tokenizer (from huggingface)
    outputs:
        model: ESM model
        tokenizer: tokenizer for the model 
    """
    # model_path = "./pepstructure/esmfold_v1/"
    if not os.path.isdir( os.path.join(os.getcwd(),model_path) ):
        raise Exception("""Model path does not exist: {}\n
                        current directory is: {}""".format(model_path, os.getcwd())
        )
    logging.info("Loading model...")
    model = EsmForProteinFolding.from_pretrained(
                        model_path,
                        local_files_only=True
    )
    tokenizer = AutoTokenizer.from_pretrained(
                        model_path,
                        local_files_only=True
    )
    return model, tokenizer

def sequence_to_pdb_esm_batch(sequences:dict[str], model):
    logging.info("Running batch ESMFold inference...")
    if not isinstance(sequences, dict):
        raise TypeError("sequences must be a dict[starpep_id]=sequence")

    pdb_outputs = {}
    count = 0
    total = len(sequences)
    t00 = time.time()
    t0 = t00
    with torch.no_grad():
        for starpep_id, sequence in sequences.items():
            if count % 50 ==0:
                logging.info(f"Processed {count} sequences of {total}. Time per 50 = {time.time()-t0}s")
                t0 = time.time()
            outputs = model.infer_pdb(sequence)
            pdb_outputs[starpep_id] = outputs
    tf = time.time()
    logging.info(f"Total inference time: {tf-t00}s for {len(sequences)} sequences")
    
    return pdb_outputs

def write_pdb_outputs_(pdb_outputs:dict[str], model_name:str, output_dir:str):
    """
    Writes the pdb outputs to the specified output directory
    Operates 'in-place'
    """
    if not os.path.isdir(output_dir):
        raise Exception("Output directory does not exist: {}".format(output_dir))
    
    for starpep_id, pdb_output in pdb_outputs.items():
        with open(os.path.join(output_dir, starpep_id+ "_" + model_name +"_prediction.pdb"), "w") as f:
            f.write(pdb_output)

