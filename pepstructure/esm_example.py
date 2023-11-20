
from transformers import AutoTokenizer, EsmForProteinFolding
import numpy as np
import torch

import time
import logging
import os
import sys

def sequence_to_pdb_esm(sequence:str):
    model_path = "./esmfold_v1/"
    if not os.path.isdir(model_path):
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

    inputs = tokenizer([sequence], return_tensors="pt", add_special_tokens=False)

    t0 = time.time()
    with torch.no_grad():
        outputs = model.infer_pdb(sequence)
    tf = time.time()
    logging.info("Inference time: {}".format(tf-t0))
    with open("esm_output.pdb", "w") as f:
        f.write(outputs)
    

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    
    cmd_line_args = sys.argv
    if len(cmd_line_args)==1:
        raise Exception("""No sequence provided as command line argument\n
                        Usage: python esm_example.py <sequence>""")
    
    sequence = sys.argv[1]
    if "<" in sequence or ">" in sequence:
        raise ValueError("Invalid sequence provided as command line argument")

    sequence_to_pdb_esm(sequence)
    logging.info("Done! Results should be in esm_output.pdb")