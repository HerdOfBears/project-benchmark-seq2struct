import pepstructure
import os
import logging

def main():
    pwd = os.getcwd()
    if pwd.split("/")[-1] == "scripts":
        os.chdir("..")

    input_fasta = "inputs/peptides_lteq100_pdbAvailable.fasta"
    output_dir  = "outputs/esmfold/"
    model_path  = "pepstructure/esmfold_v1/"

    #################
    # Read fasta file
    #################
    starpepid_to_sequence = pepstructure.read_fasta(input_fasta)

    #################
    # Load model
    #################
    model, _ = pepstructure.get_model_esm(model_path)

    #################
    # Run inference
    #################
    pdb_outputs = pepstructure.sequence_to_pdb_esm_batch(starpepid_to_sequence, model)

    #################
    # Write outputs
    #################
    pepstructure.write_pdb_outputs(pdb_outputs, "esmfold", output_dir)
    


if __name__ == '__main__':
    logname = "logs/esmfold_predictions.log"
    logging.basicConfig(
            filename=logname,
            filemode='a',
            format='%(asctime)s,%(msecs)d %(name)s %(levelname)s %(message)s',
            datefmt='%H:%M:%S',
            level=logging.INFO
    )
    main()
