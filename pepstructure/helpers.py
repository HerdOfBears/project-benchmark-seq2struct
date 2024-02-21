import os

def make_dir_for_starpepid(starpep_id:str, input_dir:str)->str:
    """
    Create a directory for a given starpep_id in the input_dir
    Args:
        starpep_id: str: the starpep id of form starPep_XXXXX
        input_dir: str: the input directory to create the starpep_id directory in
    
    Returns:
        str: the path to the created directory
    """
    if input_dir[-1] != "/":
        input_dir+= "/"
    if not os.path.exists(input_dir):
        raise Exception(f"input directory {input_dir} does not exist")
    starpep_dir = input_dir + starpep_id

    if not os.path.exists(starpep_dir):
        os.makedirs(starpep_dir)

    if starpep_dir[-1] != "/":
        starpep_dir += "/"
    return starpep_dir

def read_fasta(fasta_file):
    # read in sequences
    starpep_ids_to_seq = {}
    is_annotation = True
    with open(fasta_file, 'r') as f:
        lines = f.readlines()
        for line in lines:
            if line[0]==">":
                starpep_id = line[1:].strip()
                is_annotation = True
                starpep_ids_to_seq[starpep_id] = ""
            else:
                starpep_ids_to_seq[starpep_id] += line.strip()

    return starpep_ids_to_seq