import os

def make_command(starpep_id,mdrun_input_dir):
    """
    Create a command to run the simulation
    """        
    simulation_script = "pepstructure/simulate_protein_in_water.py"
    model_name = mdrun_input_dir.split("/")[1].strip("/")
    
    if model_name=="experimental":
        for file in os.listdir(mdrun_input_dir):
            if file.endswith(".pdb"):
                file_starpep_id = file.split("_")[0] + "_" + file.split("_")[1]
                if starpep_id==file_starpep_id:
                    pdb_file = file
    else:
        if starpep_id[0]=="s":
            pdb_file = f"{starpep_id}_{model_name}_prediction_pdbfixed.pdb"
        else:
            pdb_file = f"starPep_{starpep_id}_{model_name}_prediction_pdbfixed.pdb"
    
    
    slurm_id = "$SLURM_JOB_ID"
    prefix = ""

    if prefix != "":
        command = f"python {simulation_script} --input_dir {mdrun_input_dir} --pdb_file {pdb_file} --slurm_id {slurm_id} --prefix {prefix}"
    else:
        command = f"python {simulation_script} --input_dir {mdrun_input_dir} --pdb_file {pdb_file} --slurm_id {slurm_id}"
    return command
#python scripts/simulate_protein_in_water.py --input_dir outputs/esmfold/ --pdb_file starPep_00218_esmfold_prediction_pdbfixed.pdb --output_dir outputs/gnu_test/ --slurm_id $SLURM_JOB_ID --prefix gnu_test1

def sort_starpep_ids_by_length(starpep_ids_to_seq):
    starpep_ids_to_length = {}
    for starpep_id, seq in starpep_ids_to_seq.items():
        starpep_ids_to_length[starpep_id] = len(seq)
    
    sorted_starpep_ids = sorted(starpep_ids_to_length, key=starpep_ids_to_length.get)
    return sorted_starpep_ids

def get_starpep_ids_in_dir(input_dir):
    starpep_ids = []
    for file in os.listdir(input_dir):
        if file.endswith(".pdb"):
            starpep_id = file.split("_")[0] + "_" + file.split("_")[1]
            starpep_ids.append(starpep_id)
    return starpep_ids

def main():
    mdrun_input_dir = "outputs/experimental/"
    model_name = mdrun_input_dir.split("/")[1].strip("/")

    
    starpep_ids_to_seq = {}
    is_annotation = True
    with open("inputs/peptides_lteq100_pdbAvailable.fasta", 'r') as f:
        lines = f.readlines()
        for line in lines:
            if line[0]==">":
                starpep_id = line[1:].strip()
                is_annotation = True
                starpep_ids_to_seq[starpep_id] = ""
            else:
                starpep_ids_to_seq[starpep_id] += line.strip()
    print(f"N of starpep_ids: {len(starpep_ids_to_seq)}")
    max_len = 0
    min_len = 1000
    for starpep_id, seq in starpep_ids_to_seq.items():
        if len(seq) > max_len:
            max_len = len(seq)
        if len(seq) < min_len:
            min_len = len(seq)
    print(f"max_len: {max_len}, min_len: {min_len}")
    
    # get all starpep_ids that are available in the directory (relevant for pepfold)
    available_starpep_ids = get_starpep_ids_in_dir(mdrun_input_dir)

    # sort all starpep_ids by length
    sorted_starpep_ids = sort_starpep_ids_by_length(starpep_ids_to_seq)
    for starpep_id in sorted_starpep_ids:
        print(f"{starpep_id}: {len(starpep_ids_to_seq[starpep_id])}")

    # use sorted starpep_ids to sort the available starpep_ids
    sorted_available_starpep_ids = []
    for starpep_id in sorted_starpep_ids:
        if starpep_id in available_starpep_ids:
            sorted_available_starpep_ids.append(starpep_id)

    # iterate over all starpep_ids 
    # grab four at a time and make a command for each one
    if model_name=="esmfold":
        farm_name = "meta_farm_cases"
    else:
        farm_name = f"{model_name}_farm_cases"
    

    for i in range(0, len(sorted_available_starpep_ids), 4):
        idx = i//4

    
        starpep_ids = sorted_available_starpep_ids[i:i+4]
        cmds = []
        for starpep_id in starpep_ids:
            cmd = make_command(starpep_id, mdrun_input_dir)
            cmds.append(cmd)
            # os.mknod(f"scripts/meta_farm_cases/parallel_cmds/case{idx}.txt")
            # with open(f"scripts/{farm_name}/parallel_cmds/case{idx}.txt", "a") as f:
            #     f.write(cmd1 + "\n")
            #     f.write(cmd2 + "\n")
            #     f.write(cmd3 + "\n")
            #     f.write(cmd4)

        with open(f"scripts/{farm_name}/parallel_cmds/case{idx}.txt", "a") as f:
                for cmd in cmds:
                    f.write(cmd + "\n")

    output_file_name = f"scripts/{farm_name}/table.dat"

    meta_farm_case = "parallel -j $JOBS_PER_NODE --joblog $LOG_FILE --workdir $WDIR <"
    with open(output_file_name, "a") as f:
        for i in range(0, len(sorted_available_starpep_ids), 4):
            idx = i//4
            f.write(f"{meta_farm_case} $WDIR/scripts/{farm_name}/parallel_cmds/case{idx}.txt\n")

if __name__=="__main__":
    main()