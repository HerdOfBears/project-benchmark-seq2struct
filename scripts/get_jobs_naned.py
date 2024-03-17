import os 
import argparse


def main(input_dir):

    did_not_finish = []
    # find starPep_id directories, then find if "simulation_complete.txt" 
    # exists in the directory
    for p_ in os.listdir(input_dir):
        if os.path.isdir(p_):
            if p_.startswith("starPep"):
                if os.path.exists(input_dir + p_ + "simulation_complete.txt"):
                    continue
                else:
                    did_not_finish.append(p_)
    return did_not_finish

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_dir',  type=str, required=True, help='directory containing pdb file(s).')
    parser.add_argument(     '--wdir',  type=str, required=False, help='working directory (optional).')    
    args = parser.parse_args()
    input_dir   = args.input_dir
    wdir        = args.wdir
    if wdir:
        os.chdir(wdir)

    did_not_finish = main(input_dir)

    inputted_dir = input_dir.split("/")[-2] if input_dir.endswith("/") else input_dir.split("/")[-1]

    with open(f"{inputted_dir}_did_not_finish.txt", 'w') as f:
        for d in did_not_finish:
            f.write(d + "\n")