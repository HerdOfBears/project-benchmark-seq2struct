import os

def count_pdb_files_in_subfolders(directory_path):
    pdb_file_count = 0

    # Walk through directory
    for root, dirs, files in os.walk(directory_path):
        # Count the pdb files in each directory
        pdb_file_count += sum(file.endswith('.pdb') for file in files)

    return pdb_file_count

# Specify the top-level directory you want to search
top_level_directory = '/media/samith/My Passport/output'  # Replace with the path to your directory
number_of_pdb_files = count_pdb_files_in_subfolders(top_level_directory)

print(f"There are {number_of_pdb_files} PDB files in the directory and its subfolders.")
