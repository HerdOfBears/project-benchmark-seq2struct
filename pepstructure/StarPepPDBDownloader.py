import csv
import requests
import os
import shutil

# Paths
fasta_directory = "/home/samith/Downloads/Fasta"  # Update this path
metadata_file = "/home/samith/Downloads/peptides_lt60_metadata.csv"  # Update this path if necessary
base_output_directory = "/media/samith/My Passport/output"  # Replace with your desired output directory

def read_metadata(csv_path):
    starpep_to_pdb = {}
    with open(csv_path, mode='r', newline='', encoding='utf-8') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            peptide = row['Peptide'].strip()
            metadata_entries = row['Metadata'].split(';')
            pdb_ids = [entry.split(':')[-1].strip() for entry in metadata_entries if 'PDB:' in entry]

            if peptide not in starpep_to_pdb:
                starpep_to_pdb[peptide] = set()
            starpep_to_pdb[peptide].update(pdb_ids)

    return starpep_to_pdb

def download_pdb_files(pdb_ids, output_path):
    for pdb_id in pdb_ids:
        url = f'https://files.rcsb.org/download/{pdb_id}.pdb'
        response = requests.get(url)
        if response.status_code == 200:
            file_path = os.path.join(output_path, f'{pdb_id}.pdb')
            with open(file_path, 'wb') as file:
                file.write(response.content)
            print(f"Downloaded {pdb_id}.pdb")
        else:
            print(f"Failed to download {pdb_id}. Status code: {response.status_code}")

def copy_fasta_file(starpep_id, fasta_directory, output_path):
    fasta_file_name = f"{starpep_id}.fasta"
    fasta_file_path = os.path.join(fasta_directory, fasta_file_name)
    if os.path.exists(fasta_file_path):
        shutil.copy(fasta_file_path, output_path)
        print(f"Copied FASTA file for {starpep_id}")
    else:
        print(f"FASTA file for {starpep_id} not found in {fasta_directory}")

def main():
    peptide_to_pdb = read_metadata(metadata_file)

    for starpep_id, pdb_ids in peptide_to_pdb.items():
        starpep_output_path = os.path.join(base_output_directory, starpep_id)
        os.makedirs(starpep_output_path, exist_ok=True)

        copy_fasta_file(starpep_id, fasta_directory, starpep_output_path)
        download_pdb_files(pdb_ids, starpep_output_path)

if __name__ == "__main__":
    main()
