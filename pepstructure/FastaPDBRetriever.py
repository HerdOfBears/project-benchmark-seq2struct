import pandas as pd
import requests
import os
from pathlib import Path
import re

def read_metadata(csv_file):
    df = pd.read_csv(csv_file)
    metadata_dict = {}
    for _, row in df.iterrows():
        starpep_id = row['Peptide']
        metadata = row['Metadata']
        pdb_matches = re.findall(r'PDB: (\w+)', metadata)
        uniprot_matches = re.findall(r'UniProtKB: (\w+)', metadata)
        metadata_dict[starpep_id] = {
            'UniProt_ID': uniprot_matches[0] if uniprot_matches else None,
            'PDB_ID': pdb_matches if pdb_matches else []
        }
    return metadata_dict

def download_pdb(pdb_id, output_directory):
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    response = requests.get(url)
    if response.status_code == 200:
        pdb_file_path = Path(output_directory) / f"{pdb_id}.pdb"
        pdb_file_path.write_text(response.text)
    else:
        print(f"Failed to download PDB file for {pdb_id}")

def copy_fasta_file(fasta_file, output_directory):
    fasta_content = Path(fasta_file).read_text()
    fasta_file_path = Path(output_directory) / Path(fasta_file).name
    fasta_file_path.write_text(fasta_content)

def process_metadata(metadata_file, fasta_directory, output_directory):
    peptide_to_ids = read_metadata(metadata_file)
    for fasta_file in Path(fasta_directory).glob("*.fasta"):
        starpep_id = fasta_file.stem
        ids_info = peptide_to_ids.get(starpep_id)
        if ids_info:
            starpep_output_directory = Path(output_directory) / starpep_id
            starpep_output_directory.mkdir(parents=True, exist_ok=True)
            copy_fasta_file(str(fasta_file), str(starpep_output_directory))
            for pdb_id in ids_info['PDB_ID']:
                download_pdb(pdb_id, str(starpep_output_directory))
        else:
            print(f"No metadata found for starPep ID: {starpep_id}")

# Paths
fasta_directory = "/home/samith/Downloads/Fasta"  # Update this path
metadata_file = "/home/samith/Downloads/peptides_lt60_metadata.csv"  # Update this path if necessary
output_directory = "/media/samith/My Passport/output"  # Replace with your desired output directory

process_metadata(metadata_file, fasta_directory, output_directory)
