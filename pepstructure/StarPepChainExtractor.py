from Bio import SeqIO
from Bio.PDB import PDBParser, PPBuilder, PDBIO, Select
from Bio import pairwise2
import os

class ChainSelect(Select):
    def __init__(self, chain_letter):
        self.chain_letter = chain_letter

    def accept_chain(self, chain):
        return chain.get_id() == self.chain_letter

def extract_chain_sequences(pdb_file):
    parser = PDBParser()
    structure = parser.get_structure("PDB_structure", pdb_file)
    ppb = PPBuilder()
    chain_sequences = {}

    for chain in structure.get_chains():
        sequence = ""
        for pp in ppb.build_peptides(chain):
            sequence += str(pp.get_sequence())
        chain_sequences[chain.id] = sequence

    return chain_sequences

def align_sequences(pdb_sequences, fasta_sequence):
    best_match = {'score': -1, 'chain': None}

    for chain_id, chain_seq in pdb_sequences.items():
        alignments = pairwise2.align.globalxx(chain_seq, str(fasta_sequence))
        for alignment in alignments:
            if alignment.score > best_match['score']:
                best_match['score'] = alignment.score
                best_match['chain'] = chain_id

    return best_match

def save_chain(pdb_file, chain_id, output_file):
    parser = PDBParser()
    structure = parser.get_structure("PDB_structure", pdb_file)
    model = structure[0]  # Assuming we want the first model
    chain = model[chain_id]

    io = PDBIO()
    io.set_structure(structure)
    io.save(output_file, ChainSelect(chain_id))

def process_folder(folder_path):
    fasta_file_path = None
    pdb_files = []

    for file in os.listdir(folder_path):
        if file.endswith('.fasta'):
            fasta_file_path = os.path.join(folder_path, file)
        elif file.endswith('.pdb'):
            pdb_files.append(os.path.join(folder_path, file))

    if not fasta_file_path or not pdb_files:
        print(f"No FASTA or PDB files found in folder: {folder_path}")
        return

    fasta_sequence = SeqIO.read(fasta_file_path, "fasta").seq

    for pdb_file in pdb_files:
        pdb_sequences = extract_chain_sequences(pdb_file)
        best_match = align_sequences(pdb_sequences, fasta_sequence)

        if best_match['chain']:
            pdb_id = os.path.splitext(os.path.basename(pdb_file))[0]
            output_file = f"{pdb_id}_${best_match['chain']}.pdb"
            output_file_path = os.path.join(folder_path, output_file)
            print(f"Saving chain {best_match['chain']} of {pdb_id} to {output_file}")
            save_chain(pdb_file, best_match['chain'], output_file_path)

# The top-level directory containing all the StarPep_* folders
top_level_directory = "/media/samith/My Passport/output"

# Iterate over each folder in the top-level directory and process it
for folder_name in os.listdir(top_level_directory):
    folder_path = os.path.join(top_level_directory, folder_name)
    if os.path.isdir(folder_path) and folder_name.startswith("starPep_"):
        process_folder(folder_path)

print("All folders processed.")
