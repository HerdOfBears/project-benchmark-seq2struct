import os
from Bio import SeqIO

def create_individual_fasta_files(fasta_file, output_directory):
    # Create the output directory if it doesn't exist
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    for record in SeqIO.parse(fasta_file, "fasta"):
        file_name = f"{output_directory}/{record.id}.fasta"
        with open(file_name, "w") as output_handle:
            SeqIO.write(record, output_handle, "fasta")

# Using the provided FASTA file
fasta_file = "/home/samith/Downloads/peptides_lengthLT60_pdbAvailable.fasta"
output_directory = "/home/samith/Downloads/Fasta"  # Replace with your desired output directory

create_individual_fasta_files(fasta_file, output_directory)
