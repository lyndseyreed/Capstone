## Extract novel and canonical sequences from an input fasta/text file and save as a file with the gene name as the file name##
##Assumes each file contains one canonical and one or more novel sequences##

#Input: take in folder with one or more files,
    #parse for gene name in the first canonical sequence
    #parse for sequence information and extract unique sequences 
    #save canonical and corresponding novel sequences in a file by canoncial sequence gene name 



#Output: file with the gene name as the file name and the sequences as a list in the file within "Sequences" folder

#pip install --upgrade biopython

import os
import Bio
from Bio import SeqIO


def get_gene_name(file_path):
    """
    This function takes in a file from the input folder and
    outputs the gene name of the first fasta sequence in the file (canoncial sequence gene name)
    if no gene name is present it will return as None and will be handles during folder processing
    """
    with open(file_path, "r") as file:
        record = next(SeqIO.parse(file, "fasta"), None)
        if record:
            return record.id.split("|")[1]
        else:
            return None


def get_sequences(file_path):
    """
    This function takes in a file from the input folder and
    outputs a list of the unique sequences converted to a string
    """
    sequences = []
    with open(file_path, "r") as file:
        for record in SeqIO.parse(file_path, "fasta"):
            sequence = str(record.seq)
            if sequence not in sequences:
                sequences.append(sequence)
    return sequences


def save_seq(sequences, gene_name, count):
    """
    This function takes in a sequence list and the gene name of the canonical sequence and
    outputs those into a file in the "Sequences" folder with the name as the gene name in FASTA format
    """
    output_folder = "Sequences"
    output_file = os.path.join(output_folder, f"{gene_name}.fasta")
    with open(output_file, "w") as file:
        for sequence in sequences:
            file.write(sequence + "\n")

def process_input_folder(input_folder):
    """
    This function iterates through each file within the input folder and 
    outputs a sequence FASTA file for each one
    if the gene name is unknown the file name will be unknown_gene_x (increments up for each unknown gene)
    """
    unknown_gene_count = 0
    for filename in os.listdir(input_folder):
        file_path = os.path.join(input_folder, filename)
        gene_name = get_gene_name(file_path)
        sequences = get_sequences(file_path)
        save_seq(sequences, gene_name, unknown_gene_count)
        if gene_name:
            print(f"Sequences from {filename} saved to sequences/{gene_name}.fasta")
        else:
            print(f"Sequences from {filename} saved to sequences/Unknown_gene_{unknown_gene_count}.fasta")
            unknown_gene_count += 1



####USER INPUT AREA####
#Will be added to main eventually#

# #input folder name:
# input_folder = "Input"

# process_input_folder(input_folder)

