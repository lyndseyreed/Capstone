## Build a sequence comparison function in order to comapare the target sequences to a databse and find a tempalte sequence that is at least 30% similar##
##Assume the input is a FASTA file with a full sequence##
import Bio
from Bio import SeqIO

def get_sequences(fasta_file):
    """
    This function takes in a FASTA file and
    outputs a list of the sequences converted to a string
    """
    sequences = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences.append(str(record.seq))
    return sequences

#Test
novel_data = get_sequences("novel_input.fasta")
cannonical_data = get_sequences("cannonical_input.fasta")

print(novel_data)
print(cannonical_data)
