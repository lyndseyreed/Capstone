## Removes the lines of binary from the predictive model generated file so that it can be read by py3Dmol/stmol.
## WARNING - is very specific to the first line and last 2 lines of the file input and does not work to identify and remove binary from everything!

##INPUT - Tokens folder of files with the model applied to the amino acid sequences.

##USE - removes the first line and last two lines

##OUTPUT - the file in full PDB format into the PDB folder sorted in files by the gene

import os
import shutil

def remove_binary(input_folder, output_folder):
    """
    This function will take in the PDB output from the prediction_generator and 
    remove the first line and the last two (2) lines that are formatted as binary.
    This will allow the file to be read by the model_generator API for 3D modelling.
    """
    for root, dirs, files in os.walk(input_folder):
        relative_path = os.path.relpath(root, input_folder)
        output_root = os.path.join(output_folder, relative_path)
        os.makedirs(output_root, exist_ok=True)
        for filename in files:
            input_file_path = os.path.join(root, filename)
            output_file_path = os.path.join(output_root, filename)
            with open(input_file_path, "rb") as input_file:
                lines = input_file.readlines()[1:-2]
            with open(output_file_path, "wb") as output_file:
                output_file.writelines(lines)




####USER INPUT AREA####
#Will be added to main eventually#

# input_folder = "Tokens"  
# output_folder = "PDB" 
# remove_binary(input_folder, output_folder)