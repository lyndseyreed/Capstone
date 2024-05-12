##Using ESMFold API, input the sequences extracted from sequence_parser and obtain a pdb file with predictied 3D structure##

##Input - folder containing multiple files (Sequences). The expected format within the file is a
            #list of aa sequences with the canonical sequence as the first sequence in the list
            #and all other sequences in the file as the novel that correspond to the same gene

##Use - ESMFold API to predict the atomic locations for the 3D structures of the sequences

##Output - pdb files in the same folder format (Tokens) - one gene per file- to input into 3D modeller API
# Note: must pass this output through remove_binary_to_pdb for functionality in model_generator



# pip install --upgrade biopython transformers accelerate
import os
import torch
from transformers import AutoTokenizer, EsmForProteinFolding
from transformers.models.esm.openfold_utils.protein import to_pdb, Protein as OFProtein
from transformers.models.esm.openfold_utils.feats import atom14_to_atom37


#Load in model and tokenizer
tokenizer = AutoTokenizer.from_pretrained("facebook/esmfold_v1")
model = EsmForProteinFolding.from_pretrained("facebook/esmfold_v1", low_cpu_mem_usage=True)
model = model.to() #transfering model to GPU
model.trunk.set_chunk_size(64) #reduce chunk size to use less memory


#convert model output to PDB file - ESMFold repo function
def convert_to_pdb(outputs):
    """
    This function is sourced fro the ESMFold repo as a function to convert the model output
    into a PDB file for 3D rendering. It is utilized in the function tokenize_sequences. 
    """
    final_atom_positions = atom14_to_atom37(outputs["positions"][-1], outputs)
    outputs = {k: v.to("cpu").numpy() for k, v in outputs.items()}
    final_atom_positions = final_atom_positions.cpu().numpy()
    final_atom_mask = outputs["atom37_atom_exists"]
    pdbs = []
    for i in range(outputs["aatype"].shape[0]):
        aa = outputs["aatype"][i]
        pred_pos = final_atom_positions[i]
        mask = final_atom_mask[i]
        resid = outputs["residue_index"][i] + 1
        pred = OFProtein(
            aatype=aa,
            atom_positions=pred_pos,
            atom_mask=mask,
            residue_index=resid,
            b_factors=outputs["plddt"][i],
            chain_index=outputs["chain_index"][i] if "chain_index" in outputs else None,
        )
        pdbs.append(to_pdb(pred))
    return pdbs

        
def tokenize_sequences(input_folder, output_folder):
    """
    This function takes in a folder of sequences (Sequences), tokenizes them, applies the predictive model and outputs 
    the tokenized predicted model files to a folder (Tokens). It is used in conjunction with the function convert_to_pdb to 
    generate the full PDB formatted files.
    """
    for filename in os.listdir(input_folder):
        input_file = os.path.join(input_folder, filename)
        gene_name = os.path.splitext(filename)[0]
        gene_folder = os.path.join(output_folder, gene_name)
        os.makedirs(gene_folder, exist_ok=True)
        with open(input_file, "r") as file:
            for i, line in enumerate(file):
                sequence = line.strip()
                tokenized_sequence = tokenizer(sequence, return_tensors="pt", add_special_tokens=False)['input_ids']
                with torch.no_grad():
                    token_model = model(tokenized_sequence)
                pdb = convert_to_pdb(token_model)
                output_file = os.path.join(gene_folder, f"{gene_name}_sequence{i+1}.pdb")
                torch.save(pdb, output_file)




####USER INPUT AREA####
#Will be added to main eventually#


# Path to the input file containing protein sequences
input_folder = "Sequences"
# Path to the folder where tokenized sequences will be saved
output_folder = "Tokens"
# Tokenize sequences from the file and save them to the output folder
tokenized_sequences = tokenize_sequences(input_folder, output_folder)