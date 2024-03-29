##Using ESMFold API, input the sequences extracted from sequence_parser and obtain a pdb file with predictied 3D structure##

##Input - list of aa sequences with canonical as the first sequence in the list
    ##Use ESMFold API to predict the structures of the sequences
##Output - pdb file to input into 3D modeller

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

##WORKS UP TO HERE###

#convert model output to PDB file - ESMFold repo function
def convert_to_pdb(outputs):
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

#tokenize input sequences from file

def tokenize_sequences(input_folder, output_folder):
    # tokenized_sequences = []
    for filename in os.listdir(input_folder):
        input_file = os.path.join(input_folder, filename)
        gene_name = os.path.splitext(filename)[0]
        gene_folder = os.path.join(output_folder, gene_name)
        os.makedirs(gene_folder, exist_ok=True)
        with open(input_file, "r") as file:
            for i, line in enumerate(file):
                sequence = line.strip()
                tokenized_sequence = tokenizer(sequence, return_tensors="pt", add_special_tokens=False)['input_ids']
                print('token', tokenized_sequence)
                with torch.no_grad():
                    token_model = model(tokenized_sequence)
                # tokenized_sequences.append(token_model)
                print('token_model', token_model)
                pdb = convert_to_pdb(token_model)
                print('pdb', pdb)
                output_file = os.path.join(gene_folder, f"{gene_name}_sequence{i+1}.pdb")
                torch.save(pdb, output_file)
    # return tokenized_sequences


# Path to the input file containing protein sequences
input_folder = "Sequences"
# Path to the folder where tokenized sequences will be saved
output_folder = "Tokens"
# Tokenize sequences from the file and save them to the output folder
tokenized_sequences = tokenize_sequences(input_folder, output_folder)
