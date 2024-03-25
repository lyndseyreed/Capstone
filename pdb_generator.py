##Using ESMFold API, input the sequences extracted from sequence_parser and obtain a pdb file with predictied 3D structure##

##Input - list of aa sequences with canonical as the first sequence in the list
    ##Use ESMFold API to predict the structures of the sequences
##Output - pdb file to input into 3D modeller

# pip install --upgrade biopython transformers accelerate

import os
import torch
from transformers import AutoTokenizer, EsmForProteinFolding

#Load in model and tokenizer
tokenizer = AutoTokenizer.from_pretrained("facebook/esmfold_v1")
model = EsmForProteinFolding.from_pretrained("facebook/esmfold_v1", low_cpu_mem_usage=True)

#transfering model to GPU
model = model.to()

#reduce chunk size to use less memory
model.trunk.set_chunk_size(64)

##WORKS UP TO HERE###

#tokenize input sequences from file

def tokenize_sequences_from_file(input_folder, output_folder):
    # Initialize a list to store tokenized sequences
    tokenized_sequences = []
    # Tokenize sequences from the file and save each sequence to a separate file
    for filename in os.listdir(input_folder):
        if filename.endswith(".txt"):
            input_file = os.path.join(input_folder, filename)
            with open(input_file, "r") as file:
                for i, line in enumerate(file):
                    sequence = line.strip()  # Remove leading/trailing whitespace
                    tokenized_sequence = tokenizer(sequence, return_tensors="pt", add_special_tokens=False)['input_ids']
                    tokenized_sequences.append(tokenized_sequence)
                    # Define the output file path (e.g., PDB/TokenizedSequences/sequence1.pt)
                    output_file = os.path.join(output_folder, f"{filename}_sequence{i+1}.pt")
                    
                    # Save the tokenized sequence tensor to the output file
                    torch.save(tokenized_sequence, output_file)

                    print(f"Tokenized sequence {i+1} from {filename} saved to: {output_file}")

    return tokenized_sequences

# Path to the input file containing protein sequences
file_path = "Sequences"

# Path to the folder where tokenized sequences will be saved
output_folder = "PDB"

# Tokenize sequences from the file and save them to the output folder
tokenized_sequences = tokenize_sequences_from_file(file_path, output_folder)


#use torch for model outputs
# with torch.no_grad():
#     output = model(tokenized_input)

#convert model output to PDB file - ESMFold repo function

# from transformers.models.esm.openfold_utils.protein import to_pdb, Protein as OFProtein

# from transformers.models.esm.openfold_utils.feats import atom14_to_atom37

# def convert_outputs_to_pdb(outputs):
#     final_atom_positions = atom14_to_atom37(outputs["positions"][-1], outputs)
#     outputs = {k: v.to("cpu").numpy() for k, v in outputs.items()}
#     final_atom_positions = final_atom_positions.cpu().numpy()
#     final_atom_mask = outputs["atom37_atom_exists"]
#     pdbs = []
#     for i in range(outputs["aatype"].shape[0]):
#         aa = outputs["aatype"][i]
#         pred_pos = final_atom_positions[i]
#         mask = final_atom_mask[i]
#         resid = outputs["residue_index"][i] + 1
#         pred = OFProtein(
#             aatype=aa,
#             atom_positions=pred_pos,
#             atom_mask=mask,
#             residue_index=resid,
#             b_factors=outputs["plddt"][i],
#             chain_index=outputs["chain_index"][i] if "chain_index" in outputs else None,
#         )
#         pdbs.append(to_pdb(pred))
#     return pdbs


