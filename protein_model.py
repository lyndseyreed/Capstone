##Using ESMFold API, input the sequences extracted from sequence_parser and obtain a 3D model with corresponding novel and canonical for comparison##

##Input - list of aa sequences with canonical as the first sequence in the list
    ##Use ESMFold API to predict the structures of the sequences
##Output - 3D model that allows visualization of canonical vs novel tertiary model predictions

# pip install --upgrade biopython transformers py3Dmol accelerate

import os
import torch
#from sequence_parser import get_sequences as get_seq


#Load in model and tokenizer
from transformers import AutoTokenizer, EsmForProteinFolding
tokenizer = AutoTokenizer.from_pretrained("facebook/esmfold_v1")
model = EsmForProteinFolding.from_pretrained("facebook/esmfold_v1", low_cpu_mem_usage=True)

#transfering model to GPU
model = model.to()

#Convert language model stem to float16 (improve mem use on GPU)
model.esm = model.esm.half()

#reduce chunk size to use less memory
model.trunk.set_chunk_size(64)

##WORKS UP TO HERE###

#tokenize input sequences from file
def tokenized_sequences(file_path):
    tokenized_sequences = []
    with open(file_path, "r") as file:
        for sequence in file:
            tokenized_sequence = tokenizer([sequence], return_tensors="pt", add_special_tokens=False)['input_ids']
            tokenized_sequences.append(tokenized_sequence)
    return tokenized_sequences

file_path = "Sequences/P40121.fasta"

tokenized_input = tokenized_sequences(file_path)

tokenized_tensor = torch.cat(tokenized_input, dim=0)

#use torch for model outputs
with torch.no_grad():
    output = model(tokenized_tensor)

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

#visualize with py3Dmol
