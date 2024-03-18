# pip install --upgrade bioconductor transformers py3Dmol accelerate

from transformers.utils import send_example_telemetry
send_example_telemetry("protein_folding_notebook", framework="pytorch")

from transformers import AutoTokenizer, EsmForProteinFolding
tokenizer = AutoTokenizer.from_pretrained("facebook/esmfold_v1")
model = EsmForProteinFolding.from_pretrained("facebook/esmfold_v1", low_cpu_mem_usage=True)

model = model.to()
model.esm = model.esm.half()

import torch
model.trunk.set_chunk_size(64)


import Bio
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Seq import Seq
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

novel_data = get_sequences("novel_input.fasta")
cannonical_data = get_sequences("cannonical_input.fasta")


novel_tokenized_input = tokenizer([novel_data], return_tensors="pt", add_special_tokens=False)['input_ids']
cannonical_tokenized_input = tokenizer([cannonical_data], return_tensors="pt", add_special_tokens=False)['input_ids']

novel_tokenized_input = novel_tokenized_input.to()
cannonical_tokenized_input = cannonical_tokenized_input.to()

# with torch.no_grad():
#     novel_output = model( novel_tokenized_input)
#     cannonical_output = model(cannonical_tokenized_input)

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



# novel_pdb = convert_outputs_to_pdb(novel_output)
# cannonical_pdb = convert_outputs_to_pdb(cannonical_output)


# import py3Dmol

# view = py3Dmol.view(js='https://3dmol.org/build/3Dmol.js', width=800, height=400)
# view.addModel("".join(pdb), 'pdb')
# view.setStyle({'model': -1}, {"cartoon": {'color': 'yellow'}})
# view.addModel("".join(pdb_2), 'pdb_2')
# view.setStyle({'model': -1}, {"cartoon": {'color': 'purple'}})

# if torch.max(output['plddt']) <= 1.0:
#     vmin = 0.5
#     vmax = 0.95
# else:
#     vmin = 50
#     vmax = 95

# view.setStyle({'cartoon': {'colorscheme': {'prop':'b','gradient': 'roygb','min': vmin,'max': vmax}}})




# with open("output_structures.pdb", "w") as f:
#     f.write("".join(novel_pdb))
#     f.write("".join(cannonical_pdb))
