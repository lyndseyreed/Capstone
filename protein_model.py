##Using ESMFold API, input the sequences extracted from sequence_parser and obtain a 3D model with corresponding novel and canonical for comparison##

##Input - list of aa sequences with canonical as the first sequence in the list
    ##Use ESMFold API to predict the structures of the sequences
##Output - 3D model that allows visualization of canonical vs novel tertiary model predictions

# pip install --upgrade biopython transformers py3Dmol accelerate

import os
import torch
from esm.model import ESM
from esm.pretrained import load_model_and_alphabet
from Bio import SeqIO
from Bio.PDB import Superimposer
from sequence_parser import get_sequences

#Load in model and tokenizer
from transformers import AutoTokenizer, EsmForProteinFolding
tokenizer = AutoTokenizer.from_pretrained("facebook/esmfold_v1")
model = EsmForProteinFolding.from_pretrained("facebook/esmfold_v1", low_cpu_mem_usage=True)

#transfering model to GPU
model = model.to()

#Convert language model stem to float16 (improve mem use on GPU)
model.esm = model.esm.half()

#Enabled to speed up computation if supported on mac
import torch
model.trunk.set_chunk_size(64)

##WORKS UP TO HERE###

