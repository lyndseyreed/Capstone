import torch
import streamlit as st
from stmol import showmol
import py3Dmol
import os

#Calling Streamlit widgets
st.sidebar.title('L.Reed Capstone 2024')
st.sidebar.write('[*ESMFold*](https://esmatlas.com/about)  is utilized to render predicted 3D models and displays each canoncial sequence alongside the novel sequence(s).')
style = st.sidebar.selectbox('Style', ['cartoon', 'stick', 'sphere'])
gene_files = os.listdir("PDB")
gene = st.sidebar.selectbox('Gene', gene_files)

pdbview = py3Dmol.view()

def render_gene_files(folder_path):
    global pdbview
    # Iterate through all files in the folder
    for filename in os.listdir(folder_path):
        file_path = os.path.join(folder_path, filename)
        if os.path.isfile(file_path):
            with open(file_path, 'r') as file:
                system = "".join([x for x in file])
                pdbview.addModel(system)
                pdbview.setStyle({style:{'color':'spectrum'}})
                pdbview.setBackgroundColor('white')
        pdbview.zoomTo()
        showmol(pdbview, height = 500,width=800)


folder_path = 'PDB/' + gene
render_gene_files(folder_path)
