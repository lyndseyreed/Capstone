## Works to take in the PDB formated file and launch an app that generates the 3D model in a way that is user friendly

##INPUT - PDB files folder that contains one file per gene, inside the files are the PDB formatted files with the x,y,z coordinates of
            # each atom in the structure

##USE - Streamlit API for app development, stmol for rendering py3Dmol in python

##OUTPUT - a link or URL that takes you to the online source that displays the interactive 3D model. It has a dropdown for the style (cartoon, stick, sphere)
            # and a dropdown for the gene based on the input folder contents.
            # this API allows for easy export or printing for a COA


import torch
import streamlit as st
from stmol import showmol
import py3Dmol
import os

#Calling Streamlit widgets
st.sidebar.title('L.Reed Capstone 2024')
st.sidebar.write('[*ESMFold*](https://esmatlas.com/about)  is utilized to render predicted 3D models and displays each canoncial sequence alongside the novel sequence(s).')
st.file_uploader("Upload the fasta file you wish to 3D model", accept_multiple_files=True)
st.button("Model!", type="primary")
style = st.sidebar.selectbox('Style', ['cartoon', 'stick', 'sphere'])
gene_files = os.listdir("PDB")
gene = st.sidebar.selectbox('Gene', gene_files)

pdbview = py3Dmol.view()

def render_gene_files(folder_path):
    """
    Takes in a folder of PDB files and generates the 3D models using stmol, py3Dmol and Streamlit
    """
    global pdbview
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


####USER INPUT AREA####
#Will be added to main eventually#

folder_path = 'PDB/' + gene
render_gene_files(folder_path)
