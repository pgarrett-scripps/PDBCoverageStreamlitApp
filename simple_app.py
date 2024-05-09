


import json
from typing import Iterator
from urllib.request import urlopen

import py3Dmol
import stmol
import streamlit as st
import filterframes

DEFAULT_PROTEIN_ID = 'P00520'
DEFAULT_COV = [0, 0, 0, 0, 0, 2, 2, 2, 2, 2, 2, 2, 2, 2, 0, 0, 0, 0, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 0, 0, 0, 0, 0, 0, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 0, 0, 0, 0, 0, 0, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 0]


dta_select_file = st.file_uploader("Choose a file", type=['txt'])

if dta_select_file is not None:
    _, peptide_df, protein_df, _ = filterframes.from_dta_select_filter(dta_select_file)

# keep only rows, which have a locus that follows this regex sp|XXXX|XXXXX
protein_df = protein_df[protein_df['Locus'].str.match(r'sp\|[A-Z0-9]+\|[A-Z0-9]+')]

# split the locus column into 3 columns
protein_df[['Database', 'Protein', 'Gene']] = protein_df['Locus'].str.split('|', expand=True)

protein_of_interest = st.selectbox('Select protein', protein_df['Protein'].unique())


st.dataframe(protein_df.head())

import requests

def get_predictions(qualifier: str) -> Iterator[dict]:
    """Get all AlphaFold predictions for a UniProt accession.

    :param qualifier: A UniProt accession, e.g. P00520
    :type qualifier: str
    :return: The AlphaFold predictions
    :rtype: Iterator[dict]
    """
    url = f"https://alphafold.com/api/prediction/{qualifier}"
    # Retrieve the AlphaFold predictions with urllib
    with urlopen(url) as response:
        yield from json.loads(response.read().decode())


#download_alphafold_pdb(protein_of_interest, 'tmp.pdb')

for j in get_predictions(protein_of_interest):
    st.write(j)

    # get pdb file "pdbUrl"
    pdbUrl = j['pdbUrl']
    uniprotSequence = j['uniprotSequence']

    # download pdb file to a tmp file
    pdb_file = requests.get(pdbUrl)
    with open('tmp.pdb', 'wb') as f:
        f.write(pdb_file.content)


    def render_mol_cartoon(pdb):
        v = py3Dmol.view()
        v.addModel(pdb, 'mol')
        v.setStyle({'cartoon': {'colorscheme': 'ssPyMol'}, 'stick': {'radius': 0.05}})
        v.zoomTo()
        stmol.showmol(v, height=500, width=500)


    def render_mol(pdb, cov_arr, chain):
        view = py3Dmol.view()
        view.addModel(pdb, 'pdb')
        view.setStyle({}, {'cartoon': {}})

        view.setBackgroundColor('white')

        for i, c in enumerate(cov_arr):
            if c == 1:
                view.addStyle({"chain": chain, "resi": i, "elem": "C"},
                              {"cartoon": {"color": 'green', "radius": 0.2}})
            else:
                view.addStyle({"chain": chain, "resi": i, "elem": "C"},
                              {"cartoon": {"color": 'red', "radius": 0.2}})

        view.zoomTo()
        stmol.showmol(view, height=500, width=500)


    render_mol_cartoon(pdb_file.content.decode("utf-8"))
    render_mol(pdb_file.content.decode("utf-8"), [1, 0] * int(len(uniprotSequence)/2), 'A')


