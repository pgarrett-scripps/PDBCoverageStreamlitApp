import numpy as np
import requests
import streamlit as st
from matplotlib.pyplot import colormaps

from constants import DEFAULT_PROTEIN_ID, DEFAULT_COVERAGE_ARRAY, AA_TO_THREE_LETTER_CODE, AAS, color_maps, \
    DEFAULT_PEPTIDES
from util import serialize_coverage_array, parse_coverage_array, get_predictions, render_mol, plot_coverage_array, \
    serialize_peptides, parse_peptides, serialize_redundant_peptides, parse_redundant_peptides
import peptacular as pt

st.set_page_config(initial_sidebar_state="collapsed")

# Set defaults
qp = st.query_params
if "protein_id" not in qp:
    qp["protein_id"] = DEFAULT_PROTEIN_ID
if "input_type" not in qp:
    qp["input_type"] = "coverage_array"

if qp["input_type"] not in ["coverage_array", "peptides", "redundant_peptides"]:
    st.error(f"Invalid input type: {qp['input_type']}")
    st.stop()

if 'input' not in qp:
    if qp["input_type"] == "coverage_array" and "coverage_array" not in qp:
        SERIALIZED_COVERAGE_ARRAY = serialize_coverage_array(DEFAULT_COVERAGE_ARRAY)
        qp["input"] = SERIALIZED_COVERAGE_ARRAY
    elif qp["input_type"] == "peptides" and "peptides" not in qp:
        SERIALIZED_PEPTIDES = serialize_peptides(DEFAULT_PEPTIDES)
        qp["input"] = SERIALIZED_PEPTIDES
    elif qp["input_type"] == "redundant_peptides" and "redundant_peptides" not in qp:
        SERIALIZED_REDUNDANT_PEPTIDES = serialize_redundant_peptides(DEFAULT_PEPTIDES)
        qp["input"] = SERIALIZED_REDUNDANT_PEPTIDES
    else:
        st.error(f"Should not happen... Big oops.")
        st.stop()


with st.sidebar:
    color_map = st.selectbox("Choose a color map", color_maps)
    bcolor = st.color_picker("Background color", "#FFFFFF")
    highlight_residues = [AA_TO_THREE_LETTER_CODE[aa].upper() for aa in st.multiselect("Highlight residues", AAS)]
    pdb_style = st.selectbox("PDB style", ['cartoon', 'stick', 'sphere', 'cross'])
    binary_coverage = st.checkbox("Binary coverage", False)

    protein_id = st.text_input("Protein ID", qp["protein_id"])

    predictions = list(get_predictions(protein_id))
    if len(predictions) == 0:
        st.error(f"No predictions found for protein {protein_id}")
    elif len(predictions) > 1:
        st.warning(f"Multiple predictions found for protein {protein_id}. Using the first one.")
    pdb_url = predictions[0]['pdbUrl']
    uniprotSequence = predictions[0]['uniprotSequence']

    INPUT_TYPES = ["coverage_array", "peptides", "redundant_peptides"]
    input_type = st.radio("Input type", INPUT_TYPES, index=INPUT_TYPES.index(qp["input_type"]), horizontal=True)
    if input_type == "coverage_array":
        coverage_array = np.array(parse_coverage_array(st.text_area("Coverage array", qp["input"])))
    elif input_type == "peptides":
        c1, c2 = st.columns(2)
        strip_mods = c1.checkbox("Strip mods", False)
        filter_unqiue = c2.checkbox("Filter unique", False)
        peptides = parse_peptides(st.text_area("Peptides", qp["input"]))
        if strip_mods:
            peptides = [pt.strip_mods(peptide) for peptide in peptides]
        if filter_unqiue:
            peptides = list(set(peptides))
        coverage_array = np.array(pt.coverage(uniprotSequence, peptides, True, True))
    elif input_type == "redundant_peptides":
        c1, c2 = st.columns(2)
        strip_mods = c1.checkbox("Strip mods", False)
        filter_unqiue = c2.checkbox("Filter unique", False)
        peptides = parse_redundant_peptides(st.text_area("Peptides", qp["input"]))
        if strip_mods:
            peptides = [pt.strip_mods(peptide) for peptide in peptides]
        if filter_unqiue:
            peptides = list(set(peptides))
        coverage_array = np.array(pt.coverage(uniprotSequence, peptides, True, True))


if binary_coverage:
    coverage_array = np.where(coverage_array > 0, 1, 0)

normalized_values = (coverage_array - coverage_array.min()) / (coverage_array.max() - coverage_array.min())
color_map_function = colormaps[color_map]
color_gradient_array = color_map_function(normalized_values)
color_gradient_hex_array = [f"#{int(r * 255):02x}{int(g * 255):02x}{int(b * 255):02x}" for r, g, b, _ in
                            color_gradient_array]


if len(uniprotSequence) != len(coverage_array):
    st.error(f"Length of coverage array ({len(coverage_array)}) does not match length of UniProt sequence "
             f"({len(uniprotSequence)})")
    st.stop()

pdb_file = requests.get(pdb_url)

description = predictions[0]['uniprotDescription']
st.title(description)
st.subheader(f"{predictions[0]['uniprotAccession']}|{predictions[0]['uniprotId']}")
st.pyplot(plot_coverage_array(coverage_array, color_map))
render_mol(pdb_file.content.decode("utf-8"), color_gradient_hex_array, pdb_style, bcolor, highlight_residues)
