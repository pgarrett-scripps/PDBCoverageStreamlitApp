from urllib.error import HTTPError
import numpy as np
import requests
import streamlit as st
import streamlit_permalink as stp
from matplotlib.pyplot import colormaps
from Bio import PDB
import peptacular as pt
from streamlit_js_eval import get_page_location
import matplotlib as mpl
import matplotlib.colors as mcolors

from constants import DEFAULT_PROTEIN_ID, COLOR_MAPS, DEFAULT_PEPTIDES, DEFAULT_COLOR_MAP
from util import (
    get_predictions,
    render_mol,
    plot_coverage_array,
    compressor,
    decompressor,
    get_query_params_url,
    shorten_url,
)


st.set_page_config(layout="centered", page_title="PdbCov",
                   page_icon=":dna:", initial_sidebar_state="auto")

top_window = st.container()
bottom_window = st.container()



def change_input():

    if input_type == "Protein ID":
        url_params = {k: st.query_params.get_all(k) for k in st.query_params.keys()}
        st.session_state["saved_url_params"] = url_params
        st.query_params.clear()

    if input_type == "PDB file":
        if "saved_url_params" in st.session_state:
            # update the saved_url_params
            st.query_params.update(st.session_state["saved_url_params"])

    st.query_params["input_type"] = st.session_state["input_type"]
    st.session_state["rerun"] = True


is_stateful = False


pdb_url = None
uniprotSequence = None
pdb_content = None
description = None
uniprotAccession = None
uniprotId = None
title = None
sub_header = None

with st.sidebar:
    st.markdown(f"""
        <div style='text-align: center; padding: 15px; margin-top: 0px;'>
            <h3 style='margin: 0; font-size: 1.5em; color: #333;'>PdbCov: 3D Protein Coverage Analyzer</h3>
            <p style='font-size: 1.1em; line-height: 1.6; color: #555;'>
                Powered by 
                <a href="https://github.com/pgarrett-scripps/peptacular" target="_blank" style='color: #007BFF; text-decoration: none;'>
                    <strong>Peptacular</strong>
                </a>. 
                See the 
                <a href="https://peptacular.readthedocs.io/en/latest/modules/getting_started.html#proforma-notation" 
                target="_blank" style='color: #007BFF; text-decoration: none;'>
                    Proforma Notation Docs
                </a> for supported peptide syntax. To report any issues or suggest improvements, please visit the 
                <a href="https://github.com/pgarrett-scripps/PDBCoverageStreamlitApp" 
                target="_blank" style='color: #007BFF; text-decoration: none;'>
                    PDB Coverage Github Repo.
                </a>
            </p>
        </div>
    """, unsafe_allow_html=True)

    input_type = stp.radio(
        "Input type",
        ("Protein ID", "PDB file", "Protein Sequence"),
        key="input_type",
        horizontal=True,
        stateful=True,
        #on_change=change_input,
    )

    pdb_content, title, sub_header, uniprotSequence = None, "", "", ""
    if input_type == "PDB file":

        # Clear input params
        st.query_params.clear()
        st.query_params["input_type"] = "PDB file"

        pdb_file = st.file_uploader("Upload PDB file", type=["pdb"])

        if pdb_file is not None:
            pdb_content = pdb_file.read().decode("utf-8")

            # save the pdb file to a temporary file
            with open("temp.pdb", "w") as f:
                f.write(pdb_content)

            # Extract the protein sequence
            parser = PDB.PDBParser(QUIET=True)
            structure = parser.get_structure("uploaded_protein", "temp.pdb")

            for model in structure:
                for chain in model:
                    for residue in chain:
                        # Skip water and hetero atoms
                        if residue.id[0] == " ":
                            resname = residue.resname.strip()
                            # Convert three-letter code to one-letter code
                            one_letter = pt.constants.THREE_LETTER_CODE_TO_AA.get(
                                resname.capitalize(), "X"
                            )
                            uniprotSequence += one_letter

            title = pdb_file.name
            sub_header = None

    elif input_type == "Protein ID":
        is_stateful = True

        protein_id = stp.text_input(
            "Protein ID",
            key="protein_id",
            value=DEFAULT_PROTEIN_ID,
            stateful=is_stateful,
        )

        if protein_id:
            predictions = None
            try:
                predictions = get_predictions(protein_id)
            except HTTPError as e:
                top_window.error(
                    f"Error fetching predictions for protein {protein_id}: {e}"
                )

            if predictions is None:
                top_window.error(f"No predictions found for protein {protein_id}")
            elif len(predictions) == 0:
                top_window.error(f"No predictions found for protein {protein_id}")
            elif len(predictions) > 1:
                top_window.warning(
                    f"Multiple predictions found for protein {protein_id}. Using the first one."
                )
            else:
                pdb_url = predictions[0]["pdbUrl"]
                uniprotSequence = predictions[0]["uniprotSequence"]
                pdb_content = requests.get(pdb_url).content.decode("utf-8")
                description = predictions[0]["uniprotDescription"]
                uniprotAccession = predictions[0]["uniprotAccession"]
                uniprotId = predictions[0]["uniprotId"]
                title = predictions[0]["uniprotDescription"]
                sub_header = f"{uniprotAccession}|{uniprotId}"

    elif input_type == "Protein Sequence":
        protein_sequence = stp.text_area(
            "Protein Sequence",
            key="protein_sequence",
            height=300,
            stateful=is_stateful,
            placeholder="Enter protein sequence here.",
        )

        uniprotSequence = protein_sequence.strip().replace(" ", "").upper()
        if not uniprotSequence:
            top_window.error("Please enter a valid protein sequence.")
        else:
            title = "Custom Protein Sequence"
            sub_header = None

    peptide_str = stp.text_area(
        "Peptides (Proforma 2.0 notation)",
        value=DEFAULT_PEPTIDES if is_stateful else None,
        key="peptides",
        height=300,
        stateful=is_stateful,
        compressor=compressor,
        decompressor=decompressor,
        compress=True,
        placeholder="Enter peptides here, one per line.",
    )

    peptides = []
    if peptide_str:
        peptides = peptide_str.split("\n")

    with st.expander("Additional Options"):
        color_map = stp.selectbox(
            "Choose a color map", options=COLOR_MAPS, key="color_map", stateful=is_stateful, index=COLOR_MAPS.index(DEFAULT_COLOR_MAP)
        )
        pdb_style = stp.selectbox(
            "PDB style",
            ["cartoon", "stick", "sphere", "cross"],
            key="pdb_style",
            stateful=is_stateful,
        )
        bcolor = stp.color_picker(
            "Background color", "#FFFFFF", key="bcolor", stateful=is_stateful
        )
        selected_residue = stp.multiselect(
            "Select residue",
            pt.AMINO_ACIDS,
            key="selected_residue",
            stateful=is_stateful,
        )
        highlight_residues = [
            pt.AA_TO_THREE_LETTER_CODE[aa].upper() for aa in selected_residue
        ]
        binary_coverage = stp.checkbox(
            "Binary coverage", False, key="binary_coverage", stateful=is_stateful
        )
        reverse_protein = stp.checkbox(
            "Reverse protein", False, key="reverse_protein", stateful=is_stateful
        )
        strip_mods = stp.checkbox("Strip mods", True, stateful=is_stateful)
        filter_unique = stp.checkbox("Filter unique", False, stateful=is_stateful)

with top_window:

    if is_stateful:
        title_c, button_c = st.columns([3, 1], vertical_alignment="top")
        with title_c:
            st.header("Coverage Results")
        with button_c:
            url_btn = st.button(
                "Generate URL",
                key="generate_url",
                type="primary",
                use_container_width=True,
            )

        st.caption(
            f"""**This pages URL automatically updates with your input, and can be shared with others. 
                   You can optionally use the Generate TinyURL button to create a shortened URL.**""",
            unsafe_allow_html=True,
        )

        st.divider()

    if not uniprotSequence or not pdb_content:
        st.error("Please provide a valid input (either a PDB file or a Protein ID).")
        st.stop()

    if reverse_protein:
        uniprotSequence = uniprotSequence[::-1]

    if strip_mods:
        peptides = [pt.strip_mods(peptide) for peptide in peptides]
    if filter_unique:
        peptides = list(set(peptides))

    # check if there are non mathcing peptides
    missing_peptides = set()
    for peptide in peptides:
        if not pt.is_subsequence(peptide, uniprotSequence):
            missing_peptides.add(peptide)

    if missing_peptides:
        st.warning(f"Some peptides do not match the protein sequence!")
        with st.expander("Show missing peptides"):
            st.write(missing_peptides)

    coverage_array = np.array(pt.coverage(uniprotSequence, peptides, True, True))

    if reverse_protein:
        coverage_array = coverage_array[::-1]

    if binary_coverage:
        coverage_array = np.where(coverage_array > 0, 1, 0)

    normalized_values = (coverage_array - coverage_array.min()) / (
        coverage_array.max() - coverage_array.min()
    )
    color_map_function = colormaps[color_map]
    color_gradient_array = color_map_function(normalized_values)
    color_gradient_hex_array = [
        f"#{int(r * 255):02x}{int(g * 255):02x}{int(b * 255):02x}"
        for r, g, b, _ in color_gradient_array
    ]

    if sum(coverage_array) == 0:
        color_gradient_hex_array = ["#FFFFFF"] * len(color_gradient_hex_array)

    if len(uniprotSequence) != len(coverage_array):
        st.error(
            f"Length of coverage array ({len(coverage_array)}) does not match length of UniProt sequence "
            f"({len(uniprotSequence)})"
        )
        st.stop()

    st.title(title, anchor='PDB_TITLE')
    if sub_header:
        st.subheader(sub_header)
    st.pyplot(plot_coverage_array(coverage_array, color_map))
    render_mol(
        pdb_content, color_gradient_hex_array, pdb_style, bcolor, highlight_residues
    )

    # 2d representation

    sites = []
    for site, aa in enumerate(uniprotSequence):
        if aa in selected_residue:
            sites.append(site)

    # Create a colormap
    cmap = mpl.colormaps.get_cmap(color_map)

    def coverage_string(protein_cov_arr, stripped_protein_sequence, cmap, sites=None):
        # Find the maximum coverage value
        max_coverage = max(protein_cov_arr)

        # Color all covered amino acids based on coverage and show index on hover, using a monospace font
        protein_cov = (
            '<span style="font-family: Courier New, monospace; font-size: 16px;">'
        )
        for i, aa in enumerate(stripped_protein_sequence):
            coverage = protein_cov_arr[i]

            # Normalize the coverage based on the maximum value
            normalized_coverage = coverage / max_coverage

            # Get color from colormap
            color = cmap(normalized_coverage)
            hex_color = mcolors.to_hex(color)

            if coverage > 0 and max_coverage > 0:
                if sites and i in sites:
                    protein_cov += (
                        f'<span title="Index: {i + 1}; Coverage: {coverage}" style="background-color:red; '
                        f"color:{hex_color}; font-weight:900; padding:3px; margin:1px; border:1px solid "
                        f'#cccccc; border-radius:3px;">{aa}</span>'
                    )
                else:
                    protein_cov += (
                        f'<span title="Index: {i + 1}; Coverage: {coverage}" style="background-color'
                        f":#e0e0ff; color:{hex_color}; font-weight:900; padding:3px; margin:1px; "
                        f'border:1px solid #a0a0ff; border-radius:3px;">{aa}</span>'
                    )
            else:
                if sites and i in sites:
                    protein_cov += (
                        f'<span title="Index: {i + 1}" style="background-color:red; color:{hex_color}; '
                        f"font-weight:900; padding:3px; margin:1px; border:1px solid #cccccc;"
                        f' border-radius:3px;">{aa}</span>'
                    )
                else:
                    # Style the non-covered amino acid for a polished look with a tooltip
                    protein_cov += (
                        f'<span title="Index: {i + 1}" style="background-color:#f0f0f0; color:{hex_color}; '
                        f"font-weight:900; padding:3px; margin:1px; border:1px solid #cccccc; "
                        f'border-radius:3px;">{aa}</span>'
                    )
        protein_cov += "</span>"

        return protein_cov

    st.markdown(
        coverage_string(coverage_array, uniprotSequence, cmap, sites),
        unsafe_allow_html=True,
    )

    if is_stateful:
        @st.fragment
        def url_fragment():

            title_c, _, button_c = st.columns([2, 1, 1])
            title_c.header("PepFrag Results")

            st.caption(
                '''**This pages URL automatically updates with your input, and can be shared with others. 
            You can also click on the 'Generate TinyURL' button to create a shortened URL.**''',
                unsafe_allow_html=True,
            )

            if button_c.button("Generate TinyURL", key="generate_tinyurl", type="primary"):
                short_url = shorten_url(stp.get_page_url())
                st.caption(f"Shortened URL: {short_url}")

        url_fragment()


    st.divider()

    st.markdown(f"""
        <div style='display: flex; justify-content: space-between; align-items: center; padding: 15px 0; border-top: 0px solid #ddd;'>
            <div style='text-align: left; font-size: 1.1em; color: #555;'>
                <a href="https://github.com/pgarrett-scripps/PDBCoverageStreamlitApp" target="_blank" 
                   style='text-decoration: none; color: #007BFF; font-weight: bold;'>
                    PDBCoverage
                </a>
                <a href="https://doi.org/10.5281/zenodo.15066418" target="_blank" style="margin-left: 12px;">
                    <img src="https://zenodo.org/badge/798509918.svg" alt="DOI" 
                         style="vertical-align: middle; height: 20px;">
                </a>
            </div>
            <div style='text-align: right; font-size: 1.1em; color: #555;'>
                <a href="https://github.com/pgarrett-scripps/peptacular" target="_blank" 
                   style='text-decoration: none; color: #007BFF; font-weight: bold;'>
                    Peptacular
                </a>
                <a href="https://doi.org/10.5281/zenodo.15054278" target="_blank" style="margin-left: 12px;">
                    <img src="https://zenodo.org/badge/591504879.svg" alt="DOI" 
                         style="vertical-align: middle; height: 20px;">
                </a>
            </div>
        </div>
    """, unsafe_allow_html=True)
