from dataclasses import dataclass
from typing import Any, List, Literal, Optional
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

from app_input import get_input
from constants import DEFAULT_PROTEIN_ID, COLOR_MAPS, DEFAULT_PEPTIDES, DEFAULT_COLOR_MAP, DEFAULT_PROTEIN_SEQUENCE
from util import (
    coverage_string,
    display_header,
    get_predictions,
    render_mol,
    plot_coverage_array,
    compressor,
    decompressor,
    get_query_params_url,
    shorten_url,
    show_footer,
)

st.set_page_config(layout="centered", page_title="PdbCov",
                   page_icon=":dna:", initial_sidebar_state="auto")


with st.sidebar:
    display_header()
    cov_input = get_input()

should_render = True
try:
    cov_input.setup()
except Exception as e:
    st.error(f"Error setting up input: {e}")

    # print trace
    import traceback
    traceback.print_exc()

    
    should_render = False
    
if should_render:
    st.title(cov_input.title)
    if cov_input.subtitle:
        st.subheader(cov_input.subtitle)
        
    st.pyplot(
        plot_coverage_array(
            cov_input.color_coverage_array, 
            cov_input.color_map,
            )
        )

    if cov_input.pdb_content is not None:
        render_mol(
            cov_input.pdb_content, 
            cov_input.color_gradient_hex_array, 
            cov_input.pdb_style, 
            cov_input.bcolor, 
            cov_input.highlight_residues,
            cov_input.auto_spin,
        )

    st.markdown(
            coverage_string(cov_input.coverage_array, 
                            cov_input.protein_sequence, 
                            cov_input.cmap, 
                            color_coverage=cov_input.color_coverage_array),
            unsafe_allow_html=True,
        )
    
show_footer()