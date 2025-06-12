import json
from collections import Counter
from typing import Iterator
from urllib.parse import quote_plus
from urllib.request import urlopen

import py3Dmol
import requests
import stmol
from matplotlib import pyplot as plt
import streamlit as st

from itertools import groupby
from typing import List


def compressor(peptide_str: str) -> str:
    """Compress consecutive duplicate peptides into a compact representation."""
    if not peptide_str.strip():
        return ""

    compressed = []

    for peptide, group in groupby(peptide_str.splitlines()):
        count = sum(1 for _ in group)
        compressed.append(f"{peptide};{count}")

    return ','.join(compressed)


def decompressor(peptide_str: str) -> str:
    """Decompress the compact peptide representation back into its original form."""
    if not peptide_str.strip():
        return ""

    try:
        peptides = []

        for pair in peptide_str.split(','):
            peptide, count = pair.split(';')
            peptides.extend([peptide] * int(count))

        return '\n'.join(peptides)

    except (ValueError, IndexError) as e:
        raise ValueError(f"Invalid compressed format: {peptide_str}") from e


def serialize_peptides(peptides: List[str]) -> str:
    """Serialize a list of peptides into a single string."""
    # counter
    peptide_counts = Counter(peptides)
    serialized = ','.join([f"{peptide};{count}" for peptide, count in peptide_counts.items()])
    return serialized


def _get_predictions(qualifier: str) -> Iterator[dict]:
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


@st.cache_data()
def get_predictions(qualifier: str) -> list:
    """Get all AlphaFold predictions for a UniProt accession.

    :param qualifier: A UniProt accession, e.g. P00520
    :type qualifier: str
    :return: The AlphaFold predictions
    :rtype: list
    """
    return list(_get_predictions(qualifier))


def render_mol(pdb, cov_arr, pdb_style, bcolor, highlight_residues):
    view = py3Dmol.view()
    view.addModel(pdb, 'pdb')
    view.setStyle({}, {pdb_style: {}})

    view.setBackgroundColor(bcolor)

    for i, c in enumerate(cov_arr):
        view.addStyle({"resi": i},
                      {pdb_style: {"color": c, "radius": 0.2}})

    view.addResLabels({'resn': highlight_residues, })
    stmol.add_hover(view)

    view.zoomTo()

    stmol.showmol(view, height=500, width=700)


def plot_coverage_array(coverage_array, color_map):
    # add a color bar to understand the
    fig, ax = plt.subplots(figsize=(10, 1))
    cbar = ax.imshow([coverage_array], aspect='auto', cmap=color_map)
    fig.colorbar(cbar, orientation='horizontal')
    # set title
    ax.set_title('Protein Coverage')
    ax.set_axis_off()
    return fig


def get_query_params_url(params_dict):
    """
    Create url params from alist of parameters and a dictionary with values.

    Args:
        params_list (str) :
            A list of parameters to get the value of from `params_dict`
        parmas_dict (dict) :
            A dict with values for the `parmas_list .
        **kwargs :
            Extra keyword args to add to the url
    """
    return "?" + "&".join(
        [
            f"{key}={quote_plus(str(value))}"
            for key, values in params_dict.items()
            for value in listify(values)
        ]
    )


def listify(o=None):
    if o is None:
        res = []
    elif isinstance(o, list):
        res = o
    elif isinstance(o, str):
        res = [o]
    else:
        res = [o]
    return res


def shorten_url(url: str) -> str:
    """Shorten a URL using TinyURL."""
    api_url = f"http://tinyurl.com/api-create.php?url={url}"

    try:
        response = requests.get(api_url)
        response.raise_for_status()
        return response.text
    except requests.RequestException as e:
        return f"Error: {e}"
    

def coverage_string(protein_cov_arr, stripped_protein_sequence, cmap, sites=None):
    import matplotlib.colors as mcolors

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


def show_footer():
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
