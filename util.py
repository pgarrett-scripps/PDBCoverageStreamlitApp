import json
from collections import Counter
from typing import Iterator
from urllib.request import urlopen

import py3Dmol
import stmol
from matplotlib import pyplot as plt


def serialize_redundant_peptides(peptides: list) -> str:
    """Serialize the peptides from a list of strings to a string.

    :param peptides: The peptides as a list of strings
    :type peptides: list
    :return: The peptides as a string
    :rtype: str
    """
    # count peptides
    peptide_counter = Counter(peptides)
    return ",".join([f"{peptide};{redundancy}" for peptide, redundancy in peptide_counter.items()])


def parse_redundant_peptides(peptides: str) -> list:
    """Parse the peptides from a string to a list of strings.

    :param peptides: The peptides as a string
    :type peptides: str
    :return: The peptides as a list of strings
    :rtype: list
    """
    peptide_redundancies = peptides.split(",")
    peptides = [pr.split(";")[0] for pr in peptide_redundancies]
    redundancies = [int(pr.split(";")[1]) for pr in peptide_redundancies]
    return [peptide for peptide, redundancy in zip(peptides, redundancies) for _ in range(redundancy)]


def serialize_peptides(peptides: list) -> str:
    """Serialize the peptides from a list of strings to a string.

    :param peptides: The peptides as a list of strings
    :type peptides: list
    :return: The peptides as a string
    :rtype: str
    """
    return ",".join(peptides)


def parse_peptides(peptides: str) -> list:
    """Parse the peptides from a string to a list of strings.

    :param peptides: The peptides as a string
    :type peptides: str
    :return: The peptides as a list of strings
    :rtype: list
    """
    return peptides.split(",")


def serialize_coverage_array(coverage_array: list) -> str:
    """Serialize the coverage array from a list of integers to a string.

    :param coverage_array: The coverage array as a list of integers
    :type coverage_array: list
    :return: The coverage array as a string
    :rtype: str
    """
    return ",".join(str(c) for c in coverage_array)


def parse_coverage_array(coverage_array: str) -> list:
    """Parse the coverage array from a string to a list of integers.

    :param coverage_array: The coverage array as a string
    :type coverage_array: str
    :return: The coverage array as a list of integers
    :rtype: list
    """
    return [int(c) for c in coverage_array.split(",")]


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

    stmol.showmol(view, height=1000, width=1000)


def plot_coverage_array(coverage_array, color_map):
    # add a color bar to understand the
    fig, ax = plt.subplots(figsize=(10, 1))
    cbar = ax.imshow([coverage_array], aspect='auto', cmap=color_map)
    fig.colorbar(cbar, orientation='horizontal')
    # set title
    ax.set_title('Protein Coverage')
    ax.set_axis_off()
    return fig
