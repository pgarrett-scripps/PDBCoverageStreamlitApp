from dataclasses import dataclass
from typing import Any, List, Optional

import matplotlib as mpl
from matplotlib.pyplot import colormaps
import numpy as np
import peptacular as pt
from requests import HTTPError
import requests
import streamlit_permalink as stp
import streamlit as st
from constants import *
import tempfile
import matplotlib.colors as mcolors
from Bio import PDB

from util import (
    get_predictions,
    compressor,
    decompressor,
)

PROTEIN_ID_TYPE = "Protein ID"
PDB_FILE_TYPE = "PDB file"
PROTEIN_SEQUENCE_TYPE = "Protein Sequence"
INPUT_TYPE_OPTIONS = [PROTEIN_ID_TYPE, PDB_FILE_TYPE, PROTEIN_SEQUENCE_TYPE]

class InputType:

    def __init__(self):
        """Initialize the InputType class."""
        self.is_setup: bool = False

    def setup(self):
        """Setup the input type, fetching necessary data."""
        raise NotImplementedError(
            "This method should be implemented in subclasses to set up the input type."
        )

    @property
    def protein_sequence(self) -> Optional[str]:
        """Return the protein sequence if available."""
        raise NotImplementedError(
            "This method should be implemented in subclasses to return the protein sequence."
        )
    
    @property
    def pdb_content(self) -> Optional[Any]:
        """Return the PDB content if available."""
        raise NotImplementedError(
            "This method should be implemented in subclasses to return the PDB content."
        )
    
    @property
    def title(self) -> Optional[str]:
        """Return a description of the input type."""
        raise NotImplementedError(
            "This method should be implemented in subclasses to return a description."
        )
    
    @property
    def subtitle(self) -> Optional[str]:
        """Return the UniProt accession if available."""
        raise NotImplementedError(
            "This method should be implemented in subclasses to return the UniProt accession."
        )
    
    

class ProteinID(InputType):

    def __init__(self, protein_id: Optional[str]):
        super().__init__()
        self.protein_id = protein_id
        self.predictions = None

    def setup(self):
        if not self.protein_id:
            raise ValueError("Protein ID cannot be empty.")
        
        try:
            predictions = get_predictions(self.protein_id)
        except HTTPError as e:
            raise ValueError(
                f"Failed to fetch predictions for Protein ID {self.protein_id}: {e}"
            )
        except Exception as e:
            raise ValueError(
                f"An unexpected error occurred while fetching predictions for Protein ID {self.protein_id}: {e}"
            )
        if not predictions:
            raise ValueError(
                f"No predictions found for Protein ID {self.protein_id}."
            )
        
        self.predictions = predictions

    @property
    def protein_sequence(self) -> Optional[str]:
        """Return the protein sequence from the predictions."""        
        return self.predictions[0].get("uniprotSequence", None)


    @property
    def pdb_content(self) -> Optional[Any]:

        if "pdbUrl" not in self.predictions[0]:
            raise ValueError(
                f"No PDB URL found for Protein ID {self.protein_id}."
            )

        pdb_url = self.predictions[0].get("pdbUrl", None)

        if not pdb_url:
            raise ValueError(
                f"No PDB URL found for Protein ID {self.protein_id}."
            )
        
        return requests.get(pdb_url).content.decode("utf-8")
    
    @property
    def title(self) -> Optional[str]:
        """Return a description of the input type."""
        return self.predictions[0].get("uniprotDescription", None)
    
    @property
    def subtitle(self) -> Optional[str]:
        uniprotAccession = self.predictions[0].get("uniprotAccession", None)
        uniprotId = self.predictions[0].get("uniprotId", None)
        return f"{uniprotAccession}|{uniprotId}"

    

class PDBFile(InputType):


    def __init__(self, pdb_file: Optional[Any] = None):
        super().__init__()
        self.pdb_file = pdb_file

        self._title = None
        self._subtitle = None
        self._protein_sequence = None
        self._pdb_content = None

    def setup(self):
        if self.pdb_file is not None:
            self._pdb_content = self.pdb_file.read()

            # Create temporary PDB file using context manager
            with tempfile.NamedTemporaryFile(suffix='.pdb', delete=True) as temp_pdb:
                # Handle both string and bytes content
                if isinstance(self._pdb_content, str):
                    temp_pdb.write(self._pdb_content.encode('utf-8'))
                else:
                    temp_pdb.write(self._pdb_content)
                temp_pdb.flush()  # Ensure all data is written to disk
                # Extract the protein sequence
                
                parser = PDB.PDBParser(QUIET=True)
                structure = parser.get_structure("uploaded_protein", temp_pdb.name)

                _protein_sequence = ""
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
                                _protein_sequence += one_letter

                self._protein_sequence = _protein_sequence
                self._title = self.pdb_file.name
                self._subtitle = None
            self._pdb_content = self._pdb_content.decode("utf-8") if isinstance(self._pdb_content, bytes) else self._pdb_content
        else:
            raise ValueError("PDB file cannot be None. Please upload a valid PDB file.")

    @property
    def protein_sequence(self) -> Optional[str]:
        """Return the protein sequence if available."""
        
        if self._protein_sequence is None:
            raise ValueError("Protein sequence is not available. Please check the PDB file.")
        return self._protein_sequence 
    
    @property
    def pdb_content(self) -> Optional[Any]:
        """Return the PDB content if available."""
        if self._pdb_content is None:
            raise ValueError("PDB content is not available. Please check the PDB file.")
        return self._pdb_content
    
    @property
    def title(self) -> Optional[str]:
        """Return a description of the input type."""
        if self._title is None:
            raise ValueError("Title is not available. Please check the PDB file.")
        return self._title
    
    @property
    def subtitle(self) -> Optional[str]:
        """Return the UniProt accession if available."""
        return self._subtitle

class ProteinSequence(InputType):
    def __init__(self, sequence: Optional[str] = None):
        super().__init__()
        self.sequence = sequence

    def setup(self):
        if not self.sequence:
            raise ValueError("Protein sequence cannot be empty.")
        
    @property
    def protein_sequence(self) -> Optional[str]:
        """Return the protein sequence."""
        if not self.sequence:
            raise ValueError("Protein sequence is not set.")
        return self.sequence
    
    @property
    def pdb_content(self) -> Optional[Any]:
        """Return None as there is no PDB content for a sequence input."""
        return None
    
    @property
    def title(self) -> Optional[str]:
        """Return a description of the input type."""
        return f"Protein Sequence"
    
    @property
    def subtitle(self) -> Optional[str]:
        """Return None as there is no UniProt accession for a sequence input."""
        return f"Sequence Length: {len(self.sequence)}"

class CoverageAppConfig:


    def __init__(self,
                 input_type: InputType,
                 peptides: list[str],
                 color_map: str,
                 reverse: bool,
                 pdb_style: str,
                 bcolor: str,
                 selected_residue: List[str],
                binary_coverage: bool,
                strip_mods: bool,
                filter_unique: bool,
                auto_spin: bool,
                user_title: Optional[str],
                user_subtitle: Optional[str],
                colorbar_min: Optional[float],
                colorbar_max: Optional[float]):
        
        self.input_type = input_type
        self.peptides = peptides
        self.color_map = color_map
        self.reverse = reverse
        self.pdb_style = pdb_style
        self.bcolor = bcolor
        self.selected_residue = selected_residue
        self.binary_coverage = binary_coverage
        self.strip_mods = strip_mods
        self.filter_unique = filter_unique
        self.auto_spin = auto_spin
        self.user_title = user_title
        self.user_subtitle = user_subtitle
        self.colorbar_min = colorbar_min
        self.colorbar_max = colorbar_max


    def setup(self):
        """Setup the input type and validate the configuration."""
        self.input_type.setup()

    @property
    def highlight_residues(self) -> List[str]:
        """Return the residues to highlight based on the selected residues."""
        highlight_residues = [
            pt.AA_TO_THREE_LETTER_CODE[aa].upper() for aa in self.selected_residue
        ]

        return highlight_residues

    @property
    def filtered_peptides(self) -> List[str]:
        """Return the filtered peptides based on the configuration."""
        if self.strip_mods:
            peptides = [pt.strip_mods(p) for p in self.peptides]
        else:
            peptides = self.peptides
        
        if self.filter_unique:
            peptides = list(set(peptides))

        # remove ambiguity
        peptides = [pt.strip_mods(p) for p in peptides]
        
        return peptides

    @property
    def protein_sequence(self) -> Optional[str]:
        """Return the protein sequence from the input type."""
        if self.reverse is True:
            return self.input_type.protein_sequence[::-1]
        return self.input_type.protein_sequence
    
    @property
    def pdb_content(self) -> Optional[Any]:
        """Return the PDB content from the input type."""
        return self.input_type.pdb_content
    
    @property
    def title(self) -> Optional[str]:
        """Return a description of the input type."""
        if self.user_title:
            return self.user_title
        return self.input_type.title
    
    @property
    def subtitle(self) -> Optional[str]:
        """Return the UniProt accession from the input type."""
        if self.user_subtitle:
            return self.user_subtitle
        return self.input_type.subtitle
    
    @property
    def coverage_array(self) -> np.ndarray:
        """Return the coverage array based on the peptides and protein sequence."""
        coverage_arr = np.array(pt.coverage(self.protein_sequence, self.filtered_peptides, True, True))

        if self.binary_coverage:
            coverage_arr = np.where(coverage_arr > 0, 1, 0)
        if len(coverage_arr) != len(self.protein_sequence):
            raise ValueError(
                f"Length of coverage array ({len(coverage_arr)}) does not match length of protein sequence ({len(self.protein_sequence)})."
            )
        return coverage_arr
    
    @property
    def color_coverage_array(self):
        
        vmin = 0
        if self.colorbar_min is not None:
            vmin = self.colorbar_min
        
        vmax = max(self.coverage_array.max(), 1)
        if self.colorbar_max is not None:
            vmax = self.colorbar_max

        # Clamp values to the specified range
        clamped_array = np.clip(self.coverage_array, vmin, vmax)

        return clamped_array
    

    @property
    def color_gradient_hex_array(self) -> list[str]:

        color_min = self.colorbar_min if self.colorbar_min is not None else self.color_coverage_array.min()
        color_max = self.colorbar_max if self.colorbar_max is not None else self.color_coverage_array.max()

        normalized_values = (self.color_coverage_array - color_min) / (color_max - color_min)

        color_map_function = colormaps[self.color_map]
        color_gradient_array = color_map_function(normalized_values)
        color_gradient_hex_array = [
            f"#{int(r * 255):02x}{int(g * 255):02x}{int(b * 255):02x}"
            for r, g, b, _ in color_gradient_array
        ]

        if sum(self.coverage_array) == 0:
            color_gradient_hex_array = ["#FFFFFF"] * len(color_gradient_hex_array)

        if len(self.protein_sequence) != len(self.coverage_array):
            raise ValueError(
                f"Length of coverage array ({len(self.coverage_array)}) does not match length of protein sequence ({len(self.protein_sequence)})."
            )
        
        return color_gradient_hex_array
    
    @property
    def cmap(self) -> mcolors.Colormap:
        """Return the color map object based on the selected color map."""
        return mpl.colormaps.get_cmap(self.color_map)

    @property
    def vmin(self) -> Optional[float]:
        """Return the minimum value for the colorbar."""
        return self.colorbar_min
    
    @property
    def vmax(self) -> Optional[float]:
        """Return the maximum value for the colorbar."""
        return self.colorbar_max

    @property
    def sites(self) -> List[int]:
        sites = []
        for site, aa in enumerate(self.protein_sequence):
            if aa in self.selected_residue:
                sites.append(site)
        return sites

def get_input() -> CoverageAppConfig:
    """Get input from the user or URL parameters."""
    # Get URL parameters
    
    input_type = stp.radio(
        "Input Type",
        options=INPUT_TYPE_OPTIONS,
        index=0,
        horizontal=True,
        help="Select the type of input for the PDB Viewer.",
        key="input_type",
    )

    cov_input = None
    if input_type == PROTEIN_ID_TYPE:
        protein_id = stp.text_input(
            "Protein ID",
            value=DEFAULT_PROTEIN_ID,
            help="Enter the Protein ID (e.g., P12345).",
            key="protein_id",
        )
        cov_input = ProteinID(protein_id=protein_id)
    elif input_type == PDB_FILE_TYPE:
        pdb_file = st.file_uploader(
            "Upload PDB File",
            type=["pdb"],
            help="Upload a PDB file to visualize its coverage.",
            key="pdb_file",
        )
        cov_input = PDBFile(pdb_file=pdb_file)
    elif input_type == PROTEIN_SEQUENCE_TYPE:
        sequence = stp.text_area(
            "Protein Sequence",
            value=DEFAULT_PROTEIN_SEQUENCE,
            help="Enter the protein sequence.",
            key="protein_sequence",
            height=150,
        )
        cov_input = ProteinSequence(sequence=sequence)
    else:
        raise ValueError("Invalid input type selected.")
    

    reverse = stp.toggle(
        "Reverse Protein Sequence",
        value=False,
        help="Reverse the protein sequence for visualization.",
        key="reverse_protein",
    )

    # Get peptides input
    peptides_input = stp.text_area(
        "Peptides (Proforma 2.0 notation)",
        value=DEFAULT_PEPTIDES,
        help="Enter the peptides to visualize coverage, separated by new lines.",
        key="peptides",
        height=200,
        compressor=compressor,
        decompressor=decompressor,
        compress=True
    )

    peptides = []
    if peptides_input:
        peptides = peptides_input.split("\n")

    peptides = [p.strip() for p in peptides]



    c1, c2 = st.columns(2, vertical_alignment="bottom")

    with c1:
        color_map = stp.selectbox(
            "Choose a color map", 
            options=COLOR_MAPS, 
            index=COLOR_MAPS.index(DEFAULT_COLOR_MAP),
            help="Select a color map for visualizing coverage.",
            key="color_map", 
        )

    with c2:
        pdb_style = stp.selectbox(
            "PDB style",
            ["cartoon", "stick", "sphere", "cross"],
            key="pdb_style",
        )

    c1, c2 = st.columns([1, 3], vertical_alignment="bottom")

    with c1:
        bcolor = stp.color_picker(
            "Background", 
            "#FFFFFF", 
            key="bcolor", 
        )
    with c2:
        selected_residue = stp.multiselect(
            "Select residue",
            pt.AMINO_ACIDS,
            key="selected_residue",
        )

    c1, c2 = st.columns(2, vertical_alignment="bottom")

    with c1:
        binary_coverage = stp.checkbox(
            "Binary coverage", 
            False, 
            key="binary_coverage",
        )

        auto_spin = stp.checkbox(
            "Auto-spin",
            value=True,
            help="Enable auto-spin for the PDB structure visualization.",
            key="auto_spin",
        )

    with c2:
        strip_mods = stp.checkbox("Strip mods", 
                                False,
                                key="strip_mods",
                                help="Strip modifications from peptides before coverage calculation."
        )
        filter_unique = stp.checkbox("Filter unique", 
                                    False,
                                    key="filter_unique",
                                    help="Filter unique peptides before coverage calculation."
        )


    with st.expander('User-defined Title and Subtitle', expanded=False):
        user_title = stp.text_input(
            label="Title",
            help="Enter a title for the coverage viewer.",
            key="title",
        )
        user_subtitle = stp.text_input(
            label="Subtitle",
            help="Enter a subtitle for the coverage viewer.",
            key="subtitle",
        )

    # Colorbar min/max controls
    with st.expander("Colorbar Settings", expanded=False):
        col1, col2 = st.columns(2)
        
        with col1:

            if st.button("Reset Colorbar Min", use_container_width=True):
                # Reset colorbar min and max to None
                stp.number_input.set_url_value(
                    url_key="colorbar_min", 
                    value=None, 
                )

            colorbar_min = stp.number_input(
                "Colorbar Min",
                value=None,
                min_value=0,
                step=1,
                help="Set minimum value for colorbar scaling. Leave empty for auto-scaling.",
                key="colorbar_min",
            )
        
        with col2:

            if st.button("Reset Colorbar Max", use_container_width=True):
                # Reset colorbar max to None
                stp.number_input.set_url_value(
                    url_key="colorbar_max", 
                    value=None, 
                )

            colorbar_max = stp.number_input(
                "Colorbar Max", 
                value=None,
                min_value=0,
                step=1,
                help="Set maximum value for colorbar scaling. Leave empty for auto-scaling.",
                key="colorbar_max",
            )

    return CoverageAppConfig(
        input_type=cov_input,
        peptides=peptides,
        color_map=color_map,
        reverse=reverse,
        pdb_style=pdb_style,
        bcolor=bcolor,
        selected_residue=selected_residue,   
        binary_coverage=binary_coverage,
        strip_mods=strip_mods,
        filter_unique=filter_unique, 
        auto_spin=auto_spin,
        user_title=user_title,
        user_subtitle=user_subtitle,
        colorbar_min=colorbar_min,
        colorbar_max=colorbar_max,
    )