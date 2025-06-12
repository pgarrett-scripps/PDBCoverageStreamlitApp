import streamlit as st
import streamlit_permalink as stp
import peptacular as pt
import pandas as pd
import numpy as np
from collections import defaultdict

from constants import PDB_APP_URL
from util import serialize_peptides

st.set_page_config(layout="wide", page_title="Sage-PdbCov", page_icon=":microscope:")

with st.sidebar:
    st.title("Sage-PdbCov :microscope:")
    st.subheader("PdbCov Link Generator for Sage Parquet Files")
    st.caption("Upload a Sage parquet file to generate links to the PDB Viewer for each protein.")
    
    sage_file = st.file_uploader("Choose a Sage parquet file", type=['parquet'])
    
    st.subheader("Filter Options")
    q_value_threshold = st.slider(
        "Maximum Q-value",
        min_value=0.0,
        max_value=1.0,
        value=0.01,
        step=0.01,
        format="%.3f",
        help="Maximum Q-value for PSMs (lower values are more stringent)"
    )
    
    q_value_type = st.radio(
        "Q-value type to use for filtering",
        options=["spectrum_q", "peptide_q", "protein_q"],
        index=1,
        help="Choose which Q-value type to use for filtering"
    )

st.title("Protein Results")
st.caption("Click on the link icons to open the PDB Viewer for each protein.")

if sage_file is not None:
    # Load the parquet file
    psm_df = pd.read_parquet(sage_file)
    
    # Filter by q-value
    filtered_psm_df = psm_df[psm_df[q_value_type] <= q_value_threshold].copy()
    
    if filtered_psm_df.empty:
        st.error(f"No PSMs passed the {q_value_type} <= {q_value_threshold} filter. Try increasing the threshold.")
        st.stop()

    # split proteins by ;
    filtered_psm_df['proteins'] = filtered_psm_df['proteins'].str.split(';')

    # explode the proteins column
    filtered_psm_df = filtered_psm_df.explode('proteins')
    
    # Show filtering stats
    st.info(f"Filtered from {len(psm_df):,} to {len(filtered_psm_df):,} PSMs using {q_value_type} ≤ {q_value_threshold}")
    
    # Create required fields for PDB viewer
    filtered_psm_df['ProformaSequence'] = filtered_psm_df['peptide']
    filtered_psm_df['StrippedProformaSequence'] = filtered_psm_df['stripped_peptide']
    filtered_psm_df['ProformaSequenceCharge'] = filtered_psm_df.apply(
        lambda x: pt.add_mods(x['ProformaSequence'], {'charge': x['charge']}), axis=1
    )

    # Group by protein and collect peptides
    protein_to_peptides = defaultdict(list)
    for _, row in filtered_psm_df.iterrows():
        protein_to_peptides[row['proteins']].append(row['ProformaSequenceCharge'])
    
    # Create dataframe for proteins
    protein_data = []
    for protein, peptides in protein_to_peptides.items():
        # Parse protein ID - assuming format similar to db|id|gene
        protein_parts = protein.split('|')
        
        protein_id = protein
        db = None
        gene = None
        
        if len(protein_parts) == 3:
            db, protein_id, gene = protein_parts
        
        # Count unique peptides and spectra
        unique_peptides = len(set([pt.strip_mods(p) for p in peptides]))
        spectrum_count = len(peptides)
        
        # Serialize peptides for URL
        serialized_peptides = serialize_peptides(peptides)

        
        protein_data.append({
            'Protein': protein,
            'ProteinID': protein_id,
            'Database': db,
            'Gene': gene,
            'Unique Peptides': unique_peptides,
            'Spectrum Count': spectrum_count,
            'SerializedPeptides': serialized_peptides,
        })
    
    # Create protein dataframe
    protein_df = pd.DataFrame(protein_data)
    
    # Function to create PDB links
    def make_link(protein_id, serialized_peptides):
        params = {
            'input_type': 'Protein ID',
            'protein_id': protein_id,
            'peptides': serialized_peptides,
        }
        return stp.create_url(PDB_APP_URL, params)
    
    # Create links for each protein
    protein_df['Link'] = protein_df.apply(
        lambda x: make_link(x['ProteinID'], x['SerializedPeptides']), 
        axis=1
    )
    
    # Display protein results
    cols_to_show = ['Protein', 'Gene', 'Unique Peptides', 'Spectrum Count', 'Link']
    
    st.dataframe(
        data=protein_df[cols_to_show],
        hide_index=True,
        column_config={
            'Protein': st.column_config.TextColumn(width="medium", pinned=False),
            'Gene': st.column_config.TextColumn(width="small"),
            'Unique Peptides': st.column_config.NumberColumn(width="small"),
            'Spectrum Count': st.column_config.NumberColumn(width="small"),
            'Link': st.column_config.LinkColumn(display_text="🔗", pinned=True, width="small")
        },
        use_container_width=True
    )
    
    # Add download button
    st.download_button(
        label="Download Protein Data as CSV",
        data=protein_df[cols_to_show].to_csv(index=False).encode('utf-8'),
        file_name='sage_proteins.csv',
        mime='text/csv',
        on_click='ignore',
        use_container_width=True,
        type="primary",
    )
else:
    st.warning("No file uploaded")
    st.stop()