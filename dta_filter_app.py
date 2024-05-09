import streamlit as st
import filterframes
import peptacular as pt

from constants import PDB_APP_URL
from util import serialize_redundant_peptides

st.set_page_config(layout="wide")

dta_select_file = st.file_uploader("Choose a file", type=['txt'])

if dta_select_file is not None:
    _, peptide_df, protein_df, _ = filterframes.from_dta_select_filter(dta_select_file)
    peptide_df['ProformaSequence'] = peptide_df['Sequence'].apply(pt.convert_ip2_sequence)
    peptide_df['ProformaSequenceCharge'] = peptide_df.apply(lambda x: pt.add_mods(x['ProformaSequence'],
                                                                                  {'charge': x['Charge']}), axis=1)
    peptide_df['StrippedProformaSequence'] = peptide_df['ProformaSequence'].apply(pt.strip_mods)
    protein_df[['Database', 'Protein', 'Gene']] = protein_df['Locus'].str.split('|', expand=True)
    protein_df['Sequence Coverage'] = protein_df['Sequence Coverage'].str.rstrip('%').astype('float')
else:
    st.write("No file uploaded")
    st.stop()

protein_group_to_peptides = {}
for protein_group in protein_df['ProteinGroup'].unique():
    peptides = []
    for i, row in peptide_df[peptide_df['ProteinGroup'] == protein_group].iterrows():
        for _ in range(row['Redundancy']):
            peptides.append(row['ProformaSequenceCharge'])

    protein_group_to_peptides[protein_group] = peptides

protein_df['SerializedPeptides'] = protein_df['ProteinGroup'].apply(
    lambda x: serialize_redundant_peptides(protein_group_to_peptides[x]))


def make_link(protein_id, serialized_peptides):
    return f'{PDB_APP_URL}?protein_id={protein_id}&input={serialized_peptides}&input_type=redundant_peptides'


protein_df['Link'] = protein_df.apply(lambda x: make_link(x['Protein'], x['SerializedPeptides']), axis=1)

cols_to_keep = ['Locus', 'Descriptive Name', 'Sequence Count', 'Spectrum Count', 'Sequence Coverage', 'Length', 'Link']

st.dataframe(data=protein_df[cols_to_keep],
             hide_index=True,
             column_config={
                 'Locus': st.column_config.TextColumn(width="medium"),
                 'Sequence Count': st.column_config.NumberColumn(width="small"),
                 'Spectrum Count': st.column_config.NumberColumn(width="small"),
                 'Sequence Coverage': st.column_config.NumberColumn(width="small"),
                 'Length': st.column_config.NumberColumn(width="small"),
                 'Descriptive Name': st.column_config.TextColumn(width="large"),
                 'Link': st.column_config.LinkColumn(display_text="PDB Viewer")
             },
             use_container_width=True)
