from collections import Counter

import streamlit as st
import pandas as pd
import urllib.parse

from util import serialize_redundant_peptides



def main():
    st.title("Coverage App URL Generator")

    # Upload the data file
    uploaded_file = st.file_uploader("Upload your data file", type=['csv', 'tsv', 'txt'])

    if uploaded_file is not None:
        # Determine the separator based on file extension
        if uploaded_file.name.endswith('.csv'):
            df = pd.read_csv(uploaded_file)
        elif uploaded_file.name.endswith('.tsv') or uploaded_file.name.endswith('.txt'):
            df = pd.read_csv(uploaded_file, sep='\t')
        else:
            st.error("Unsupported file type.")
            return

        # Ensure required columns exist
        required_columns = ['Protein.Ids', 'Stripped.Sequence']
        if not all(column in df.columns for column in required_columns):
            st.error(f"The uploaded file must contain the following columns: {', '.join(required_columns)}")
            return

        # Handle cases where 'Protein.Ids' contains multiple IDs separated by ';'
        df['Protein.Ids'] = df['Protein.Ids'].astype(str)  # Ensure the column is string type
        df['Protein.Ids'] = df['Protein.Ids'].str.split(';')  # Split the IDs into lists
        df = df.explode('Protein.Ids')  # Expand the lists into rows

        # Remove any leading/trailing whitespace in 'Protein.Ids'
        df['Protein.Ids'] = df['Protein.Ids'].str.strip()

        sample_names = df['Sample.Name'].unique()

        #sort the sample names
        sample_names = sorted(sample_names)

        selected_samples = st.multiselect("Select samples", sample_names, default=sample_names)

        df = df[df['Sample.Name'].isin(selected_samples)]


        # Replace with the actual URL of your coverage app
        coverage_app_base_url = 'https://pdb-coverage.streamlit.app/'  # Update this URL

        # Group the data by Protein ID
        grouped = df.groupby('Protein.Ids')

        url_list = []
        for protein_id, group in grouped:
            peptides = group['Stripped.Sequence'].dropna().astype(str).tolist()
            if not peptides:
                continue  # Skip if no peptides are available

            # Use the redundant peptide serializer
            serialized_peptides = serialize_redundant_peptides(peptides)

            # Construct query parameters
            params = {
                'protein_id': protein_id,
                'input_type': 'redundant_peptides',
                'input': serialized_peptides
            }

            # Encode the parameters
            query_string = urllib.parse.urlencode(params, safe=';/')

            # Construct the full URL
            url = f"{coverage_app_base_url}?{query_string}"

            url_list.append({
                'Protein ID': protein_id,
                'Peptides': serialized_peptides,
                'URL': url,
                'Peptide Count': len(peptides)
            })

        # Create a dataframe with URLs
        url_df = pd.DataFrame(url_list)

        # Display the dataframe
        st.subheader("Generated URLs")
        st.dataframe(url_df,
                     column_config={
            "URL": st.column_config.LinkColumn(
                "URL",
                width='small',
            ),
            "Peptide Count": st.column_config.NumberColumn(
                "Peptide Count",
            ),
        },
                     column_order=["Protein ID", "Peptide Count", "URL"],
                     use_container_width=True, hide_index=True)

    else:
        st.info("Please upload a data file to proceed.")

if __name__ == '__main__':
    main()
