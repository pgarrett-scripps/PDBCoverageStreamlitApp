import streamlit as st


from app_input import get_input
from util import (
    apply_expanded_sidebar,
    coverage_string,
    display_header,
    render_mol,
    plot_coverage_array,
    show_footer,
)

st.set_page_config(layout="centered", page_title="PdbCov",
                   page_icon=":dna:", initial_sidebar_state="auto")


with st.sidebar:
    apply_expanded_sidebar()
    display_header()
    cov_input = get_input()

should_render = True
try:
    cov_input.setup()
except Exception as e:
    st.error(f"Error setting up input: {e}")
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
            vmin=cov_input.vmin,
            vmax=cov_input.vmax,
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
                            color_coverage=cov_input.color_coverage_array,
                            vmin=cov_input.vmin,
                            vmax=cov_input.vmax),
            unsafe_allow_html=True,
        )
    
show_footer()