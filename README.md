# PDB Coverage App

Made using [Peptacular](https://github.com/pgarrett-scripps/peptacular): [![DOI](https://zenodo.org/badge/591504879.svg)](https://doi.org/10.5281/zenodo.15054278)

This repository contains the **PDB Coverage Streamlit app**, available online at:  
🔗 [PDB Coverage App](https://pdb-cov.streamlit.app/)

Additionally, it includes two result viewers for **IP2** and **DIANN** files, which parse identified proteins and peptides to generate PDB coverage-compatible links:  
- 🔗 [DIANN Coverage Viewer](https://diann-coverage.streamlit.app/)  
- 🔗 [DTAFilter PDB Coverage Viewer](https://dtafilter-pdb-coverage.streamlit.app/)  

## Features
The PDB Coverage app supports two input modes:  
1. **UniProt ID Mode:**  
   - Retrieves the latest **AlphaFold-predicted structure** for the associated PDB file.  
   - The app is stateful in this mode, meaning the URL reflects all input values, allowing you to share the current view with others.  
2. **PDB File Upload:**  
   - Allows you to upload custom PDB files for coverage analysis.  

## Running Locally
To run the app on your local machine:  
```bash
git clone https://github.com/pgarrett-scripps/PDBCoverageStreamlitApp.git
cd PDBCoverageStreamlitApp
pip install -r requirements.txt
streamlit run app.py
```

## Citation

If you use [PDBCoverage](https://github.com/pgarrett-scripps/PDBCoverageStreamlitApp) in a publication, please cite: [![DOI](https://zenodo.org/badge/798509918.svg)](https://doi.org/10.5281/zenodo.15066418)