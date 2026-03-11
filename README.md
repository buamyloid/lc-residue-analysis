# Light Chain Sequence Analysis

A Streamlit web application for analyzing antibody light chain variable domain protein sequences with IMGT numbering and V(D)J gene assignment. 

Implements the residue frequency analysis described in Morgan and Prokaeva (2025) Frontiers in Immunology https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2025.1622207/full

Please cite that paper if you use the app.

## Features

- **IMGT Alignment**: Automatic sequence alignment to standardized IMGT numbering
- **V/D/J Assignment**: Identification of variable, diversity, and joining genes
- **N-Glycosylation Detection**: Identifies potential glycosylation sites (N[^P][ST])
- **Residue Frequency Analysis**: Compare against healthy antibody repertoires
- **Multiple Exports**: CSV, FASTA, JSON, Excel formats
- **Interactive Editing**: Manually review and edit alignments

## Quick Start

### Web Version (Recommended)

Visit the live app: [pending]

No installation required - just paste your sequence and click "Run ANARCI"!

**Note:** First alignment may take 1-2 minutes (one-time only). Subsequent submissions are fast.

### Local Installation

This allows you to use custom germline sequences and frequency tables, e.g. for heavy chains.

ANARCI requires HMMER. The easiest way to install this is into a dedicated Conda environment.

1. **Clone the repository:**
   ```bash
   git clone https://github.com/buamyloid/lc-residue-analysis.git
   cd lc-residue-analysis

2. **Build a Conda environment**
    ```bash
    conda create -n light-chain python=3.10

3. **Activate environment**
    ```bash
    conda activate light-chain

4. **Install dependencies (including HMMER)**
    ```bash
    conda install -c bioconda anarci
    conda install -c conda-forge streamlit pandas numpy matplotlib biopython openpyxl requests
    
5. **Run the app**
    ```bash
    streamlit run app.py

6. **Opens in web browser**
Navigate to http://localhost:8501 if Streamlit doesn't send you there automatically.