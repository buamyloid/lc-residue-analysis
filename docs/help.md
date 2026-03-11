# Help & Citation

## Welcome to Light Chain Sequence Analysis

This web application analyzes antibody light chain sequences with IMGT numbering and VJ gene assignment.
Written mostly by Anthropic Claude with help from Gareth Morgan, Boston University Amyloidosis Center. 
Implements the residue frequency analysis described in Morgan and Prokaeva (2025) Frontiers in Immunology https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2025.1622207/full 
Please cite this paper if you publish anything that uses this analysis. 
The goal is to identify residues that are infrequently observed among healthy light chain sequences, since these may contribute to amyloid propensity.

---

## Quick Start

1. **Enter your sequence** in the "Input & Align" tab
   - Paste a protein sequence (single letter amino acid code)
   - Example: `DIQMTQSPSVSVAPGKTARISCSGDGSYNN...`

2. **Click "Run ANARCI"**
   - First time may take 1-2 minutes (downloads alignment models)
   - Subsequent submissions are instant

3. **Review results**
   - View V/D/J gene assignments
   - Check alignment metrics
   - Identify N-glycosylation sites

4. **Analyze residues** (optional)
   - Compare against healthy repertoire frequencies
   - Identify rare or unusual residues

5. **Export results**
   - Download as CSV, FASTA, JSON, or Excel

---

## Features

### Tab 1: Input & Alignment

**What it does:**
- Aligns your sequence to IMGT numbering scheme
- Identifies V, D, and J germline genes
- Detects potential N-glycosylation sites

**Metrics explained:**
- **Query Length**: Total amino acids in your sequence
- **Chain**: Antibody chain type (H or L)
- **Unmatched Indels**: Gaps that don't align with germline
- **Percent Identity**: Similarity to germline sequence
- **Residue Changes**: Number of mutations from germline

**N-Glycosylation Sites:**
- Shows the N[not P][S/T] motif
- Potential sites for protein glycosylation
- Position in IMGT numbering scheme

### Tab 2: Edit Alignment & Export

**What it does:**
- Review the aligned sequence
- Make manual edits if needed
- Export in multiple formats

**Export formats:**
- **CSV**: Spreadsheet format (Excel, Google Sheets)
- **FASTA**: Standard sequence format
- **JSON**: Structured data format
- **Excel**: Native Excel workbook

### Tab 3: Residue Frequency Analysis

**What it does:**
- Compares your residues to healthy antibody repertoires
- Shows how common each residue is
- Identifies unusual or rare positions

**Understanding the chart:**
- **Grey bars** (>10%): Common residues (normal)
- **Light blue bars** (0.1-10%): Rare residues (unusual)
- **Navy bars** (<0.1%): Very rare residues (check for errors)

**Exports:**
- High-resolution PNG image
- PDF document
- CSV data file

---

## Glossary

| Term | Meaning |
|------|---------|
| **IMGT** | International ImMunoGeneTics - standardized numbering system |
| **V gene** | Variable gene segment (most variable region) |
| **D gene** | Diversity gene segment (heavy chain only) |
| **J gene** | Joining gene segment |
| **Germline** | Original, unmutated reference sequence |
| **Alignment** | Positioning sequence to standard numbering |
| **Frequency** | How often a residue appears in healthy antibodies |

---

## Frequently Asked Questions

### Q: Why does the first alignment (sometimes) take 1-2 minutes?

**A:** The app downloads alignment models (~150 MB) on first use. These are cached afterwards, so all subsequent alignments are instant.

### Q: What sequence format should I use?

**A:** 
- Single letter amino acid codes (ACDEFGHIKLMNPQRSTVWY)
- Can include gaps (-) but they will be ignored
- Copy directly from FASTA files or sequence databases
- No headers or special characters needed

### Q: What if I get an error?

**A:** Common solutions:
- **"ANARCI initialization failed"** → Refresh the page and try again
- **"No sequence returned"** → Check your sequence is at least 50 amino acids
- **"Export not working"** → Try a different browser or disable pop-up blockers

### Q: Can I use this for heavy chains?

**A:** The aligner will work but the residue frequency data is currently only available for light chains.

### Q: What if my gene isn't recognized?

**A:** The app searches for matches in the IMGT database:
- Human and mouse genes supported
- Other species: install locally for custom databases
- Partial matches may appear if exact gene not found

---

## Limitations of the Web Version

**Available:**
- ✅ IMGT numbering alignment
- ✅ V/D/J gene assignment
- ✅ Human and mouse species
- ✅ Residue frequency analysis
- ✅ Multiple export formats
- ✅ N-glycosylation detection

**Not available in web version:**
- ❌ Custom germline databases
- ❌ Custom species/organisms
- ❌ Advanced ANARCI options
- ❌ Batch processing
- ❌ Custom frequency matrices
- ❌ Local data storage

---

## For Advanced Users

### Install Locally from GitHub

For more options and advanced features:

1. **Visit GitHub:** https://github.com/buamyloid/lc-residue-analysis

2. **Clone the repository:**
   ```bash
   git clone https://github.com/buamyloid/lc-residue-analysis.git
   cd lc-residue-analysis

3. **Install dependencies:**
   pip install -r requirements.txt

4. **Run locally:**
streamlit run app.py

Local Installation Advantages
🔧 Custom germline and frequency databases
🌍 Support for any organism
💾 Keep data private (no upload to servers)
See the GitHub repository for detailed setup instructions.

## Citation
If you use this tool in your research, please cite:

Original Manuscript:

Morgan and Prokaeva
Clone-specific residue changes at multiple positions are associated with amyloid formation by antibody light chains
Frontiers in Immunology (2025)
https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2025.1622207/full

## Key Methods Used:

ANARCI: Dunbar J, et al. (2016) SAbDab: the structural antibody database. Nucleic Acids Research
IMGT Numbering: Lefranc MP, et al. (2015) IMGT, the International ImMunoGeneTics information system. Nucleic Acids Research

## Data Privacy
### Web Version:

Your sequences are processed in real-time
Not permanently stored on servers
Results downloaded to your computer
Session data cleared when you close the app

### Local Version:

All data remains on your computer
No internet connection required (after setup)
Full privacy and control

## Getting Help
For the Web Version:
Check the FAQ above
Try refreshing the page
Clear your browser cache
Try a different browser

## Acknowledgments
This tool uses:

ANARCI for antibody numbering
IMGT reference databases
OAS (Observed Antibody Space) for frequency data
Streamlit for the web interface

Last updated: 2026
Version: 0.1