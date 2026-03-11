# ============================================================================
# LIGHT CHAIN SEQUENCE ANALYSIS - INITIALIZATION
# Complete startup code with lazy-loaded ANARCI (no pre-caching)
# ============================================================================

# ============================================================================
# IMPORTS
# ============================================================================

# Data handling
import pandas as pd
import numpy as np
from io import BytesIO
from pathlib import Path

# UI
import streamlit as st

# Visualization
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')

# Biology
from Bio import SeqIO

# System utilities
import json
import warnings
import re
from typing import Tuple, Optional, Dict, List

# Suppress warnings
warnings.filterwarnings("ignore")

# ============================================================================
# STREAMLIT CONFIGURATION
# ============================================================================

st.set_page_config(
    page_title="Light chain sequence analysis",
    layout="wide",
    initial_sidebar_state="expanded"
)

# ============================================================================
# CONSTANTS
# ============================================================================

MAX_SEQUENCE_LENGTH = 1000
MIN_SEQUENCE_LENGTH = 50
N_GLYCOSYLATION_PATTERN = r'N[^P][ST]'
RARE_RESIDUE_THRESHOLD_1 = 0.001
RARE_RESIDUE_THRESHOLD_2 = 0.1
FREQUENCY_LOOKUP_CACHE_SIZE = 1000

# ============================================================================
# SESSION STATE INITIALIZATION
# ============================================================================

if "initialized" not in st.session_state:
    st.session_state.show_alignment = False
    st.session_state.alignment_df = None
    st.session_state.alignment_result = None
    st.session_state.seq_name = ""
    st.session_state.input_seq = ""
    st.session_state.initialized = True

# ============================================================================
# LAZY ANARCI LOADER
# ============================================================================

@st.cache_resource
def get_anarci():
    """
    Lazy load ANARCI only when needed.
    
    This delays ANARCI import until first alignment request.
    HMMER models download on first use (~1-2 minutes).
    Subsequent uses are instant (cached in memory).
    
    Returns:
        anarci module or None if import fails
    """
    try:
        import anarci
        return anarci
    except ImportError as e:
        st.error(f"❌ Failed to import ANARCI: {e}")
        return None

# ============================================================================
# DATABASE LOADING
# ============================================================================

def load_germline_database(fasta_file: Optional[str] = None) -> Dict[str, str]:
    """
    Load germline sequences from FASTA file.
    
    Args:
        fasta_file: Path to FASTA file (auto-searches if None)
        
    Returns:
        Dictionary mapping germline IDs to sequences
    """
    germlines = {}
    
    if fasta_file is None:
        possible_paths = [
            Path("data/imgt_iglv_iglj.fas"),
            Path("./data/imgt_iglv_iglj.fas"),
            Path("data/imgt_iglv_iglj.fasta"),
        ]
        
        for path in possible_paths:
            if path.exists():
                fasta_file = str(path)
                break
    
    if fasta_file is None or not Path(fasta_file).exists():
        return {}
    
    try:
        for record in SeqIO.parse(fasta_file, "fasta"):
            germline_id = record.id
            sequence = str(record.seq).upper()
            germlines[germline_id] = sequence
        
        return germlines
    
    except Exception as e:
        return {}


def load_frequency_matrix(matrix_file: Optional[str] = None) -> pd.DataFrame:
    """
    Load OAS residue frequency matrix.
    
    Args:
        matrix_file: Path to frequency matrix (auto-searches if None)
        
    Returns:
        DataFrame with frequency data
    """
    if matrix_file is None:
        possible_paths = [
            Path("data/oas_matrices_dash.txt"),
            Path("./data/oas_matrices_dash.txt"),
        ]
        
        for path in possible_paths:
            if path.exists():
                matrix_file = str(path)
                break
    
    if matrix_file is None or not Path(matrix_file).exists():
        return pd.DataFrame()
    
    try:
        df = pd.read_csv(matrix_file, sep="\t")
        return df
    except Exception as e:
        return pd.DataFrame()


@st.cache_resource(show_spinner=False)
def get_germline_db():
    """Load and cache germline database."""
    return load_germline_database()


@st.cache_resource(show_spinner=False)
def get_freq_matrix():
    """Load and cache frequency matrix."""
    return load_frequency_matrix()


# Initialize databases
germline_db = get_germline_db()
freq_matrix = get_freq_matrix()

# ============================================================================
# GERMLINE LOOKUP
# ============================================================================

def get_germline_sequence(germline_id: str, germline_db: Dict[str, str]) -> Optional[str]:
    """
    Get germline sequence from database by ID.
    
    Args:
        germline_id: Germline identifier
        germline_db: Germline database dictionary
        
    Returns:
        Germline sequence or None
    """
    if not germline_id:
        return None
    
    if germline_id in germline_db:
        return germline_db[germline_id]
    
    for db_id, seq in germline_db.items():
        if db_id.upper() == germline_id.upper():
            return seq
    
    base_id = germline_id.split("*")[0]
    for db_id, seq in germline_db.items():
        if base_id in db_id or db_id in base_id:
            return seq
    
    return None

# ============================================================================
# V GENE NORMALIZATION
# ============================================================================

def normalize_v_gene(v_gene: str) -> str:
    """
    Normalize V gene name by removing alleles and distal letters.
    
    Args:
        v_gene: V gene name (e.g., IGLV1D-40*01)
        
    Returns:
        Normalized name (e.g., IGLV1-40)
    """
    if not v_gene:
        return None
    
    # Remove allele (*01, *02, etc.)
    v_gene = v_gene.split("*")[0]
    
    # Remove distal letter (D, O, etc.)
    v_gene = re.sub(r'([A-Z]+\d+)[A-Z](-\d+)', r'\1\2', v_gene)
    
    return v_gene

# ============================================================================
# ANARCI ALIGNMENT
# ============================================================================

def run_anarci_alignment(sequence: str, seq_id: str = "query") -> Optional[Dict]:
    """
    Run ANARCI on input sequence with germline assignment.
    ANARCI is lazy-loaded on first call.
    
    Args:
        sequence: Protein sequence
        seq_id: Sequence identifier
        
    Returns:
        Dictionary with alignment results
    """
    try:
        # Lazy load ANARCI
        anarci = get_anarci()
        if anarci is None:
            return {"status": "error", "message": "ANARCI not available"}
        
        sequence = sequence.replace("\n", "").replace(" ", "").upper()
        
        if len(sequence) < 50:
            return {"status": "error", "message": "Sequence too short (need ~100+ AA)"}
        
        # Run alignment using lazy-loaded module
        results = anarci.anarci(
            sequences=[(seq_id, sequence)],
            scheme="imgt",
            assign_germline=True,
            allowed_species=['human', 'mouse']
        )
        
        if not results or len(results) < 2:
            return {"status": "failed", "message": "ANARCI returned unexpected structure"}
        
        alignments = results[0]
        
        if not alignments or len(alignments) == 0:
            return {"status": "failed", "message": "No alignments returned"}
        
        seq_alignment = alignments[0]
        alignment_tuple = seq_alignment[0]
        alignment_data = alignment_tuple[0]
        
        if not alignment_data:
            return {"status": "failed", "message": "No alignment data"}
        
        numbering = []
        query_seq = []
        
        for item in alignment_data:
            position_info, amino_acid = item
            numbering.append(position_info)
            query_seq.append(amino_acid)
        
        chain = "H"
        v_gene = None
        j_gene = None
        d_gene = None
        germline_species = None
        evalue = None
        bitscore = None
        query_range = None
        
        if len(results) > 1 and results[1]:
            matches = results[1]
            if matches and len(matches) > 0:
                match_info = matches[0]
                if match_info and len(match_info) > 0:
                    first_match = match_info[0]
                    if isinstance(first_match, dict):
                        chain = first_match.get('chain_type', 'H')
                        evalue = first_match.get('evalue', None)
                        bitscore = first_match.get('bitscore', None)
                        query_start = first_match.get('query_start', 0)
                        query_end = first_match.get('query_end', len(query_seq))
                        query_range = (query_start, query_end)
                        
                        germlines = first_match.get('germlines', {})
                        if germlines:
                            if 'v_gene' in germlines:
                                v_info = germlines['v_gene']
                                if v_info and len(v_info) > 0:
                                    germline_species = v_info[0][0]
                                    v_gene = v_info[0][1]
                            if 'd_gene' in germlines:
                                d_info = germlines['d_gene']
                                if d_info and len(d_info) > 0 and d_info[0]:
                                    d_gene = d_info[0][1]
                            if 'j_gene' in germlines:
                                j_info = germlines['j_gene']
                                if j_info and len(j_info) > 0 and j_info[0]:
                                    j_gene = j_info[0][1]
        
        v_sequence = get_germline_sequence(v_gene, germline_db) if v_gene else None
        j_sequence = get_germline_sequence(j_gene, germline_db) if j_gene else None
        
        germline_seq = None
        if v_sequence and j_sequence:
            germline_seq = v_sequence + j_sequence
        elif v_sequence:
            germline_seq = v_sequence
        elif j_sequence:
            germline_seq = j_sequence
        
        return {
            "query_seq": "".join(query_seq),
            "germline_seq": germline_seq,
            "numbering": numbering,
            "v_gene": v_gene,
            "d_gene": d_gene,
            "j_gene": j_gene,
            "germline_species": germline_species,
            "evalue": evalue,
            "bitscore": bitscore,
            "query_range": query_range,
            "chain": chain,
            "status": "success"
        }
    
    except Exception as e:
        import traceback
        traceback.print_exc()
        return {"status": "error", "message": f"{type(e).__name__}: {str(e)}"}

# ============================================================================
# ALIGNMENT BUILDING
# ============================================================================

def build_alignment_df(query_seq: str, germline_seq: Optional[str], numbering: list) -> pd.DataFrame:
    """
    Build alignment DataFrame for visualization.
    
    Args:
        query_seq: Query sequence
        germline_seq: Germline sequence (optional)
        numbering: IMGT numbering list
        
    Returns:
        DataFrame with alignment
    """
    data = []
    
    if germline_seq is None:
        germline_seq = "-" * len(query_seq)
    
    max_len = max(len(query_seq), len(germline_seq))
    query_seq = query_seq.ljust(max_len, "-")
    germline_seq = germline_seq.ljust(max_len, "-")
    
    for i, (pos_info, query_aa) in enumerate(zip(numbering, query_seq)):
        if pos_info and len(pos_info) >= 2:
            pos_num, pos_letter = pos_info
            imgt_pos = f"{pos_num}{pos_letter}"
        else:
            imgt_pos = "?"
        
        germ_aa = germline_seq[i] if i < len(germline_seq) else "-"
        
        data.append({
            "IMGT_Pos": imgt_pos,
            "Query": query_aa,
            "Germline": germ_aa,
            "Notes": ""
        })
    
    return pd.DataFrame(data)


def alignment_to_dict(df: pd.DataFrame) -> Dict:
    """Convert alignment DataFrame to dictionary."""
    return df.to_dict(orient="records")


def export_alignment_fasta(df: pd.DataFrame, name_prefix: str = "sequence") -> str:
    """Export alignment as FASTA format."""
    query_seq = "".join(df["Query"].astype(str))
    germline_seq = "".join(df["Germline"].astype(str))
    
    fasta = f">{name_prefix}_query\n{query_seq}\n>{name_prefix}_germline\n{germline_seq}\n"
    return fasta

# ============================================================================
# ALIGNMENT METRICS
# ============================================================================

@st.cache_data
def calculate_alignment_metrics(alignment_df: pd.DataFrame) -> Dict:
    """
    Calculate alignment quality metrics.
    
    Args:
        alignment_df: Alignment DataFrame
        
    Returns:
        Dictionary with metrics
    """
    query_seq = alignment_df["Query"].astype(str)
    germline_seq = alignment_df["Germline"].astype(str)
    
    unmatched_indels = sum(
        (query_seq == "-") != (germline_seq == "-")
    )
    
    matches = 0
    total_non_gap_matches = 0
    
    for q, g in zip(query_seq, germline_seq):
        if q == "-" and g == "-":
            continue
        total_non_gap_matches += 1
        if q == g:
            matches += 1
    
    percent_identity = (matches / total_non_gap_matches * 100) if total_non_gap_matches > 0 else 0
    
    residue_changes = sum(
        (query_seq != "-") & (germline_seq != "-") & (query_seq != germline_seq)
    )
    
    return {
        "unmatched_indels": int(unmatched_indels),
        "percent_identity": float(percent_identity),
        "residue_changes": int(residue_changes)
    }

# ============================================================================
# POST-TRANSLATIONAL MODIFICATIONS
# ============================================================================

def find_n_glycosylation_sites(sequence: str, numbering: list) -> List[Dict]:
    """
    Find potential N-glycosylation sites (N[^P][ST]).
    
    Args:
        sequence: Protein sequence
        numbering: IMGT numbering list
        
    Returns:
        List of dicts with glycosylation site information
    """
    sites = []
    
    for match in re.finditer(N_GLYCOSYLATION_PATTERN, sequence):
        pos_in_seq = match.start()
        motif = match.group()
        
        if pos_in_seq < len(numbering):
            imgt_pos = numbering[pos_in_seq]
            if imgt_pos and len(imgt_pos) >= 2:
                imgt_str = f"{imgt_pos[0]}{imgt_pos[1]}"
            else:
                imgt_str = "?"
        else:
            imgt_str = "?"
        
        sites.append({
            "Seq_Position": pos_in_seq + 1,
            "IMGT_Position": imgt_str,
            "Motif": motif,
            "X_residue": motif[1],
            "S/T_residue": motif[2]
        })
    
    return sites

# ============================================================================
# RESIDUE FREQUENCY ANALYSIS
# ============================================================================

@st.cache_data(max_entries=FREQUENCY_LOOKUP_CACHE_SIZE)
def get_residue_frequencies_cached(v_gene: str, position: int, residue: str) -> Optional[float]:
    """
    Look up residue frequency in matrix (cached).
    
    Args:
        v_gene: V gene name
        position: IMGT position number
        residue: Amino acid
        
    Returns:
        Normalized frequency (0-1) or None
    """
    if freq_matrix.empty:
        return None
    
    v_gene = normalize_v_gene(v_gene)
    
    mask = (freq_matrix['v_gene'] == v_gene) & (freq_matrix['imgt'] == position)
    
    if mask.sum() == 0:
        return None
    
    row = freq_matrix[mask].iloc[0]
    
    if residue in row.index:
        count = row[residue]
        total = row.get('total', None)
        
        if pd.isna(count) or pd.isna(total) or total == 0:
            return None
        
        return float(count) / float(total)
    
    return None


def build_frequency_analysis(alignment_df: pd.DataFrame, v_gene: str) -> pd.DataFrame:
    """
    Build frequency analysis table.
    
    Args:
        alignment_df: Alignment DataFrame
        v_gene: V gene name
        
    Returns:
        DataFrame with frequency data
    """
    data = []
    
    for idx, row in alignment_df.iterrows():
        imgt_pos_str = row["IMGT_Pos"]
        query_residue = row["Query"]
        
        try:
            imgt_pos = int(''.join(c for c in imgt_pos_str if c.isdigit()))
        except:
            imgt_pos = None
        
        frequency = None
        if imgt_pos and v_gene:
            frequency = get_residue_frequencies_cached(v_gene, imgt_pos, query_residue)
        
        data.append({
            "Position": imgt_pos_str,
            "Residue": query_residue,
            "Frequency": frequency if frequency is not None else "N/A",
            "Frequency_num": frequency
        })
    
    return pd.DataFrame(data)

# ============================================================================
# UI HELPER FUNCTIONS
# ============================================================================

def create_excel_buffer(df: pd.DataFrame, sheet_name: str = "Data") -> BytesIO:
    """Create Excel file in memory."""
    buffer = BytesIO()
    df.to_excel(buffer, index=False, sheet_name=sheet_name)
    buffer.seek(0)
    return buffer

# ============================================================================
# APP TITLE
# ============================================================================

st.title("🧪 Light chain sequence analysis")
st.markdown("**Analyze protein sequences to IMGT numbering with V/D/J gene assignment**")

# Info banner on first visit
if "info_shown" not in st.session_state:
    st.info("""
    ⏱️ **First Time?** When you submit your first sequence, ANARCI models will download (~1-2 minutes).
    This is a one-time process. All subsequent submissions will be instant!
    """)
    st.session_state.info_shown = True

# ============================================================================
# SIDEBAR
# ============================================================================

with st.sidebar:
    st.header("⚙️ Settings")
    
    st.markdown("**Instructions:**")
    st.markdown("""
    1. Paste your query sequence
    2. Click "Run ANARCI"
    3. V/D/J genes and germline will be automatically assigned
    4. Review and edit alignment if needed
    5. Check residue frequency analysis
    6. Export when satisfied
    """)
    
    st.markdown("---")
    st.markdown("**Database Status:**")
    
    if germline_db and len(germline_db) > 0:
        st.success(f"✅ {len(germline_db)} germline sequences loaded")
    else:
        st.warning("⚠️ No germline database found")
        st.caption("Expected: data/imgt_iglv_iglj.fas")
    
    if not freq_matrix.empty:
        st.success(f"✅ Residue frequency matrix loaded")
    else:
        st.warning("⚠️ Frequency matrix not loaded")
        st.caption("Expected: data/oas_matrices_dash.txt")
    
    st.markdown("---")
    st.markdown("**⚠️ Note on First Run**")
    st.markdown("""
    ANARCI models (~150 MB) download when you submit your first sequence. 
    This takes 1-2 minutes. Subsequent submissions are instant.
    """)

    
# ============================================================================
# TABS 
# ============================================================================

tab1, tab2, tab3 = st.tabs(["Input & Align", "Edit Alignment & Export", "Residue Frequency"])

# ============================================================================
# TAB 1: INPUT & ALIGNMENT
# ============================================================================

with tab1:
    st.subheader("📥 Input Sequence")
    
    # Sequence name input
    seq_name = st.text_input(
        "Sequence name (optional)",
        value=st.session_state.get("seq_name", ""),
        placeholder="e.g., Patient_001, Sample_A1",
        key="seq_name_input"
    )
    st.session_state.seq_name = seq_name
    
    col1, col2 = st.columns([0.95, 0.05])
    
    with col2:
        reset_clicked = st.button("🔄", help="Reset input", key="reset_btn")
        if reset_clicked:
            st.session_state.input_seq = ""
            st.session_state.seq_name = ""
            st.session_state.show_alignment = False
            st.session_state.alignment_df = None
            st.session_state.alignment_result = None
            st.rerun()
    
    with col1:
        with st.form("input_form"):
            input_seq = st.text_area(
                "Paste protein sequence (single letter code)",
                height=150,
                placeholder="DIQMTQSPSVSVAPGKTARISCSGDGSYNN...",
                key="input_seq",
                label_visibility="visible"
            )
            
            query_len = len(input_seq.replace("\n", "").replace(" ", ""))
            st.caption(f"Length: {query_len} AA")
            
            run_clicked = st.form_submit_button("🔄 Run ANARCI Alignment", width='stretch')
    
    if run_clicked:
        clean_query = "".join(input_seq.split()).upper()
        
        with st.spinner("Running ANARCI with germline assignment..."):
            result = run_anarci_alignment(clean_query, seq_id="query")
        
        if result["status"] == "success":
            st.success("✅ Alignment successful!")
            if result.get("v_gene"):
                genes = f"V: {result['v_gene']}"
                if result.get("j_gene"):
                    genes += f", J: {result['j_gene']}"
                if result.get("d_gene"):
                    genes += f", D: {result['d_gene']}"
                st.info(f"✅ {genes}")
            if result.get("germline_seq"):
                st.success("✅ Germline sequence found (V+J concatenated)")
            else:
                st.warning(f"⚠️ Germline sequence not found in database")
            
            st.session_state.alignment_result = result
            st.session_state.show_alignment = True
        else:
            st.error(f"❌ Alignment failed: {result.get('message', 'Unknown error')}")
    
    if st.session_state.get("show_alignment"):
        result = st.session_state.alignment_result
        seq_name_display = st.session_state.seq_name if st.session_state.seq_name else "Query"
        
        st.subheader(f"📊 Alignment Result - {seq_name_display}")
        
        col1, col2, col3, col4 = st.columns(4)
        with col1:
            st.metric("Query Length", len(result["query_seq"]))
        with col2:
            st.metric("Chain", result["chain"])
        with col3:
            st.metric("V Gene", result.get("v_gene", "N/A"))
        with col4:
            st.metric("J Gene", result.get("j_gene", "N/A"))
        
        col1, col2, col3 = st.columns(3)
        with col1:
            if result.get("d_gene"):
                st.metric("D Gene", result.get("d_gene"))
        with col2:
            if result.get("bitscore"):
                st.metric("Bit Score", f"{result['bitscore']:.1f}")
        with col3:
            if result.get("evalue"):
                st.metric("E-value", f"{result['evalue']:.2e}")
        
        st.markdown("---")
        
        alignment_df = build_alignment_df(
            result["query_seq"],
            result.get("germline_seq"),
            result["numbering"]
        )
        
        st.session_state.alignment_df = alignment_df
        
        # Calculate and display alignment metrics
        metrics = calculate_alignment_metrics(alignment_df)
        
        col1, col2, col3 = st.columns(3)
        with col1:
            st.metric("Unmatched Indels", metrics["unmatched_indels"])
        with col2:
            st.metric("Percent Identity", f"{metrics['percent_identity']:.1f}%")
        with col3:
            st.metric("Residue Changes", metrics["residue_changes"])
        
        st.markdown("---")
        
        # N-glycosylation site analysis
        st.subheader("🔬 Post-Translational Modifications")
        
        ngly_sites = find_n_glycosylation_sites(result["query_seq"], result["numbering"])
        
        if ngly_sites:
            st.success(f"🧬 Found {len(ngly_sites)} potential N-glycosylation site(s) (N[^P][ST])")
            
            ngly_df = pd.DataFrame(ngly_sites)
            st.dataframe(
                ngly_df,
                width='stretch',
                hide_index=True
            )
            
            # Show in context
            with st.expander("📍 View sites in sequence context"):
                for site in ngly_sites:
                    start = max(0, site["Seq_Position"] - 4)
                    end = min(len(result["query_seq"]), site["Seq_Position"] + 3)
                    
                    context_seq = result["query_seq"][start:end]
                    highlight_pos = site["Seq_Position"] - start - 1
                    
                    # Display with highlighting
                    before = context_seq[:highlight_pos]
                    motif = context_seq[highlight_pos:highlight_pos+3]
                    after = context_seq[highlight_pos+3:]
                    
                    st.markdown(
                        f"**Position {site['IMGT_Position']}**: "
                        f"`{before}`**`{motif}`**`{after}`"
                    )
        else:
            st.info("ℹ️ No N-glycosylation sites (N[^P][ST]) found")
        
        st.markdown("---")
        
        st.write("**Full Aligned Sequence:**")
        
        st.dataframe(
            alignment_df,
            width='stretch',
            height=500
        )
        
        st.markdown("---")
        with st.expander("📋 View as Text (FASTA)"):
            query_text = "".join(alignment_df["Query"].astype(str))
            germline_text = "".join(alignment_df["Germline"].astype(str))
            
            fasta_text = f">{seq_name_display}\n{query_text}\n>germline_V+J\n{germline_text}"
            
            st.code(fasta_text, language="fasta")

# ============================================================================
# TAB 2: MANUAL EDITING & EXPORT
# ============================================================================

with tab2:
    st.subheader("✏️ Manual Alignment Review")
    
    if st.session_state.get("alignment_df") is not None:
        df = st.session_state.alignment_df.copy()
        seq_name = st.session_state.get("seq_name", "alignment")
        
        st.info("⚠️ Please verify that any gaps in the sequence are positioned correctly")
        st.markdown("Gaps are highlighted in light red below. Click the editable table to make changes.")
        
        st.markdown("---")
        st.markdown("**Alignment with Gap Highlighting (Reference):**")
        
        # Show highlighted version as reference
        st.dataframe(
            df.style.map(lambda val: 'background-color: #ffcccc' if val == "-" else '', subset=["Query", "Germline"]),
            width='stretch',
            height=300
        )
        
        st.markdown("---")
        st.markdown("**Editable Alignment Table:**")
        
        # Editable data editor
        edited_df = st.data_editor(
            df,
            width='stretch',
            height=300,
            key="alignment_editor"
        )
        
        # Update session state with edited dataframe
        st.session_state.alignment_df = edited_df
        
        st.markdown("---")
        st.subheader("🔧 Quick Actions")
        
        col1, col2, col3 = st.columns(3)
        with col1:
            if st.button("Reset to original alignment", width='stretch'):
                if "alignment_result" in st.session_state:
                    result = st.session_state.alignment_result
                    df = build_alignment_df(
                        result["query_seq"],
                        result.get("germline_seq"),
                        result["numbering"]
                    )
                    st.session_state.alignment_df = df
                    st.info("✅ Reset to original alignment")
                    st.rerun()
        
        with col2:
            if st.button("📊 Rerun Frequency Analysis", width='stretch', type="primary"):
                st.info("✅ Frequency analysis will update with your edited alignment")
                st.rerun()
        
        with col3:
            if st.button("View statistics", width='stretch'):
                query_seq = "".join(edited_df["Query"].astype(str))
                germline_seq = "".join(edited_df["Germline"].astype(str))
                
                gaps_query = (query_seq.count("-") / len(query_seq)) * 100 if len(query_seq) > 0 else 0
                gaps_germ = (germline_seq.count("-") / len(germline_seq)) * 100 if len(germline_seq) > 0 else 0
                
                col1, col2, col3 = st.columns(3)
                with col1:
                    st.metric("Query Gap %", f"{gaps_query:.1f}%")
                with col2:
                    st.metric("Germline Gap %", f"{gaps_germ:.1f}%")
                with col3:
                    metrics = calculate_alignment_metrics(edited_df)
                    st.metric("Percent Identity", f"{metrics['percent_identity']:.1f}%")
        
        st.markdown("---")
        st.subheader("📤 Export Alignment")
        
        result = st.session_state.alignment_result
        
        col1, col2, col3, col4 = st.columns(4)
        
        with col1:
            st.download_button(
                "📋 CSV",
                data=edited_df.to_csv(index=False),
                file_name=f"{seq_name}_alignment.csv",
                mime="text/csv",
                width='stretch'
            )
        
        with col2:
            st.download_button(
                "🧬 FASTA",
                data=export_alignment_fasta(edited_df),
                file_name=f"{seq_name}_alignment.fasta",
                mime="text/plain",
                width='stretch'
            )
        
        with col3:
            st.download_button(
                "📄 JSON",
                data=json.dumps(alignment_to_dict(edited_df), indent=2),
                file_name=f"{seq_name}_alignment.json",
                mime="application/json",
                width='stretch'
            )
        
        with col4:
            excel_buffer = BytesIO()
            edited_df.to_excel(excel_buffer, index=False, sheet_name="Alignment")
            excel_buffer.seek(0)
            
            st.download_button(
                "📊 Excel",
                data=excel_buffer.getvalue(),
                file_name=f"{seq_name}_alignment.xlsx",
                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                width='stretch'
            )
    
    else:
        st.info("⏳ Run ANARCI alignment first (Tab 1)")

# ============================================================================
# TAB 3: RESIDUE FREQUENCY ANALYSIS
# ============================================================================

with tab3:
    st.subheader("📊 Residue Frequency Analysis")
    
    if st.session_state.get("alignment_df") is not None:
        df = st.session_state.alignment_df
        result = st.session_state.alignment_result
        seq_name = st.session_state.get("seq_name", "sequence")
        
        v_gene = result.get("v_gene")
        v_gene_normalized = normalize_v_gene(v_gene) if v_gene else None
        
        if not freq_matrix.empty:
            # Find first and last non-gap positions in query
            query_seq = df["Query"].astype(str)
            first_non_gap = None
            last_non_gap = None
            
            for i, res in enumerate(query_seq):
                if res != "-":
                    if first_non_gap is None:
                        first_non_gap = i
                    last_non_gap = i
            
            # Build frequency analysis (freq_matrix is global)
            freq_analysis_full = build_frequency_analysis(df, v_gene)
            
            # Filter to include only positions within the non-gap range
            if first_non_gap is not None and last_non_gap is not None:
                freq_analysis = freq_analysis_full.iloc[first_non_gap:last_non_gap+1].reset_index(drop=True)
            else:
                freq_analysis = freq_analysis_full
            
            st.markdown(f"**Sequence:** {seq_name}")
            st.markdown(f"**V Gene:** {v_gene} → {v_gene_normalized}")
            st.markdown("Showing **frequency** of residue at each IMGT position from healthy repertoire (excluding terminal gaps)")
            
            st.markdown("---")
            st.subheader("📈 Residue Frequency Distribution")
            
            # Get numeric frequencies for plotting
            freqs = freq_analysis["Frequency_num"].dropna()
            positions = freq_analysis.loc[freqs.index, "Position"]
            residues = freq_analysis.loc[freqs.index, "Residue"]
            
            if len(freqs) > 0:
                # Create bar chart
                fig, ax = plt.subplots(figsize=(16, 6))
                
                # Determine colors based on thresholds
                colors = []
                for freq in freqs:
                    if freq < RARE_RESIDUE_THRESHOLD_1:
                        colors.append('navy')
                    elif freq < RARE_RESIDUE_THRESHOLD_2:
                        colors.append('cornflowerblue')
                    else:
                        colors.append('grey')
                
                x_pos = np.arange(len(freqs))
                
                # Create bars that drop from top
                bars = ax.bar(x_pos, 1.0 - freqs, bottom=freqs, 
                             color=colors, edgecolor='black', linewidth=0.5, width=0.8)
                
                # Add query residue labels above bars
                for i, (x, residue, freq) in enumerate(zip(x_pos, residues, freqs)):
                    if freq < RARE_RESIDUE_THRESHOLD_1:
                        color = 'navy'
                    elif freq < RARE_RESIDUE_THRESHOLD_2:
                        color = 'cornflowerblue'
                    else:
                        color = 'black'
                    ax.text(x, 1.05, residue, ha='center', va='bottom', 
                           fontsize=9, fontweight='bold', color=color)
                
                # Set y-axis to log scale
                ax.set_yscale('log')
                ax.set_ylim([0.001, 1.5])
                
                ax.set_yticks([1.0, 0.1, 0.01, 0.001])
                ax.set_yticklabels(['100%', '10%', '1%', '0.1%'])
                ax.set_ylabel("Residue frequency", fontsize=12, fontweight='bold')
                
                # Add horizontal dashed lines at thresholds
                ax.axhline(y=0.1, linestyle='--', color='black', linewidth=0.8, alpha=0.5)
                ax.axhline(y=0.01, linestyle='--', color='black', linewidth=0.8, alpha=0.5)
                ax.axhline(y=0.001, linestyle='--', color='black', linewidth=0.8, alpha=0.5)
                
                # Configure x-axis - only label every 10th
                ax.set_xticks(x_pos)
                tick_labels = [positions.iloc[i] if i % 10 == 0 else '' for i in range(len(positions))]
                ax.set_xticklabels(tick_labels, fontsize=10)
                ax.set_xlabel("IMGT position", fontsize=12, fontweight='bold')
                
                ax.set_xlim([-0.5, len(freqs) - 0.5])
                
                # Title
                ax.set_title(f"Residue Frequency Profile - {seq_name} ({v_gene_normalized})", fontsize=14, fontweight='bold')
                
                # Grid
                ax.grid(axis='y', alpha=0.3)
                
                plt.tight_layout()
                st.pyplot(fig)
                
                # Export buttons for frequency profile
                st.markdown("---")
                st.subheader("📤 Export Frequency Profile")
                
                col1, col2, col3 = st.columns(3)
                
                with col1:
                    # Export as PNG
                    png_buffer = BytesIO()
                    fig.savefig(png_buffer, format='png', dpi=300, bbox_inches='tight')
                    png_buffer.seek(0)
                    st.download_button(
                        "🖼️ PNG (High-res)",
                        data=png_buffer.getvalue(),
                        file_name=f"{seq_name}_frequency_profile.png",
                        mime="image/png",
                        width='stretch'
                    )
                
                with col2:
                    # Export as PDF
                    pdf_buffer = BytesIO()
                    fig.savefig(pdf_buffer, format='pdf', bbox_inches='tight')
                    pdf_buffer.seek(0)
                    st.download_button(
                        "📄 PDF",
                        data=pdf_buffer.getvalue(),
                        file_name=f"{seq_name}_frequency_profile.pdf",
                        mime="application/pdf",
                        width='stretch'
                    )
                
                with col3:
                    # Export frequency data as CSV
                    csv_data = freq_analysis.copy()
                    csv_data["Frequency"] = csv_data["Frequency_num"].apply(lambda x: f"{x:.4f}" if isinstance(x, float) else x)
                    csv_data = csv_data[["Position", "Residue", "Frequency"]]
                    
                    st.download_button(
                        "📋 CSV Data",
                        data=csv_data.to_csv(index=False),
                        file_name=f"{seq_name}_frequency_profile.csv",
                        mime="text/csv",
                        width='stretch'
                    )
                
                st.markdown("---")
                st.subheader("📊 Summary Statistics")
                
                col1, col2, col3, col4 = st.columns(4)
                with col1:
                    st.metric("Mean Frequency", f"{freqs.mean():.2%}")
                with col2:
                    st.metric("Median Frequency", f"{freqs.median():.2%}")
                with col3:
                    st.metric("Min Frequency", f"{freqs.min():.2%}")
                with col4:
                    st.metric("Max Frequency", f"{freqs.max():.2%}")
                
                # Show positions with very rare residues
                very_rare = freq_analysis[freq_analysis["Frequency_num"] < RARE_RESIDUE_THRESHOLD_1]
                if len(very_rare) > 0:
                    st.markdown("---")
                    
                    # Proper grammar for singular/plural
                    count_text = f"{len(very_rare)} position" if len(very_rare) == 1 else f"{len(very_rare)} positions"
                    has_text = "has" if len(very_rare) == 1 else "have"
                    residue_text = "residue" if len(very_rare) == 1 else "residues"
                    
                    st.warning(f"🔴 {count_text} {has_text} very rare {residue_text} (frequency < 0.1%)")
                    st.dataframe(
                        very_rare[["Position", "Residue", "Frequency"]].assign(
                            Frequency=very_rare["Frequency_num"].apply(lambda x: f"{x:.2%}")
                        )[["Position", "Residue", "Frequency"]],
                        width='stretch'
                    )
                
                # Show positions with rare residues
                rare_positions = freq_analysis[freq_analysis["Frequency_num"] < RARE_RESIDUE_THRESHOLD_2]
                if len(rare_positions) > 0:
                    st.markdown("---")
                    
                    # Proper grammar for singular/plural
                    count_text = f"{len(rare_positions)} position" if len(rare_positions) == 1 else f"{len(rare_positions)} positions"
                    has_text = "has" if len(rare_positions) == 1 else "have"
                    residue_text = "residue" if len(rare_positions) == 1 else "residues"
                    
                    st.info(f"ℹ️ {count_text} {has_text} rare {residue_text} (frequency < 10%)")
                    st.dataframe(
                        rare_positions[["Position", "Residue", "Frequency"]].assign(
                            Frequency=rare_positions["Frequency_num"].apply(lambda x: f"{x:.2%}")
                        )[["Position", "Residue", "Frequency"]],
                        width='stretch'
                    )
                
                # Hidden frequency table
                st.markdown("---")
                with st.expander("📋 View Frequency Table"):
                    display_df = freq_analysis.copy()
                    display_df["Frequency"] = display_df["Frequency"].apply(
                        lambda x: f"{x:.2%}" if isinstance(x, float) else x
                    )
                    
                    st.dataframe(
                        display_df[["Position", "Residue", "Frequency"]],
                        width='stretch',
                        height=400
                    )
            else:
                st.warning("No frequency data available for this sequence")
        else:
            st.error("❌ Residue frequency matrix not loaded")
            st.caption("Expected: data/oas_matrices_dash.txt")
    
    else:
        st.info("⏳ Run ANARCI alignment first (Tab 1)")