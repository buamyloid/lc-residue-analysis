# ============================================================================
# LIGHT CHAIN SEQUENCE ANALYSIS - INITIALIZATION
# Pairwise light-chain germline alignment workflow
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
from Bio.Align import PairwiseAligner

# System utilities
import html
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
GERMLINE_BEST_OPTION = "Choose best within selected group"

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
# PAIRWISE GERMLINE ALIGNMENT
# ============================================================================

def infer_chain_type(germline_id: str) -> Optional[str]:
    """Infer chain type from germline ID."""
    if germline_id.startswith("IGK"):
        return "K"
    if germline_id.startswith("IGL"):
        return "L"
    return None


def infer_gene_segment(germline_id: str) -> Optional[str]:
    """Infer gene segment (V/J) from germline ID."""
    if "KV" in germline_id or "LV" in germline_id:
        return "V"
    if "KJ" in germline_id or "LJ" in germline_id:
        return "J"
    return None


def load_imgt_v_template_table(tsv_file: Optional[str] = None) -> pd.DataFrame:
    """Load VJ IMGT-gapped templates from TSV and annotate V-gene metadata."""
    if tsv_file is None:
        possible_paths = [
            Path("data/imgt.all.lc.alleles.txt"),
            Path("./data/imgt.all.lc.alleles.txt"),
        ]
        for path in possible_paths:
            if path.exists():
                tsv_file = str(path)
                break

    if tsv_file is None or not Path(tsv_file).exists():
        return pd.DataFrame()

    try:
        df = pd.read_csv(tsv_file, sep="\t")
    except Exception:
        return pd.DataFrame()

    df.columns = [str(col).strip().strip('"') for col in df.columns]
    for col in ["v_allele", "j_allele", "vj_sequence_imgt"]:
        if col not in df.columns:
            return pd.DataFrame()

    df["v_allele"] = df["v_allele"].astype(str).str.strip().str.strip('"')
    df["j_allele"] = df["j_allele"].astype(str).str.strip().str.strip('"')
    df["vj_sequence_imgt"] = df["vj_sequence_imgt"].astype(str).str.strip().str.strip('"')
    df = df[(df["v_allele"] != "") & (df["vj_sequence_imgt"] != "")]

    # Keep one record per unique V/J combination + IMGT-gapped template.
    df = df.drop_duplicates(subset=["v_allele", "j_allele", "vj_sequence_imgt"]).copy()
    df["v_gene"] = df["v_allele"].apply(normalize_v_gene)
    df["chain"] = df["v_allele"].apply(infer_chain_type)

    return df


@st.cache_resource(show_spinner=False)
def get_imgt_v_template_table() -> pd.DataFrame:
    """Load and cache IMGT VJ-gapped template table."""
    return load_imgt_v_template_table()


@st.cache_resource(show_spinner=False)
def get_germline_records() -> List[Dict]:
    """Build structured VJ-template records while tracking collapsed V-gene IDs."""
    records = []
    template_df = get_imgt_v_template_table()
    for _, row in template_df.iterrows():
        v_allele = str(row["v_allele"])
        j_allele = str(row["j_allele"])
        vj_imgt = str(row["vj_sequence_imgt"])
        records.append({
            "gene_id": f"{v_allele}|{j_allele}",
            "v_gene": normalize_v_gene(v_allele),
            "v_allele": v_allele,
            "j_allele": j_allele,
            "aligned_seq": vj_imgt,
            "ungapped_seq": vj_imgt.replace("-", ""),
            "chain": infer_chain_type(v_allele),
            "segment": "V",
        })
    return records


def get_gene_choices(chain_filter: Optional[str]) -> List[str]:
    """Return selectable collapsed V-gene IDs for a chain filter."""
    choices = []
    for record in get_germline_records():
        if chain_filter and record["chain"] != chain_filter:
            continue
        if record.get("v_gene"):
            choices.append(record["v_gene"])
    return sorted(set(choices))


def _alignment_to_strings(target: str, query: str, alignment) -> Tuple[str, str]:
    """Convert PairwiseAligner coordinates into aligned target/query strings."""
    target_parts: List[str] = []
    query_parts: List[str] = []

    coords = alignment.coordinates
    for i in range(coords.shape[1] - 1):
        t0, t1 = int(coords[0, i]), int(coords[0, i + 1])
        q0, q1 = int(coords[1, i]), int(coords[1, i + 1])
        t_len = t1 - t0
        q_len = q1 - q0

        if t_len > 0 and q_len > 0:
            target_parts.append(target[t0:t1])
            query_parts.append(query[q0:q1])
        elif t_len > 0:
            target_parts.append(target[t0:t1])
            query_parts.append("-" * t_len)
        elif q_len > 0:
            target_parts.append("-" * q_len)
            query_parts.append(query[q0:q1])

    return "".join(target_parts), "".join(query_parts)


def _project_alignment_to_gapped_template(
    template_gapped: str,
    aligned_template_ungapped: str,
    aligned_query: str,
) -> Tuple[str, str, List[Optional[int]]]:
    """Project ungapped alignment onto IMGT-gapped template columns.

    Returns:
        tuple(aligned_template, aligned_query, template_column_map)
        template_column_map entries are 1-based template column indices or None
        for true query insertions that fall outside the template columns.
    """
    template_positions: List[str] = []
    query_positions: List[str] = []
    template_column_map: List[Optional[int]] = []

    aligned_pairs = list(zip(aligned_template_ungapped, aligned_query))
    pair_idx = 0

    def flush_inserts(before_ungapped_index: int) -> None:
        nonlocal pair_idx
        seen_ungapped = 0
        for i in range(pair_idx):
            if aligned_pairs[i][0] != "-":
                seen_ungapped += 1
        while pair_idx < len(aligned_pairs):
            t_char, q_char = aligned_pairs[pair_idx]
            if t_char != "-":
                break
            template_positions.append("-")
            query_positions.append(q_char)
            template_column_map.append(None)
            pair_idx += 1

    ungapped_seen = 0
    flush_inserts(ungapped_seen)

    for col_idx, t_char in enumerate(template_gapped, start=1):
        if t_char == "-":
            template_positions.append("-")
            query_positions.append("-")
            template_column_map.append(col_idx)
            continue

        flush_inserts(ungapped_seen)

        if pair_idx >= len(aligned_pairs):
            template_positions.append(t_char)
            query_positions.append("-")
            template_column_map.append(col_idx)
        else:
            aligned_t_char, aligned_q_char = aligned_pairs[pair_idx]
            if aligned_t_char == "-":
                template_positions.append(t_char)
                query_positions.append("-")
                template_column_map.append(col_idx)
            else:
                template_positions.append(t_char)
                query_positions.append(aligned_q_char)
                template_column_map.append(col_idx)
                pair_idx += 1
        ungapped_seen += 1

    while pair_idx < len(aligned_pairs):
        aligned_t_char, aligned_q_char = aligned_pairs[pair_idx]
        if aligned_t_char == "-":
            template_positions.append("-")
            query_positions.append(aligned_q_char)
            template_column_map.append(None)
        pair_idx += 1

    return "".join(template_positions), "".join(query_positions), template_column_map


def run_global_pairwise_alignment(query_seq: str, germline_seq_gapped: str) -> Dict:
    """Align query to a VJ template and project onto IMGT gapped columns."""
    aligner = PairwiseAligner(mode="global")
    aligner.match_score = 2.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -8.0
    aligner.extend_gap_score = -0.5

    germline_seq_ungapped = germline_seq_gapped.replace("-", "")
    alignment = aligner.align(germline_seq_ungapped, query_seq)[0]
    aligned_template_ungapped, aligned_query_ungapped = _alignment_to_strings(
        germline_seq_ungapped,
        query_seq,
        alignment,
    )
    aligned_germline, aligned_query, template_column_map = _project_alignment_to_gapped_template(
        germline_seq_gapped,
        aligned_template_ungapped,
        aligned_query_ungapped,
    )

    non_gap_pairs = 0
    matches = 0
    for g_res, q_res in zip(aligned_germline, aligned_query):
        if g_res == "-" or q_res == "-":
            continue
        non_gap_pairs += 1
        if g_res == q_res:
            matches += 1

    identity = (matches / non_gap_pairs) if non_gap_pairs else 0.0
    gaps = aligned_germline.count("-") + aligned_query.count("-")

    return {
        "aligned_germline": aligned_germline,
        "aligned_query": aligned_query,
        "template_column_map": template_column_map,
        "score": float(alignment.score),
        "identity": identity,
        "gaps": gaps,
    }


def build_numbering_from_alignment(
    aligned_query: str,
    aligned_germline: str,
    template_column_map: Optional[List[Optional[int]]] = None,
) -> List[Tuple[int, str]]:
    """Derive numbering directly from template columns.

    Template columns (including '-' columns in vj_sequence_imgt) keep their own
    sequential numeric position. Only true extra query insertions (no template
    column) receive insertion suffixes.
    """
    numbering = []
    current_pos = 0
    insert_counts: Dict[int, int] = {}

    if template_column_map is None:
        template_column_map = [i + 1 for i in range(len(aligned_germline))]

    for q_res, col_idx in zip(aligned_query, template_column_map):
        if col_idx is not None:
            current_pos = col_idx
            numbering.append((current_pos, ""))
            insert_counts[current_pos] = 0
        else:
            anchor_pos = current_pos if current_pos > 0 else 1
            insert_counts[anchor_pos] = insert_counts.get(anchor_pos, 0) + 1
            suffix_idx = insert_counts[anchor_pos]
            suffix = chr(ord("A") + suffix_idx - 1) if suffix_idx <= 26 else f"I{suffix_idx}"
            numbering.append((anchor_pos, suffix))

    return numbering


def _select_germline_candidates(chain_filter: Optional[str], selected_gene: Optional[str]) -> List[Dict]:
    """Select germline candidates based on chain filter and optional explicit gene."""
    records = get_germline_records()
    if selected_gene:
        selected = [r for r in records if r.get("v_gene") == selected_gene]
        if chain_filter:
            selected = [r for r in selected if r["chain"] == chain_filter]
        return selected

    candidates = []
    for record in records:
        if chain_filter and record["chain"] != chain_filter:
            continue
        if record["segment"] != "V":
            continue
        candidates.append(record)
    return candidates


def run_germline_pairwise_alignment(
    sequence: str,
    chain_filter: Optional[str] = None,
    selected_gene: Optional[str] = None,
) -> Dict:
    """Align query sequence against best-matching V-gene IMGT template."""
    try:
        sequence = sequence.replace("\n", "").replace(" ", "").upper()
        if len(sequence) < MIN_SEQUENCE_LENGTH:
            return {
                "status": "error",
                "message": f"Sequence too short (need at least {MIN_SEQUENCE_LENGTH} AA)",
            }

        candidates = _select_germline_candidates(chain_filter, selected_gene)
        if not candidates:
            return {
                "status": "error",
                "message": "No compatible germline candidates found for selected options",
            }

        best_record = None
        best_alignment = None
        best_rank = None

        for record in candidates:
            alignment = run_global_pairwise_alignment(sequence, record["aligned_seq"])
            rank = (
                alignment["score"],
                alignment["identity"],
                -alignment["gaps"],
                record["gene_id"],
            )
            if best_rank is None or rank > best_rank:
                best_rank = rank
                best_record = record
                best_alignment = alignment

        if best_record is None or best_alignment is None:
            return {"status": "failed", "message": "No valid alignment generated"}

        v_gene = best_record.get("v_gene")
        j_gene = None
        chain = best_record["chain"] or "?"
        numbering = build_numbering_from_alignment(
            best_alignment["aligned_query"],
            best_alignment["aligned_germline"],
            best_alignment.get("template_column_map"),
        )

        return {
            "query_seq": best_alignment["aligned_query"],
            "germline_seq": best_alignment["aligned_germline"],
            "numbering": numbering,
            "v_gene": v_gene,
            "j_gene": j_gene,
            "selected_germline": v_gene,
            "selected_allele": best_record.get("v_allele"),
            "matched_j_allele": best_record.get("j_allele"),
            "chain": chain,
            "alignment_score": best_alignment["score"],
            "alignment_identity": best_alignment["identity"],
            "status": "success",
        }

    except Exception as e:
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


def load_help_markdown() -> str:
    """Load the help content from the docs markdown file."""
    candidates = [
        Path(__file__).resolve().parent / "docs" / "help.md",
        Path("docs/help.md"),
    ]

    for path in candidates:
        if path.exists():
            return path.read_text(encoding="utf-8")

    return "Help content not found."


def format_highlighted_fasta_html(
    query_seq: str,
    germline_seq: str,
    query_header: str,
    germline_header: str,
) -> str:
    """Render FASTA-style alignment where residue substitutions are bold red."""
    max_len = max(len(query_seq), len(germline_seq))
    query_seq = query_seq.ljust(max_len, "-")
    germline_seq = germline_seq.ljust(max_len, "-")

    query_tokens: List[str] = []
    germline_tokens: List[str] = []

    for q_res, g_res in zip(query_seq, germline_seq):
        is_change = (q_res != g_res) and (q_res != "-") and (g_res != "-")
        q_char = html.escape(q_res)
        g_char = html.escape(g_res)
        if is_change:
            query_tokens.append(f"<span style='color:#b00020; font-weight:700;'>{q_char}</span>")
            germline_tokens.append(f"<span style='color:#b00020; font-weight:700;'>{g_char}</span>")
        else:
            query_tokens.append(q_char)
            germline_tokens.append(g_char)

    query_line = "".join(query_tokens)
    germline_line = "".join(germline_tokens)

    blocks = [f"&gt;{html.escape(query_header)}"]
    blocks.append(query_line)
    blocks.append(f"&gt;{html.escape(germline_header)}")
    blocks.append(germline_line)

    html_lines = "<br/>".join(blocks)
    return (
        "<div style='font-family: monospace; white-space: nowrap; overflow-x: auto; line-height: 1.4;'>"
        + html_lines
        + "</div>"
    )

# ============================================================================
# APP TITLE
# ============================================================================

st.title("🧪 Light chain sequence analysis")
st.markdown("**Analyze protein sequences with pairwise germline alignment and IMGT-like numbering**")

# Info banner on first visit
if "info_shown" not in st.session_state:
    st.info("""
    Use **Choose best** to auto-select the highest-scoring germline,
    or select **Kappa/Lambda** and optionally lock to a specific germline gene.
    """)
    st.session_state.info_shown = True

# ============================================================================
# SIDEBAR
# ============================================================================

with st.sidebar:
    st.header("Instructions")
    template_df = get_imgt_v_template_table()
    
    st.markdown("""
    1. Paste your query sequence
    2. Select chain mode: Choose best, Kappa, or Lambda
    3. Optionally pre-specify a V gene
    4. Click "Run Pairwise Alignment"
    5. Alignment to IMGT-gapped VJ template is generated
    6. Review and edit alignment if needed
    7. Check residue frequency analysis
    8. Export when satisfied
    """)
    
    st.markdown("---")
    st.markdown("**Database Status:**")
    
    if not template_df.empty:
        n_alleles = template_df["v_allele"].nunique()
        n_genes = template_df["v_gene"].nunique()
        n_templates = len(template_df)
        st.success(f"✅ {n_templates} IMGT-gapped VJ templates loaded ({n_alleles} V alleles / {n_genes} V genes)")
    else:
        st.warning("⚠️ No IMGT V-template table found")
        st.caption("Expected: data/imgt.all.lc.alleles.txt")
    
    if not freq_matrix.empty:
        st.success(f"✅ Residue frequency matrix loaded")
    else:
        st.warning("⚠️ Frequency matrix not loaded")
        st.caption("Expected: data/oas_matrices_dash.txt")

    st.markdown("---")
    st.caption("Pairwise alignment runs against IMGT-gapped VJ templates from the bundled TSV table.")

    
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

            group_label_to_code = {
                "Choose best": None,
                "Kappa (IGK)": "K",
                "Lambda (IGL)": "L",
            }
            selected_group_label = st.selectbox(
                "Germline group",
                options=list(group_label_to_code.keys()),
                index=0,
                key="germline_group_selector",
                help="Choose best searches all light-chain germlines. Kappa/Lambda constrains candidates.",
            )
            selected_chain_filter = group_label_to_code[selected_group_label]

            available_genes = get_gene_choices(selected_chain_filter)
            gene_options = [GERMLINE_BEST_OPTION] + available_genes
            selected_gene_choice = st.selectbox(
                "Optional V-gene override",
                options=gene_options,
                index=0,
                key="germline_gene_selector",
                help="Pick a specific variable gene to constrain candidates, or keep choose best.",
            )
            
            run_clicked = st.form_submit_button("🔄 Run Pairwise Alignment", width='stretch')
    
    if run_clicked:
        clean_query = "".join(input_seq.split()).upper()
        selected_gene = None if selected_gene_choice == GERMLINE_BEST_OPTION else selected_gene_choice
        
        with st.spinner("Running pairwise alignment against germline database..."):
            result = run_germline_pairwise_alignment(
                clean_query,
                chain_filter=selected_chain_filter,
                selected_gene=selected_gene,
            )
        
        if result["status"] == "success":
            st.success("✅ Alignment successful!")
            if result.get("v_gene") or result.get("j_gene"):
                genes = []
                if result.get("v_gene"):
                    genes.append(f"V: {result['v_gene']}")
                if result.get("j_gene"):
                    genes.append(f"J: {result['j_gene']}")
                st.info(f"✅ {', '.join(genes)}")
            if result.get("selected_germline"):
                selected_v_gene = str(result["selected_germline"]).replace("*", "\\*")
                selected_text = f"Selected V gene: {selected_v_gene}"
                if result.get("selected_allele"):
                    selected_allele = str(result["selected_allele"]).replace("*", "\\*")
                    selected_text += f" (best allele: {selected_allele})"
                if result.get("matched_j_allele"):
                    matched_j_allele = str(result["matched_j_allele"]).replace("*", "\\*")
                    selected_text += f", matched J: {matched_j_allele}"
                st.info(selected_text)
            if result.get("germline_seq"):
                st.success("✅ Germline alignment generated")
            else:
                st.warning("⚠️ Germline sequence not found in database")
            
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
            st.metric("Best Allele", result.get("selected_allele", "N/A"))
        
        col1, col2, col3 = st.columns(3)
        with col1:
            st.metric("Selected Germline", result.get("selected_germline", "N/A"))
        with col2:
            if result.get("alignment_score") is not None:
                st.metric("Alignment Score", f"{result['alignment_score']:.1f}")
        with col3:
            if result.get("alignment_identity") is not None:
                st.metric("Pairwise Identity", f"{result['alignment_identity'] * 100:.1f}%")

        st.info(
            "If the alignment is not giving the expected V gene, try identifying the V gene with "
            "IMGT's DomainGapAlign tool: https://www.imgt.org/3Dstructure-DB/cgi/DomainGapAlign.cgi"
        )
        
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
        
        st.write("**FASTA Alignment (Residue Changes in Bold):**")
        query_text = "".join(alignment_df["Query"].astype(str))
        germline_text = "".join(alignment_df["Germline"].astype(str))
        fasta_html = format_highlighted_fasta_html(
            query_text,
            germline_text,
            seq_name_display,
            "germline_VJ_template",
        )
        st.markdown(fasta_html, unsafe_allow_html=True)

        st.markdown("---")
        with st.expander("📋 View Full Aligned Table"):
            st.dataframe(
                alignment_df,
                width='stretch',
                height=500
            )

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
        st.info("⏳ Run pairwise alignment first (Tab 1)")

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
            if not v_gene:
                st.warning("Residue frequency analysis needs a V-gene alignment. Choose best or select a V germline in Tab 1.")
            else:
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
                    ax.bar(
                        x_pos,
                        1.0 - freqs,
                        bottom=freqs,
                        color=colors,
                        edgecolor='black',
                        linewidth=0.5,
                        width=0.8,
                    )

                    # Add query residue labels above bars
                    for x, residue, freq in zip(x_pos, residues, freqs):
                        if freq < RARE_RESIDUE_THRESHOLD_1:
                            color = 'navy'
                        elif freq < RARE_RESIDUE_THRESHOLD_2:
                            color = 'cornflowerblue'
                        else:
                            color = 'black'
                        ax.text(x, 1.05, residue, ha='center', va='bottom', fontsize=9, fontweight='bold', color=color)

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
        st.info("⏳ Run pairwise alignment first (Tab 1)")