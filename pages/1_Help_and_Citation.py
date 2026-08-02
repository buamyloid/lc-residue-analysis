import streamlit as st
from pathlib import Path


def load_help_markdown() -> str:
    """Load the help content from the docs markdown file."""
    candidates = [
        Path(__file__).resolve().parents[1] / "docs" / "help.md",
        Path("docs/help.md"),
    ]

    for path in candidates:
        if path.exists():
            return path.read_text(encoding="utf-8")

    return "Help content not found."


st.set_page_config(page_title="Help & Citation", page_icon="📘", layout="wide")
st.title("Help & Citation")
st.markdown(load_help_markdown())
