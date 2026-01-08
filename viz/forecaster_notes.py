# viz/forecaster_notes.py
import streamlit as st

def display_notes(notes=""):
    """
    Show a simple forecaster notes panel.
    """
    st.subheader("📝 Forecaster Notes")
    notes = st.text_area("Enter notes here:", value=notes, height=200)
    return notes

