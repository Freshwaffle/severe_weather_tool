import streamlit as st
import requests

st.set_page_config(page_title="Severe Weather Tool", layout="wide")

st.title("🌩️ Severe Weather Tool — Base Build")

st.markdown("""
This is a clean baseline build.
If this loads, dependencies are working.
""")

if st.button("Test NWS API"):
    try:
        r = requests.get(
            "https://api.weather.gov/points/35.23,-97.46",
            headers={"User-Agent": "(severe-weather-tool)"}
        )
        if r.status_code == 200:
            st.success("NWS API reachable ✅")
            st.json(r.json()["properties"])
        else:
            st.error(f"NWS API error: {r.status_code}")
    except Exception as e:
        st.error(str(e))
