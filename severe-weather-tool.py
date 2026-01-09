import requests
import folium
import streamlit as st
from streamlit_folium import st_folium

@st.cache_data(ttl=900)
def fetch_spc_day1_cat():
    url = "https://www.spc.noaa.gov/products/outlook/day1otlk_cat.geojson"
    r = requests.get(url, timeout=10)
    r.raise_for_status()
    return r.json()

def build_spc_map(geojson):
    m = folium.Map(location=[39, -98], zoom_start=4, tiles="cartodbpositron")

    folium.GeoJson(
        geojson,
        name="SPC Day 1 Categorical Outlook",
        style_function=lambda f: {
            "fillColor": "none",
            "color": "red",
            "weight": 2,
        },
        tooltip=folium.GeoJsonTooltip(
            fields=["risk"],
            aliases=["Risk:"],
        ),
    ).add_to(m)

    folium.LayerControl().add_to(m)
    return m

st.set_page_config(layout="wide")
st.title("Severe Weather Forecast Tool")

geojson = fetch_spc_day1_cat()
m = build_spc_map(geojson)

st_folium(m, width=1000, height=600)
