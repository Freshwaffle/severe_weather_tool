import streamlit as st
import requests
import json
import matplotlib.pyplot as plt


@st.cache_data(ttl=1800)
def get_spc_day1_cat():
    url = "https://www.spc.noaa.gov/products/outlook/day1otlk_cat.lyr.geojson"
    r = requests.get(url)
    r.raise_for_status()
    return r.json()



def plot_spc_outlook(geojson):
    fig, ax = plt.subplots(figsize=(10, 6))

    for feature in geojson["features"]:
        coords = feature["geometry"]["coordinates"]

        if feature["geometry"]["type"] == "Polygon":
            for ring in coords:
                xs = [pt[0] for pt in ring]
                ys = [pt[1] for pt in ring]
                ax.plot(xs, ys)

        elif feature["geometry"]["type"] == "MultiPolygon":
            for poly in coords:
                for ring in poly:
                    xs = [pt[0] for pt in ring]
                    ys = [pt[1] for pt in ring]
                    ax.plot(xs, ys)

    ax.set_title("SPC Day 1 Categorical Outlook")
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.grid(True)

    return fig



st.subheader("🟨 SPC Day 1 Convective Outlook")

try:
    geojson = get_spc_day1_cat()
    fig = plot_spc_outlook(geojson)
    st.pyplot(fig)
except Exception as e:
    st.error(f"SPC load failed: {e}")
