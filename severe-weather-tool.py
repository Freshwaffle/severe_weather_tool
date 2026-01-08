# --- Standard Library ---
import os
import sys
from datetime import datetime

# --- Third-Party Libraries ---
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import requests
from plyer import notification
import folium
from streamlit_folium import st_folium
import streamlit as st
from metpy.units import units
from metpy.calc import cape_cin
from siphon.simplewebservice.ndbc import NDBC

# --- Your Local Modules ---
from data import radar, stations
from analysis import forecasting, fronts
from viz import forecaster_notes, map as viz_map, radar_overlay


st.set_page_config(page_title="Severe Weather Tool", layout="wide")

st.title("🌩 Severe Weather Analysis Tool")

# --- Sidebar: Station Selection ---
st.sidebar.header("Station Selection")

station_type = st.sidebar.radio("Station Type", ["Surface", "Upper-Air"])

if station_type == "Surface":
    df_stations = stations.load_surface_stations()
else:
    df_stations = stations.load_upper_air_stations()

station_name = st.sidebar.selectbox(
    "Select Station",
    df_stations["name"].tolist()
)
station_row = df_stations[df_stations["name"] == station_name].iloc[0]
lat, lon = station_row["lat"], station_row["lon"]

st.sidebar.write(f"Coordinates: {lat:.2f}, {lon:.2f}")

# --- SPC Risk ---
st.subheader("⚠️ SPC Risk")
risk = forecasting.get_spc_risk(lat, lon)
st.write(f"Current risk at {station_name}: **{risk}**")

# --- Map Display ---
st.subheader("🗺 Map Overview")
m = viz_map.create_base_map(lat, lon, zoom_start=6)

# Add radar overlays for selected station (optional)
if st.checkbox("Show Radar Overlay"):
    radar_url = radar.radar_loop_url("KTLX")  # Example: radar near Oklahoma City
    m = radar_overlay.add_radar_overlay(m, radar_url)

# Add placeholder for SPC polygons (requires geojson)
# Example: user can upload a geojson
geojson_file = st.file_uploader("Upload SPC GeoJSON", type=["geojson"])
if geojson_file:
    geojson_data = pd.read_json(geojson_file)
    m = viz_map.add_spc_polygon(m, geojson_data)

# Display map
st_data = st_folium(m, width=800, height=500)

# --- Forecaster Notes ---
st.sidebar.subheader("Forecaster Notes")
notes = forecaster_notes.display_notes()

# --- Dryline / Front Proximity ---
st.subheader("🌡 Fronts / Drylines")
dryline_dist = fronts.dryline_proximity(lat, lon)
if dryline_dist:
    st.write(f"Distance to nearest dryline/front: {dryline_dist:.1f} km")
else:
    st.write("Front/dryline proximity data not available.")

# --- Optional: Show Station Data ---
if st.checkbox("Show Station Data"):
    st.write(df_stations)
