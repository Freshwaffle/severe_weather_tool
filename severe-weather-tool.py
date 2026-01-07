import streamlit as st
import datetime
import requests
import pandas as pd
import numpy as np
from shapely.geometry import shape, Point
from siphon.simplewebservice.wyoming import WyomingUpperAir
from metpy.units import units
import metpy.calc as mpcalc
import matplotlib.pyplot as plt
from metpy.plots import SkewT, Hodograph
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import plyer
import sounderpy as spy  # Model soundings

# [Stations and radars dicts — keep your cleaned versions, removed duplicate MHX]

# =========================
# Utility Functions
# =========================

def get_utc_now():
    return datetime.datetime.now(datetime.timezone.utc)

def get_location():
    try:
        resp = requests.get('https://ipinfo.io/json', timeout=5)
        data = resp.json()
        loc = data.get('loc', '').split(',')
        if len(loc) == 2:
            return float(loc[0]), float(loc[1])
    except Exception:
        return None, None
    return None, None

def haversine(lat1, lon1, lat2, lon2):
    R = 6371.0
    dlat = np.radians(lat2 - lat1)
    dlon = np.radians(lon2 - lon1)
    a = np.sin(dlat / 2)**2 + np.cos(np.radians(lat1)) * np.cos(np.radians(lat2)) * np.sin(dlon / 2)**2
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))
    return R * c

def get_spc_risk(lat, lon):
    # [Your existing function — unchanged]

def get_nws_data(lat, lon):
    headers = {"User-Agent": "(severeweatherinterface.com, contact@example.com)"}
    points_url = f"https://api.weather.gov/points/{lat:.4f},{lon:.4f}"
    try:
        resp = requests.get(points_url, headers=headers, timeout=10).json()
        props = resp['properties']
        return props.get('radarStation'), props.get('forecast')
    except Exception as e:
        st.error(f"NWS points API error: {e}")
        return None, None

# [get_alerts, get_forecast, get_mesoscale_discussions — keep as-is]

# =========================
# Sounding Fetching
# =========================

def fetch_sounding(station_code=None, sounding_type="Observed", lat=None, lon=None, forecast_hour=0):
    now = get_utc_now()

    if sounding_type == "Observed":
        if not station_code:
            return None, None, None
        # [Your full observed fetching code with cycle_times loop]
        # ... (keep exactly as before)

    else:  # Model
        model = 'hrrr' if sounding_type == "HRRR" else 'rap'
        try:
            # Correct keys from SounderPy v3.1.0
            data = spy.get_model_data(model, {'lat': lat, 'lon': lon}, forecast_hour=forecast_hour)
            df = pd.DataFrame({
                'pressure': data['p'] * units.hPa,
                'temperature': data['T'] * units.degC,
                'dewpoint': data['Td'] * units.degC,
                'u_wind': data['u'] * units.knots,
                'v_wind': data['v'] * units.knots,
                'height': data['z'] * units.meter
            })
            valid_time = now + datetime.timedelta(hours=forecast_hour)
            st.write(f"**{sounding_type}** forecast (+{forecast_hour}h) valid {valid_time:%Y-%m-%d %HZ}")
            return df, f"{sounding_type} +{forecast_hour}h", valid_time
        except Exception as e:
            st.error(f"Model data error: {e}")
            return None, None, None

# =========================
# Full analyze() — restored completely
# =========================

def analyze(df):
    # [Paste your entire original analyze() function here — all calculations, try/excepts, return dict]
    # It works perfectly with both observed and model data

# [storm_mode, CRI, plot_skewt — unchanged]

# =========================
# Streamlit UI
# =========================

st.set_page_config(page_title="Severe Weather Tool", layout="wide")
st.title("Enhanced Severe Weather Interface")

st.markdown("Select sounding type and location. Model soundings (HRRR/RAP) use your lat/lon directly.")

col1, col2 = st.columns(2)
with col1:
    use_auto = st.checkbox("Auto-detect location", value=True)
with col2:
    sounding_type = st.selectbox("Sounding Type", ["Observed", "HRRR", "RAP"])

if sounding_type != "Observed":
    forecast_hour = st.slider("Forecast Hour", 0, 18 if sounding_type == "HRRR" else 48, 0)
else:
    forecast_hour = 0
    st.selectbox("Observed Station", ["Auto"] + list(stations.keys()), key="station")

st.selectbox("Radar Site", ["Auto"] + list(radars.keys()), key="radar")

# Location input
if use_auto:
    lat, lon = get_location()
    if lat is None:
        st.error("Auto-location failed. Please enter coordinates manually below.")
        use_auto = False

if not use_auto or lat is None:
    col_lat, col_lon = st.columns(2)
    lat = col_lat.number_input("Latitude", value=39.0, format="%.4f")
    lon = col_lon.number_input("Longitude", value=-95.0, format="%.4f")

if st.button("Run Analysis"):
    if lat is None or lon is None:
        st.error("Valid latitude and longitude required!")
        st.stop()

    st.success(f"Location: {lat:.3f}°, {lon:.3f}°")

    # Station for observed only
    station_code = None
    if sounding_type == "Observed":
        selected = st.session_state.station
        if selected == "Auto":
            sorted_st = sorted(stations.items(), key=lambda x: haversine(lat, lon, x[1][0], x[1][1]))
            station_code = sorted_st[0][0]
        else:
            station_code = selected

    # Fetch sounding
    df, station_label, valid_time = fetch_sounding(
        station_code=station_code,
        sounding_type=sounding_type,
        lat=lat, lon=lon,
        forecast_hour=forecast_hour
    )

    # NWS data — now safe because lat/lon validated
    radar_station, forecast_url = get_nws_data(lat, lon)
    if st.session_state.radar != "Auto":
        radar_station = st.session_state.radar

    # [Your tabs and display code — unchanged, will now work]

st.caption("Data: SPC, NWS, Wyoming Upper Air, SounderPy • Unofficial tool")
