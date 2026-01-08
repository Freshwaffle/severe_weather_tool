import streamlit as st
import datetime
import requests
import pandas as pd
import numpy as np
from shapely.geometry import shape, Point
from siphon.simplewebservice.wyoming import WyomingUpperAir
from metpy.units import units
import metpy.calc as mpcalc
import sounderpy as spy

# =========================
# Stations Dictionary (defined early)
# =========================
stations = {
    'ABQ': (35.04, -106.60),
    'ABR': (45.45, -98.40),
    'ABX': (35.15, -106.82),
    'AKQ': (36.97, -76.30),
    'ALY': (42.66, -73.80),
    'AMA': (35.22, -101.70),
    'APX': (43.75, -86.25),
    'BMX': (33.33, -86.75),
    'BOU': (39.75, -105.00),
    'BIS': (46.80, -100.78),
    'BRO': (27.97, -97.38),
    'BUF': (42.17, -78.88),
    'BTV': (44.50, -73.15),
    'CAE': (33.98, -81.12),
    'CAR': (35.22, -80.84),
    'CHS': (32.90, -80.03),
    'CLE': (41.42, -81.86),
    'CWA': (44.78, -91.50),
    'DDC': (37.74, -98.03),
    'DLH': (46.84, -92.20),
    'DMX': (41.53, -93.59),
    'DVN': (41.60, -90.58),
    'FSD': (43.55, -96.73),
    'FGF': (44.09, -96.21),
    'GID': (43.21, -96.97),
    'GRR': (42.88, -85.52),
    'GSP': (35.22, -82.43),
    'HUN': (34.73, -92.20),
    'ICT': (37.65, -97.42),
    'ILX': (40.10, -88.24),
    'IND': (39.71, -86.29),
    'IWX': (41.55, -86.38),
    'JAN': (32.31, -90.11),
    'JKL': (37.30, -83.08),
    'KEY': (25.70, -80.20),
    'LUB': (33.55, -101.88),
    'LWX': (38.84, -77.02),
    'MEG': (44.33, -70.72),
    'MKX': (43.05, -88.55),
    'MPX': (45.03, -93.25),
    'MQT': (46.54, -87.55),
    'OHX': (36.00, -86.75),
    'OKX': (40.88, -73.88),
    'OUN': (35.22, -97.47),
    'PAH': (38.63, -90.15),
    'PBZ': (39.20, -78.88),
    'PHI': (39.88, -75.23),
    'PIH': (43.52, -112.04),
    'PUB': (38.28, -104.49),
    'RAH': (33.94, -84.52),
    'REV': (49.05, -122.28),
    'RIW': (43.33, -108.23),
    'RLX': (38.30, -81.61),
    'RNK': (37.44, -80.88),
    'SAW': (38.85, -77.05),
    'SGF': (37.22, -93.38),
    'SHV': (32.38, -94.78),
    'SJT': (31.36, -100.49),
    'SLC': (40.78, -111.98),
    'TOP': (36.36, -97.39),
    'TSA': (36.12, -95.94),
    'UNR': (39.52, -119.81),
    'VTX': (39.28, -119.01),
    'WFO': (46.59, -112.02),
    'YKN': (42.00, -97.40),
    'YUM': (40.70, -122.88),
}

# =========================
# Utility & Helper Functions
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
    except:
        return None, None

def haversine(lat1, lon1, lat2, lon2):
    R = 6371.0
    dlat = np.radians(lat2 - lat1)
    dlon = np.radians(lon2 - lon1)
    a = np.sin(dlat/2)**2 + np.cos(np.radians(lat1)) * np.cos(np.radians(lat2)) * np.sin(dlon/2)**2
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))
    return R * c

def get_spc_risk(lat, lon):
    url = "https://mapservices.weather.noaa.gov/vector/rest/services/outlooks/SPC_wx_outlks/MapServer/1/query"
    params = {"where": "1=1", "outFields": "*", "f": "geojson", "outSR": "4326"}
    try:
        data = requests.get(url, params=params, timeout=10).json()
    except:
        return "Error fetching SPC outlook"

    pt = Point(lon, lat)
    hierarchy = ["TSTM", "MRGL", "SLGT", "ENH", "MDT", "HIGH"]
    found = []
    for f in data.get("features", []):
        try:
            poly = shape(f["geometry"])
            if poly.contains(pt):
                label = f["properties"].get("LABEL", "TSTM")
                found.append(label)
        except:
            continue
    if not found:
        return "No Risk"
    found.sort(key=lambda x: hierarchy.index(x))
    return found[-1]

def get_nws_data(lat, lon):
    headers = {"User-Agent": "(severeweatherinterface.com, contact@example.com)"}
    points_url = f"https://api.weather.gov/points/{lat:.4f},{lon:.4f}"
    try:
        resp = requests.get(points_url, headers=headers, timeout=10).json()
        return resp['properties'].get('radarStation'), resp['properties'].get('forecast')
    except:
        return None, None

def get_alerts(lat, lon):
    headers = {"User-Agent": "(severeweatherinterface.com, contact@example.com)"}
    alerts_url = f"https://api.weather.gov/alerts/active?point={lat:.4f},{lon:.4f}"
    try:
        resp = requests.get(alerts_url, headers=headers, timeout=10).json()
        features = resp.get('features', [])
        severe_alerts = [f['properties'] for f in features if 'Severe' in f['properties'].get('severity', '') or 'Tornado' in f['properties'].get('event', '') or 'Thunderstorm' in f['properties'].get('event', '')]
        return severe_alerts
    except:
        return []

def get_forecast(forecast_url):
    if not forecast_url:
        return []
    headers = {"User-Agent": "(severeweatherinterface.com, contact@example.com)"}
    try:
        resp = requests.get(forecast_url, headers=headers, timeout=10).json()
        return resp['properties']['periods']
    except:
        return []

# =========================
# Fetch Sounding Function (using sounderpy)
# =========================
def fetch_sounding(station_code=None, sounding_type="Observed", lat=None, lon=None, forecast_hour=0):
    now = get_utc_now()

    if sounding_type == "Observed":
        if not station_code:
            return None, None, None
        cycle_hour = 12 if now.hour >= 12 else 0
        cycle_times = [
            now.replace(hour=cycle_hour, minute=0, second=0, microsecond=0),
            (now - datetime.timedelta(hours=12)).replace(hour=(0 if cycle_hour == 12 else 12), minute=0, second=0, microsecond=0)
        ]
        for cycle in cycle_times:
            try:
                df = WyomingUpperAir.request_data(cycle, station_code)
                df = df.dropna(subset=["pressure", "temperature", "dewpoint", "u_wind", "v_wind", "height"])
                if len(df) >= 8:
                    st.write(f"Using observed sounding from **{station_code}** at {cycle:%Y-%m-%d %HZ}")
                    return df, station_code, cycle
            except:
                continue
        st.warning("No recent observed sounding available.")
        return None, None, None
    else:
        # Use correct model names supported by sounderpy
        model = 'rap-ruc' if sounding_type == "HRRR" else 'ncep'
        forecast_time = now + datetime.timedelta(hours=forecast_hour)

        # NOAA data typically available within last 48 hours
        max_hours = 48
        earliest_time = now - datetime.timedelta(hours=max_hours)

        if forecast_time > now:
            st.error("Forecast time is in the future. Please select a more recent forecast.")
            return None, None, None
        if forecast_time < earliest_time:
            st.error(f"Forecast time {forecast_time} is too far in the past. Data available for last {max_hours} hours.")
            return None, None, None

        # Extract date parts
        year = forecast_time.year
        month = forecast_time.month
        day = forecast_time.day
        hour = forecast_time.hour

        try:
            data = spy.get_model_data(model, [lat, lon], year=year, month=month, day=day, hour=hour)
            df = pd.DataFrame({
                'pressure': data['p'] * units.hPa,
                'temperature': data['T'] * units.degC,
                'dewpoint': data['Td'] * units.degC,
                'u_wind': data['u'] * units.knots,
                'v_wind': data['v'] * units.knots,
                'height': data['z'] * units.meter
            })
            st.write(f"Using **{sounding_type}** forecast (+{forecast_hour}h) valid {forecast_time:%Y-%m-%d %HZ}")
            return df, f"{sounding_type} +{forecast_hour}h", forecast_time
        except Exception as e:
            st.error(f"Model error: {e}")
            return None, None, None

# =========================
# Main Streamlit App
# =========================
st.set_page_config(page_title="Severe Weather Tool", layout="wide")
st.title("🌪️ Enhanced Severe Weather Interface")

st.markdown("Auto-detect or enter location • Choose observed or model sounding (HRRR/RAP)")

# Location input and detection
use_auto = st.checkbox("Auto-detect my location", value=True)
lat, lon = None, None

if use_auto:
    lat, lon = get_location()
    if lat is None:
        st.warning("Auto-location failed — please enter manually.")
        use_auto = False

if not use_auto:
    col_lat, col_lon = st.columns(2)
    lat = col_lat.number_input("Latitude", value=35.23, format="%.4f")
    lon = col_lon.number_input("Longitude", value=-97.46, format="%.4f")

# Sounding type selection
sounding_type = st.selectbox("Sounding Type", ["Observed", "HRRR", "RAP"])
if sounding_type != "Observed":
    forecast_hour = st.slider("Forecast Hour", 0, 48, 0)
else:
    forecast_hour = 0

# Station and Radar selection
selected_station = st.selectbox("Observed Station (Auto = nearest)", ["Auto"] + list(stations.keys()))
selected_radar = st.selectbox("Radar Site (Auto = NWS default)", ["Auto"])

# Run analysis button
if st.button("🚀 Run Analysis", type="primary"):
    if lat is None or lon is None:
        st.error("Valid coordinates required!")
    else:
        st.success(f"Analyzing: {lat:.3f}°, {lon:.3f}")

        # Determine station code if needed
        station_code = None
        if sounding_type == "Observed":
            if selected_station == "Auto":
                distances = {code: haversine(lat, lon, coords[0], coords[1]) for code, coords in stations.items()}
                station_code = min(distances, key=distances.get)
            else:
                station_code = selected_station

        # Fetch sounding data
        df, station_label, valid_time = fetch_sounding(station_code, sounding_type, lat, lon, forecast_hour)

        # Get NWS radar and forecast info
        radar_station, forecast_url = get_nws_data(lat, lon)
        if selected_radar != "Auto":
            radar_station = selected_radar

        # Tabs for different views
        tab1, tab2, tab3, tab4, tab5, tab6 = st.tabs(["Overview", "Parameters", "Skew-T", "Radar & Sat", "Alerts & MDs", "Forecast"])

        with tab1:
            risk = get_spc_risk(lat, lon)
            st.write(f"**SPC Day 1 Risk:** {risk}")
            st.image("https://www.spc.noaa.gov/products/outlook/day1otlk.gif")
            if df is not None:
                # Your analyze() function should be here; assuming you have it.
                # For demo, display info
                st.write("Sounding loaded successfully.")
                # You can call your analyze() function here if defined
            else:
                st.info("Run analysis to view parameters.")

        with tab2:
            if df is not None:
                # Your detailed parameters display
                st.write("Detailed parameters display here.")
            else:
                st.info("No sounding data to analyze.")

        with tab3:
            if df is not None:
                # Plot skew-t or display info
                st.write("Skew-T plot placeholder.")
            else:
                st.info("No sounding to plot.")

        with tab4:
            if radar_station:
                st.image(f"https://radar.weather.gov/ridge/standard/{radar_station}_loop.gif", caption=f"{radar_station} Radar")
            st.image("https://cdn.star.nesdis.noaa.gov/GOES19/ABI/CONUS/GEOCOLOR/GOES19-CONUS-GEOCOLOR-625x375.gif", caption="GOES-19 GeoColor")

        with tab5:
            alerts = get_alerts(lat, lon)
            if alerts:
                for a in alerts:
                    st.write(f"**{a['event']}** — {a['headline']}")
            else:
                st.info("No active severe alerts")
            st.subheader("SPC Mesoscale Discussions")
            # Placeholder for MDs
            st.write("No MDs loaded.")

        with tab6:
            periods = get_forecast(forecast_url)
            for p in periods[:7]:
                st.write(f"**{p['name']}**: {p['detailedForecast']}")

st.caption("Unofficial tool • Data: NWS, SPC, Wyoming, SounderPy • Built for storm season 2026")
