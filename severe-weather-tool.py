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
import sounderpy as spy

# =========================
# Utility Function
# =========================
def safe_float(val):
    """Safely convert a masked array or value to a float, returning NaN if invalid."""
    if val is None:
        return np.nan
    if np.ma.is_masked(val):
        return np.nan
    try:
        return float(val.magnitude if hasattr(val, 'magnitude') else val)
    except:
        return np.nan

# =========================
# Stations & Radars
# =========================
stations = {
    'ABQ': (35.04, -106.60), 'ABR': (45.45, -98.40), 'ABX': (35.15, -106.82), 'AFC': (61.27, -149.99),
    'AKQ': (37.08, -76.63), 'ALB': (42.75, -73.80), 'AMA': (35.22, -101.70), 'AMC': (38.53, -121.30),
    'AMX': (25.62, -80.42), 'APX': (44.90, -84.72), 'BIS': (46.77, -100.75), 'BMX': (33.17, -86.77),
    'BOI': (43.57, -116.22), 'BRO': (25.90, -97.43), 'BUF': (42.93, -78.73), 'CAR': (46.87, -68.02),
    'CHS': (32.90, -80.03), 'CHX': (41.67, -69.97), 'CKL': (32.90, -87.25), 'CRP': (27.77, -97.50),
    'DDC': (37.77, -99.97), 'DEN': (39.75, -104.87), 'DNR': (39.75, -104.87), 'DRT': (29.37, -100.92),
    'DTX': (42.70, -83.47), 'DVN': (41.61, -90.58), 'DYX': (32.54, -99.25), 'EKA': (40.80, -124.16),
    'EPZ': (31.90, -106.70), 'EWX': (29.70, -98.03), 'FFC': (33.36, -84.56), 'FGZ': (35.23, -111.82),
    'FWD': (32.79, -97.30), 'GGW': (48.22, -106.62), 'GJT': (39.11, -108.53), 'GRB': (44.48, -88.13),
    'GRR': (42.89, -85.54), 'GSO': (36.08, -79.95), 'GYX': (43.89, -70.26), 'HNX': (36.31, -119.63),
    'IAD': (38.95, -77.45), 'ILN': (39.42, -83.82), 'ILX': (40.15, -89.34), 'INL': (48.57, -93.38),
    'JAN': (32.32, -90.08), 'JAX': (30.48, -81.70), 'KEY': (24.55, -81.75), 'LBF': (41.13, -100.68),
    'LCH': (30.12, -93.22), 'LIX': (30.34, -89.83), 'LKN': (40.87, -115.73), 'LOT': (41.60, -88.08),
    'LZK': (34.84, -92.26), 'MAF': (31.95, -102.18), 'MFL': (25.75, -80.38), 'MHX': (34.78, -76.88),
    'MOB': (30.68, -88.24), 'MPX': (44.85, -93.57), 'MTR': (37.73, -122.28), 'OAK': (37.72, -122.22),
    'OKX': (40.87, -72.86), 'OTX': (47.68, -117.63), 'OUN': (35.23, -97.46), 'PDT': (45.70, -118.97),
    'PIT': (40.50, -80.22), 'PSR': (33.43, -112.02), 'REV': (39.57, -119.79), 'RIW': (43.07, -108.47),
    'RNK': (37.20, -80.41), 'SGX': (32.87, -117.13), 'SHV': (32.45, -93.84), 'SJT': (31.37, -100.49),
    'SLE': (44.92, -123.00), 'SLC': (40.78, -111.97), 'SGF': (37.23, -93.38), 'TBW': (27.70, -82.40),
    'TFX': (47.46, -111.38), 'TLH': (30.40, -84.35), 'TOP': (39.07, -95.62), 'TWC': (32.22, -110.96),
    'UIL': (47.95, -124.55), 'UNR': (44.07, -103.12), 'VBG': (34.73, -120.58), 'VEF': (36.05, -115.18),
    'WAL': (37.85, -75.48), 'YAK': (59.51, -139.67), 'YMO': (46.48, -84.51),
}

radars = {
    'KABR': (45.45, -98.41), 'KABX': (35.15, -106.82), 'KAKQ': (36.98, -77.00), 'KAMA': (35.22, -101.71),
    # ... (rest of your radars)
}

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
    except:
        return None, None

def haversine(lat1, lon1, lat2, lon2):
    R = 6371.0
    dlat = np.radians(lat2 - lat1)
    dlon = np.radians(lon2 - lon1)
    a = np.sin(dlat / 2)**2 + np.cos(np.radians(lat1)) * np.cos(np.radians(lat2)) * np.sin(dlon / 2)**2
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

# =========================
# NWS API Functions
# =========================
def get_nws_data(lat, lon):
    headers = {"User-Agent": "(severeweatherinterface.com, contact@example.com)"}
    points_url = f"https://api.weather.gov/points/{lat:.4f},{lon:.4f}"
    try:
        resp = requests.get(points_url, headers=headers, timeout=10).json()
        props = resp['properties']
        return props.get('radarStation'), props.get('forecast')
    except Exception as e:
        st.error(f"NWS points error: {e}")
        return None, None

def get_alerts(lat, lon):
    headers = {"User-Agent": "(severeweatherinterface.com, contact@example.com)"}
    alerts_url = f"https://api.weather.gov/alerts/active?point={lat:.4f},{lon:.4f}"
    try:
        resp = requests.get(alerts_url, headers=headers, timeout=10).json()
        features = resp.get('features', [])
        severe_alerts = [f['properties'] for f in features if 'Severe' in f['properties'].get('severity', '') or 'Tornado' in f['properties'].get('event', '') or 'Thunderstorm' in f['properties'].get('event', '')]
        return severe_alerts
    except Exception as e:
        st.error(f"Alerts error: {e}")
        return []

def get_forecast(forecast_url):
    if not forecast_url:
        return []
    headers = {"User-Agent": "(severeweatherinterface.com, contact@example.com)"}
    try:
        resp = requests.get(forecast_url, headers=headers, timeout=10).json()
        return resp['properties']['periods']
    except Exception as e:
        st.error(f"Forecast error: {e}")
        return []

def get_mesoscale_discussions():
    url = "https://www.spc.noaa.gov/products/md/"
    try:
        resp = requests.get(url, timeout=10).text
        mds = [line.strip() for line in resp.splitlines() if "Mesoscale Discussion" in line][:5]
        return mds or ["No current Mesoscale Discussions."]
    except:
        return ["Unable to load MDs."]

# =========================
# Sounding Fetching
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
        model = 'hrrr' if sounding_type == "HRRR" else 'rap'
        try:
            data = spy.get_model_data(model, [lat, lon], forecast_hour=forecast_hour)
            df = pd.DataFrame({
                'pressure': data['p'] * units.hPa,
                'temperature': data['T'] * units.degC,
                'dewpoint': data['Td'] * units.degC,
                'u_wind': data['u'] * units.knots,
                'v_wind': data['v'] * units.knots,
                'height': data['z'] * units.meter
            })
            valid_time = now + datetime.timedelta(hours=forecast_hour)
            st.write(f"Using **{sounding_type}** forecast (+{forecast_hour}h) valid {valid_time:%Y-%m-%d %HZ}")
            return df, f"{sounding_type} +{forecast_hour}h", valid_time
        except Exception as e:
            st.error(f"Model error: {e}")
            return None, None, None

# =========================
# Analysis (with fixed safe_float)
# =========================
def analyze(df):
    p = df['pressure'].values * units.hPa
    T = df['temperature'].values * units.degC
    Td = df['dewpoint'].values * units.degC
    u = (df['u_wind'].values * units.knots).to('m/s')
    v = (df['v_wind'].values * units.knots).to('m/s')
    z = df['height'].values * units.meter

    # Calculate CAPE and CIN
    mlcape, mlcin = mpcalc.mixed_layer_cape_cin(p, T, Td)
    mucape, mucin = mpcalc.most_unstable_cape_cin(p, T, Td)
    sbcape, sbcin = mpcalc.surface_based_cape_cin(p, T, Td)

    # LCL and LFC
    lcl_p, _ = mpcalc.lcl(p[0], T[0], Td[0])
    lcl_z = mpcalc.pressure_to_height_std(lcl_p)
    lfc_p, _ = mpcalc.lfc(p, T, Td)
    el_p, _ = mpcalc.el(p, T, Td)

    # Lapse rates
    def lapse(z_slice, T_slice):
        dz = np.diff(z_slice.magnitude)
        dT = np.diff(T_slice.magnitude)
        lr = -1000 * dT / dz
        return np.nanmean(lr) if len(lr) > 0 else np.nan

    lr_0_3 = lapse(z[z < 3000*units.meter], T[z < 3000*units.meter])
    lr_700_500 = lapse(z[(p <= 700*units.hPa) & (p >= 500*units.hPa)], T[(p <= 700*units.hPa) & (p >= 500*units.hPa)])
    lr_850_500 = lapse(z[(p <= 850*units.hPa) & (p >= 500*units.hPa)], T[(p <= 850*units.hPa) & (p >= 500*units.hPa)])

    # Downward CAPE
    try:
        dcape_val = float(mpcalc.downdraft_cape(p, T, Td).magnitude)
    except:
        dcape_val = np.nan

    # Storm motion and storm-relative helicity
    try:
        storm_motion, _, _ = mpcalc.bunkers_storm_motion(p, u, v, z)
        rm_speed = float(mpcalc.wind_speed(storm_motion[0], storm_motion[1]).magnitude)
    except:
        storm_motion = (np.nan, np.nan)
        rm_speed = np.nan

    try:
        srh_1 = float(mpcalc.storm_relative_helicity(z, u, v, depth=1*units.km, storm_u=storm_motion[0], storm_v=storm_motion[1])[0].magnitude)
    except:
        srh_1 = np.nan
    try:
        srh_3 = float(mpcalc.storm_relative_helicity(z, u, v, depth=3*units.km, storm_u=storm_motion[0], storm_v=storm_motion[1])[0].magnitude)
    except:
        srh_3 = np.nan
    try:
        srh_eff = float(mpcalc.storm_relative_helicity(z, u, v, depth=mpcalc.effective_layer(p, T, Td, z)[0])[0].magnitude)
    except:
        srh_eff = np.nan

    # Shear magnitudes
    try:
        shear_1 = mpcalc.bulk_shear(p, u, v, height=z, depth=1*units.km)
        shear_1_mag = safe_float(shear_1)
    except:
        shear_1_mag = np.nan

    try:
        shear_3 = mpcalc.bulk_shear(p, u, v, height=z, depth=3*units.km)
        shear_3_mag = safe_float(shear_3)
    except:
        shear_3_mag = np.nan

    try:
        shear_6 = mpcalc.bulk_shear(p, u, v, height=z, depth=6*units.km)
        shear_6_mag = safe_float(shear_6)
    except:
        shear_6_mag = np.nan

    # Wind parameters
    wind_speed = mpcalc.wind_speed(u, v).to('knots')
    wind_dir = mpcalc.wind_direction(u, v)

    # Wrap all mpcalc functions with safe_float
    mlcape = safe_float(mlcape)
    mlcin = safe_float(mlcin)
    mucape = safe_float(mucape)
    mucin = safe_float(mucin)
    sbcape = safe_float(sbcape)
    sbcin = safe_float(sbcin)
    lcl_z = safe_float(lcl_z)
    lfc_p = safe_float(lfc_p)
    el_p = safe_float(el_p)
    dcape_val = safe_float(dcape_val)
    srh_1 = safe_float(srh_1)
    srh_3 = safe_float(srh_3)
    srh_eff = safe_float(srh_eff)
    shear_1_mag = safe_float(shear_1_mag)
    shear_3_mag = safe_float(shear_3_mag)
    shear_6_mag = safe_float(shear_6_mag)

    # Calculate EHI
    try:
        ehi = mpcalc.energy_helicity_index(mlcape * units('J/kg'), srh_3 * units('m^2/s^2'))
        ehi = safe_float(ehi)
    except:
        ehi = np.nan

    # Storm motion and shear magnitude
    try:
        storm_motion, _, _ = mpcalc.bunkers_storm_motion(p, u, v, z)
        rm_speed = float(mpcalc.wind_speed(storm_motion[0], storm_motion[1]).magnitude)
    except:
        rm_speed = np.nan

    # Significant Hail Parameter (SHIP)
    try:
        freezing_level = safe_float(mpcalc.freezing_level(z, T))
        ship_param = mpcalc.significant_hail(mlcape, freezing_level * units('m'), T[p.argmin()], lr_700_500)
        ship_value = safe_float(ship_param)
    except:
        ship_value = np.nan

    # Supercell and tornado indices
    try:
        sbc = mpcalc.significant_tornado(sbcape, sbcin, lcl_z, mpcalc.bulk_shear(p, u, v, height=z, depth=6*units.km))
        sbc = safe_float(sbc)
    except:
        sbc = np.nan

    try:
        scp = mpcalc.supercell_composite(mucape, rm_speed * units('m/s'), srh_eff)
        scp = safe_float(scp)
    except:
        scp = np.nan

    # Final dictionary
    return {
        "SBCAPE": float(sbcape),
        "SBCIN": float(sbcin),
        "MLCAPE": float(mlcape),
        "MLCIN": float(mlcin),
        "MUCAPE": float(mucape),
        "MUCIN": float(mucin),
        "LCL": float(lcl_z),
        "LFC": float(lfc_p) if lfc_p is not None else np.nan,
        "EL": float(el_p) if el_p is not None else np.nan,
        "DCAPE": dcape_val,
        "LR_0_3": float(lr_0_3) if not np.isnan(lr_0_3) else np.nan,
        "LR_700_500": float(lr_700_500) if not np.isnan(lr_700_500) else np.nan,
        "LR_850_500": float(lr_850_500) if not np.isnan(lr_850_500) else np.nan,
        "SRH_1": srh_1,
        "SRH_3": srh_3,
        "SRH_EFF": srh_eff,
        "SHEAR_1": shear_1_mag,
        "SHEAR_3": shear_3_mag,
        "SHEAR_6": shear_6_mag,
        "SWEAT": safe_float(mpcalc.sweat_index(p, T, Td, wind_speed, wind_dir)),
        "K_INDEX": safe_float(mpcalc.k_index(p, T, Td)),
        "TT_INDEX": safe_float(mpcalc.total_totals_index(p, T, Td)),
        "SHOWALTER": safe_float(mpcalc.showalter_index(p, T, Td)),
        "LIFTED_INDEX": safe_float(mpcalc.lifted_index(p, T, Td)),
        "EHI": ehi,
        "SHIP": ship_value,
        "STP": sbc,
        "SCP": scp,
        "RM_SPEED": rm_speed
    }

# --- Your existing other functions and main Streamlit code ---

# Start of your main Streamlit app
st.set_page_config(page_title="Severe Weather Tool", layout="wide")
st.title("🌪️ Enhanced Severe Weather Interface")

st.markdown("Auto-detect or enter location • Choose observed or model sounding (HRRR/RAP)")

col1, col2 = st.columns(2)
with col1:
    use_auto = st.checkbox("Auto-detect my location", value=True)
with col2:
    sounding_type = st.selectbox("Sounding Type", ["Observed", "HRRR", "RAP"])

if sounding_type != "Observed":
    forecast_hour = st.slider("Forecast Hour", 0, 18 if sounding_type == "HRRR" else 48, 0)
else:
    forecast_hour = 0

selected_station = st.selectbox("Observed Station (Auto = nearest)", ["Auto"] + list(stations.keys()))
selected_radar = st.selectbox("Radar Site (Auto = NWS default)", ["Auto"] + list(radars.keys()))

if use_auto:
    lat, lon = get_location()
    if lat is None:
        st.warning("Auto-location failed — enter manually below")
        use_auto = False

if not use_auto or lat is None:
    col_lat, col_lon = st.columns(2)
    lat = col_lat.number_input("Latitude", value=35.23, format="%.4f")
    lon = col_lon.number_input("Longitude", value=-97.46, format="%.4f")

if st.button("🚀 Run Analysis", type="primary"):
    if lat is None or lon is None:
        st.error("Valid coordinates required!")
    else:
        st.success(f"Analyzing: {lat:.3f}°, {lon:.3f}")

        station_code = None
        if sounding_type == "Observed":
            if selected_station == "Auto":
                distances = {code: haversine(lat, lon, coords[0], coords[1]) for code, coords in stations.items()}
                station_code = min(distances, key=distances.get)
            else:
                station_code = selected_station

        df, station_label, _ = fetch_sounding(station_code, sounding_type, lat, lon, forecast_hour)

        radar_station, forecast_url = get_nws_data(lat, lon)
        if selected_radar != "Auto":
            radar_station = selected_radar

        tab1, tab2, tab3, tab4, tab5, tab6 = st.tabs(["Overview", "Parameters", "Skew-T", "Radar & Sat", "Alerts & MDs", "Forecast"])

        with tab1:
            risk = get_spc_risk(lat, lon)
            st.write(f"**SPC Day 1 Risk:** {risk}")
            st.image("https://www.spc.noaa.gov/products/outlook/day1otlk.gif")

            if df is not None:
                p = analyze(df)
                st.write(f"**Storm Mode:** {storm_mode(p)}")
                st.write(f"**CRI Score:** {CRI(p)} / 10")

                key_df = pd.DataFrame({
                    "Parameter": ["MLCAPE", "SRH 0-1km", "Shear 0-6km", "DCAPE", "EHI", "SHIP", "STP", "SCP"],
                    "Value": [p[k] for k in ["MLCAPE", "SRH_1", "SHEAR_6", "DCAPE", "EHI", "SHIP", "STP", "SCP"]],
                    "Units": ["J/kg", "m²/s²", "kt", "J/kg", "", "", "", ""]
                }).round(1)
                st.table(key_df)

                if CRI(p) >= 6 or risk in ["ENH", "MDT", "HIGH"]:
                    st.error("⚠️ Elevated severe weather risk — stay alert!")
                    try:
                        plyer.notification.notify(title="Severe Weather", message=f"{risk} | CRI {CRI(p)}", timeout=10)
                    except:
                        pass

        with tab2:
            st.subheader("🔍 Detailed Parameter Breakdown")
            if df is not None:
                p = analyze(df)

                st.markdown("#### All Calculated Parameters")
                df_params = pd.DataFrame.from_dict(p, orient='index', columns=['Value'])
                df_params = df_params.round(1)
                df_params.index.name = "Parameter"
                df_params = df_params.reset_index()

                def get_color(param, val):
                    if pd.isna(val):
                        return ''
                    if param in ["MLCAPE", "MUCAPE", "SBCAPE"] and val > 2000:
                        return 'background-color: #ff9999'
                    if param == "SRH_1" and val > 150:
                        return 'background-color: #ffcc99'
                    if param == "SHEAR_6" and val > 40:
                        return 'background-color: #ffeb99'
                    if param == "DCAPE" and val > 1000:
                        return 'background-color: #ffeb99'
                    if param == "STP" and val > 1:
                        return 'background-color: #ff9999'
                    if param == "SCP" and val > 2:
                        return 'background-color: #ffcc99'
                    if param == "SHIP" and val > 1:
                        return 'background-color: #ffff99'
                    return ''

                styled = df_params.style.apply(lambda row: [get_color(row['Parameter'], row['Value']) if i == 1 else '' for i in range(len(row))], axis=1)
                st.dataframe(styled, width='stretch')

                st.markdown("---")
                st.markdown("#### Parameter Interpretations")
                explanations = {
                    # your explanations here...
                }
                categories = {
                    # your categories here...
                }
                for category_name, params in categories.items():
                    st.markdown(f"**{category_name}**")
                    cols = st.columns(3)
                    for i, param in enumerate(params):
                        with cols[i % 3]:
                            val = p.get(param, np.nan)
                            display_val = "N/A" if val is None or np.isnan(val) else f"{float(val):.1f}"
                            with st.expander(f"{param}: {display_val}"):
                                st.write(explanations.get(param, "No detailed explanation available."))

            else:
                st.info("Run analysis to view detailed parameters.")

        with tab3:
            if df is not None:
                plot_skewt(df, station_label)
                st.image("skewt_hodograph.png")
            else:
                st.info("No sounding to display")

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
            for md in get_mesoscale_discussions():
                st.write(md)

        with tab6:
            periods = get_forecast(forecast_url)
            for p in periods[:7]:
                st.write(f"**{p['name']}**: {p['detailedForecast']}")

st.caption("Unofficial tool • Data: NWS, SPC, Wyoming, SounderPy • Built for storm season 2026")
