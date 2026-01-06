import streamlit as st
import datetime
import requests
import pandas as pd
import numpy as np
import geocoder
from shapely.geometry import shape, Point
from siphon.simplewebservice.wyoming import WyomingUpperAir
from metpy.units import units
import metpy.calc as mpcalc
import matplotlib.pyplot as plt
from metpy.plots import SkewT, Hodograph
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import plyer

# Comprehensive list of upper air stations (expanded from searches)
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
    'MHX': (34.78, -76.88), 'MOB': (30.68, -88.24), 'MPX': (44.85, -93.57), 'MTR': (37.73, -122.28),
    'OAK': (37.72, -122.22), 'OKX': (40.87, -72.86), 'OTX': (47.68, -117.63), 'OUN': (35.23, -97.46),
    'PDT': (45.70, -118.97), 'PIT': (40.50, -80.22), 'PSR': (33.43, -112.02), 'REV': (39.57, -119.79),
    'RIW': (43.07, -108.47), 'RNK': (37.20, -80.41), 'SGX': (32.87, -117.13), 'SHV': (32.45, -93.84),
    'SJT': (31.37, -100.49), 'SLE': (44.92, -123.00), 'SLC': (40.78, -111.97), 'SGF': (37.23, -93.38),
    'TBW': (27.70, -82.40), 'TFX': (47.46, -111.38), 'TLH': (30.40, -84.35), 'TOP': (39.07, -95.62),
    'TWC': (32.22, -110.96), 'UIL': (47.95, -124.55), 'UNR': (44.07, -103.12), 'VBG': (34.73, -120.58),
    'VEF': (36.05, -115.18), 'WAL': (37.85, -75.48), 'YAK': (59.51, -139.67), 'YMO': (46.48, -84.51),
    # Add more if needed from searches
}

# List of NWS radar sites (from searches, e.g., WSR-88D)
radars = {
    'KABR': (45.45, -98.41), 'KABX': (35.15, -106.82), 'KAKQ': (36.98, -77.00), 'KAMA': (35.22, -101.71),
    'KAMX': (25.61, -80.41), 'KAPX': (44.91, -84.72), 'KARX': (43.82, -91.19), 'KATX': (47.68, -122.50),
    'KBBX': (38.82, -121.63), 'KBGM': (42.20, -75.98), 'KBHX': (40.50, -124.29), 'KBIS': (46.77, -100.76),
    'KBLX': (45.85, -108.61), 'KBMX': (33.17, -86.77), 'KBOX': (41.96, -71.14), 'KBRO': (25.92, -97.42),
    'KBUF': (42.95, -78.74), 'KBYX': (24.60, -81.70), 'KCAE': (33.95, -81.12), 'KCBW': (46.03, -67.81),
    'KCBX': (43.49, -116.24), 'KCCX': (40.92, -78.00), 'KCLE': (41.41, -81.86), 'KCLX': (32.66, -81.04),
    'KCRP': (27.78, -97.51), 'KCXX': (44.51, -73.17), 'KDAX': (38.50, -121.68), 'KDDC': (37.76, -99.97),
    'KDGX': (32.28, -89.98), 'KDLH': (46.84, -92.21), 'KDMX': (41.73, -93.72), 'KDOX': (38.83, -75.44),
    'KDTX': (42.70, -83.47), 'KDVN': (41.61, -90.58), 'KDYX': (32.54, -99.25), 'KEAX': (38.81, -94.26),
    'KEMX': (31.89, -110.63), 'KENX': (42.59, -74.06), 'KEOX': (31.46, -85.46), 'KEPZ': (31.87, -106.70),
    'KEVX': (30.56, -85.92), 'KEWX': (29.70, -98.03), 'KEYX': (35.10, -117.56), 'KFCX': (37.02, -80.27),
    'KFDX': (34.63, -103.63), 'KFFC': (33.36, -84.57), 'KFDR': (34.36, -98.98), 'KFSX': (35.29, -111.82),
    'KFTG': (39.79, -104.55), 'KFWS': (32.57, -97.30), 'KGGW': (48.21, -106.63), 'KGJX': (39.06, -108.13),
    'KGLD': (39.07, -101.70), 'KGRB': (44.50, -88.11), 'KGRK': (30.72, -97.38), 'KGRR': (42.89, -85.55),
    'KGSP': (34.88, -82.22), 'KGUAM': (13.46, 144.81), 'KGWX': (32.67, -90.90), 'KGYX': (43.89, -70.26),
    'KHDX': (33.08, -106.12), 'KHGX': (29.47, -95.08), 'KHNX': (36.31, -119.63), 'KHPX': (36.74, -87.29),
    'KHTX': (34.93, -86.08), 'KICT': (37.65, -97.44), 'KICX': (37.59, -112.86), 'KILN': (39.42, -83.82),
    'KILX': (40.15, -89.34), 'KIND': (39.71, -86.28), 'KINX': (36.18, -95.56), 'KIWA': (33.29, -111.67),
    'KIWX': (41.41, -85.70), 'KJAN': (32.32, -90.08), 'KJAX': (30.48, -81.70), 'KJGX': (32.67, -83.35),
    'KJKL': (37.59, -83.31), 'KLBB': (33.65, -101.81), 'KLCH': (30.13, -93.22), 'KLGX': (47.12, -124.11),
    'KLIX': (30.34, -89.83), 'KLNX': (41.96, -100.58), 'KLOT': (41.60, -88.08), 'KLOX': (34.20, -119.13),
    'KLRX': (40.74, -116.80), 'KLSX': (38.70, -90.68), 'KLTX': (33.99, -78.43), 'KLVX': (37.98, -85.94),
    'KLWX': (38.98, -77.48), 'KLZK': (34.84, -92.26), 'KMAF': (31.94, -102.19), 'KMAX': (42.08, -122.72),
    'KMBX': (48.39, -100.86), 'KMHX': (34.78, -76.88), 'KMKX': (42.97, -88.55), 'KMLB': (28.11, -80.65),
    'KMOB': (30.68, -88.24), 'KMPX': (44.85, -93.57), 'KMQT': (46.53, -87.55), 'KMSX': (47.04, -113.99),
    'KMTX': (41.26, -112.45), 'KMUX': (37.16, -121.90), 'KMVX': (47.53, -97.53), 'KMXX': (32.54, -85.79),
    'KNKX': (32.92, -117.04), 'KNQA': (35.34, -89.87), 'KOAX': (41.32, -96.37), 'KOHX': (36.25, -86.56),
    'KOKX': (40.87, -72.86), 'KOTX': (47.68, -117.63), 'KOUN': (35.24, -97.46), 'KPAH': (37.07, -88.77),
    'KPBZ': (40.53, -80.22), 'KPDT': (45.69, -118.99), 'KPDX': (45.60, -122.60), 'KPOE': (31.16, -92.98),
    'KPSR': (33.43, -112.02), 'KPUX': (38.46, -104.18), 'KRAX': (35.67, -78.49), 'KRIW': (43.07, -108.48),
    'KRLX': (38.31, -81.72), 'KRTX': (45.72, -122.97), 'KSFX': (43.11, -112.69), 'KSGF': (37.24, -93.40),
    'KSHV': (32.45, -93.84), 'KSJT': (31.37, -100.49), 'KSOX': (33.82, -117.64), 'KSRX': (35.29, -94.36),
    'KTLH': (30.40, -84.33), 'KTLX': (35.33, -97.28), 'KTWX': (38.99, -96.23), 'KTYX': (43.76, -75.68),
    'KUDX': (44.12, -102.83), 'KUEX': (40.32, -98.44), 'KVAX': (30.90, -83.00), 'KVBX': (34.84, -120.40),
    'KVNX': (36.74, -98.13), 'KVTX': (34.41, -119.18), 'KVWX': (37.59, -87.72), 'KYUX': (43.89, -75.03),
    # Add more if needed
}

# -------------------------
# Utility functions
# -------------------------

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
    pass
return None, None

def get_spc_risk(lat, lon):
    url = "https://mapservices.weather.noaa.gov/vector/rest/services/outlooks/SPC_wx_outlks/MapServer/1/query"
    params = {"where": "1=1", "outFields": "*", "f": "geojson", "outSR": "4326"}
    try:
        data = requests.get(url, params=params, timeout=10).json()
    except Exception as e:
        return f"Error fetching SPC: {str(e)}"

    pt = Point(lon, lat)
    hierarchy = ["TSTM", "MRGL", "SLGT", "ENH", "MDT", "HIGH"]
    found = []

    for f in data.get("features", []):
        try:
            poly = shape(f["geometry"])
            if poly.contains(pt):
                found.append(f["properties"].get("LABEL", "TSTM"))
        except:
            continue

    if not found:
        return "No Risk"
    found.sort(key=lambda x: hierarchy.index(x))
    return found[-1]

def haversine(lat1, lon1, lat2, lon2):
    R = 6371.0
    dlat = np.radians(lat2 - lat1)
    dlon = np.radians(lon2 - lon1)
    a = np.sin(dlat / 2)**2 + np.cos(np.radians(lat1)) * np.cos(np.radians(lat2)) * np.sin(dlon / 2)**2
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))
    return R * c

def fetch_sounding(station_code):
    now = get_utc_now()
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
                st.write(f"Using sounding from {station_code} at {cycle:%Y-%m-%d %HZ}")
                return df, station_code, cycle
        except:
            continue
    return None, None, None

def analyze(df):
    p = df.pressure.values * units.hPa
    T = df.temperature.values * units.degC
    Td = df.dewpoint.values * units.degC
    u = (df.u_wind.values * units.knots).to(units.meter / units.second)
    v = (df.v_wind.values * units.knots).to(units.meter / units.second)
    z = df.height.values * units.meter
    # Basic parameters
    mlcape, mlcin = mpcalc.mixed_layer_cape_cin(p, T, Td)
    mucape, mucin = mpcalc.most_unstable_cape_cin(p, T, Td)
    sbcape, sbcin = mpcalc.surface_based_cape_cin(p, T, Td)
    lcl_p, lcl_T = mpcalc.lcl(p[0], T[0], Td[0])
    lcl_z = mpcalc.pressure_to_height_std(lcl_p)
    lfc_p, lfc_T = mpcalc.lfc(p, T, Td)
    el_p, el_T = mpcalc.el(p, T, Td)

    # Lapse rates
    def lapse(z_slice, T_slice):
        dz = np.diff(z_slice.magnitude)
        dT = np.diff(T_slice.magnitude)
        lr = -1000 * dT / dz
        return np.nanmean(lr)
    lr_0_3 = lapse(z[z < 3000 * units.meter], T[z < 3000 * units.meter])
    lr_700_500 = lapse(z[(p <= 700 * units.hPa) & (p >= 500 * units.hPa)], T[(p <= 700 * units.hPa) & (p >= 500 * units.hPa)])
    lr_850_500 = lapse(z[(p <= 850 * units.hPa) & (p >= 500 * units.hPa)], T[(p <= 850 * units.hPa) & (p >= 500 * units.hPa)])

    # DCAPE
    try:
        dcape_val = mpcalc.downdraft_cape(p, T, Td).magnitude
    except:
        dcape_val = 0

    # Storm motion and helicity
    try:
        rm, lm, mean = mpcalc.bunkers_storm_motion(p, u, v, z)
        rm_speed = mpcalc.wind_speed(rm[0], rm[1]).magnitude
    except:
        rm = (0 * units.knots, 0 * units.knots)
        rm_speed = 0
    try:
        srh_1 = mpcalc.storm_relative_helicity(z, u, v, depth=1 * units.km, storm_u=rm[0], storm_v=rm[1])[0].magnitude
        srh_3 = mpcalc.storm_relative_helicity(z, u, v, depth=3 * units.km, storm_u=rm[0], storm_v=rm[1])[0].magnitude
        srh_eff = mpcalc.storm_relative_helicity(z, u, v, depth=mpcalc.effective_layer(p, T, Td, z))[0].magnitude
    except:
        srh_1 = srh_3 = srh_eff = 0

    # Shear
    try:
        shear_1 = mpcalc.bulk_shear(p, u, v, height=z, depth=1 * units.km)
        shear_1_mag = mpcalc.wind_speed(*shear_1).magnitude
        shear_3 = mpcalc.bulk_shear(p, u, v, height=z, depth=3 * units.km)
        shear_3_mag = mpcalc.wind_speed(*shear_3).magnitude
        shear_6 = mpcalc.bulk_shear(p, u, v, height=z, depth=6 * units.km)
        shear_6_mag = mpcalc.wind_speed(*shear_6).magnitude
    except:
        shear_1_mag = shear_3_mag = shear_6_mag = 0

    # Additional indices
    wind_speed = mpcalc.wind_speed(u, v).to(units.knots)
    wind_dir = mpcalc.wind_direction(u, v)
    sweat = mpcalc.sweat_index(p, T, Td, wind_speed, wind_dir).magnitude
    k_index = mpcalc.k_index(p, T, Td).magnitude
    tt_index = mpcalc.total_totals_index(p, T, Td).magnitude
    showalter = mpcalc.showalter_index(p, T, Td).magnitude
    lifted_index = mpcalc.lifted_index(p, T, Td).magnitude
    try:
        ehi = mpcalc.energy_helicity_index(mlcape, srh_3).magnitude
    except:
        ehi = 0
    try:
        ship = mpcalc.significant_hail(mlcape, mpcalc.freezing_level(z, T), T[p.argmin()], lr_700_500).magnitude
    except:
        ship = 0
    try:
        stp = mpcalc.significant_tornado(sbcape, sbcin, lcl_z, mpcalc.bulk_shear(p, u, v, height=z, bottom=0*units.meter, top=6000*units.meter)).magnitude
    except:
        stp = 0
    try:
        scp = mpcalc.supercell_composite(mucape, rm_speed, srh_eff).magnitude
    except:
        scp = 0

    return {
        "SBCAPE": sbcape.magnitude, "SBCIN": sbcin.magnitude,
        "MLCAPE": mlcape.magnitude, "MLCIN": mlcin.magnitude,
        "MUCAPE": mucape.magnitude, "MUCIN": mucin.magnitude,
        "LCL": lcl_z.magnitude, "LFC": lfc_p.magnitude if lfc_p else np.nan, "EL": el_p.magnitude if el_p else np.nan,
        "DCAPE": dcape_val,
        "LR_0_3": lr_0_3, "LR_700_500": lr_700_500, "LR_850_500": lr_850_500,
        "SRH_1": srh_1, "SRH_3": srh_3, "SRH_EFF": srh_eff,
        "SHEAR_1": shear_1_mag, "SHEAR_3": shear_3_mag, "SHEAR_6": shear_6_mag,
        "SWEAT": sweat, "K_INDEX": k_index, "TT_INDEX": tt_index,
        "SHOWALTER": showalter, "LIFTED_INDEX": lifted_index,
        "EHI": ehi, "SHIP": ship, "STP": stp, "SCP": scp,
        "RM_SPEED": rm_speed
    }

def storm_mode(p):
    if p["MLCAPE"] < 250: return "Weak / Elevated"
    if p["SHEAR_6"] > 40 and p["SRH_1"] > 150: return "Discrete Supercells"
    if p["SHEAR_6"] > 35 and p["DCAPE"] > 1000: return "QLCS / Derecho"
    if p["SHEAR_6"] > 30: return "Linear / Embedded Supercells"
    return "Pulse / Multicell"

def CRI(p):
    score = 0
    score += min(p["MLCAPE"]/1000, 3)
    score += min(p["SHEAR_6"]/20, 3)
    score += min(p["SRH_1"]/150, 3)
    if p["DCAPE"] > 1000: score += 1
    return round(score, 1)

def plot_skewt(df, station):
    fig = plt.figure(figsize=(12, 12))
    skew = SkewT(fig, rotation=45)
    p = df.pressure.values * units.hPa
    T = df.temperature.values * units.degC
    Td = df.dewpoint.values * units.degC
    u = df.u_wind.values * units.knots
    v = df.v_wind.values * units.knots

    skew.plot(p, T, 'r', label='Temp')
    skew.plot(p, Td, 'g', label='Dewpoint')
    skew.plot_barbs(p, u, v)
    skew.ax.set_ylim(1050, 100)
    skew.ax.set_xlim(-50, 60)
    skew.ax.legend()
    ax_hod = inset_axes(skew.ax, '40%', '40%', loc='upper right')
    h = Hodograph(ax_hod, component_range=80)
    h.add_grid(increment=20)
    wind_speed = mpcalc.wind_speed(u, v)
    h.plot_colormapped(u, v, wind_speed)
    plt.title(f"Skew-T & Hodograph: {station}")
    plt.savefig('skewt_hodograph.png')
    plt.close()

# NWS API functions
def get_nws_data(lat, lon):
    headers = {"User-Agent": "(severeweatherinterface.com, contact@example.com)"}
    points_url = f"https://api.weather.gov/points/{lat},{lon}"
    try:
        resp = requests.get(points_url, headers=headers).json()
        props = resp['properties']
        radar_station = props['radarStation']
        forecast_url = props['forecast']
        return radar_station, forecast_url
    except Exception as e:
        st.error(f"Error fetching NWS points data: {str(e)}")
        return None, None

def get_alerts(lat, lon):
    headers = {"User-Agent": "(severeweatherinterface.com, contact@example.com)"}
    alerts_url = f"https://api.weather.gov/alerts/active?point={lat},{lon}"
    try:
        resp = requests.get(alerts_url, headers=headers).json()
        features = resp.get('features', [])
        severe_alerts = [f['properties'] for f in features if 'Severe' in f['properties'].get('severity', '') or 'Tornado' in f['properties'].get('event', '') or 'Thunderstorm' in f['properties'].get('event', '')]
        return severe_alerts
    except Exception as e:
        st.error(f"Error fetching alerts: {str(e)}")
        return []

def get_forecast(forecast_url):
    headers = {"User-Agent": "(severeweatherinterface.com, contact@example.com)"}
    try:
        resp = requests.get(forecast_url, headers=headers).json()
        periods = resp['properties']['periods']
        return periods
    except Exception as e:
        st.error(f"Error fetching forecast: {str(e)}")
        return []

def get_mesoscale_discussions():
    url = "https://www.spc.noaa.gov/products/md/"
    try:
        resp = requests.get(url).text
        # Simple parse for latest MDs
        mds = [line.strip() for line in resp.splitlines() if "Mesoscale Discussion" in line]
        return mds[:5]  # Top 5 recent
    except:
        return []

# Parameter explanations
param_explanations = {
    "SBCAPE": "Surface-Based CAPE: Energy available for convection from surface parcel.",
    "SBCIN": "Surface-Based CIN: Inhibition to convection from surface.",
    "MLCAPE": "Mixed-Layer CAPE: Average parcel from lowest 100 mb.",
    "MLCIN": "Mixed-Layer CIN.",
    "MUCAPE": "Most Unstable CAPE: Parcel with highest theta-e in lowest 300 mb.",
    "MUCIN": "Most Unstable CIN.",
    "LCL": "Lifted Condensation Level: Height where cloud base forms.",
    "LFC": "Level of Free Convection: Height where parcel becomes buoyant.",
    "EL": "Equilibrium Level: Height where parcel becomes stable again.",
    "DCAPE": "Downdraft CAPE: Potential for strong downdrafts/gusty winds.",
    "LR_0_3": "0-3 km Lapse Rate: Steep rates favor strong updrafts.",
    "LR_700_500": "700-500 mb Lapse Rate: Mid-level cooling for instability.",
    "LR_850_500": "850-500 mb Lapse Rate: Overall mid-low instability.",
    "SRH_1": "0-1 km Storm-Relative Helicity: Rotation potential for low-level mesocyclones/tornadoes.",
    "SRH_3": "0-3 km SRH: For supercell rotation.",
    "SRH_EFF": "Effective SRH: Over effective inflow layer.",
    "SHEAR_1": "0-1 km Bulk Shear: Low-level shear for tornado potential.",
    "SHEAR_3": "0-3 km Bulk Shear.",
    "SHEAR_6": "0-6 km Bulk Shear: Deep-layer shear for organized convection.",
    "SWEAT": "Severe Weather Threat Index: Combines instability, moisture, shear for severe potential.",
    "K_INDEX": "K-Index: Measures thunderstorm potential from moisture and lapse rates.",
    "TT_INDEX": "Total Totals Index: Instability from 850-500 mb temp/dewpoint.",
    "SHOWALTER": "Showalter Index: Stability for elevated convection.",
    "LIFTED_INDEX": "Lifted Index: Stability at 500 mb.",
    "EHI": "Energy Helicity Index: Combines CAPE and SRH for tornado potential.",
    "SHIP": "Significant Hail Parameter: For large hail (>2\").",
    "STP": "Significant Tornado Parameter: For strong tornadoes.",
    "SCP": "Supercell Composite Parameter: For supercell thunderstorms.",
    "RM_SPEED": "Right-Moving Storm Speed: Estimated supercell motion."
}

# -------------------------
# Streamlit app
# -------------------------
st.title("Enhanced Severe Weather Interface")

st.write("Enter location or auto-detect. Select custom sounding station or radar if desired.")

use_auto_location = st.checkbox("Use my current location", value=True)
selected_station = st.selectbox("Select Sounding Station (optional, defaults to nearest)", ["Auto"] + list(stations.keys()))
selected_radar = st.selectbox("Select Radar Site (optional, defaults to nearest)", ["Auto"] + list(radars.keys()))

if use_auto_location:
    lat, lon = get_location()
    if lat is None:
        st.error("Could not detect location. Enter manually.")
        use_auto_location = False

if not use_auto_location:
    lat = st.number_input("Latitude", value=38.0, format="%f")
    lon = st.number_input("Longitude", value=-77.0, format="%f")

if st.button("Run Analysis"):
    if lat is not None and lon is not None:
        st.write(f"Analyzing location: {lat:.2f}, {lon:.2f}")

        # Determine station
        if selected_station == "Auto":
            sorted_stations = sorted(stations.items(), key=lambda x: haversine(lat, lon, x[1][0], x[1][1]))
            station_code = sorted_stations[0][0]
        else:
            station_code = selected_station
        df, station, cycle = fetch_sounding(station_code)
        
        # Determine radar
        radar_station, forecast_url = get_nws_data(lat, lon)
        if selected_radar != "Auto":
            radar_station = selected_radar

        # Tabs
        tab1, tab2, tab3, tab4, tab5, tab6 = st.tabs(["Overview", "Detailed Breakdown", "Sounding", "Radar & Satellite", "Alerts & MDs", "Forecast"])

        with tab1:
            st.subheader("Overview")
            spc = get_spc_risk(lat, lon)
            st.write(f"SPC Day 1 Risk Level: {spc}")
            st.image("https://www.spc.noaa.gov/products/outlook/day1otlk.gif", caption="SPC Day 1 Convective Outlook")

            if df is None:
                st.warning("No valid sounding data available.")
            else:
                p = analyze(df)
                mode = storm_mode(p)
                cri = CRI(p)
                st.write(f"Storm Mode: {mode}")
                st.write(f"CRI Score: {cri}")
                st.write(f"STP: {p['STP']:.1f} | SCP: {p['SCP']:.1f}")
                st.subheader("Key Parameters")
                key_params = {k: v for k, v in p.items() if k in ["MLCAPE", "SRH_1", "SHEAR_6", "DCAPE", "EHI", "SHIP"]}
                params_df = pd.DataFrame.from_dict(key_params, orient='index', columns=['Value']).round(0)
                st.table(params_df)

                # Alert
                elevated = ["SLGT", "ENH", "MDT", "HIGH"]
                if cri >= 6 or spc in elevated:
                    message = f"SPC: {spc}\nCRI: {cri}\nMode: {mode}\nMLCAPE: {p['MLCAPE']:.0f}\nSRH(1km): {p['SRH_1']:.0f}\nShear: {p['SHEAR_6']:.0f} kt"
                    try:
                        plyer.notification.notify(title="⚠️ Severe Weather Alert", message=message, timeout=10)
                    except:
                        st.warning("Notifications not supported.")
                    st.error("⚠️ Severe Weather Possible!")

        with tab2:
            st.subheader("Detailed Environmental Breakdown")
            if 'p' in locals():
                params_df = pd.DataFrame.from_dict(p, orient='index', columns=['Value']).round(1)
                for param, value in params_df.iterrows():
                    with st.expander(f"{param}: {value['Value']} - Explanation"):
                        st.write(param_explanations.get(param, "No explanation available."))
                st.subheader("Interpretation")
                st.write("High CAPE (>2000 J/kg) indicates strong instability. High shear (>40 kt) favors organized storms. High SRH (>150 m²/s²) increases tornado risk. Steep lapse rates (>7 C/km) support hail and strong updrafts.")
            else:
                st.info("Run analysis first.")

        with tab3:
            st.subheader("Sounding Analysis")
            if 'df' in locals() and df is not None:
                plot_skewt(df, station)
                st.image("skewt_hodograph.png", caption="Skew-T & Hodograph")
            else:
                st.info("Run analysis to view.")

        with tab4:
            st.subheader("Radar & Satellite Imagery")
            if radar_station:
                radar_url = f"https://radar.weather.gov/ridge/standard/{radar_station}_loop.gif"
                st.image(radar_url, caption=f"Animated Radar Loop for {radar_station}")
            st.subheader("Satellite")
            st.image("https://cdn.star.nesdis.noaa.gov/GOES16/ABI/CONUS/GEOCOLOR/GOES16-CONUS-GEOCOLOR-625x375.gif", caption="GOES-16 CONUS GeoColor")

        with tab5:
            st.subheader("Active Severe Alerts")
            alerts = get_alerts(lat, lon)
            if alerts:
                for alert in alerts:
                    st.write(f"**Event:** {alert['event']}")
                    st.write(f"**Headline:** {alert['headline']}")
                    st.write(f"**Description:** {alert['description']}")
                    st.write(f"**Effective:** {alert['effective']} to {alert['expires']}")
                    st.write("---")
            else:
                st.info("No active severe alerts.")
            st.subheader("SPC Mesoscale Discussions")
            mds = get_mesoscale_discussions()
            for md in mds:
                st.write(md)

        with tab6:
            st.subheader("Weather Forecast")
            if forecast_url:
                periods = get_forecast(forecast_url)
                if periods:
                    for period in periods:
                        st.write(f"**{period['name']}:** {period['detailedForecast']}")
                else:
                    st.info("No forecast data.")
            else:
                st.info("Unable to fetch forecast.")
