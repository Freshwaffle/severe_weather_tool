import pandas as pd

def load_radars():
    """
    Load all NEXRAD radar sites.
    Returns DataFrame with columns: id, lat, lon, name
    """
    try:
        url = "https://www.roc.noaa.gov/WSR88D/DownloadableSites.aspx"  # or use a static CSV if preferred
        # Example: if CSV already exists:
        # df = pd.read_csv(url)
        # Placeholder static list for now:
        df = pd.DataFrame([
            {"id": "KABR", "lat": 45.45, "lon": -98.41, "name": "Aberdeen"},
            {"id": "KAMA", "lat": 35.22, "lon": -101.71, "name": "Amarillo"},
            {"id": "KTLX", "lat": 35.35, "lon": -97.44, "name": "Oklahoma City"},
            {"id": "KFWS", "lat": 32.57, "lon": -97.30, "name": "Fort Worth"},
            {"id": "KBOX", "lat": 41.96, "lon": -71.14, "name": "Boston"},
            {"id": "KBUF", "lat": 42.95, "lon": -78.74, "name": "Buffalo"},
            # Add remaining NEXRAD sites or load dynamically from official CSV
        ])
        return df
    except Exception:
        # fallback minimal
        return pd.DataFrame([
            {"id": "KTLX", "lat": 35.35, "lon": -97.44, "name": "Oklahoma City"},
            {"id": "KFWS", "lat": 32.57, "lon": -97.30, "name": "Fort Worth"},
        ])


def radar_loop_url(radar_id):
    """
    Returns the NWS radar loop GIF URL for a given radar site.
    """
    return f"https://radar.weather.gov/ridge/standard/{radar_id}_loop.gif"
