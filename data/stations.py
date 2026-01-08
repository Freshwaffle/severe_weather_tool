import pandas as pd

def load_surface_stations():
    """
    Load US surface observation stations (METAR/ISD).
    """
    try:
        url = "https://www.ncei.noaa.gov/data/global-hourly/access/isd-history.csv"
        df = pd.read_csv(url)

        df = df.rename(
            columns={
                "USAF": "id",
                "LAT": "lat",
                "LON": "lon",
                "STATION NAME": "name",
            }
        )

        df = df[["id", "lat", "lon", "name"]].dropna()
        return df

    except Exception:
        # fallback static
        return pd.DataFrame([
            {"id": "72295", "lat": 38.889, "lon": -77.035, "name": "Washington DC"},
            {"id": "72469", "lat": 39.872, "lon": -75.243, "name": "Philadelphia"},
        ])


def load_upper_air_stations():
    """
    Load US upper-air (RAOB) stations.
    """
    try:
        url = "https://www.spc.noaa.gov/exper/soundings/raob_sites.csv"
        df = pd.read_csv(url)

        df = df.rename(
            columns={
                "ID": "id",
                "LAT": "lat",
                "LON": "lon",
                "NAME": "name",
            }
        )

        return df[["id", "lat", "lon", "name"]]

    except Exception:
        return pd.DataFrame([
            {"id": "72295", "lat": 38.889, "lon": -77.035, "name": "Washington DC"},
            {"id": "72469", "lat": 39.872, "lon": -75.243, "name": "Philadelphia"},
        ])
