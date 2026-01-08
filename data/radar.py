import pandas as pd

# NOAA-maintained NEXRAD site table
NEXRAD_TABLE_URL = "https://www.ncdc.noaa.gov/nexradinv/choosesite.jsp"


def load_nexrad_sites() -> pd.DataFrame:
    """
    Load nationwide NEXRAD radar metadata.
    Returns DataFrame with: id, lat, lon, name
    """
    tables = pd.read_html(NEXRAD_TABLE_URL)
    df = tables[0]

    df = df.rename(
        columns={
            "ICAO": "id",
            "Lat": "lat",
            "Lon": "lon",
            "Location": "name",
        }
    )

    df = df[["id", "lat", "lon", "name"]]
    df["lat"] = pd.to_numeric(df["lat"])
    df["lon"] = pd.to_numeric(df["lon"])

    return df


def find_nearest_radar(lat: float, lon: float, radars: pd.DataFrame):
    """
    Find nearest radar to a given lat/lon
    """
    radars = radars.copy()
    radars["dist"] = (radars["lat"] - lat) ** 2 + (radars["lon"] - lon) ** 2
    return radars.sort_values("dist").iloc[0]
