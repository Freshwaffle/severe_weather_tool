# analysis/forecasting.py
import requests
from shapely.geometry import shape, Point

SPC_HIERARCHY = ["TSTM", "MRGL", "SLGT", "ENH", "MDT", "HIGH"]

def get_spc_risk(lat, lon):
    url = "https://mapservices.weather.noaa.gov/vector/rest/services/outlooks/SPC_wx_outlks/MapServer/1/query"
    params = {"where": "1=1", "outFields": "*", "f": "geojson", "outSR": "4326"}
    try:
        data = requests.get(url, params=params, timeout=10).json()
    except:
        return "Error fetching SPC outlook"
    pt = Point(lon, lat)
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
    found.sort(key=lambda x: SPC_HIERARCHY.index(x))
    return found[-1]

def get_forecast(forecast_url):
    if not forecast_url:
        return []
    headers = {"User-Agent": "(severeweatherinterface.com, contact@example.com)"}
    try:
        resp = requests.get(forecast_url, headers=headers, timeout=10).json()
        return resp['properties']['periods']
    except:
        return []
