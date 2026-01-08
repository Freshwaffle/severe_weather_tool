# viz/map.py
import folium

def create_base_map(lat=39.8283, lon=-98.5795, zoom_start=4):
    """
    Create a Folium base map of the US
    """
    m = folium.Map(location=[lat, lon], zoom_start=zoom_start)
    return m

def add_spc_polygon(m, geojson):
    """
    Add SPC outlook polygons to the map
    """
    folium.GeoJson(geojson, name="SPC Outlook").add_to(m)
    return m

