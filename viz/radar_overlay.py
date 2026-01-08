# viz/radar_overlay.py
import folium

def add_radar_overlay(m, radar_url, bounds=[[24.396308,-125.0],[49.384358,-66.93457]]):
    """
    Overlay radar GIF on map
    """
    folium.raster_layers.ImageOverlay(
        image=radar_url,
        bounds=bounds,
        opacity=0.6,
        interactive=True,
        cross_origin=False,
        zindex=10
    ).add_to(m)
    return m

