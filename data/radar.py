# data/radar.py

radars = {
    'KABR': (45.45, -98.41),
    'KABX': (35.15, -106.82),
    'KAKQ': (36.98, -77.00),
    'KAMA': (35.22, -101.71),
    'KAMX': (25.61, -80.41),
    # ... all other radar sites
}

def get_radar_url(radar_code):
    """Return the NWS radar loop GIF URL."""
    return f"https://radar.weather.gov/ridge/standard/{radar_code}_loop.gif"
