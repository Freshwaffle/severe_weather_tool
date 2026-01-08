# data/radar.py
radars = {
    'KABR': (45.45, -98.41), 'KABX': (35.15, -106.82), 'KAKQ': (36.98, -77.00),
    'KAMA': (35.22, -101.71), 'KAMX': (25.61, -80.41), 'KAPX': (44.91, -84.72),
    'KARX': (43.82, -91.19), 'KATX': (47.68, -122.50)
    # ... add all others or dynamically load if desired
}

def radar_loop_url(station_code):
    """
    Returns the standard NWS radar loop URL for a given station code.
    """
    return f"https://radar.weather.gov/ridge/standard/{station_code}_loop.gif"
