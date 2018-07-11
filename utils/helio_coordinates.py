import numpy as np

sin = np.sin
cos = np.cos


def hg2hcc(r, lon, lat, B0, L0):
    lon_L0 = lon - L0
    x = r * cos(lat) * sin(lon)
    y = r * (sin(lat) * cos(B0) - cos(lat) * cos(lon_L0) * sin(B0))
    z = r * (sin(lat) * sin(B0) + cos(lat) * cos(lon_L0) * cos(B0))
    return x, y, z


def hcc2hg(x, y, z, B0, L0):
    r = np.sqrt(x ** 2 + y ** 2 + z ** 2)
    lat = np.arcsin((y * cos(B0) + z * sin(B0)) / r)
    lon = L0 + np.arctan2(x, z * cos(B0) - y * sin(B0))
    return r, lon, lat
