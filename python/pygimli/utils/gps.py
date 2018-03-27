# -*- coding: utf-8 -*-
"""Utility methods to read GPS data and convert them via pyproj."""

import sys
import os

from math import floor

import numpy as np

from pygimli.io import opt_import


def handleWPTS(wpts):
    """Handler for Waypoints in gpx xml-dom."""
    w = []

    for wpt in wpts:
        if wpt.hasAttribute('lat'):
            lat = float(wpt.getAttribute('lat'))
        else:
            continue
        if wpt.hasAttribute('lon'):
            lon = float(wpt.getAttribute('lon'))
        else:
            continue

        name = 'None'
        if wpt.getElementsByTagName('name'):
            name = wpt.getElementsByTagName('name')[0].childNodes[0].data

        time = 0
        if wpt.getElementsByTagName('time'):
            time = wpt.getElementsByTagName('time')[0].childNodes[0].data

        w.append((lon, lat, name, time))
    return w


def readGPX(fileName):
    """Extract GPS Waypoint from GPS Exchange Format (GPX).

    Currently only simple waypoint extraction is supported.


    <gpx version="1.0" creator="creator">

    <metadata>
    <name>Name</name>
    </metadata>
    <wpt lat="50.0" lon="10.">
    <name>WPT1 Name</name>
    </wpt>

    <wpt lat="51." lon="11.">
    <name>WPT2 Name</name>
    </wpt>

    </gpx>

    """
    from xml.dom.minidom import parse

    dom = parse(fileName)

    name = fileName
    meta = dom.getElementsByTagName("metadata")
    if len(meta) > 0:
        if meta[0].getElementsByTagName('name'):
            name = meta[0].getElementsByTagName('name')[0].childNodes[0].data

    wpts = dom.getElementsByTagName("wpt")

    return handleWPTS(wpts), name


def readSimpleLatLon(filename, verbose=False):
    """Read Lat Lon coordinates from file.

    Try converting automatic.
    To be sure, provide format without "d" to ensure floating point format:

    lon lat

    or

    marker lat lon

    Returns
    -------
    list: []
        lon lat name time
    """
    def conv_(deg):
        """Convert degree into floating vpoint."""
        ret = 0.0
        if 'd' in deg:
            # 10d???
            d = deg.split('d')
            ret = float(d[0])

            if "'" in d[1]:
                # 10d23'2323.44''
                m = d[1].split("'")
                if len(m) > 1:
                    ret += float(m[0]) / 60.
                    ret += float(m[1]) / 3600.
            else:
                # 10d23.2323
                ret += float(d[1]) / 60.
        else:
            # 10.4432323
            ret = float(deg)

        return ret
    # def conv_(...):

    w = []

    with open(filename, 'r') as fi:
        content = fi.readlines()

    for line in content:
        if line[0] == '#':
            continue

        vals = line.split()

        if len(vals) == 2:  # lon lat
            w.append((conv_(vals[1]), conv_(vals[0]), '', 'time'))
        if len(vals) == 3:
            # marker lat lon
            w.append((conv_(vals[2]), conv_(vals[1]), vals[0], 'time'))

        if verbose:
            print(w[-1])

    return w


def GK2toUTM(ea, no=None, zone=32):
    """Transform Gauss-Krueger zone 2 into UTM (for backward compatibility)."""
    return GKtoUTM(ea, no, zone, gkzone=2)


def GK3toUTM(ea, no=None, zone=32):
    """Transform Gauss-Krueger zone 3 into UTM (for backward compatibility)."""
    return GKtoUTM(ea, no, zone, gkzone=3)


def GK4toUTM(ea, no=None, zone=32):
    """Transform Gauss-Krueger zone 4 into UTM (for backward compatibility)."""
    return GKtoUTM(ea, no, zone, gkzone=4)


def GKtoUTM(ea, no=None, zone=32, gk=None, gkzone=None):
    """Transform any Gauss-Krueger to UTM autodetect GK zone from offset."""
    if gk is None and gkzone is None:
        if no is None:
            rr = ea[0][0]
        else:
            if isinstance(ea, list) or isinstance(ea, tuple):
                rr = ea[0]
            else:
                rr = ea

        gkzone = int(floor(rr * 1e-6))
        print(gkzone)

        if gkzone <= 0 or gkzone >= 5:
            print("cannot detect valid GK zone")

    pyproj = opt_import('pyproj', 'coordinate transformations')
    if pyproj is None:
        return None

    gk = pyproj.Proj(init="epsg:"+str(31464+gkzone))
    wgs84 = pyproj.Proj(init="epsg:4326")  # pure ellipsoid to doubel transform
    utm = pyproj.Proj(proj='utm', zone=zone, ellps='WGS84')  # UTM
    if no is None:  # two-column matrix
        lon, lat = pyproj.transform(gk, wgs84, ea[0], ea[1])
    else:
        lon, lat = pyproj.transform(gk, wgs84, ea, no)

    return utm(lon, lat)


def convddmm(num):
    """Convert numeric position into degree and minute."""
    dd = np.floor(num / 100.)
    r1 = num - dd * 100.
    return dd + r1 / 60.


def readGeoRefTIF(file_name):
    """Read geo-referenced TIFF file and return image and bbox.

    plt.imshow(im, ext = bbox.ravel()), bbox might need transform.
    """
    try:
        import gdal
        from osgeo.gdalconst import GA_ReadOnly
    except ImportError:
        sys.stderr.write("no module osgeo\n")

    dataset = gdal.Open(file_name, GA_ReadOnly)
    geotr = dataset.GetGeoTransform()
    projection = dataset.GetProjection()

    im = np.flipud(mpimg.imread(file_name))

    tifx, tify, dx = geotr[0], geotr[3], geotr[1]

    bbox = [[tifx, tifx + im.shape[1] * dx],
            [tify - im.shape[0] * dx, tify]]

    return im, bbox, projection


