# -*- coding: utf-8 -*-
"""
    Utility methods to read GPS data and convert them via pyproj.
"""

import sys

import matplotlib.image as mpimg
from math import floor
import numpy as np
import urllib
import os

try:
    from pyproj import Proj, transform
    gk2 = Proj(init="epsg:31466")  # GK zone 2
    gk3 = Proj(init="epsg:31467")  # GK zone 3
    gk4 = Proj(init="epsg:31468")  # GK zone 3
    wgs84 = Proj(init="epsg:4326")  # pure ellipsoid for step-wise change
except ImportError:
    sys.stderr.write("no module pyproj\n")


def handleWPTS(wpts):
    """ Handler for Waypoints in gpx xml-dom """
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

        name = wpt.getElementsByTagName('name')[0].childNodes[0].data
        time = wpt.getElementsByTagName('time')[0].childNodes[0].data

        w.append((lon, lat, name, time))
    return w
# def findWPTS( ... )


def readGPX(filename):
    """
    Extract GPS Waypoint from GPS Exchange Format (GPX).

    Currently only simple waypoint extraction is supported.
    """
    from xml.dom.minidom import parse
    
    dom = parse(filename)
    wpts = dom.getElementsByTagName("wpt")

    return handleWPTS(wpts)
# def readGPX( ... )


def readSimpleLatLon(filename, verbose=False):
    """
        Read a list of the following formats. Try converting automatically
        To be sure, provide format without "d" to ensure floating point format:

        lon lat

        or

        marker lat lon

        return list:
            lon lat name time
    """
    def conv_(deg):
        """convert degree into floating vpoint."""
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


def GK2toUTM(R, H=None, zone=32):
    """ Transform Gauss-Krueger zone 2 into UTM 
        
        Note the double transformation (1-ellipsoid, 2-projection)
        default zone is 32 
    """
    return GKtoUTM(R, H=None, zone=32, gk=gk2)
    #utm = Proj(proj='utm', zone=zone, ellps='WGS84')  # UTM

    #if H is None:  # two-column matrix
        #lon, lat = transform(gk2, wgs84, R[0], R[1])
    #else:
        #lon, lat = transform(gk2, wgs84, R, H)

    #return utm(lon, lat)


def GK3toUTM(R, H=None, zone=32):
    """ Transform Gauss-Krueger zone 3 into UTM 
        
        Note the double transformation (1-ellipsoid, 2-projection)
        default zone is 32 
    """
    return GKtoUTM(R, H=None, zone=32, gk=gk3)
    #utm = Proj(proj='utm', zone=zone, ellps='WGS84')  # UTM

    #if H is None:  # two-column matrix
        #lon, lat = transform(gk3, wgs84, R[0], R[1])
    #else:
        #lon, lat = transform(gk3, wgs84, R, H)

    #return utm(lon, lat)


def GK4toUTM(R, H=None, zone=32):
    """ Transform Gauss-Krueger zone 4 into UTM 
        
        Note the double transformation (1-ellipsoid, 2-projection)
        default zone is 32.
    """
    return GKtoUTM(R, H=None, zone=32, gk=gk4)
    #utm = Proj(proj='utm', zone=zone, ellps='WGS84')  # UTM

    #if H is None:  # two-column matrix
        #lon, lat = transform(gk4, wgs84, R[0], R[1])
    #else:
        #lon, lat = transform(gk4, wgs84, R, H)

    #return utm(lon, lat)


def GKtoUTM(R, H=None, zone=32, gk=None):
    """ Transforms any Gauss-Krueger to UTM autodetect GK zone from offset. """
    if gk is None:
        
        if H is None:
            rr = R[0][0]
        else:
            if isinstance(R, list) or isinstance(R, tuple):
                rr = R[0]
            else:
                rr = R

        gkzone = int(floor(rr * 1e-6))
        print(gkzone)
    
        if gkzone <= 0 or gkzone >= 5:
            print("cannot detect valid GK zone")
    
        gk = Proj(init="epsg:"+str(31464+gkzone))
        
    if H is None:  # two-column matrix
        lon, lat = transform(gk, wgs84, R[0], R[1])
    else:
        lon, lat = transform(gk, wgs84, R, H)
        
    utm = Proj(proj='utm', zone=zone, ellps='WGS84')  # UTM
        
    return utm(lon, lat)


def convddmm(num):
    """ Convert numeric position into degree and minute. """
    dd = np.floor(num / 100.)
    r1 = num - dd * 100.
    return dd + r1 / 60.


def readGeoRefTIF(file_name):
    """ Read geo-referenced TIFF file and return image and bbox.

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


def getBKGaddress(xlim, ylim, imsize=1000, zone=32, service='dop40',
                  usetls=False, uuid='', fmt='image/jpeg'):
    """ Generate address for rendering web service image from BKG.
    
        Assumes UTM in given zone.
    """
    url = 'http://sg.geodatenzentrum.de/wms_' + service
    if usetls:
        url = 'https://sgtls12.geodatenzentrum.de/wms_' + service  # new
    stdarg = '&SERVICE=WMS&VERSION=1.1.0&LAYERS=0&STYLES=default&FORMAT=' + fmt
    srsstr = 'SRS=EPSG:' + str(25800 + zone)  # EPSG definition of UTM

    if imsize is None or imsize <= 1:
        imsize = int((xlim[1] - xlim[0])/0.4) + 1  # take 40cm DOP resolution
        print('choose image size ', imsize)
    box = ','.join(str(int(v)) for v in [xlim[0], ylim[0], xlim[1], ylim[1]])
    ysize = int((imsize - 1.) * (ylim[1] - ylim[0]) / (xlim[1] - xlim[0])) + 1
    sizestr = 'WIDTH=' + str(imsize) + '&HEIGHT=' + '%d' % ysize
    addr = url + '__' + uuid + '?REQUEST=GetMap' + stdarg + '&' + srsstr + \
        '&' + 'BBOX=' + box + '&' + sizestr

    return addr, box


def underlayBKGMap(ax, mode='DOP', utmzone=32, imsize=2500, uuid='',
                   usetls=False):
    """ Underlay digital orthophoto or topographic (mode='DTK') map under axes.
    
        At first access, the image is retrieved from BKG and saved then loaded.
    """
    ext = {'DOP': '.jpg', 'DTK': '.png'}  # extensions for different map types
    wms = {'DOP': 'dop40', 'DTK': 'dtk25'}  # wms service name for map types
    fmt = {'DOP': 'image/jpeg', 'DTK': 'image/png'}  # format
    ad, box = getBKGaddress(ax.get_xlim(), ax.get_ylim(), imsize, zone=utmzone,
                            service=wms[mode], usetls=usetls, uuid=uuid,
                            fmt=fmt[mode])
    imname = mode + box + ext[mode]
    if not os.path.isfile(imname):  # not already existing
        print('Retrieving file from geodatenzentrum.de using URL:')
        print(ad)
        url = urllib.URLopener()
        url.retrieve(ad, imname)

    im = mpimg.imread(imname)
    bb = [int(bi) for bi in box.split(',')]  # bounding box
    ax.imshow(im, extent=[bb[0], bb[2], bb[1], bb[3]],
              interpolation='nearest')
