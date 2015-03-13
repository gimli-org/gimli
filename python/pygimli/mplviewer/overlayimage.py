# -*- coding: utf-8 -*-

""" Overlay / Underlay an image or a geo referenced map to mpl.axes."""
import os

import math
import random

import matplotlib.image as mpimg
import numpy as np

import urllib


class OverlayImageMPL(object):

    """ What is this? """

    def __init__(self, imageFileName, axes):
        self.axes = axes
        self.imAxes = None
        self.image = mpimg.open(imageFileName)
        self.figure = self.axes.get_figure()

    def clear(self):
        """ What is this? """
        if self.imAxes in self.figure.axes:
            self.figure.delaxes(self.imAxes)

    def setPosition(self, posX, posY, axes=None):
        """ What is this? """
        if axes is not None:
            self.axes = axes
        self.dx = float(self.image.size[0]) / \
            self.figure.get_dpi() / self.figure.get_size_inches()[0]
        self.dy = float(self.image.size[0]) / \
            self.figure.get_dpi() / self.figure.get_size_inches()[1]

        xRange = self.axes.get_xlim()[1] - self.axes.get_xlim()[0]
        yRange = self.axes.get_ylim()[1] - self.axes.get_ylim()[0]

        x = (posX - self.axes.get_xlim()[0]) / xRange
        y = (posY - self.axes.get_ylim()[0]) / yRange

        x *= (self.axes.get_position().x1 - self.axes.get_position().x0)
        y *= (self.axes.get_position().y1 - self.axes.get_position().y0)

        # print self.imAxes
        # print self.figure.axes
        if self.imAxes not in self.figure.axes:
            if (x + self.axes.get_position().x0) > 10:
                print(("overlay size out of range",
                       (x + self.axes.get_position().x0)))
                print((posX, posY))
                print((xRange, yRange))
                print((x, y))
                print((self.axes.get_position().x0,
                       self.axes.get_position().x1))
                print((self.figure.get_size_inches()))
                print(("add axes",
                       [x + self.axes.get_position().x0 - self.dx / 6.0,
                        y + self.axes.get_position().y0, self.dx, self.dy]))
                # hackish
                return

            self.imAxes = self.figure.add_axes([
                x + self.axes.get_position().x0 - self.dx / 6.0,
                y + self.axes.get_position().y0,
                self.dx, self.dy],
                frameon=False, axisbg='y')
        else:
            self.imAxes.set_position([
                x + self.axes.get_position().x0 - self.dx / 6.0,
                y + self.axes.get_position().y0,
                self.dx, self.dy])

        if (len(self.imAxes.get_xticks()) > 0):
            print("overlay imshow")
            self.imAxes.imshow(self.image, origin='lower')
            self.imAxes.set_xticks([])
            self.imAxes.set_yticks([])


def deg2MapTile(lon_deg, lat_deg, zoom):
    """ What is this? """
    lat_rad = math.radians(lat_deg)
    n = 2.0 ** zoom
    xtile = int((lon_deg + 180.0) / 360.0 * n)
    ytile = int((1.0 - math.log(math.tan(lat_rad) +
                                (1 / math.cos(lat_rad))) / math.pi) / 2.0 * n)
    return (xtile, ytile)


def mapTile2deg(xtile, ytile, zoom):
    """
    Returns the NW-corner of the square.

    Use the function with xtile+1 and/or ytile+1 to get the other corners.
    With xtile+0.5  ytile+0.5 it will return the center of the tile.
    """
    n = 2.0 ** zoom
    lon_deg = xtile / n * 360.0 - 180.0
    lat_rad = math.atan(math.sinh(math.pi * (1 - 2 * ytile / n)))
    lat_deg = math.degrees(lat_rad)
    return (lon_deg, lat_deg)


def cacheFileName(fullname, vendor):
    """ Utility. Create filename and path to cache download data. """
    (dirName, fileName) = os.path.split(fullname)

    path = './' + vendor + '/' + dirName

    try:
        os.makedirs(path)
    except OSError:
        pass

    return path + '/' + fileName


def getMapTile(xtile, ytile, zoom, vendor='OSM', verbose=False):
    """
    Get a map tile from public mapping server.

    Parameters
    ----------
    xtile : int

    ytile : int

    zoom : int

    vendor : str
        . 'OSM' or 'Open Street Map' (tile.openstreetmap.org)
        . 'GM' or 'Google Maps' (mt.google.com)

    verbose : bool [false]
        be verbose
    """
    imagename = str(zoom) + '/' + str(xtile) + '/' + str(ytile)

    if vendor == 'OSM' or vendor == 'Open Street Map':
        # http://[abc].tile.openstreetmap.org
        serverName = 'tile.openstreetmap.org'
        url = 'http://a.' + serverName + '/' + imagename + '.png'
        imFormat = '.png'
    elif vendor == 'GM' or vendor == 'Google Maps':
        # http://mt1.google.com/vt/x=70389&s=&y=43016&z=17
        # http://mt.google.com/vt/x=70389&s=&y=43016&z
        serverName = 'mt.google.com'
        nr = random.randint(0, 3)
        url = 'http://mt' + str(nr) + '.google.com/vt/x=' + \
            str(xtile) + '&y=' + str(ytile) + '&z=' + str(zoom)
        imFormat = '.png'
    elif vendor == 'GMS' or vendor == 'Google Maps Satellite':
        serverName = 'khm.google.com'
        nr = random.randint(0, 3)
        url = 'http://khm' + str(nr) + '.google.com/kh/v=60&x=' + \
            str(xtile) + '&y=' + str(ytile) + '&z=' + str(zoom)
        imFormat = '.jpeg'
#        http://khm0.google.com/kh/v=60&x=2197&y=1346&z=12
    else:
        raise "Vendor: " + vendor + \
            " not supported (currently only OSM (Open Street Map))"

    filename = cacheFileName(imagename, serverName) + imFormat

    if os.path.exists(filename):
        if verbose:
            print(("Read image from disk", filename))
        image = mpimg.imread(filename)
    else:
        if verbose:
            print(("Get map from url maps", url))

        opener1 = urllib.request.build_opener()
        filehandle = opener1.open(url, timeout=15)
        data = filehandle.read()
        opener1.close()

        if verbose:
            print(imagename)

        fi = open(filename, 'w')
        fi.write(data)
        fi.close()
        image = mpimg.imread(filename)

    if imFormat == '.jpeg':
        image = image[::-1, ...] / 256.
    return image
# def getMapTile(...)


def underlayMap(axes, proj, vendor='OSM', zoom=-1, pixelLimit=None,
                verbose=False, fitMap=False):
    """
    Get a map from public mapping server and underlay it on the given axes

    Parameters
    ----------
    axes : matplotlib.axes

    proj : pyproy

        Proj Projection

    vendor : str

        . 'OSM' or 'Open Street Map' (tile.openstreetmap.org)
        . 'GM' or 'Google Maps' (mt.google.com)

    zoom : int [-1]

        Zoom level. If zoom is set to -1, the pixel size of the resulting
        image is lower than pixelLimit.

    pixelLimit : [int,int]

    verbose : bool [false]

        be verbose

    fitMap : bool

        The axes is resized to fit the whole map.
    """
    if pixelLimit is None:
        pixelLimit = [1024, 1024]

    origXLimits = axes.get_xlim()
    origYLimits = axes.get_ylim()

    ul = proj(axes.get_xlim()[0], axes.get_ylim()[1], inverse=True)
    lr = proj(axes.get_xlim()[1], axes.get_ylim()[0], inverse=True)

    if zoom == -1:

        nXtiles = 1e99
        nYtiles = 1e99
        zoom = 19

        while ((nYtiles * 256) > pixelLimit[0] or
               (nXtiles * 256) > pixelLimit[1]):
            zoom = zoom - 1
            startTile = deg2MapTile(ul[0], ul[1], zoom)
            endTile = deg2MapTile(lr[0], lr[1], zoom)

            nXtiles = (endTile[0] - startTile[0]) + 1
            nYtiles = (endTile[1] - startTile[1]) + 1
            if verbose:
                print(("tiles: ", zoom, nYtiles, nXtiles))
            if nXtiles == 1 and nYtiles == 1:
                break

        if verbose:
            print(("zoom set to ", zoom))

    startTile = deg2MapTile(ul[0], ul[1], zoom)
    endTile = deg2MapTile(lr[0], lr[1], zoom)

    nXtiles = (endTile[0] - startTile[0]) + 1
    nYtiles = (endTile[1] - startTile[1]) + 1

    image = np.zeros(shape=(256 * nYtiles, 256 * nXtiles, 3))

    if verbose:
        print(("Mapimage size:", image.shape))

    for i in range(nXtiles):
        for j in range(nYtiles):
            im = getMapTile(startTile[0] + i, startTile[1] + j,
                            zoom, vendor, verbose=verbose)

            image[(j * 256): ((j + 1) * 256),
                  (i * 256): ((i + 1) * 256)] = im

    lonLatStart = mapTile2deg(startTile[0], startTile[1], zoom)
    lonLatEnd = mapTile2deg(endTile[0] + 1, endTile[1] + 1, zoom)

    imUL = proj(lonLatStart[0], lonLatStart[1])
    imLR = proj(lonLatEnd[0], lonLatEnd[1])

    extent = np.asarray([imUL[0], imLR[0], imLR[1], imUL[1]])

    axes.imshow(image, extent=extent)

    if not fitMap:
        axes.set_xlim(origXLimits)
        axes.set_ylim(origYLimits)
