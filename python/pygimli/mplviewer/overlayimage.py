# -*- coding: utf-8 -*-
"""Overlay / Underlay an image or a geo referenced map to mpl.ax."""

import os

import math
import random

import numpy as np
import matplotlib.image as mpimg


class OverlayImageMPL(object):
    """TODO Documentme."""

    def __init__(self, imageFileName, ax):
        """Constructor."""
        self.ax = ax
        self.imAxes = None
        self.image = mpimg.open(imageFileName)
        self.figure = self.ax.get_figure()
        self.dx = 0
        self.dy = 0

    def clear(self):
        """TODO Documentme."""
        if self.imAxes in self.figure.ax:
            self.figure.delax(self.imAxes)

    def setPosition(self, posX, posY, ax=None):
        """TODO Documentme."""
        if ax is not None:
            self.ax = ax
        self.dx = float(self.image.size[0]) / \
            self.figure.get_dpi() / self.figure.get_size_inches()[0]
        self.dy = float(self.image.size[0]) / \
            self.figure.get_dpi() / self.figure.get_size_inches()[1]

        xRange = self.ax.get_xlim()[1] - self.ax.get_xlim()[0]
        yRange = self.ax.get_ylim()[1] - self.ax.get_ylim()[0]

        x = (posX - self.ax.get_xlim()[0]) / xRange
        y = (posY - self.ax.get_ylim()[0]) / yRange

        x *= (self.ax.get_position().x1 - self.ax.get_position().x0)
        y *= (self.ax.get_position().y1 - self.ax.get_position().y0)

        # print self.imAxes
        # print self.figure.ax
        if self.imAxes not in self.figure.ax:
            if (x + self.ax.get_position().x0) > 10:
                print(("overlay size out of range",
                       (x + self.ax.get_position().x0)))
                print((posX, posY))
                print((xRange, yRange))
                print((x, y))
                print((self.ax.get_position().x0,
                       self.ax.get_position().x1))
                print((self.figure.get_size_inches()))
                print(("add ax",
                       [x + self.ax.get_position().x0 - self.dx / 6.0,
                        y + self.ax.get_position().y0, self.dx, self.dy]))
                # hackish
                return

            self.imAxes = self.figure.add_ax([
                x + self.ax.get_position().x0 - self.dx / 6.0,
                y + self.ax.get_position().y0,
                self.dx, self.dy], frameon=False, axisbg='y')
        else:
            self.imAxes.set_position([
                x + self.ax.get_position().x0 - self.dx / 6.0,
                y + self.ax.get_position().y0,
                self.dx, self.dy])

        if len(self.imAxes.get_xticks()) > 0:
            print("overlay imshow")
            self.imAxes.imshow(self.image, origin='lower')
            self.imAxes.set_xticks([])
            self.imAxes.set_yticks([])


def deg2MapTile(lon_deg, lat_deg, zoom):
    """TODO Documentme."""
    lat_rad = math.radians(lat_deg)
    n = 2.0 ** zoom
    xtile = int((lon_deg + 180.0) / 360.0 * n)
    ytile = int((1.0 - math.log(math.tan(lat_rad) +
                                (1 / math.cos(lat_rad))) / math.pi) / 2.0 * n)
    return (xtile, ytile)


def mapTile2deg(xtile, ytile, zoom):
    """Calculate the NW-corner of the square.

    Use the function with xtile+1 and/or ytile+1 to get the other corners.
    With xtile+0.5  ytile+0.5 it will return the center of the tile.
    """
    n = 2.0 ** zoom
    lon_deg = xtile / n * 360.0 - 180.0
    lat_rad = math.atan(math.sinh(math.pi * (1 - 2 * ytile / n)))
    lat_deg = math.degrees(lat_rad)
    return (lon_deg, lat_deg)


def cacheFileName(fullname, vendor):
    """Createfilename and path to cache download data."""
    (dirName, fileName) = os.path.split(fullname)

    path = './' + vendor + '/' + dirName

    try:
        os.makedirs(path)
    except OSError:
        pass

    return path + '/' + fileName


def getMapTile(xtile, ytile, zoom, vendor='OSM', verbose=False):
    """Get a map tile from public mapping server.

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
        image = image[:,:,0:3]
    else:
        if verbose:
            print(("Get map from url maps", url))

        image = mpimg.imread(url)
        if verbose:
            print(imagename)
        mpimg.imsave(filename, image)

    if imFormat == '.jpeg':
        image = image[::-1, ...] / 256.
    return image
# def getMapTile(...)


def underlayMap(ax, proj, vendor='OSM', zoom=-1, pixelLimit=None,
                verbose=False, fitMap=False):
    """Get a map from public mapping server and underlay it on the given ax.

    Parameters
    ----------
    ax : matplotlib.ax

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

        The ax is resized to fit the whole map.
    """
    if pixelLimit is None:
        pixelLimit = [1024, 1024]

    origXLimits = ax.get_xlim()
    origYLimits = ax.get_ylim()

    ul = proj(ax.get_xlim()[0], ax.get_ylim()[1], inverse=True)
    lr = proj(ax.get_xlim()[1], ax.get_ylim()[0], inverse=True)

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

    ax.imshow(image, extent=extent)

    if not fitMap:
        ax.set_xlim(origXLimits)
        ax.set_ylim(origYLimits)
