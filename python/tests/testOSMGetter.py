#!/usr/bin/env python

import pygimli as g
import urllib

import pylab as P

import math

def deg2tile(lon_deg, lat_deg, zoom):
    lat_rad = math.radians(lat_deg)
    n = 2.0 ** zoom
    xtile = int((lon_deg + 180.0) / 360.0 * n)
    ytile = int((1.0 - math.log(math.tan(lat_rad) + (1 / math.cos(lat_rad))) / math.pi) / 2.0 * n)
    return (xtile, ytile)
    
    
def tile2deg(xtile, ytile, zoom):
    n = 2.0 ** zoom
    lon_deg = xtile / n * 360.0 - 180.0
    lat_rad = math.atan(math.sinh(math.pi * (1 - 2 * ytile / n)))
    lat_deg = math.degrees(lat_rad)
    return (lat_deg, lon_deg)
    
zoom = 12
    
xt,yt = deg2tile( 6.80898269772, 51.5228851788, zoom )
#print deg2tile( 6.86827602838, 51.4880055697, zoom )

#http://tile.openstreetmap.org/17/70389/43016.png

imagename = str(zoom) + '/' + str( yt )+ '/' + str( xt )
url='http://tile.openstreetmap.org/' + imagename + '.png'

def replaceFilename( fname ): return fname.replace('/', '_' )

try:
    print "Read image from disk"
    image = P.imread( replaceFilename( imagename ) + ".png" )
except:
    print "Get map from url maps", url
    filehandle = urllib.urlopen( url, proxies={})
    data = filehandle.read()
    filehandle.close()
    
    print imagename
    fi = open( replaceFilename( imagename ) + ".png", 'w')
    fi.write( data )
    fi.close()
    image = P.imread( replaceFilename( imagename ) + ".png" )


P.imshow( image )

P.show()





#6.80898269772,51.5228851788
#6.86827602838,51.4880055697
