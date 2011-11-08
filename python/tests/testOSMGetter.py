#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pygimli as g
from pygimli.mplviewer import underlayMap
from pygimli.importexport import readGPX

import pyproj

import pylab as P

pnts = readGPX( 'rieseltag.gpx' ); zone = 32
#pnts = readGPX( 'gps.gpx' ); zone = 29

print pnts
zoom = -1

proj = pyproj.Proj( proj = 'utm', zone = zone, ellps = 'WGS84' )

#pnts = [ ( 7.34737639312, 51.5834143599 ),
#           ( 7.43579248219, 51.5306718274 ) ]

for p in pnts:
    x,y = proj( p[0], p[1] )
    P.plot( x, y, '.', markersize = 18 )

axes = P.gca()

underlayMap( axes, proj, vendor = 'GMS', zoom = 17, pixelLimit = [1024,1024], verbose = True, fitMap = True )
    
axes.grid()

P.show()
