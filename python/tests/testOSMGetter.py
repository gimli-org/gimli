#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pygimli as g
from pygimli.mplviewer import underlayMap
from pygimli.importexport import readGPX

import pyproj

import pylab as P

pnts = readGPX( 'rieseltag.gpx' )

print pnts
zoom = 12

proj = pyproj.Proj( proj = 'utm', zone = 32, ellps = 'WGS84' )

#pnts = [ ( 7.34737639312, 51.5834143599 ),
#           ( 7.43579248219, 51.5306718274 ) ]

for p in pnts:
    x,y = proj( p[0], p[1] )
    P.plot( x, y, '.', markersize = 18 )

axes = P.gca()

underlayMap( axes, proj, vendor = 'OSM', zoom = zoom, verbose = True )
    
axes.grid()

P.show()
