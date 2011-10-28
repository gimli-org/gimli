#!/usr/bin/env python

import pygimli as g
import urllib
import pyproj

import pylab as P

import math

def deg2tile(lon_deg, lat_deg, zoom):
    lat_rad = math.radians(lat_deg)
    n = 2.0 ** zoom
    xtile = int((lon_deg + 180.0) / 360.0 * n)
    ytile = int((1.0 - math.log(math.tan(lat_rad) + (1 / math.cos(lat_rad))) / math.pi) / 2.0 * n)
    return (xtile, ytile)
    
    
def tile2deg(xtile, ytile, zoom):
    ''' This returns the NW-corner of the square. Use the function with xtile+1 and/or ytile+1 to get the other corners.
        With xtile+0.5  ytile+0.5 it will return the center of the tile.
    '''
    n = 2.0 ** zoom
    lon_deg = xtile / n * 360.0 - 180.0
    lat_rad = math.atan(math.sinh(math.pi * (1 - 2 * ytile / n)))
    lat_deg = math.degrees(lat_rad)
    return (lat_deg, lon_deg)
    
def replaceFilename( fname ): return fname.replace('/', '_' )

def getTile( xtile, ytile, zoom, server, verbose = False ):
    imagename = str(zoom) + '/' + str( xtile )+ '/' + str( ytile )
    url = server + imagename + '.png'

    try:
        if verbose: print "Read image from disk"
        image = P.imread( replaceFilename( imagename ) + ".png" )
    except:
        if verbose: print "Get map from url maps", url
        filehandle = urllib.urlopen( url, proxies={})
        data = filehandle.read()
        filehandle.close()

        if verbose: print imagename
        fi = open( replaceFilename( imagename ) + ".png", 'w')
        fi.write( data )
        fi.close()
        image = P.imread( replaceFilename( imagename ) + ".png" )

    return image

zoom = 15

proj = pyproj.Proj( proj = 'utm', zone = 32, ellps = 'WGS84' )

points = [ ( 6.80898269772, 51.5228851788 ),
         ( 6.86827602838, 51.4880055697 ) ]

for p in points:
    x,y = proj( p[0], p[1] )
    P.plot( x, y, '.', markersize = 18 )

axes = P.gca()    

upperleft  = proj( axes.get_xlim( )[0], axes.get_ylim( )[1], inverse = True )
lowerright = proj( axes.get_xlim( )[1], axes.get_ylim( )[0], inverse = True )


#xStart, yStart = deg2tile( 6.80898269772, 51.5228851788, zoom )
#xEnd, yEnd= deg2tile( 6.86827602838, 51.4880055697, zoom )


#server = 'http://tile.openstreetmap.org/'
#nXtiles = xEnd - xStart
#nYtiles = yEnd - yStart


#image = P.zeros( shape = (256 * nYtiles, 256 * nXtiles, 3) )

#print image[0,0,:]
#for i in range( nXtiles ):
    #for j in range( nYtiles ):
        #print i, j
        #im = getTile( xStart + i, yStart + j, zoom, server, verbose = True )
        #image[ ( j * 256 ) : ( ( j + 1 ) * 256 ),
               #( i * 256 ) : ( ( i + 1 ) * 256 ) ] = im
        ##age[ i *256, j * 256, : ]

#print image.shape


#P.imshow( image )

P.show()





#6.80898269772,51.5228851788
#6.86827602838,51.4880055697
