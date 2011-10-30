##try:
##    import Image
##except ImportError, exc:
##    raise SystemExit("PIL must be installed to run this example")

import urllib
import math
import pylab as P
import numpy as np

class OverlayImageMPL( ):
        ''
        ' What is this? '
        ''
   
        def __init__( self, imageFileName, axes ):
            self.axes = axes
            self.imAxes  = None
            self.image = Image.open( imageFileName )
            self.figure = self.axes.get_figure()
            
        def clear( self ):
            if self.imAxes in self.figure.axes:
                self.figure.delaxes( self.imAxes )
        
        def setPosition( self, posX, posY, axes = None ):
            if axes is not None:
                self.axes = axes
            self.dx = float( self.image.size[0])/self.figure.get_dpi()/self.figure.get_size_inches()[0]
            self.dy = float( self.image.size[0])/self.figure.get_dpi()/self.figure.get_size_inches()[1]
            
            xRange = self.axes.get_xlim( )[1]-self.axes.get_xlim( )[0]
            yRange = self.axes.get_ylim( )[1]-self.axes.get_ylim( )[0]
            
            x = ( posX - self.axes.get_xlim( )[0] ) / xRange
            y = ( posY - self.axes.get_ylim( )[0] ) / yRange
        
            x *= ( self.axes.get_position().x1 - self.axes.get_position().x0  )
            y *= ( self.axes.get_position().y1 - self.axes.get_position().y0  )
            
            #print self.imAxes
            #print self.figure.axes
            if self.imAxes not in self.figure.axes:
                if ( x + self.axes.get_position().x0 ) > 10:
                    print "overlay size out of range", (x + self.axes.get_position().x0 )
                    print posX, posY
                    print xRange, yRange
                    print x, y
                    print self.axes.get_position().x0, self.axes.get_position().x1
                    print self.figure.get_size_inches()
                    print "add axes", [ x + self.axes.get_position().x0 - self.dx/6.0
                                                    , y + self.axes.get_position().y0, self.dx, self.dy
                                                    ]
                    #hackish
                    return 
                
                
                self.imAxes = self.figure.add_axes([ x + self.axes.get_position().x0 - self.dx/6.0
                                                    , y + self.axes.get_position().y0, self.dx, self.dy
                                                    ]
                                                , frameon=False, axisbg='y')
            else:
                self.imAxes.set_position( [ x + self.axes.get_position().x0 - self.dx/6.0
                                            , y + self.axes.get_position().y0, self.dx, self.dy ] )
                                            
            if ( len( self.imAxes.get_xticks() ) > 0 ):
                print "overlay imshow"
                self.imAxes.imshow( self.image, origin='lower' )
                self.imAxes.set_xticks( [] )
                self.imAxes.set_yticks( [] )


def deg2MapTile(lon_deg, lat_deg, zoom):
    ''
    ''
    ''
    lat_rad = math.radians(lat_deg)
    n = 2.0 ** zoom
    xtile = int((lon_deg + 180.0) / 360.0 * n)
    ytile = int((1.0 - math.log(math.tan(lat_rad) + (1 / math.cos(lat_rad))) / math.pi) / 2.0 * n)
    return (xtile, ytile)


def mapTile2deg(xtile, ytile, zoom):
    ''
    ' This returns the NW-corner of the square. Use the function with xtile+1 and/or ytile+1 to get the other corners.'
    '    With xtile+0.5  ytile+0.5 it will return the center of the tile.'
    ''
    n = 2.0 ** zoom
    lon_deg = xtile / n * 360.0 - 180.0
    lat_rad = math.atan(math.sinh(math.pi * (1 - 2 * ytile / n)))
    lat_deg = math.degrees(lat_rad)
    return (lon_deg, lat_deg)
    

def getMapTile( xtile, ytile, zoom, vendor = 'OSM', verbose = False ):
    ''
    ' vendor = OSM for Open Street Map'
    ''

    def filenameProxi( fname ): return vendor + '-' + fname.replace('/', '_' )

    if vendor != 'OSM': raise "Vendor: "+ vendor + " not supported (currently only OSM (Open Street Map) )"
    
    server = 'http://tile.openstreetmap.org/'

    imagename = str(zoom) + '/' + str( xtile )+ '/' + str( ytile )
    url = server + imagename + '.png'

    try:
        if verbose: print "Read image from disk"
        
        image = P.imread( filenameProxi( imagename ) + ".png" )
    except:
        if verbose: print "Get map from url maps", url
        
        filehandle = urllib.urlopen( url, proxies = {} )
        data = filehandle.read()
        filehandle.close()

        if verbose: print imagename
        
        fi = open( filenameProxi( imagename ) + ".png", 'w')
        fi.write( data )
        fi.close()
        image = P.imread( filenameProxi( imagename ) + ".png" )

    return image
            
def underlayMap( axes, proj, vendor = 'OSM-', zoom = -1, verbose = False ):
    ''
    ' vendor = OSM for Open Street Map'
    ''
    if zoom == -1:
        raise "automatically zoom not yet supported"
    
    ul = proj( axes.get_xlim( )[0], axes.get_ylim( )[1], inverse = True )
    lr = proj( axes.get_xlim( )[1], axes.get_ylim( )[0], inverse = True )

    startTile = deg2MapTile( ul[ 0 ], ul[ 1 ], zoom )
    endTile   = deg2MapTile( lr[ 0 ], lr[ 1 ], zoom )

    nXtiles = ( endTile[0] - startTile[0] ) + 1
    nYtiles = ( endTile[1] - startTile[1] ) + 1

    image = np.zeros( shape = (256 * nYtiles, 256 * nXtiles, 3) )

    print "Mapimage size:", image.shape


    for i in range( nXtiles ):
        for j in range( nYtiles ):
            im = getMapTile( startTile[0] + i, startTile[1] + j, zoom, vendor, verbose = verbose )

            image[ ( j * 256 ) : ( ( j + 1 ) * 256 ),
                   ( i * 256 ) : ( ( i + 1 ) * 256 ) ] = im

    lonLatStart = mapTile2deg( startTile[0], startTile[1], zoom)
    lonLatEnd = mapTile2deg( endTile[0] + 1, endTile[1] + 1, zoom)

    imUL = proj( lonLatStart[0], lonLatStart[1] )
    imLR = proj( lonLatEnd[0], lonLatEnd[1] )
    
    extent = np.asarray( [imUL[0], imLR[0], imLR[1], imUL[1] ])

    axes.imshow( image, extent = extent )
#def underlayMap(  )

