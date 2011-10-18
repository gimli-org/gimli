##try:
##    import Image
##except ImportError, exc:
##    raise SystemExit("PIL must be installed to run this example")

class OverlayImageMPL( ):
   
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
            

def underlayGoogleMaps( axes, proj, proxies = {} ):
    '''
        Get map image from googlemaps and underlay the current axes
        The image is proxied on disk
        axes must contain coordinates based on a valid projection
    '''   
    from PIL import Image
    import urllib
    import numpy as np
    
    upperleft = proj( axes.get_xlim( )[0], axes.get_ylim( )[1], inverse = True )
    lowerright= proj( axes.get_xlim( )[1], axes.get_ylim( )[0], inverse = True )

    imagename = 'size=640x640&sensor=false'
    imagename += "&markers="
    imagename += str(upperleft[1])+ "," + str(upperleft[0])
    imagename += "|"
    imagename += str(lowerright[1])+ "," + str(lowerright[0])
            
    image = None
    try:
        print "Read image from disk"
        image = Image.open( imagename + ".png" )
    except:
        print "Get map from google maps"
        googlemap="http://maps.google.com/maps/api/staticmap?"+imagename
        print googlemap
        filehandle = urllib.urlopen(googlemap, proxies=proxies)
        data = filehandle.read()
        filehandle.close()
        imagename1 = "bla"
        fi = open( imagename1 + ".png", 'w')
        fi.write( data )
        image = Image.open( imagename1 + ".png" )
        fi.close()

    extent = np.asarray( [axes.get_xlim( )[0], axes.get_xlim( )[1], axes.get_ylim( )[0], axes.get_ylim( )[1] ] )
    x= extent[1]-extent[0]
    y= extent[3]-extent[2]
    scale = 0.1
    extent += np.asarray( [ -x*scale, x*scale, -y*scale, y*scale ] )
    print extent
    axes.imshow( image, origin='lower'
                , extent = extent )
    
# def underlayGoogleMaps( axes, proj )            

