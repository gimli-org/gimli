#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This program is part of pybert
Visit http://www.resistivity.net for further information or the latest version.
"""
import sys

from os import system, path

try:
    import pygimli as g
    import pygimli.mplviewer
    
except ImportError:
    sys.stderr.write('''ERROR: cannot import the library 'pygimli'. Ensure that pygimli is in your PYTHONPATH ''')
    sys.exit( 1 )

import pylab


def main( argv ):
    
    from optparse import OptionParser
    parser = OptionParser( "usage: %prog [options] *.dat|*.shm|*.ohm" )
    
    parser.add_option("-v", "--verbose", dest="verbose", action="store_true"
                            , help="be verbose", default=False )
    parser.add_option("-S", "--silent", dest="silent", action="store_true"
                            , help="set viewer silent mode", default=False)
    parser.add_option("-o", "--output", dest="outFileName",
                            help="filename for the resulting picture. suffix define fileformat (.pdf|.png|.svg)", metavar="File" )
    parser.add_option("-d", "--data", dest="token"
                            , help="data token or vector to be displayed", metavar="String"
                            , default = "rhoa" )
    
    parser.add_option("", "--cMin", dest="cMin",
                            help="minimum colour", type="float" )
    parser.add_option("", "--cMax", dest="cMax",
                            help="maximum colour", type="float" )
                            
    parser.add_option("", "--cLin", dest="cLin", default=False
                            , help="scale colorbar linear", action="store_true" )
                            
    parser.add_option("", "--label", dest="label"
                            , help="label to be displayed", metavar="String"
                            , default = "Apparent resistivity [$\Omega $m]" )
                            
    parser.add_option("", "--title", dest="title", default ="", help="draw title" )

    parser.add_option("", "--sheme", dest="sheme", default ="A_M", help="tmp hack until autodetect draw pseudosection sheme" )
    
    parser.add_option("", "--patchData", dest = "patchData", default = False, action="store_true"
                        , help = "Draw patch grafic instead of simple matrix." )

    (options, args) = parser.parse_args()

    if options.verbose:
        print options, args
        
    g.checkAndFixLocaleDecimal_point( options.verbose )   

    if len( args ) == 0:
        parser.print_help()
        print "Please add a data file"
        sys.exit( 2 )
    else:
        dataFileName = args[ 0 ];
        
    data = g.DataContainer( dataFileName )
    
    try:
        vals = data( options.token )

        #if options.token == 'rhoa':
        #if min( vals ) == max( vals ):

        if options.token == 'err':
            vals *= 100
            options.label = "Error [%]"
    except:
        vals = g.RVector( options.token )
    
    if options.verbose:
        print "Token: ", options.token
        print "min max ", g.min( vals ), g.max( vals )
    
    data.showInfos()
    fig = pylab.figure()
    axes = fig.add_subplot( 111 )
    
    try:
        pseudosheme = getattr( g.mplviewer.Pseudotype, options.sheme)
    except:
        print "no pseudosheme ", options.sheme, " found. Falling back to A_M"
        print dir( g.mplviewer.Pseudotype )
        pseudosheme = g.mplviewer.Pseudotype.A_M

    if g.min( vals ) == g.max( vals ):
        g.mplviewer.drawDataAsMarker( axes, data, shemetype = pseudosheme )
    else:
        if options.patchData :
            gci = g.mplviewer.drawDataAsPatches( axes, data, vals, shemetype = pseudosheme )
        else:
            gci = g.mplviewer.drawDataAsMatrix( axes, data, vals, pseudosheme, logScale = not options.cLin )
        #cmap = pylab.cm.get_cmap( 'jet', 12 )
        #gci.set_cmap( cmap )
   
        g.mplviewer.colorbar.createColorbar( gci, nLevs = 5
                                            , cMin = options.cMin, cMax = options.cMax
                                            , label = options.label )

    if options.title:
        axes.set_title( options.title )

    if options.outFileName:
        if options.verbose:
            print "writing: ", options.outFileName
        fig = axes.figure
        scale = 0.5
        fig.set_size_inches( float( 1024 *scale)/fig.get_dpi()
                            , float( 1453 *scale)/fig.get_dpi() )
        
	(fileBaseName, fileExtension) = path.splitext( options.outFileName )

	if ( fileExtension == '.svg' ):
            pylab.savefig( options.outFileName )
        elif ( fileExtension == '.pdf' ):
            pylab.savefig( options.outFileName, bbox_inches='tight' )
        elif ( fileExtension == '.png' ):
            pylab.savefig( options.outFileName )
        else:
            assert False, ('format %s unknown. (available( svg, png, pdf) )' % fileExtension )
    
    if ( not options.silent ):
        pylab.show()

    #pylab.show()

if __name__ == "__main__":
    main( sys.argv[ 1: ] )