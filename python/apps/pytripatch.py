#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This program is part of pygimli
Visit http://www.resistivity.net for further information or the latest version.
"""

import sys, math

__verbose__ = False

def sign( i ): return abs( i ) / i

# for system call
from os import system, path

try:
    import pygimli as g
except ImportError:
    sys.stderr.write('''ERROR: cannot import the library 'pygimli'. Ensure that pygimli is in your PYTHONPATH ''')
    sys.exit( 1 )

import pygimli.mplviewer
from pygimli.viewer import *

import matplotlib

# use this only if silent does not work
# matplotlib.use('agg')
#

import matplotlib as mpl
from matplotlib import colors, ticker
from matplotlib import pyplot as pylab
import pylab

import numpy as np
from numpy import arange, sin, pi, linspace, meshgrid, array, zeros

def on_draw( event = None):
    def _getBB( items, fig ):
        bboxes = []
        #for label in labels:
        for item in items:
            bbox = item.get_window_extent()
            bboxi = bbox.inverse_transformed(fig.transFigure)
            bboxes.append(bboxi)

        bbox = mtransforms.Bbox.union( bboxes )
    return bbox

    ybbox = _getBB( ax.get_yticklabels() + [ax.yaxis.label], fig )
    fig.subplots_adjust( left = 1.05 * ybbox.width )
    xbbox = _getBB( ax.get_xticklabels() + [ax.xaxis.label], fig )
    fig.subplots_adjust( bottom = 1.05 * xbbox.height )

    print ybbox.width, xbbox.width
    fig.subplots_adjust( right = 0.99 ) # pad a little
    fig.subplots_adjust( top = 0.99 ) # pad a little

#    print bbox
#    print bbox.width
#   if fig.subplotpars.left < bbox.width:
#       # we need to move it over
#    fig.subplots_adjust(left=bbox.width) # pad a little
    #fig.canvas.draw()
    #fig.update(  )
    return False

def applyPublishStyle( style ):
    vals = style.split(':')

    if vals[0] == 'w':
        "not yet done"
    elif vals[0] == 'h':
        "not yet done"
    else:
        print "publish dominant dimension not known", vals[ 0 ]

    paper = vals[1]
    margin = float( vals[2] )
    widthScale = float( vals[3] )
    heightScale = float( vals[4] )
    fontsize = int( vals[5] )
    scale = float( vals[6] )
    wOffset = float( vals[7] )
    hOffset = float( vals[8] )

    pygimli.mplviewer.setOutputStyle( dim = 'w', paperMargin = margin
                                   , xScale = widthScale, yScale = heightScale,
                                   fontsize = fontsize, scale = scale )

    return wOffset, hOffset
#    fig.canvas.mpl_connect('draw_event', on_draw)
#def applyPublishStyle( ... )


#from MatplotPanel import MatplotPanel
#from MiscUtils import *
#import wx

class MyLinearSegmentedColormapAlpha( mpl.colors.LinearSegmentedColormap ):
    def __init__(self, name, segmentdata, N=256):
        mpl.colors.LinearSegmentedColormap.__init__( self, name, segmentdata, N )

    def _init( self ):
        self._lut = np.ones((self.N + 3, 4), np.float)
        self._lut[:-3, 0] = mpl.colors.makeMappingArray(self.N, self._segmentdata['red'])
        self._lut[:-3, 1] = mpl.colors.makeMappingArray(self.N, self._segmentdata['green'])
        self._lut[:-3, 2] = mpl.colors.makeMappingArray(self.N, self._segmentdata['blue'])
        if self._segmentdata.has_key('alpha'):
            self._lut[:-3, 3] = mpl.colors.makeMappingArray( self.N, self._segmentdata['alpha'] )
            #print "found alpha"

        self._isinit = True
        self._set_extremes()

def showTriMesh( meshname, modelname, contour = False, constraintMat = None, cWeight = None, cMin = None, drawEdges = False,
                    cMax = None, coverage = None, showCbar = True, label = "", linear = False , offset = g.RVector3( 0.0, 0.0 ) ):
    mesh = g.Mesh( )
    if ( meshname.rfind( '.vtk' ) != -1 ):
        mesh.importVTK(meshname)
    else:
        mesh.load( meshname )
    
    if __verbose__:
        print mesh
        print "mesh data are:"
    for key, val in mesh.exportDataMap():
        print key, val

    mesh.translate( offset )
    data = g.RVector();

    fig = pylab.figure()
    
    if __verbose__:
        print "create figure"
    
    axis = fig.add_subplot(111)

    # overide default draw function that ignore zorder of images, what we need for coverage image overlay
    #mypa = MPLAxesPatch( axis )
    #axis.draw = mypa.draw

    axis.set_aspect( 'equal' )

    if ( constraintMat ):
        showConstrainMat( axis, mesh, constraintMat, cWeight )
    elif ( modelname ):

        if modelname == 'marker':
            data = pylab.asarray( mesh.cellMarker() )
        elif modelname == 'attribute':
            data = pylab.asarray( mesh.cellAttributes() )
        elif modelname in mesh.exportDataMap().keys():
            data = mesh.exportData( modelname )
        else:
            g.load( data, modelname, g.Ascii )

        print "data min/max:", min(data), max(data)
        cov = None
        
        if coverage:
            cov = g.RVector()
            try:
                if coverage in mesh.exportDataMap().keys():
                    cov = mesh.exportData( coverage )
#                    if coverage.find( 'log10' ) > 0:
#                        cov = g.exp10( cov )
                else:
                    g.load( cov, coverage )
                print "coverage min/max:", min(cov), max(cov)
            except Exception as e:
                print e
                "coverage not found, ignoring"
                cov = None

        if ( contour ):
            showMeshInterpolated( axis, mesh, data, cov = cov, cMin = cMin, cMax = cMax, linear = linear )
        else:
            showMeshPatch( axis, mesh, data, cov = cov, cMin = cMin, cMax = cMax, showCbar = showCbar, label = label, linear = linear )
    else:
        g.mplviewer.drawMeshBoundaries( axis, mesh )
        pass

    #ol = g.Mesh()
    #ol.createEdge( ol.createNode( 148.5,   0.0, 0.0 ),
                   #ol.createNode( 0.0, -100.0, 0.0 ) )
    #overlaysLines = ol.boundaries()
    ##overlaysLines.append( g.Edge(  ) )

    #g.mplviewer.drawSelectedMeshBoundaries( axis, overlaysLines, color = ( 0.0, 0.0, 0.0, 1.0 ),
                                                    #linewidth = 2.0 )

    m = mesh.findBoundaryByMarker(2)
    print "boundary > 2 " , len(m)
    if len(m)>0:
        g.mplviewer.drawSelectedMeshBoundaries( axis, filter( lambda b: b.marker() > 1, mesh.boundaries() )
                                , color = ( 0.0, 0.0, 0.0, 1.0 )
                                , linewidth = 1.0 )
    elif drawEdges:
        g.mplviewer.drawMeshBoundaries( axis, mesh )

    return axis

def showDC2DInvResMod( modfile, contour, cMin = None, cMax = None, label = ""):
    mesh = g.Mesh();
    mesh.importMod( modfile );
    mesh.showInfos();

    fig = pylab.figure()
    axis = fig.add_subplot(111)
    axis.set_aspect( 'equal' )

    data = mesh.exportData( "rho/Ohmm" );
    cov  = mesh.exportData( "coverage" );
    if ( contour ):
        showMeshInterpolated( axis, mesh, data, cov = cov )
    else:
        showMeshPatch( axis, mesh, data, cov = cov, cMin = cMin, cMax = cMax, label = label)

    return axis

def showMeshPatch( axis, mesh, data, cov = None, cMin = None, cMax = None, showCbar = True, 
                   label = "", linear = False, nLevs = 5, orientation = 'horizontal' ):

    patches = pygimli.mplviewer.drawModel( axis, mesh, data, cMin = cMin, cMax = cMax
               , showCbar = showCbar, linear = linear, label = label
               , nLevs = nLevs, orientation = orientation )
    
    #from pygimli.mplviewer import colorbar
    #cmap = colorbar.blueRedCMap
    #patches.set_cmap( cmap )
    
    patches.set_edgecolor( 'face' )
    patches.set_antialiased( False )
    patches.set_linewidth( 0.001 )
    
    alphaPatch = True

    if cov is not None:
        if alphaPatch:
            patches.set_antialiaseds( True )
            patches.set_linewidth( 0.000 )

            # generate individual color values here
            patches.update_scalarmappable()
    
            cols = patches.get_facecolor( )
            
            C = np.asarray( cov )
            print np.min( C ), np.max( C )
            if ( np.min( C ) < 0. ) | ( np.max ( C ) > 1. ) | ( np.max( C ) < 0.5 ): # not already alpha map
                ( nn, hh ) = np.histogram( C, 50 )
                nnn = nn.cumsum( axis = 0 ) / float( len( C ) )
                print "min-max nnn ", min(nnn), max(nnn)
                mi = hh[ min( np.where( nnn > 0.02 )[0] ) ]
                ma = hh[ max( np.where( nnn < 0.4 )[0] ) ]
                #mi = hh[ min( np.where( nnn > 0.2 )[0] ) ]
                #ma = hh[ max( np.where( nnn < 0.7 )[0] ) ]
                C = ( C - mi ) / ( ma - mi )
                C[ np.where( C < 0. ) ] = 0.0
                C[ np.where( C > 1. ) ] = 1.0

            # add alpha value to the color values
            cols[:,3 ] = C

            patches._facecolors = cols
            patches._edgecolor = 'None'

            # delete patch data to avoid automatically rewrite of _facecolors
            patches._A = None
        else:
            addCoverageImageOverlay( axis, mesh, cov )
            
    axis.set_aspect('equal')


def addCoverageImageOverlay( axis, mesh, cov ):
    Nx = 200
    Ny = 100
    tix = linspace( mesh.xmin(), mesh.xmax(), Nx )
    tiy = linspace( mesh.ymin(), mesh.ymax(), Ny )
    (X,Y) = meshgrid( tix, tiy )
    extent = X.min(), X.max(), Y.min(), Y.max()
#           print "interpolate prep t = ", swatch.duration( True )

    c = arange( 0, Nx * Ny ); c[ : ] = 0.0
    c = g.interpolate( mesh, cov, g.ListToRVector( X.flat[:].tolist() )
                                  , g.ListToRVector( Y.flat[:].tolist() )
                                  , g.RVector( len( Y.flat[:] ), 0.0 ) )
    c = asarray( c )

    print "coverage min: ", min(c), "max: ", max( c )

    ( nn, hh ) = np.histogram( c, bins = 50 )
    nnn = nn.cumsum( axis = 0 ) / float( len( c ) )
    mi = hh[ min( np.where( nnn > 0.02 )[0] ) ]
    ma = hh[ max( np.where( nnn < 0.5 )[0] ) ]
    C = array( c ).reshape( Ny, Nx )
    C = ( C - mi ) / ( ma - mi )
    C[ np.where( C < 0 ) ] = 0.0
    C[ np.where( C > 1 ) ] = 1.0

            #(Nhist , xHist) = np.histogram( c, bins = 100 );
    #covMin = xHist[ 0 ];
    #covMax = xHist[ 80 ];
    ##covMin = -3.51
    ##covMax = -3.5

    covColors = np.ones( (Ny, Nx, 4) )

    for i, row in enumerate( C ):
        for j, val in enumerate( row ):
            covColors[ i, j ] = [ 1.0, 1.0, 1.0, 1.0 - C[ i, j ] ]

    #zorder default for patches = 1, axis = 2.5
    cso = axis.imshow( covColors, extent = extent
                                , origin = 'lower'
                                , zorder = 1.5
                                , alpha = 1.0
                     )
    return cso

def showConstrainMat( axes, mesh, constraintMat, cWeight = None  ):
    g.mplviewer.drawMeshBoundaries( axes, mesh )

    C = g.RMatrix();
    g.loadMatrixCol( C, constraintMat );
    cw = None

    if cWeight is not None:
        cw = g.RVector( cWeight );

    drawParameterConstraints( axes, mesh, C, cw)



def showMeshInterpolated( axis, mesh, data, cov = None, cMin = None, cMax = None, linear = False ):
    withCoverageOverlayImage = True
    swatch = g.Stopwatch( True )
    Nx = 200
    Ny = 100
    nLevels = 12

    dropColLimitsPerc = 5
    axis.set_xlim( mesh.xmin(), mesh.xmax() );
    axis.set_ylim( mesh.ymin(), mesh.ymax() );

    tix = linspace( mesh.xmin(), mesh.xmax(), Nx )
    tiy = linspace( mesh.ymin(), mesh.ymax(), Ny )
    (X,Y) = meshgrid( tix, tiy )
    extent = X.min(), X.max(), Y.min(), Y.max()
    print "interpolate prep t = ", swatch.duration( True )

    z = arange( 0, Nx * Ny )

    if ( data.size() > 0 ):
        z = g.interpolate( mesh, data, g.ListToRVector( X.flat[:].tolist() )
                                       , g.ListToRVector( Y.flat[:].tolist() )
                                       , g.RVector( len( Y.flat ), 0.0 )
                        )
        z = asarray ( z )

    print "interpolate t = ", swatch.duration( True )

    Z, cMin, cMax = g.mplviewer.findAndMaskBestClim( z, cMin, cMax, not(linear) );
    Z = ma.masked_where( z <= 0.0, Z )
    Z = Z.reshape( Ny, Nx )
    #print cMin, cMax

    if ( linear ):
	    levs = g.mplviewer.createLinLevs( cMin, cMax, nLevels + 1 )
    else:
	    levs = g.mplviewer.createLogLevs( cMin, cMax, nLevels + 1 )

    #print np.min( Z )
    print levs

    levs[0] = levs[0]* 0.999
    levs[ len(levs)-1] = levs[len(levs)-1]* 1.001

    cmap = matplotlib.cm.get_cmap( name='jet' );

    cmap.set_over( color='#001111', alpha=0.5 )
    cmap.set_under( color=(1.0, 1.0, 0.0 ), alpha=0.5 )

    cs = axis.contourf( X, Y, Z
                                , levs
                                , cmap = cmap
                                #, norm   = mpl.colors.LogNorm()
                        )
    cs.set_clim( cMin, cMax )

    if __verbose__:
        print "plotting t = ", swatch.duration( True )

    if cov:
        addCoverageImageOverlay( axis, mesh, cov )

    g.mplviewer.createColorbar( cs, cMin = cMin, cMax = cMax, nLevs = 5 )
    return cs


def main( argv ):
    global __verbose__
    
    from optparse import OptionParser

    parser = OptionParser( "usage: %prog [options] mesh|mod", version="%prog: " + g.versionStr()  )
    parser.add_option("-v", "--verbose", dest="verbose", action="store_true"
                            , help="be verbose", default=False )
    parser.add_option("-S", "--silent", dest="silent", action="store_true"
                            , help="set viewer silent mode", default=False)
    parser.add_option("-L", "--linear", dest="linear", action="store_true"
                            , help="set linear color scale", default=False)
    parser.add_option("-C", "--contour", dest="contourplot",action="store_true"
                            , help="Create contourplot instead of simple patch", default=False)
    parser.add_option("-E", "--drawEdges", dest="drawEdges",action="store_true"
                            , help="Force drawing all edges independent on marker", default=False)
    parser.add_option("-d", "--data", dest="datafile",
                            help="data file", metavar="File" )
    parser.add_option("-o", "--output", dest="outFileName",
                            help="filename for the resulting picture. suffix define fileformat (.pdf|.png|.svg)", metavar="File" )
    parser.add_option("-c", "--constraintMat", dest="constraintMat",
                            help="show mesh constraints connections", metavar="File", default='' )
    parser.add_option("-e", "--electrodes", dest="electrodes",
                            help="Show electrode positions as black dots. Give datafile.", metavar="File" )
    parser.add_option("", "--cMin", dest="cMin",
                            help="minimum colour", type="float" )
    parser.add_option("", "--cMax", dest="cMax",
                            help="maximum colour", type="float" )
    parser.add_option("", "--coverage", dest="coverage",
                            help="coverage vector", metavar="File" )
    parser.add_option("", "--cWeight", dest="cWeight",
                            help="constraint weight vector", metavar="File" )
    parser.add_option("", "--no-cbar", dest="showCbar", default =True
                            , help="show colorbar", action="store_false" )
    parser.add_option("", "--cbar-only", dest="cbarOnly", default =False
                            , help="show colorbar only", action="store_true" )
    parser.add_option("", "--label", dest="label"
                            , help="label to be displayed on colorbar", metavar="String"
                            , default = "Resistivity [$\Omega $m]" )
    parser.add_option("", "--xLabel", dest="xlabel"
                            , help="label to x-axes", metavar="String"
                            , default = None )
    parser.add_option("", "--yLabel", dest="ylabel"
                            , help="label to y-axes", metavar="String"
                            , default = None )
    parser.add_option("", "--title", dest="title", default ="", help="draw title" )
    parser.add_option("", "--aspect", dest="aspect" , default ="equal"
                            , help="set aspect ratio 'auto' or 'equal' [equal] " )
    parser.add_option("", "--xOffset", dest="xoffset" , default = 0.0, type = "float"
                            , help="set x-coordinate offset for first electrode [0.0] " )
    parser.add_option("", "--zOffset", dest="zoffset" , default = 0.0, type = "float"
                            , help="set z-coordinate offset for first electrode [0.0] " )
    parser.add_option("", "--outSize", dest="outSize", default = None, type = "string"
                            , help="set the x:y pixel for the resulting figure file" )
    parser.add_option("", "--dpi", dest="dpi", default = 600, type = "float"
                            , help="set the dpi for pixel output or preview" )
    parser.add_option("", "--publish", dest="publish", default = None, type = "string"
                            , help="set output style for publishing " +
                             " dominant-dim:paperSize:margin:wScale:hScale:Fontsize:scale" +
                             " e.g., w:a4:5:0.5:0.2:9:2 (size:width of a4-5cm*0.5 plus 9pt Font scal everything by 2)"
                            )

    (options, args) = parser.parse_args()

    if options.verbose:
        __verbose__ = True 
        print "matplotlib-", mpl.__version__
        print options, args

    wOffset = 0.05
    hOffset = 0.05
    if options.publish:
        wOffset,hOffset= applyPublishStyle( options.publish )

    axes = None

    if options.cbarOnly:
        print "cbar only"

        fig = pylab.figure()
        #axes = fig.add_axes([0.023, 0.25, 0.967, 0.1])

        #Horizontal
        axes = fig.add_axes([0.035, 0.6, 0.93, 0.05]); orientation='horizontal'

        # Vertical
        # axes = fig.add_axes([0.30, 0.02, 0.22, 0.96]); orientation='vertical'

        cmin = 1
        if options.cMin:
            cmin = options.cMin
        cmax = 100
        if options.cMax:
            cmax = options.cMax

        norm = None
        if cmin > 0:
            norm = mpl.colors.LogNorm( vmin = cmin, vmax = cmax)
        else:
            norm = mpl.colors.NoNorm( vmin = cmin, vmax = cmax)

        cmap = mpl.cm.get_cmap( 'jet', 256 )
        cmap.set_bad( [1.0, 1.0, 1.0, 0.0 ] )

        #from pygimli.mplviewer import colorbar
        #cmap = colorbar.blueRedCMap

        cbar = mpl.colorbar.ColorbarBase(axes, norm = norm, cmap = cmap
                                        ,orientation = orientation
        #                               ,drawedges='True'
                                       )
        g.mplviewer.setCbarLevels( cbar, cMin = None, cMax = None, nLevs = 5 )

        #cbar.labelpad = -20
        #cbar.ax.yaxis.set_label_position('left')
        cbar.set_label( options.label )

    else:
        if len( args ) == 0:
            parser.print_help()
            print "Please add a mesh or model name."
            sys.exit( 2 )
        else:
            meshname = args[ 0 ];

        if ( options.verbose ):
            print "verbose =", options.verbose
            print "silent =", options.silent
            print "linear =", options.linear
            print "drawEdges =", options.drawEdges
            print "meshname =", meshname
            print "output =", options.outFileName
            print "data =", options.datafile
            print "coverage =", options.coverage
            print "cMin =", options.cMin, type( options.cMin )
            print "cMax =", options.cMax

        axes = None
        try:
            if ( meshname.rfind( '.mod' ) != -1 ):
                axes = showDC2DInvResMod( meshname, options.contourplot, cMin = options.cMin, cMax = options.cMax
                                        , label = options.label)
            elif ( ( meshname.rfind( '.bms' ) != -1 ) | ( meshname.rfind( '.vtk' ) != -1 ) ):
                axes = showTriMesh( meshname, options.datafile, options.contourplot
                        , options.constraintMat, cWeight = options.cWeight
                        , cMin = options.cMin, cMax = options.cMax
                        , coverage = options.coverage
                        , showCbar = options.showCbar
                        , label = options.label
                        , linear = options.linear
                        , drawEdges = options.drawEdges
                        , offset = g.RVector3( options.xoffset, options.zoffset )
                        )
            else:
                print "Cannot determine format for input mesh. Available are *.bms, *.mod"
                exit(2)
        except RuntimeError,err:
            print err
            print "something goes wrong while drawing mesh"
            exit(2)

        if options.electrodes:
            try:
                d = g.DataContainer( options.electrodes )
                elPos = d.sensorPositions()
                pygimli.mplviewer.drawSensors( axes, elPos, diam = None )
            except Exception as e:
                print ( e + "Cannot determine electrode informations from file:" + str( options.electrodes ) )

            axes.figure.canvas.draw()

        if options.title:
            axes.set_title( options.title )

        if options.aspect:
            axes.set_aspect( options.aspect )

        if options.xlabel:
            axes.set_xlabel( options.xlabel )

        if options.ylabel:
            axes.set_ylabel( options.ylabel )

            if 'Depth' in options.ylabel or 'Tiefe' in options.ylabel:
                print "fixing 'Depth' to be positive values"
                ticks = axes.yaxis.get_majorticklocs()
                tickLabels=[]
                for t in ticks:
                    print t
                    tickLabels.append( str( int( abs( t ) ) ) )

                axes.set_yticklabels( tickLabels )
                print tickLabels

        
    # else not cbar only

    if options.outFileName:
        print "writing: ", options.outFileName
        fig = axes.figure

        if options.publish:
            axes.figure.show()
            axes.figure.canvas.draw()

            def _getBB( items, fig ):
                import matplotlib.transforms as mtransforms
                bboxes = []
                #for label in labels:
                for item in items:
                    print item._renderer, axes.figure.canvas.renderer
                    bbox = item.get_window_extent(axes.figure.canvas.renderer)
                    bboxi = bbox.inverse_transformed(fig.transFigure)
                    bboxes.append(bboxi)

                bbox = mtransforms.Bbox.union( bboxes )
                return bbox

            #ybbox = _getBB( axes.get_yticklabels() + [axes.yaxis.label], fig )
            #print wOffset, hOffset
            #fig.subplots_adjust( left = wOffset + ybbox.width )

            #xbbox = _getBB( axes.get_xticklabels() + [axes.xaxis.label], fig )
            #print xbbox
            #print xbbox.height

            fig.subplots_adjust( bottom = hOffset, left = wOffset )

            #print ybbox.width, xbbox.width
            fig.subplots_adjust( right = 0.99 ) # pad a little
            fig.subplots_adjust( top = 0.99 ) # pad a little

            fig.patch.set_alpha( 0.0 )
            axes.patch.set_alpha( 1.0 )
        # if options.publish

	(fileBaseName, fileExtension) = path.splitext( options.outFileName )

	if ( fileExtension == '.svg' ):
            pylab.savefig( options.outFileName )
        elif ( fileExtension == '.pdf' ):
            pylab.savefig( options.outFileName, bbox_inches='tight' )
        elif ( fileExtension == '.png' ):
            pylab.savefig( options.outFileName, dpi=options.dpi, bbox_inches='tight' )
        elif ( fileExtension == '.ps' ):
            pylab.savefig( options.outFileName, dpi=(600) )
        elif ( fileExtension == '.eps' ):
            pylab.savefig( options.outFileName, dpi=(1200), bbox_inches='tight' )
        else:
            assert False, ('format %s unknown. (available( svg, png, pdf) )' % fileExtension )

    if not options.silent:
        pylab.show()

if __name__ == "__main__":
    main( sys.argv[ 1: ] )
