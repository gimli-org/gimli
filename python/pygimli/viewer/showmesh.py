# -*- coding: utf-8 -*-

try:
    import pygimli as g
    import pygimli.mplviewer
except ImportError:
    sys.stderr.write('''ERROR: cannot import the library 'pygimli'. Ensure that pygimli is in your PYTHONPATH ''')
    sys.exit( 1 )

from pygimli.misc import streamline

import pylab
import matplotlib as mpl
import numpy as np

def showMesh( mesh, data = None, *args, **kwargs):
    '''
        syntactic sugar, short-cut to create axes and plot node or cell values
    '''
    fig = pylab.figure()
    a = fig.add_subplot( 111 )
    if data is None:
        drawMesh( a, mesh )
    else:
        if len( data ) == mesh.cellCount():
            drawModel( a, mesh, data, *args, **kwargs )
        elif len( data ) == mesh.nodeCount():
            drawField( a, mesh, data, *args, **kwargs )

    pylab.show()

    #fig.show()
    #fig.canvas.draw()

    return a

def drawMesh( axes, mesh ):

    g.mplviewer.drawMeshBoundaries( axes, mesh )
    axes.set_aspect( 'equal')
    axes.set_xlim( mesh.xmin(), mesh.xmax() )
    axes.set_ylim( mesh.ymin(), mesh.ymax() )

def drawModel( axes, mesh, data = None, cMin = None, cMax = None
               , showCbar = True , linear = False, label = "" ):
    gci = g.mplviewer.createMeshPatches( axes, mesh, alpha = 1.0 )
    axes.set_aspect( 'equal')

    gci.set_antialiased( True )
    gci.set_linewidth( None )

    if data is None:
        data = g.RVector( mesh.cellCount() )

    if len( data ) != mesh.cellCount():
        viewdata = data( mesh.cellMarker() )
    else:
        viewdata = data

    if min( data ) < 0:
        linear = True

    g.mplviewer.setMappableData( gci, viewdata, cMin = cMin, cMax = cMax, logScale = not(linear)  )

    if showCbar and (cMin is not cMax):
        g.mplviewer.createColorbar( gci, cMin = cMin, cMax = cMax, nLevs = 5, label = label )
    return gci


def drawField( axes, mesh, data = None, filled = False, *args, **kwargs ):
    '''
    '''
    print kwargs
    import matplotlib.tri as tri

    x = np.zeros( mesh.nodeCount() )
    y = np.zeros( mesh.nodeCount() )

    for i, p in enumerate( mesh.positions() ):
        x[i] = p[0]
        y[i] = p[1]

    triangles = np.zeros( (mesh.cellCount(), 3 ) )
    for i, c in enumerate( mesh.cells() ):
        triangles[ i, 0 ] = c.node(0).id()
        triangles[ i, 1 ] = c.node(1).id()
        triangles[ i, 2 ] = c.node(2).id()

    #triang = tri.Triangulation(x, y)

    #mask = np.where( data == 0, 1, 0)
    #triang.set_mask(mask)

    #axes.tricontour( triang, z )
    if filled:
        axes.tricontourf(x, y, triangles, data, *args, **kwargs)
    else:
        axes.tricontour(x, y, triangles, data, *args, **kwargs)

# def drawField(...)

def drawStreamCircular( a, mesh, u, pos, rad, nLines = 20, step = 0.1, showStartPos = False ):
    '''
        draw nLines streamlines for u circular around pos staring at radius rad
    '''
    for i in np.linspace( 0, 2. * np.pi, nLines ):
        start = pos + g.RVector3( 1.0, 0.0, 0.0 ) * rad * np.cos( i ) + \
                g.RVector3( 0.0, 1.0, 0.0 ) * rad * np.sin( i )
        x,y = streamline( mesh, u, start, step, maxSteps=50000, koords=[0,1] )
        a.plot( x,y, color = 'black', linewidth = 0.6, linestyle = 'solid' )

        if showStartPos:
            a.plot( [start[0], start[0]], [start[1], start[1]], color = 'blue', linewidth = 2, linestyle = 'solid' )
#def drawStreamCircular( ... )


def drawStreamLinear( a, mesh, u, start, end, nLines = 50, step = 0.01, showStartPos = True, color = 'black' ):
    '''
        draw nLines streamlines for u linear from start to end
    '''
    for i in range( nLines ):
        s = start + (end-start)/float( (nLines-1)) * float(i)

        x,y = streamline( mesh, u, s, step, maxSteps=50000, koords=[0,1] )
        a.plot( x,y, color = color, linewidth = 0.6, linestyle = 'solid' )

        if showStartPos:
            a.plot( [s[0], s[0]], [s[1], s[1]], color = 'blue', linewidth = 2, linestyle = 'solid' )

#def drawStreamLinear( ... )


def drawElectrodes( axes, elecs, diam = None ):
    eCircles = []
    eSpacing = elecs[ 0 ].distance( elecs[ 1 ] )
    if diam is None:
        diam = eSpacing / 10.0

    for e in elecs:
        eCircles.append( mpl.patches.Circle( (e[0], e[1]), diam ) )

    p = mpl.collections.PatchCollection( eCircles, color=(0.0, 0.0, 0.0) )
    axes.add_collection( p )

def createParameterContraintsLines( mesh, cMat, cWeight = None ):
    C = g.RMatrix()
    if type( cMat ) == g.DSparseMapMatrix:
        cMat.save( 'tmpC.matrix' )
        g.loadMatrixCol( C, 'tmpC.matrix');
    else:
        C = cMat

    paraMarker = mesh.cellMarker()
    cellList = dict()
    for cId, marker in enumerate( paraMarker ):
        if not cId in cellList:
            cellList[ cId ] = []
        cellList[ cId ].append( mesh.cell( cId ) )

    paraCenter = dict()
    for id, vals in cellList.items():
        p = g.RVector3( 0.0, 0.0, 0.0);
        for c in vals:
            p += c.center()
        p /= float( len( vals ) )
        paraCenter[ id ] = p

    nConstraints = C[ 0 ].size()
    start      = []
    end        = []
    swatch = g.Stopwatch( True )
    for i in range( 0, nConstraints / 2 ):
        #print i
        #if i == 1000: break;
        idL = int( C[ 1 ][ i * 2 ] )
        idR = int( C[ 1 ][ i * 2 + 1] )
        #leftCells = []
        #rightCells = []
        #for c, index in enumerate( paraMarker ):
            #if idL == index:
                #leftCells.append( mesh.cell( c ) )
            #if idR == index:
                #rightCells.append( mesh.cell( c ) )

        #p1 = g.RVector3(0.0,0.0);
        #for c in leftCells:
            #p1 += c.center()
        #p1 /= float( len( leftCells) )

        #p2 = g.RVector3(0.0,0.0);
        #for c in rightCells:
            #p2 += c.center()
        ##print cWeight[ i ]
        #p2 /= float( len( rightCells) )
        p1 = paraCenter[ idL ]; p2 = paraCenter[ idR ]

        if cWeight is not None:
            pa = g.RVector3( p1 + (p2-p1)/2.0 * ( 1.0 - cWeight[ i ] ) )
            pb = g.RVector3( p2 + (p1-p2)/2.0 * ( 1.0 - cWeight[ i ] ) )
        else:
            pa = p1; pb = p2

        start.append( pa )
        end.append( pb )

    print "createParameterContraintsLines t = ", swatch.duration( True )
    return start, end

def drawParameterConstraints( axes, mesh, cMat, cWeight = None ):

    start, end = createParameterContraintsLines( mesh, cMat, cWeight )

    lines       = []
    colors      = []
    linewidths  = []
    for i, v in enumerate( start ):
        lines.append( zip( [ start[i].x(), end[i].x() ], [ start[i].y(), end[i].y() ] ) )

        linewidth = 0.5
        col = ( 0.0, 0.0, 1.0, 1.0 )
        colors.append( col )
        linewidths.append( linewidth );

    linCol = mpl.collections.LineCollection( lines, antialiaseds = True )

    linCol.set_color( colors )
    linCol.set_linewidth( linewidths )
    axes.add_collection( linCol )

