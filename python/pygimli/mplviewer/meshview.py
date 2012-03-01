# -*- coding: utf-8 -*-
import pygimli as g

from numpy import arange, array, ma
import matplotlib as mpl

from colorbar import *

def drawMesh( axes, mesh ):
    ''
    ' Draw a 2d mesh into a given axes. Set the limits of the axes tor the mesh extent. '
    ''
    g.mplviewer.drawMeshBoundaries( axes, mesh )
    axes.set_aspect( 'equal')
    axes.set_xlim( mesh.xmin(), mesh.xmax() )
    axes.set_ylim( mesh.ymin(), mesh.ymax() )
# def drawMesh( ... )

def drawModel( axes, mesh, data = None, cMin = None, cMax = None
               , showCbar = True , linear = False, label = "", cmap=None
               , nLevs = 5, orientation = 'horizontal', xlab=None, ylab=None):
    ''
    ' Draw a 2d mesh and color the cell by the data '
    ''
        
    gci = g.mplviewer.createMeshPatches( axes, mesh, alpha = 1.0 )

    if cmap is not None: eval('mpl.pyplot.'+cmap+'()')

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

    if showCbar:
        
        patches = g.mplviewer.createColorbar( gci, cMin = cMin, cMax = cMax,
                                    nLevs = nLevs, label = label, orientation = orientation 
                                    )
    if xlab is not None: axes.set_xlabel( xlab )
    if ylab is not None: axes.set_ylabel( ylab )
    
    return gci
# def drawModel( ... )

def drawSelectedMeshBoundaries( axes, boundaries, color = ( 0.0, 0.0, 0.0, 1.0 ), linewidth = 1.0 ):
    ''
    ' Draw selected mesh boundaries into a given axes'
    ''
    #print "drawSelectedMeshBoundaries", boundaries

    drawAA = True
    lines  = []

    for bound in boundaries:
        lines.append( zip( [ bound.node( 0 ).x(), bound.node( 1 ).x() ],
                           [ bound.node( 0 ).y(), bound.node( 1 ).y() ]) )

    lineCollection = mpl.collections.LineCollection( lines, antialiaseds = drawAA )

    lineCollection.set_color( color )
    lineCollection.set_linewidth( linewidth )
    axes.add_collection( lineCollection )

    return lineCollection

def drawSelectedMeshBoundariesShadow( axes, boundaries, first='x', second='y', color=( 0.5, 0.5, 0.5, 1.0 ) ):
    ''' 
        
    '''
    polys = []
    print len( boundaries )
    for cell in boundaries:
        polys.append( zip( [  getattr( cell.node( 0 ), first )( )
                            , getattr( cell.node( 1 ), first )( )
                            , getattr( cell.node( 2 ), first )( ) ],
                           [  getattr( cell.node( 0 ), second )( )
                            , getattr( cell.node( 1 ), second )( )
                            , getattr( cell.node( 2 ), second )( ) ] ) )

    collection = mpl.collections.PolyCollection( polys, antialiaseds = True )

    collection.set_color( color )
    collection.set_edgecolor( color )
    collection.set_linewidth( 0.2 )
    axes.add_collection( collection )

def drawMeshBoundaries( axes, mesh, fitView = True):
    ''
    ' Draw all mesh boundaries '
    ''
    if not mesh:
        print "drawMeshBoundaries( axes, mesh ): invalid mesh"
        return

    if mesh.nodeCount() < 2:
        print "drawMeshBoundaries( axes, mesh ): to few nodes"
        return

    if fitView:
        axes.set_xlim( mesh.xmin() - 0.05, mesh.xmax() + 0.05 )
        axes.set_ylim( mesh.ymin() - 0.05, mesh.ymax() + 0.05 )

    drawAA = True;
    swatch = g.Stopwatch( True )
    mesh.createNeighbourInfos()

    drawSelectedMeshBoundaries( axes, mesh.findBoundaryByMarker( 0 )
                                , color = ( 0.0, 0.0, 0.0, 1.0 )
                                , linewidth = 0.3 )
    #return
    drawSelectedMeshBoundaries( axes, mesh.findBoundaryByMarker( g.MARKER_BOUND_HOMOGEN_NEUMANN )
                                , color = ( 0.0, 1.0, 0.0, 1.0 )
                                , linewidth = 1.0 )
    drawSelectedMeshBoundaries( axes, mesh.findBoundaryByMarker( g.MARKER_BOUND_MIXED )
                                , color = ( 1.0, 0.0, 0.0, 1.0 )
                                , linewidth = 1.0 )
    drawSelectedMeshBoundaries( axes, filter( lambda b: b.marker() > 0, mesh.boundaries() )
                                , color = ( 0.0, 0.0, 0.0, 1.0 )
                                , linewidth = 1.0 )
    drawSelectedMeshBoundaries( axes, filter( lambda b: b.marker() < -4, mesh.boundaries() )
                                , color = ( 0.0, 0.0, 0.0, 1.0 )
                                , linewidth = 1.0 )

    #drawSelectedMeshBoundaries( axes, [mesh.boundary( 344 )]
                                #, color = ( 1.0, 0.0, 0.0, 1.0 )
                                #, linewidth = 5.5 )

def createMeshPatches( axes, mesh, **kwarg ):
    ''
    ' Utility function to create 2d mesh patches in a axes'
    ''
    if not mesh:
        print "drawMeshBoundaries( axes, mesh ): invalid mesh"
        return

    if mesh.nodeCount() < 2:
        print "drawMeshBoundaries( axes, mesh ): to few nodes"
        return

    swatch = g.Stopwatch( True )

    axes.set_xlim( mesh.xmin(), mesh.xmax() )
    axes.set_ylim( mesh.ymin(), mesh.ymax() )

    polys = []

    for cell in mesh.cells():
        if ( cell.shape().nodeCount() == 3 ):
            polys.append( zip( [ cell.node( 0 ).x(), cell.node( 1 ).x(), cell.node( 2 ).x() ],
                               [ cell.node( 0 ).y(), cell.node( 1 ).y(), cell.node( 2 ).y() ] ) )
        elif ( cell.shape().nodeCount() == 4 ):
            polys.append( zip( [ cell.node( 0 ).x(), cell.node( 1 ).x(), cell.node( 2 ).x(),
                                    cell.node( 3 ).x() ],
                               [ cell.node( 0 ).y(), cell.node( 1 ).y(), cell.node( 2 ).y(),
                                    cell.node( 3 ).y() ] ) )
        else:
            print "unknown shape to patch: " , cell.shape(), cell.shape().nodeCount()

    patches = mpl.collections.PolyCollection( polys, antialiaseds = False, lod = True, **kwarg)

    #patches.set_edgecolor( None )
    patches.set_edgecolor( 'face' )
    #patches.set_linewidth( 1.001 )
    axes.add_collection( patches )

    print "plotting time = ", swatch.duration( True )
    return patches
# def createMeshPatches( ... )

def drawMeshPotential( ax, mesh, u, x=[-10.0, 50.0], z=[-50.0, 0.0]
                    , dx = 1, nLevs = 20, title = None, verbose = False, maskZero = False ):
    ''
    ' Draw the potential that is associated to a mesh '
    ''
    
    swatch = g.Stopwatch( True )
    if ( verbose ):
        print "start interpolation:", swatch.duration( True )
        
    xg = createLinLevs( x[ 0 ], x[ 1 ], int( ( x[1] - x[0] ) / dx ) )
    yg = createLinLevs( z[ 0 ], z[ 1 ], int( ( z[1] - z[0] ) / dx ) )
    X,Y = np.meshgrid( xg, yg )
    
    uI = g.interpolate( mesh, u
                    , g.ListToRVector( list( X.flat ) )
                    , g.RVector( len( Y.flat ), 0.0 )
                    , g.ListToRVector( list( Y.flat ) ), verbose )
    
    if ( verbose ):
        print "interpolation:", swatch.duration( True )

    zi = np.asarray( uI )
    if maskZero:
        zi = np.ma.masked_where( zi <= 0.0, zi )
    Z = zi.reshape( X.shape )

    maxZ = max( min( zi ), max( zi ) )
    epsZ = min( abs( zi ) )
    
    if min( zi ) < 0:
        potLevs = np.linspace( -maxZ, -epsZ, nLevs/2.)
        print potLevs
        potLevs = np.hstack( ( potLevs, potLevs[::-1] * -1. ) )
    else:
        potLevs = np.linespace( 0, maxZ, nLevs )
        
    print potLevs
    linestyles = ['solid'] * len( potLevs )
    
    gci = ax.contourf( X, Y, Z, potLevs )
    ax.contour( X, Y, Z, potLevs, colors = 'white', linewidths = 0.3, linestyles = linestyles )
    ax.set_aspect('equal')
    
    ax.set_xlim( x )
    ax.set_ylim( z )
    
    ax.set_ylabel( 'Depth [m]')
    ax.set_xlabel( '$x$ [m]')
    
    if title is not None:
        ax.set_title( title )
    
    if ( verbose ):
        print "time:", swatch.duration( True )
        
    print "fixing 'Depth' to be positive values"
    ticks = ax.yaxis.get_majorticklocs()
    tickLabels=[]
    
    for t in ticks:
        tickLabels.append( str( int( abs( t ) ) ) )
        ax.set_yticklabels( tickLabels )
        
    return gci
#def drawMeshPotential( ... )

def drawField( axes, mesh, data = None, filled = False, *args, **kwargs ):
    ''
    ' What is this? '
    ''
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
    ''
    ' Draw nLines streamlines for u circular around pos staring at radius rad '
    ''
    for i in np.linspace( 0, 2. * np.pi, nLines ):
        start = pos + g.RVector3( 1.0, 0.0, 0.0 ) * rad * np.cos( i ) + \
                g.RVector3( 0.0, 1.0, 0.0 ) * rad * np.sin( i )
        x,y = streamline( mesh, u, start, step, maxSteps=50000, koords=[0,1] )
        a.plot( x,y, color = 'black', linewidth = 0.6, linestyle = 'solid' )

        if showStartPos:
            a.plot( [start[0], start[0]], [start[1], start[1]], color = 'blue', linewidth = 2, linestyle = 'solid' )
#def drawStreamCircular( ... )


def drawStreamLinear( a, mesh, u, start, end, nLines = 50, step = 0.01, showStartPos = True, color = 'black' ):
    ''
    '  draw nLines streamlines for u linear from start to end '
    ''
    for i in range( nLines ):
        s = start + (end-start)/float( (nLines-1)) * float(i)

        x,y = streamline( mesh, u, s, step, maxSteps=50000, koords=[0,1] )
        a.plot( x,y, color = color, linewidth = 0.6, linestyle = 'solid' )

        if showStartPos:
            a.plot( [s[0], s[0]], [s[1], s[1]], color = 'blue', linewidth = 2, linestyle = 'solid' )

#def drawStreamLinear( ... )

def drawSensors( axes, sensors, diam = None ):
    ''
    ''
    ''
    eCircles = []
    eSpacing = sensors[ 0 ].distance( sensors[ 1 ] )
    
    if diam is None:
        diam = eSpacing / 8.0

    for e in sensors:
        eCircles.append( mpl.patches.Circle( (e[0], e[2]), diam ) )

    p = mpl.collections.PatchCollection( eCircles, color=(0.0, 0.0, 0.0) )
    axes.add_collection( p )
# def drawSensors( ... )


def createParameterContraintsLines( mesh, cMat, cWeight = None ):
    ''
    ''
    ''
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
    ''
    ''
    ''
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



