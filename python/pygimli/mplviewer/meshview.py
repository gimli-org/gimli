# -*- coding: utf-8 -*-
import pygimli as g

from numpy import arange, array, ma
import matplotlib as mpl

from colorbar import *

def drawSelectedMeshBoundaries( axes, boundaries, color = ( 0.0, 0.0, 0.0, 1.0 ), linewidth = 1.0 ):
    """Draw selected mesh boundaries"""
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
    '''
        Draw all mesh boundaries
    '''
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
    #patches.set_linewidth( 0.001 )
    axes.add_collection( patches )

    print "plotting time = ", swatch.duration( True )
    return patches
# def createMeshPatches( ... )

def drawMeshPotential( ax, mesh, u, x=[-10.0, 50.0], z=[-50.0, 0.0]
                    , dx = 1, nLevs = 20, title = None, verbose = False, maskZero = False ):
    '''
        Draw the potential that is associated to a mesh
    '''
    
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




