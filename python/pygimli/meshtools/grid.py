# -*- coding: utf-8 -*-

import pygimli as g
from os import system

from pygimli.polytools import *

def appendTriangleBoundary( mesh, xbound = 10, ybound = 10, marker = 1
                            , quality = 34.0, isSubSurface = False, verbose = False ):
    ''' adds a box around an existing mesh (for forward calculation)
        outputMesh = appendTriangleBoundary( inputMesh <parameters> )        
        xbound/ybound defines the x/y distance to the box
        marker defines the marker of the new elements
        quality is the triangle quality
        if isSubsurface also distance positive y is added 
    '''         
           
    ''' function comparing x for using sort '''
    def sortNodeY( n1, n2 ):
        return cmp( n1.pos().y(), n2.pos().y() )
    ''' function comparing y for using sort '''
    def sortNodeX( n1, n2 ):
        return cmp( n1.pos().x(), n2.pos().x() )

    surface = 0.0
    ''' find boundaries on left/right/bottom/top side  '''
    le = filter( lambda b: b.center().x() == mesh.xmin(), mesh.boundaries() )
    bo = filter( lambda b: b.center().y() == mesh.ymin(), mesh.boundaries() )
    ri = filter( lambda b: b.center().x() == mesh.xmax(), mesh.boundaries() )
    top = filter( lambda b: b.center().y() == mesh.ymax(), mesh.boundaries() )

    ''' gather all right boundary nodes after sorting in boundaryNodes  '''
    tmp = []
    for b in ri:
        if b.node( 0 ) not in tmp:
            tmp.append( b.node( 0 ) );
        if b.node( 1 ) not in tmp:
            tmp.append( b.node( 1 ) );
    tmp.sort( sortNodeY ); tmp.reverse()
    boundaryNodes = tmp;

    ''' gather all bottom boundary nodes and add them to boundaryNodes '''
    boNode = []
    for b in bo:
        if b.node( 0 ) not in boNode + boundaryNodes:
            boNode.append( b.node( 0 ) );
        if b.node( 1 ) not in boNode + boundaryNodes:
            boNode.append( b.node( 1 ) );

    boNode.sort( sortNodeX ); boNode.reverse()
    boundaryNodes = boundaryNodes + boNode;

    ''' gather all left boundary nodes and add them to boundaryNodes '''
    tmp = []
    for b in le:
        if b.node( 0 ) not in tmp + boundaryNodes:
            tmp.append( b.node( 0 ) );
        if b.node( 1 ) not in tmp + boundaryNodes:
            tmp.append( b.node( 1 ) );

    tmp.sort( sortNodeY );
    boundaryNodes = boundaryNodes + tmp;

    if isSubSurface:
        ''' gather all top boundary nodes and add them to boundaryNodes '''
        topNodes = []
        for boundary in top:
            if boundary.node( 0 ) not in topNodes + boundaryNodes:
                topNodes.append( boundary.node( 0 ) )
            if boundary.node( 1 ) not in topNodes + boundaryNodes:
                topNodes.append( boundary.node( 1 ) )
        topNodes.sort( sortNodeX )
        boundaryNodes = boundaryNodes + topNodes


    poly = g.Mesh();

    if isSubSurface:
        ''' add all boundary nodes '''
        for n in boundaryNodes:
            poly.createNode( n.pos() );
        ''' and connect them by a closed polygon '''
        for id in range( 0, poly.nodeCount() ):
            poly.createEdge( poly.node( id ), poly.node( (id + 1)%poly.nodeCount() ), marker )
        ''' add four corners of the world box '''
        n1 = poly.createNode( g.RVector3( mesh.xmax() + xbound, surface, 0.0 ) );
        n2 = poly.createNode( g.RVector3( mesh.xmin() - xbound, surface, 0.0 ) );
        n3 = poly.createNode( g.RVector3( mesh.xmin() - xbound, mesh.ymin() - ybound, 0.0 ) );
        n4 = poly.createNode( g.RVector3( mesh.xmax() + xbound, mesh.ymin() - ybound, 0.0 ) );
        ''' and connect them by a closed polygon '''
        poly.createEdge( n1, n2, g.MARKER_BOUND_HOMOGEN_NEUMANN )
        poly.createEdge( n2, n3, g.MARKER_BOUND_MIXED )
        poly.createEdge( n3, n4, g.MARKER_BOUND_MIXED )
        poly.createEdge( n4, n1, g.MARKER_BOUND_MIXED )

    else:
        ''' add top right node and boundary nodes '''
        poly.createNode( g.RVector3( mesh.xmax() + xbound, mesh.ymax(), 0.0 ) );

        for n in boundaryNodes:
            poly.createNode( n.pos() );

        ''' add top left, bottom left and bottom right node '''
        poly.createNode( g.RVector3( mesh.xmin() - xbound, mesh.ymax(), 0.0 ) );
        poly.createNode( g.RVector3( mesh.xmin() - xbound, mesh.ymin() - ybound, 0.0 ) );
        poly.createNode( g.RVector3( mesh.xmax() + xbound, mesh.ymin() - ybound, 0.0 ) );

        ''' create a closed polygone through all new nodes '''
        for id in range( 0, poly.nodeCount() ):
            poly.createEdge( poly.node( id ), poly.node( (id + 1)%poly.nodeCount() ), marker )
    
    poly.exportVTK('out.poly')
    mesh2 = g.Mesh()
    ''' call triangle mesh generation '''
    triswitches= '-pzeAfaq' + str( quality )
    if not verbose:
        triswitches += 'Q'
    
    if isSubSurface:
        tri = g.TriangleWrapper( poly );
        tri.addRegionMarkerTmp( 0, g.RVector3( mesh.xmin() +1 , mesh.ymax() - 1 ), -1 );
        tri.setSwitches( triswitches );
        tri.generate( mesh2 );
    else:
        g.TriangleWrapper( poly, mesh2, triswitches );

    #mesh3.smooth( 1, 4, 1, 1 )
    map( lambda cell: cell.setMarker( marker ), mesh2.cells() )

    for cell in mesh.cells():
        mesh2.copyCell( cell );

    mesh2.createNeighbourInfos()
    for b in mesh2.boundaries():
        if b.leftCell() is None or b.rightCell() is None:
            if b.center().y() == mesh2.ymax():
                b.setMarker( g.MARKER_BOUND_HOMOGEN_NEUMANN )
            else:
                b.setMarker( g.MARKER_BOUND_MIXED )

    return mesh2


def appendTetrahedronBoundary( mesh, xbound = 100, ybound = 100, zbound = 100, marker = 1
                                , quality = 2, isSubSurface = False, verbose = False ):
    '''
        boundary of mesh need marker == 1
    '''                                
    
    ### create boundary for mesh from boundary marker == 1
    meshBoundary = g.Mesh()
    meshBoundary.createH2Mesh( mesh )

    bounds = []
    for b in meshBoundary.boundaries():
        if b.marker() == 1:
            bounds.append( b )

    meshBoundaryPoly = g.Mesh()
    meshBoundaryPoly.createMeshByBoundaries( meshBoundary, meshBoundary.findBoundaryByMarker( 1 ) );
    
    meshBoundaryPoly.exportAsTetgenPolyFile( "paraBoundary.poly" )
    #system( 'polyConvert -V paraBoundary' )

    ### create worldSurface.poly including boundary mesh for a nice surface mesh
    ### it will be later the tetgen input with preserve boundary
    polyCreateWorld( 'worldSurface', x = xbound, y = ybound, depth = zbound, marker = 1, verbose = verbose )
    system( 'polyMerge -N worldSurface paraBoundary worldSurface' )
    polyAddVIP( 'worldSurface', mesh.cell( 0 ).center(), isHoleMarker=True, verbose = verbose )
    worldBoundary = tetgen( 'worldSurface', quality = 1.12, verbose = verbose )
    #worldBoundary.exportBoundaryVTU('worldSurface')
    
    worldPoly = g.Mesh()
    worldPoly.createMeshByBoundaries( worldBoundary, worldBoundary.findBoundaryByMarker( -2, 0 ) );
    worldPoly.exportAsTetgenPolyFile( "worldSurface.poly" )
    
    system( 'polyMerge -N worldSurface paraBoundary boundaryWorld' )
    
    ### mesh should have to be a hole
    polyAddVIP( 'boundaryWorld', mesh.cell( 0 ).center(), isHoleMarker = True, verbose = verbose )
    #system( 'polyConvert -o world-poly -V boundaryWorld' )  

    boundMesh = tetgen( 'boundaryWorld', quality = quality, preserveBoundary = True, verbose = verbose )
    #boundMesh.exportVTK( 'boundaryWorld' )
    
    ### merge mesh and worldBoundary
    for c in boundMesh.cells():
        c.setMarker( marker )
        
    if verbose: 
        print "merge grid and boundary"
        
    swatch = g.Stopwatch( True )
    for c in meshBoundary.cells():
        nodes = g.stdVectorNodes( )
        for i, n in enumerate( c.nodes() ):
            nodes.append( boundMesh.createNodeWithCheck( n.pos() ) )
                
        boundMesh.createCell( nodes, c.marker() );
    if verbose:
        print " done.", swatch.duration( True )
    
    try:
        os.remove( 'boundaryWorld.bms' )
        os.remove( 'worldSurface.bms' )
        os.remove( 'boundaryWorld.poly')
        os.remove( 'paraBoundary.poly')
        os.remove( 'worldSurface.poly')
    except:
        None
    

    return boundMesh