# -*- coding: utf-8 -*-

import pygimli as g
from os import system

from pygimli.polytools import *

import numpy as np


def appendTriangleBoundary(mesh, xbound=10, ybound=10, marker=1, quality=34.0,
                           smooth=False, isSubSurface=False, verbose=False):
    """
    Returns a new mesh that contains a triangulated box around a given mesh
    suitable for geo-simulation (surface boundary at top).

    Parameters
    ----------
    mesh : mesh object
        Mesh to which the triangle boundary should be appended.
    xbound : float, optional
        Horizontal prolongation distance.
    ybound : float, optional
        Vertical prolongation distance.
    marker : int, optional
        Marker of new cells.
    quality : float, optional
        Triangle quality.
    smooth : boolean, optional
        Apply mesh smoothing.
    isSubSurface : boolean, optional
        Apply boundary conditions suitable for geo-simulaion and prolongate
        mesh to the surface if necessary.
    verbose : boolean, optional
        Be verbose.

    See Also
    --------
    appendTetrahedronBoundary
    """

    def sortNodeY(n1, n2):
        """ function comparing x for using sort """
        return cmp(n1.pos().y(), n2.pos().y())
    def sortNodeX(n1, n2):
        """ function comparing y for using sort """
        return cmp(n1.pos().x(), n2.pos().x())

    surface = 0.0
    # find boundaries on left/right/bottom/top side
    le = filter(lambda b: b.center().x() == mesh.xmin(), mesh.boundaries())
    bo = filter(lambda b: b.center().y() == mesh.ymin(), mesh.boundaries())
    ri = filter(lambda b: b.center().x() == mesh.xmax(), mesh.boundaries())
    top = filter(lambda b: b.center().y() == mesh.ymax(), mesh.boundaries())

    # gather all right boundary nodes after sorting in boundaryNodes
    tmp = []
    for b in ri:
        if b.node(0) not in tmp:
            tmp.append(b.node(0))
        if b.node(1) not in tmp:
            tmp.append(b.node(1))

    tmp.sort(sortNodeY)
    tmp.reverse()
    boundaryNodes = tmp

    # gather all bottom boundary nodes and add them to boundaryNodes
    boNode = []
    for b in bo:
        if b.node(0) not in boNode + boundaryNodes:
            boNode.append(b.node(0))
        if b.node(1) not in boNode + boundaryNodes:
            boNode.append(b.node(1))

    boNode.sort(sortNodeX)
    boNode.reverse()
    boundaryNodes = boundaryNodes + boNode

    # gather all left boundary nodes and add them to boundaryNodes
    tmp = []
    for b in le:
        if b.node(0) not in tmp + boundaryNodes:
            tmp.append(b.node(0))
        if b.node(1) not in tmp + boundaryNodes:
            tmp.append(b.node(1))

    tmp.sort(sortNodeY)
    boundaryNodes = boundaryNodes + tmp

    if isSubSurface:
        # gather all top boundary nodes and add them to boundaryNodes
        topNodes = []
        for boundary in top:
            if boundary.node(0) not in topNodes + boundaryNodes:
                topNodes.append(boundary.node(0))
            if boundary.node(1) not in topNodes + boundaryNodes:
                topNodes.append(boundary.node(1))
        topNodes.sort(sortNodeX)
        boundaryNodes = boundaryNodes + topNodes

    poly = g.Mesh()

    preserveSwitch = ''

    if isSubSurface:
        # add all boundary nodes
        for n in boundaryNodes:
            poly.createNode(n.pos())

        # and connect them by a closed polygon
        for i in range(0, poly.nodeCount()):
            poly.createEdge(
                poly.node(i), poly.node((i + 1) %
                                        poly.nodeCount()), marker)

        # add four corners of the world box
        xtLen = 12
        # x bottom boundary sampling points
        #xBottom = g.asvector( np.linspace( mesh.xmin() - xbound, mesh.xmax() + xbound, xtLen ) )

        n1 = poly.createNode(g.RVector3(mesh.xmax() + xbound, surface, 0.0))
        n2 = poly.createNode(g.RVector3(mesh.xmin() - xbound, surface, 0.0))
        n3 = poly.createNode(
            g.RVector3(mesh.xmin() - xbound,
                       mesh.ymin() - ybound,
                       0.0))
        n4 = poly.createNode(
            g.RVector3(mesh.xmax() + xbound,
                       mesh.ymin() - ybound,
                       0.0))
        # and connect them by a closed polygon
        poly.createEdge(n1, n2, g.MARKER_BOUND_HOMOGEN_NEUMANN)
        poly.createEdge(n2, n3, g.MARKER_BOUND_MIXED)
        poly.createEdge(n3, n4, g.MARKER_BOUND_MIXED)
        poly.createEdge(n4, n1, g.MARKER_BOUND_MIXED)

    else:
        # add top right node and boundary nodes
        xtLen = 12

        dxMin = boNode[0].pos().distance(boNode[1].pos()) * 1.5
        # x top boundary sampling points
        xTop = g.increasingRange(dxMin, xbound, xtLen)
        # y boundary sampling points
        yLeft = g.increasingRange(
            xTop[xtLen - 1] - xTop[xtLen - 2],
            abs(mesh.ymin() - ybound),
            xtLen)

        # x bottom boundary sampling points
        xBottom = g.asvector(
            np.linspace(mesh.xmin() - xbound,
                        mesh.xmax() + xbound,
                        xtLen))

        for t in g.fliplr(xTop)(0, len(xTop) - 1):
            poly.createNode(g.RVector3(mesh.xmax() + t, mesh.ymax(), 0.0))

        for n in boundaryNodes:
            poly.createNode(n.pos())

        # add top left, bottom left and bottom right node

        for t in xTop(1, len(xTop)):
            poly.createNode(g.RVector3(mesh.xmin() - t, mesh.ymax(), 0.0))

        for t in yLeft(1, len(yLeft)):
            poly.createNode(
                g.RVector3(mesh.xmin() - xbound,
                           mesh.ymax() - t,
                           0.0))

        for t in xBottom(1, len(xBottom) - 1):
            poly.createNode(g.RVector3(t, mesh.ymin() - ybound, 0.0))

        for t in g.fliplr(yLeft)(0, len(yLeft) - 1):
            poly.createNode(
                g.RVector3(mesh.xmax() + xbound,
                           mesh.ymax() - t,
                           0.0))

        # create a closed polygon through all new nodes
        for i in range(0, poly.nodeCount()):
            poly.createEdge(
                poly.node(i), poly.node((i + 1) %
                                        poly.nodeCount()), marker)

        preserveSwitch = 'Y'
    # poly.exportVTK('out.poly')

    mesh2 = g.Mesh()

    # call triangle mesh generation
    triswitches = '-pzeAfa' + preserveSwitch + 'q' + str(quality)

    if not verbose:
        triswitches += 'Q'

    if isSubSurface:
        tri = g.TriangleWrapper(poly)
        # area -1.0 means this is a hole
        tri.addRegionMarkerTmp(
            0,
            g.RVector3(mesh.xmin() + 0.0001,
                       mesh.ymax() - 0.0001),
            -1.0)
        tri.setSwitches(triswitches)
        tri.generate(mesh2)
    else:
        g.TriangleWrapper(poly, mesh2, triswitches)

    if smooth:
        mesh2.smooth(
            nodeMoving=True,
            edgeSwapping=True,
            smoothFunction=1,
            smoothIteration=2)

    map(lambda cell: cell.setMarker(marker), mesh2.cells())

    #! map copy the cell not the reference, this should not happen
    #! map( lambda cell: mesh2.copyCell( cell ), mesh2.cells() )
    for cell in mesh.cells():
        mesh2.copyCell(cell)

    #! old neighbor infos need to be cleaned since the new cells are added
    mesh2.createNeighbourInfos(force=True)

    for b in mesh2.boundaries():
        if b.leftCell() is None or b.rightCell() is None:
            if b.center().y() == mesh2.ymax():
                b.setMarker(g.MARKER_BOUND_HOMOGEN_NEUMANN)
            else:
                b.setMarker(g.MARKER_BOUND_MIXED)

    return mesh2


def appendTetrahedronBoundary(mesh, xbound=100, ybound=100, zbound=100,
                              marker=1, quality=2, isSubSurface=False,
                              verbose=False):
    """
    Returns a new mesh that contains a tetrahedral box around a given mesh
    suitable for geo-simulation (surface boundary at top).

    Parameters
    ----------
    mesh : mesh object
        Mesh to which the tetrahedron boundary should be appended.
    xbound : float, optional
        Horizontal prolongation distance in x-direction.
    ybound : float, optional
        Horizonal prolongation distance in y-direction.
    zbound : float, optional
        Vertical prolongation distance.
    marker : int, optional
        Marker of new cells.
    quality : float, optional
        Triangle quality.
    isSubSurface : boolean, optional
        Apply boundary conditions suitable for geo-simulaion and prolongate
        mesh to the surface if necessary.
    verbose : boolean, optional
        Be verbose.

    See Also
    --------
    appendTriangleBoundary

    Notes
    -----
    Boundaries of mesh need marker 1.
    """
    # create boundary for mesh from boundary marker == 1
    meshBoundary = g.Mesh()
    meshBoundary.createH2Mesh(mesh)

    bounds = []
    for b in meshBoundary.boundaries():
        if b.marker() == 1:
            bounds.append(b)

    meshBoundaryPoly = g.Mesh()
    meshBoundaryPoly.createMeshByBoundaries(
        meshBoundary,
        meshBoundary.findBoundaryByMarker(1))

    meshBoundaryPoly.exportAsTetgenPolyFile("paraBoundary.poly")
    #system( 'polyConvert -V paraBoundary' )

    # create worldSurface.poly including boundary mesh for a nice surface mesh
    # it will be later the tetgen input with preserve boundary
    polyCreateWorld(
        'worldSurface',
        x=xbound,
        y=ybound,
        depth=zbound,
        marker=1,
        verbose=verbose)
    system('polyMerge -N worldSurface paraBoundary worldSurface')
    polyAddVIP(
        'worldSurface',
        mesh.cell(0).center(),
        isHoleMarker=True,
        verbose=verbose)
    worldBoundary = tetgen('worldSurface', quality=1.12, verbose=verbose)
    # worldBoundary.exportBoundaryVTU('worldSurface')

    worldPoly = g.Mesh()
    worldPoly.createMeshByBoundaries(
        worldBoundary,
        worldBoundary.findBoundaryByMarker(
            -2,
            0))
    worldPoly.exportAsTetgenPolyFile("worldSurface.poly")

    system('polyMerge -N worldSurface paraBoundary boundaryWorld')

    # mesh should have to be a hole
    polyAddVIP(
        'boundaryWorld',
        mesh.cell(0).center(),
        isHoleMarker=True,
        verbose=verbose)
    #system( 'polyConvert -o world-poly -V boundaryWorld' )

    boundMesh = tetgen(
        'boundaryWorld',
        quality=quality,
        preserveBoundary=True,
        verbose=verbose)
    #boundMesh.exportVTK( 'boundaryWorld' )

    # merge mesh and worldBoundary
    for c in boundMesh.cells():
        c.setMarker(marker)

    if verbose:
        print "merge grid and boundary"

    swatch = g.Stopwatch(True)
    for c in meshBoundary.cells():
        nodes = g.stdVectorNodes()
        for i, n in enumerate(c.nodes()):
            nodes.append(boundMesh.createNodeWithCheck(n.pos()))

        boundMesh.createCell(nodes, c.marker())
    if verbose:
        print " done.", swatch.duration(True)

    try:
        os.remove('boundaryWorld.bms')
        os.remove('worldSurface.bms')
        os.remove('boundaryWorld.poly')
        os.remove('paraBoundary.poly')
        os.remove('worldSurface.poly')
    except:
        None

    return boundMesh
