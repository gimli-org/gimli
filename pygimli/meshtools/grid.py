# -*- coding: utf-8 -*-
"""General grid generation and maintenance."""

import os

import numpy as np

import pygimli as pg

from .polytools import polyCreateWorld, syscallTetgen


def createGrid(x=None, y=None, z=None, **kwargs):
    """Create grid style mesh.

    Generate simple grid with defined node positions for each dimension.
    The resulting grid depends on the amount of given coordinate arguments
    and consists out of edges (1D - x), quads (2D- x and y), or
    hexahedrons(3D- x, y, and z).

    Parameters
    ----------
    kwargs:
        x: array
            x-coordinates for all Nodes (1D, 2D, 3D)
        y: array
            y-coordinates for all Nodes (2D, 3D)
        z: array
            z-coordinates for all Nodes (3D)
        marker : int = 0
            Marker for resulting cells.
        worldBoundaryMarker : bool = False
            Boundaries are enumerated with world marker, i.e., Top = -1
            All remaining = -2.
            Default marker are left=1, right=2, top=3, bottom=4, front=5, back=6

    Examples
    --------
    >>> import pygimli as pg
    >>> mesh = pg.meshtools.createGrid(x=[0, 1, 1.5, 2], y=[-1, -0.5, -0.25, 0], 
    ...                                marker=2)
    >>> print(mesh)
    Mesh: Nodes: 16 Cells: 9 Boundaries: 24
    >>> fig, axs = plt.subplots(1, 2)
    >>> _ = pg.show(mesh, markers=True, showMesh=True, ax=axs[0])
    >>> mesh = pg.meshtools.createGrid(x=[0, 1, 1.5, 2], 
    ...                                y=[-1, -0.5, -0.25, 0], 
    ...                                worldBoundaryMarker=True, marker=2)
    >>> print(mesh)
    Mesh: Nodes: 16 Cells: 9 Boundaries: 24
    >>> _ = pg.show(mesh, markers=True, showBoundaries=True, 
    ...             showMesh=True, ax=axs[1])
    """
    if 'degree' in kwargs:
        return createGridPieShaped(x, **kwargs)

    if x is not None:
        if isinstance(x, int):
            x = list(range(x))
        kwargs['x'] = x
    if y is not None:
        if isinstance(y, int):
            y = list(range(y))
        kwargs['y'] = y
    if z is not None:
        if isinstance(z, int):
            z = list(range(z))
        kwargs['z'] = z

    return pg.core._pygimli_.createGrid(**kwargs)


def createGridPieShaped(x, degree=10.0, h=2, marker=0):
    """Create a 2D pie shaped grid (segment from annulus or cirlce).

    TODO:
    ----
    * degree: > 90 .. 360

    Arguments
    ---------
    x: array
        x-coordinates for all Nodes (2D). If you need it 3D, you can apply :py:mod:`pygimli.meshtools.extrudeMesh` on it.

    degree: float [None]
        Create a pie shaped grid for a value between 0 and 90.
        Creates an optional inner boundary (marker=2) for a annulus with x[0] > 0. Outer boundary marker is 1. Optional h refinement. Center node is the first for circle segment.

    h: int [2]
        H-Refinement for degree option.

    marker: int = 0
        Marker for resulting cells.

    Returns
    -------
    mesh: :gimliapi:`GIMLI::Mesh`

    Examples
    --------
    >>> import pygimli as pg
    >>> mesh = pg.meshtools.createGridPieShaped(x=[0, 1, 3], degree=45, h=3)
    >>> print(mesh)
    Mesh: Nodes: 117 Cells: 128 Boundaries: 244
    >>> _ = pg.show(mesh)
    >>> mesh = pg.meshtools.createGridPieShaped(x=[1, 2, 3], degree=45, h=3)
    >>> print(mesh)
    Mesh: Nodes: 153 Cells: 128 Boundaries: 280
    >>> _ = pg.show(mesh)
    """
    mesh = pg.Mesh(dim=2)

    for i in range(0, len(x)):
        mesh.createNodeWithCheck([x[i], 0.0])

        mesh.createNodeWithCheck([x[i]*np.cos(degree*np.pi/180),
                                  x[i]*np.sin(degree*np.pi/180)])

    if abs(x[0]) < 1e-6:
        mesh.createCell([0, 1, 2])
        for i in range(0, (len(x)-2)*2-1, 2):
            c = mesh.createCell([i+1, i+3, i+4, i+2])
    else:
        for i in range(0, len(x)*2-2, 2):
            c = mesh.createCell([i, i+2, i+3, i+1])
        mesh.createBoundary([0, 1], marker=1)

    mesh.createBoundary([mesh.nodeCount()-2, mesh.nodeCount()-1], marker=2)

    for i in range(h):
        mesh = mesh.createH2()
    mesh.createNeighbourInfos()

    for b in mesh.boundaries():
        if b.outside() and b.marker() == 0:
            if b.norm()[1] == 0.0:
                b.setMarker(4) # bottom
            else:
                b.setMarker(3)

    meshR = pg.Mesh(mesh)

    ## move all nodes on the inner boundary to rw
    for b in mesh.boundaries():
        line = pg.Line(b.node(0).pos(), b.node(1).pos())

        rSoll = line.intersect([0.0, 0.0], [1.0, 0.0])[0]
        if rSoll > 1e-4:
            for n in b.nodes():
                scale = rSoll/n.pos().abs()
                if scale > 1:
                    meshR.node(n.id()).setPos(pg.Line([0.0, 0.0], n.pos()).at(scale))

    if marker != 0:
        for c in meshR.cells():
            c.setMarker(marker)

    return meshR


def appendTriangleBoundary(mesh, xbound=10, ybound=10, marker=1,        
                           isSubSurface=False, **kwargs):
    """Add a triangle mesh boundary to a given mesh.

    Returns a new mesh that contains a triangulated box around a given mesh
    suitable for geo-simulation (surface boundary with marker = -1  at top and marker = -2 in the inner subsurface).
    The old boundary marker from mesh will be preserved, except for marker == -2 which will be switched to 2 since we assume
    -2 is the world marker for outer boundaries in the subsurface.

    Note, this all will only work stable if the mesh generator (triangle) preserve all input boundaries. 
    This will lead to bad quality meshes for the boundary region so its a good idea to play with the addNodes keword argument 
    to manually refine the newly created outer boundaries.

    Parameters
    ----------
    mesh : mesh object
        Mesh to which the triangle boundary should be appended.
    xbound : float, optional
        Absolute horizontal prolongation distance.
    ybound : float, optional
        Absolute vertical prolongation distance. 
    marker : int, optional
        Marker of new cells.
    isSubSurface : boolean, optional
        Apply boundary conditions suitable for geo-simulation and prolongate
        mesh to the surface if necessary.
    
    Additional Args
    ---------------
    ** kargs forwarded to pg.createMesh
    
    quality : float, optional
        Triangle quality.
    area: float, optional
        Triangle max size within the boundary.
    smooth : boolean, optional
        Apply mesh smoothing.
        
    addNodes : int[5], iterable
        Add aditional nodes on the outer boundaries.
    
    See Also
    --------
    appendTetrahedronBoundary
    
    Examples
    --------
    >>> import matplotlib.pyplot as plt
    >>> import pygimli as pg
    >>> from pygimli.viewer.mpl import drawMesh, drawModel
    >>> import pygimli.meshtools as mt
    >>> inner = pg.createGrid(range(5), range(5), marker=1)
    >>> fig, axs = plt.subplots(2,3)
    >>> ax, _ = pg.show(inner, markers=True, showBoundaries=True, showMesh=True, ax=axs[0][0])
    >>> m = mt.appendTriangleBoundary(inner, xbound=3, ybound=3, marker=2, addNodes=0, isSubSurface=False)
    >>> ax, _ = pg.show(m, markers=True, showBoundaries=True, showMesh=True, ax=axs[0][1])
    >>> m = mt.appendTriangleBoundary(inner, xbound=4, ybound=1, marker=2, addNodes=5, isSubSurface=False)
    >>> ax, _ = pg.show(m, markers=True, showBoundaries=True, showMesh=True, ax=axs[0][2])
    >>> m = mt.appendTriangleBoundary(inner, xbound=4, ybound=4, marker=2, addNodes=0, isSubSurface=True)
    >>> ax, _ = pg.show(m, markers=True, showBoundaries=True, showMesh=True, ax=axs[1][0])
    >>> m = mt.appendTriangleBoundary(inner, xbound=4, ybound=4, marker=2, addNodes=5, isSubSurface=True)
    >>> ax, _ = pg.show(m, markers=True, showBoundaries=True, showMesh=True, ax=axs[1][1])
    >>> surf = mt.createPolygon([[0, 0],[5, 3], [10, 1]], boundaryMarker=-1, addNodes=5, interpolate='spline')
    >>> m = mt.appendTriangleBoundary(surf, xbound=4, ybound=4, marker=2, addNodes=5, isSubSurface=True)
    >>> ax, _ = pg.show(m, markers=True, showBoundaries=True, showMesh=True, ax=axs[1][2])
    """
    poly = pg.Mesh(isGeometry=True)

    if isSubSurface == True:

        bs = mesh.findBoundaryByMarker(pg.core.MARKER_BOUND_HOMOGEN_NEUMANN)

        if len(bs) == 0:
            for b in mesh.boundaries():
                if b.outside() and b.norm()[1] == 1.0:
                    bs.append(b)
        
        paths = mesh.findPaths(bs)

        if len(paths) > 0:
            startPoint = mesh.node(paths[0][0]).pos()
            endPoint = mesh.node(paths[0][-1]).pos()

            if startPoint[0] > endPoint[0]:
                startPoint = endPoint
                endPoint = mesh.node(paths[0][0]).pos()
        else:
            pg.critical("Can't identify upper part of the mesh to be moved to the surface.",
                        "Maybe you can define them with Marker==-1")

        boundPoly = [pg.Pos(startPoint)]
        boundPoly.append(boundPoly[-1] - pg.Pos(xbound, 0))
        boundPoly.append(boundPoly[-1] - pg.Pos(0, mesh.ymax()- mesh.ymin() + ybound))
        boundPoly.append(boundPoly[-1] + pg.Pos((endPoint-startPoint)[0] + 2*xbound, 0))
        boundPoly.append(pg.Pos(endPoint) + pg.Pos(xbound, 0))
        boundPoly.append(pg.Pos(endPoint))

        poly = pg.meshtools.createPolygon(boundPoly, isClosed=False, 
                                          addNodes=kwargs.pop('addNodes', 5))
        poly.addRegionMarker(boundPoly[1] + [xbound/10, -ybound/10], marker=marker)
        
        if mesh.cellCount() > 0:
            poly.addHoleMarker(boundPoly[0] + [xbound/10, -ybound/10])

    else:  # no isSubSurface

        boundPoly = [
                        [mesh.xmin() - xbound, mesh.ymin() - ybound],
                        [mesh.xmin() - xbound, mesh.ymax() + ybound],
                        [mesh.xmax() + xbound, mesh.ymax() + ybound],
                        [mesh.xmax() + xbound, mesh.ymin() - ybound],
                    ]

        poly = pg.meshtools.createPolygon(boundPoly, isClosed=True, marker=marker,
                                          addNodes=kwargs.pop('addNodes', 5))

    for b in mesh.boundaries():
        if b.outside() or b.marker() == -1:
            poly.copyBoundary(b)

    preserveSwitch = 'Y'
    # pg.show(poly, boundaryMarkers=True, showNodes=True)
    # pg.wait()

    mesh2 = pg.meshtools.createMesh(poly, preserveBoundary=preserveSwitch, **kwargs)
    
    # pg.show(mesh2, boundaryMarkers=True, showNodes=True)

    ## start extracting all cells with marker from mesh2 and all orginal cells from mesh
    mesh3 = pg.Mesh(2)

    for c in mesh2.cells():
        if c.marker() == marker:
            mesh3.copyCell(c)
        
    ##! map does copies the cell not the reference, this should not happen **TODO check 20210305
    # map(lambda cell: mesh2.copyCell(cell), mesh2.cells())
    for c in mesh.cells():
        mesh3.copyCell(c)
    
    ## we need to delete the old boundary markers or the new neighbour infos will fail for old outside boundaries
    mesh3.setBoundaryMarkers(np.zeros(mesh3.boundaryCount()))
    mesh3.createNeighborInfos(force=True)

    for b in mesh.boundaries():
        if b.marker() != 0:
            b2 = mesh3.copyBoundary(b)
            
            # some automagic .. original mesh contains bmarker == -2 which means mixed condition
            # this special marker will be switched to 2
            if b.marker() == -2:
                b2.setMarker(2)
                
    for b in mesh3.boundaries():
        if b.outside() and b.marker() > -1:
            if b.norm().x() != 0 or b.norm().y() == -1.0:
                b.setMarker(pg.core.MARKER_BOUND_MIXED)
            else:
                b.setMarker(pg.core.MARKER_BOUND_HOMOGEN_NEUMANN)

    return mesh3


def appendTetrahedronBoundary(mesh, xbound=100, ybound=100, zbound=100,
                              marker=1, quality=2, isSubSurface=False,
                              verbose=False):
    """    Return new mesh surrounded by tetrahedron boundary box.

    Creates a tetrahedral box around a given mesh
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
    DeprecationWarning()
    # create boundary for mesh from boundary marker == 1

    if isSubSurface:
        raise Exception('Implement me')

    meshBoundary = pg.Mesh()
    meshBoundary.createH2()

    bounds = []
    for b in meshBoundary.boundaries():
        if b.marker() == 1:
            bounds.append(b)

    meshBoundaryPoly = pg.Mesh()
    meshBoundaryPoly.createMeshByBoundaries(
        meshBoundary, meshBoundary.findBoundaryByMarker(1))

    meshBoundaryPoly.exportAsTetgenPolyFile("paraBoundary.poly")
    # system( 'polyConvert -V paraBoundary' )

    # create worldSurface.poly including boundary mesh for a nice surface mesh
    # it will be later the tetgen input with preserve boundary
    polyCreateWorld('worldSurface', x=xbound, y=ybound, depth=zbound, marker=1,
                    verbose=verbose)
    os.system('polyMerge -N worldSurface paraBoundary worldSurface')
    # There syscalls are not working!
    polyAddVIP('worldSurface', mesh.cell(0).center(), isHoleMarker=True,
               verbose=verbose)
    worldBoundary = tetgen('worldSurface', quality=1.12, verbose=verbose)
    # worldBoundary.exportBoundaryVTU('worldSurface')

    worldPoly = pg.Mesh()
    worldPoly.createMeshByBoundaries(worldBoundary,
                                     worldBoundary.findBoundaryByMarker(-2, 0))
    worldPoly.exportAsTetgenPolyFile("worldSurface.poly")

    os.system('polyMerge -N worldSurface paraBoundary boundaryWorld')

    # mesh should have to be a hole
    polyAddVIP('boundaryWorld', mesh.cell(0).center(), isHoleMarker=True,
               verbose=verbose)
    # system( 'polyConvert -o world-poly -V boundaryWorld' )

    boundMesh = tetgen('boundaryWorld', quality=quality, preserveBoundary=True,
                       verbose=verbose)
    # boundMesh.exportVTK( 'boundaryWorld' )

    # merge mesh and worldBoundary
    for c in boundMesh.cells():
        c.setMarker(marker)

    if verbose:
        print("merge grid and boundary")

    swatch = pg.core.Stopwatch(True)
    for c in meshBoundary.cells():
        nodes = pg.stdVectorNodes()
        for n in c.nodes():
            nodes.append(boundMesh.createNodeWithCheck(n.pos()))

        boundMesh.createCell(nodes, c.marker())
    if verbose:
        print(" done.", swatch.duration(True))

    try:
        os.remove('boundaryWorld.bms')
        os.remove('worldSurface.bms')
        os.remove('boundaryWorld.poly')
        os.remove('paraBoundary.poly')
        os.remove('worldSurface.poly')
    except BaseException as e:
        print(e)

    return boundMesh


if __name__ == "__main__":
    pass
