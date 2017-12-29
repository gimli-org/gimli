# -*- coding: utf-8 -*-
"""Tools to create or manage PLC.

Please note there is currently no collision or intersection check at all.

Volunteers welcome to help creating, adapting or interfacing a basic
geometry system. A lot of thinks are needed:

    * 2D
    * 3D
    * More geometric primitives
    * Boolean operations (union, intersection, difference)
    * Collision recognizing
    * Cubic spline interpolation for polygons (partly done)
    * GUI .. interactive creation
    *
"""

import math
import os
from os import system

import numpy as np

import pygimli as pg


def polyCreateDefaultEdges_(poly, boundaryMarker=1,
                            isClosed=True, isHole=False):
    """INTERNAL."""
    nEdges = poly.nodeCount() - 1 + isClosed
    bm = None
    if hasattr(boundaryMarker, '__len__'):
        if len(boundaryMarker) == nEdges:
            bm = boundaryMarker
        else:
            raise Exception("marker length != nEdges", len(boundaryMarker),
                            nEdges)
    else:
        bm = [boundaryMarker] * nEdges

    for i in range(poly.nodeCount() - 1):
        poly.createEdge(poly.node(i), poly.node(i + 1), bm[i])

    if isClosed:
        poly.createEdge(poly.node(poly.nodeCount() - 1), poly.node(0), bm[-1])


def createRectangle(start=None, end=None, pos=None, size=None, **kwargs):
    """Create rectangle polygon.

    Create rectangle with start position and a given size.
    Give either start and end OR pos and size.

    Parameters
    ----------
    start : [x, y]
        Left upper corner. Default [-0.5, 0.5]
    end : [x, y]
        Right lower corner. Default [0.5, -0.5]
    pos : [x, y]
        Center position. The rectangle will be moved.
    size : [x, y]
        width and height. The rectangle will be scaled.

    **kwargs:

        marker : int [1]
            Marker for the resulting triangle cells after mesh generation
        area : float [0]
            Maximum cell size for resulting triangles after mesh generation
        boundaryMarker : int [1]
            Marker for the resulting boundary edges
        leftDirection : bool [True]
            TODO Rotational direction
        isHole : bool [False]
            The Polygone will become a hole instead of a triangulation
        isClosed : bool [True]
            Add closing edge between last and first node.

    Returns
    -------
    poly : :gimliapi:`GIMLI::Mesh`
        The resulting polygon is a :gimliapi:`GIMLI::Mesh`.

    Examples
    --------
    >>> from pygimli.meshtools import createRectangle
    >>> from pygimli.mplviewer import drawMesh
    >>> import matplotlib.pyplot as plt
    >>> rectangle = createRectangle(start=[4, -4], end=[6, -6],
    ...                             marker=4, area=0.1)
    >>>
    >>> fig, ax = plt.subplots()
    >>> drawMesh(ax, rectangle)
    >>> plt.show()
    """
    if not ((start and end) or (pos and size)):
        raise BaseException("createRectangle pls. give either start and end"
                            "OR pos and size.")
    if start is None:
        start = [-0.5, 0.5]
    if end is None:
        end = [0.5, -0.5]
    poly = pg.Mesh(2)

    sPos = pg.RVector3(start)
    ePos = pg.RVector3(end)

    poly.createNode(sPos)
    poly.createNode([sPos[0], ePos[1]])
    poly.createNode(ePos)
    poly.createNode([ePos[0], sPos[1]])

    if kwargs.pop('isHole', False):
        poly.addHoleMarker(sPos + (ePos - sPos) * 0.2)
    else:
        poly.addRegionMarker(sPos + (ePos - sPos) * 0.2,
                             marker=kwargs.pop('marker', 1),
                             area=kwargs.pop('area', 0))

    if size is not None:
        poly.scale(size)
    if pos is not None:
        poly.translate(pos)

    polyCreateDefaultEdges_(poly, **kwargs)

    return poly


def createWorld(start, end, marker=1, area=0., layers=None, worldMarker=True):
    """Create simple rectangular world.

    Create simple rectangular world with appropriate boundary conditions.
    Surface boundary is set do pg.MARKER_BOUND_HOMOGEN_NEUMANN, i.e, -1 and
    inner subsurface is set to pg.MARKER_BOUND_MIXED, i.e., -2 or
    Numbered in ascending order in left direction starting upper left if
    worldMarker is set to false.

    Parameters
    ----------
    start : [x, y]
        Upper/Left Corner
    end : [x, y]
        Lower/Right Corner
    marker : int
        Marker for the resulting triangle cells after mesh generation.
    area : float | list
        Maximum cell size for resulting triangles after mesh generation.
        If area is a float set it global, if area is a list set it per layer.

    layers : [float]
        Add some layers to the world.
    worldMarker : [bool]
        Specify kind of preset boundary marker [-1, -2] or [1, 2, 3, 4 ..]

    Returns
    -------
    poly : :gimliapi:`GIMLI::Mesh`
        The resulting polygon is a :gimliapi:`GIMLI::Mesh`.

    Examples
    --------
    >>> from pygimli.meshtools import createWorld
    >>> from pygimli.mplviewer import drawMesh
    >>> import matplotlib.pyplot as plt
    >>> world = createWorld(start=[-5, 0], end=[5, -5], layers=[-1,-2,-3])
    >>>
    >>> fig, ax = plt.subplots()
    >>> drawMesh(ax, world)
    >>> plt.show()
    """
    z = [start[1]]

    if layers is not None:
        z = z + list(layers)

    z.append(end[1])

    # ensure - decreasing order if layers are out of bounding box
    z = pg.sort(z)[::-1]

    poly = pg.Mesh(2)

    if type(area) == float or type(area) == int:
        area = np.ones(len(z)) * float(area)

    for i, depth in enumerate(z):
        n = poly.createNode([start[0], depth])

        if i > 0:
            if len(z) == 2:
                poly.addRegionMarker(n.pos() + [0.2, 0.2], marker=marker,
                                     area=area[0])
            else:
                poly.addRegionMarker(n.pos() + [0.2, 0.2], marker=i,
                                     area=area[i-1])

    for i, depth in enumerate(z[::-1]):
        poly.createNode([end[0], depth])

    polyCreateDefaultEdges_(poly,
                            boundaryMarker=range(1, poly.nodeCount() + 1))

    if worldMarker:
        for b in poly.boundaries():
            if b.norm()[1] == 1.0:
                b.setMarker(pg.MARKER_BOUND_HOMOGEN_NEUMANN)
            else:
                b.setMarker(pg.MARKER_BOUND_MIXED)

    if layers is not None:
        for i in range(len(layers)):
            poly.createEdge(
                poly.node(i + 1), poly.node(poly.nodeCount() - i - 2),
                poly.boundaryCount() + 1)

    # pg.warnNonEmptyArgs(kwargs)
    return poly


def createCircle(pos=None, radius=1, segments=12, start=0, end=2. * math.pi,
                 **kwargs):
    """Create simple circle polygon.

    Create simple circle polygon with given attributes.

    Parameters
    ----------
    pos : [x, y] [[0.0, 0.0]]
        Center position
    radius : float | [a,b] [1]
        radius or halfaxes of the circle
    segments : int
        Discrete amount of segments for the circle
    start : double [0]
        Starting angle in radians
    end : double [2*pi]
        Ending angle in radians

    **kwargs:

        marker : int [1]
            Marker for the resulting triangle cells after mesh generation
        area : float [0]
            Maximum cell size for resulting triangles after mesh generation
        boundaryMarker : int [1]
            Marker for the resulting boundary edges
        leftDirection : bool [True]
            Rotational direction
        isHole : bool [False]
            The Polygone will become a hole instead of a triangulation
        isClosed : bool [True]
            Add closing edge between last and first node.

    Returns
    -------
    poly : :gimliapi:`GIMLI::Mesh`
        The resulting polygon is a :gimliapi:`GIMLI::Mesh`.

    Examples
    --------
    >>> import matplotlib.pyplot as plt
    >>> import math
    >>> from pygimli.mplviewer import drawMesh
    >>> from pygimli.meshtools import polytools as plc
    >>> c0 = plc.createCircle(pos=(-5.0, 0.0), radius=2, segments=6)
    >>> c1 = plc.createCircle(pos=(0.0, 0.0), segments=5, start=0, end=math.pi)
    >>> c2 = plc.createCircle(pos=(5.0, 0.0), segments=3, start=math.pi,
    ...                       end=1.5*math.pi, isClosed=False)
    >>> plc = plc.mergePLC([c0, c1, c2])
    >>> fig, ax = plt.subplots()
    >>> drawMesh(ax, plc)
    >>> plt.show()
    """
    if pos is None:
        pos = [0.0, 0.0]

    poly = pg.Mesh(2)

    dPhi = (end - start) / (segments)
    nPhi = segments + 1

    if abs((end % (2. * math.pi) - start)) < 1e-6:
        nPhi = segments

    for i in range(0, nPhi):
        if kwargs.pop('leftDirection', True):
            phi = start + i * dPhi
        else:
            phi = start - i * dPhi

        xp = np.cos(phi)
        yp = np.sin(phi)
        poly.createNode([xp, yp])

    if kwargs.pop('isHole', False):
        poly.addHoleMarker([0.0, 0.0])
    else:
        poly.addRegionMarker([0.0, 0.0], marker=kwargs.pop('marker', 1),
                             area=kwargs.pop('area', 0))

    if hasattr(radius, '__len__'):
        poly.scale(radius)
    else:
        poly.scale([radius, radius])

    poly.translate(pos)

    polyCreateDefaultEdges_(poly, **kwargs)

    # need a better way mess with these or wrong kwargs
    # pg.warnNonEmptyArgs(kwargs)

    return poly


def createLine(start, end, segments, **kwargs):
    """Create simple line polygon.

    Create simple line polygon from start to end.

    Parameters
    ----------
    start : [x, y]
        start position
    end : [x, y]
        end position
    segments : int
        Discrete amount of segments for the line

    **kwargs:

        boundaryMarker : int [1]
            Marker for the resulting boundary edges
        leftDirection : bool [True]
            Rotational direction

    Returns
    -------
    poly : :gimliapi:`GIMLI::Mesh`
        The resulting polygon is a :gimliapi:`GIMLI::Mesh`.

    Examples
    --------
    >>>  # no need to import matplotlib. pygimli's show does
    >>> import pygimli as pg
    >>> import pygimli.meshtools as mt
    >>>
    >>> w = mt.createWorld(start=[0, 0], end=[3, 3])
    >>> l1 = mt.createLine(start=[1, 1], end=[1, 2], segments=1,
    ...                    leftDirection=False)
    >>> l2 = mt.createLine(start=[1, 1], end=[2, 1], segments=20,
    ...                    leftDirection=True)
    >>>
    >>> ax, _ = pg.show(mt.createMesh([w, l1, l2,]))
    >>> ax, _ = pg.show([w, l1, l2,], ax=ax)
    >>> pg.wait()
    """
    poly = pg.Mesh(2)
    startPos = pg.RVector3(start)
    endPos = pg.RVector3(end)
    a = endPos - startPos

    dt = 1. / segments
    left = kwargs.pop('leftDirection', True)

    for i in range(0, segments + 1):
        if left:
            p = startPos + a * (dt * i)
        else:
            p = endPos - a * (dt * i)

        poly.createNode(p)

    polyCreateDefaultEdges_(poly, isClosed=False, **kwargs)
    return poly


def createPolygon(verts, isClosed=False, isHole=False, **kwargs):
    """Create a polygon from a list of vertices.

    All vertices needs to be unique and duplicate vertices will be ignored.
    If you want the polygon be a closed region you can set the 'isCloses' flag.
    Closed region can be attributed by assigning a region marker.
    The automatic region marker is set in the center of all vertices.

    Parameters
    ----------
    verts : []
        * List of x y pairs [[x0, y0], ... ,[xN, yN]]

    **kwargs:

        * boundaryMarker : int [1]
            Marker for the resulting boundary edges
        * leftDirection : bool [True]
            Rotational direction
        * marker : int [None]
            Marker for the resulting triangle cells after mesh generation.
        * area : float [0]
            Maximum cell size for resulting triangles after mesh generation
        * isHole : bool [False]
            The Polygone will become a hole instead of a triangulation

    isClosed : bool [True]
        Add closing edge between last and first node.

    Returns
    -------
    poly : :gimliapi:`GIMLI::Mesh`
        The resulting polygon is a :gimliapi:`GIMLI::Mesh`.

    Examples
    --------
    >>>  # no need to import matplotlib. pygimli's show does
    >>> import pygimli as pg
    >>> import pygimli.meshtools as mt
    >>> p1 = mt.createPolygon([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0]],
    ...                       isClosed=True, marker=3, area=0.1)
    >>> p2 = mt.createPolygon([[0.3, 0.15], [0.85, 0.15], [0.85, 0.7]],
    ...                       isClosed=True, marker=3, area=0.1, isHole=True)
    >>> ax, _ = pg.show(mt.mergePLC([p1,p2]))
    """
    poly = pg.Mesh(2)

    for v in verts:
        poly.createNodeWithCheck(v, warn=True)

    marker = kwargs.pop('marker', None)
    area = kwargs.pop('area', 0)
    isHole = kwargs.pop('isHole', False)

    polyCreateDefaultEdges_(poly, isClosed=isClosed, isHole=False)

    if isClosed and marker is not None or area > 0:
        if isHole:
            poly.addHoleMarker(pg.center(poly.positions()))
        else:
            if marker is None:  # in case marker is None but area is given
                marker = 0
            poly.addRegionMarker(pg.center(poly.positions()),
                                 marker=marker, area=area)

    # set a regionmarker here .. somewhere

    return poly


def mergePLC(pols):
    """Merge multiply polygons.

    Merge multiply polygons into a single polygon.
    Common nodes and common edges will be checked and removed.

    Crossing or touching edges or Node/Edge intersections will NOT be
    recognized yet. -> TODO

    Parameters
    ----------
    pols: [:gimliapi:`GIMLI::Mesh`]
        List of polygons that need to be merged

    Returns
    -------
    poly : :gimliapi:`GIMLI::Mesh`
        The resulting polygon is a :gimliapi:`GIMLI::Mesh`.

    Examples
    --------
    >>> from pygimli.meshtools import polytools as plc
    >>> from pygimli.meshtools import createMesh
    >>> from pygimli.mplviewer import drawMesh
    >>> import matplotlib.pyplot as plt
    >>> world = plc.createWorld(start=[-10, 0], end=[10, -10], marker=1)
    >>> c1 = plc.createCircle([-1, -4], radius=1.5, area=0.1,
    ...                       marker=2, segments=4)
    >>> c2 = plc.createCircle([-6, -5], radius=[1.5, 3.5], isHole=1)
    >>> r1 = plc.createRectangle(pos=[3, -5], size=[2, 2], marker=3)
    >>> r2 = plc.createRectangle(start=[4, -4], end=[6, -6],
    ...                          marker=4, area=0.1)
    >>> plc = plc.mergePLC([world, c1, c2, r1, r2])
    >>> fig, ax = plt.subplots()
    >>> drawMesh(ax, plc)
    >>> drawMesh(ax, createMesh(plc))
    >>> plt.show()
    """
    poly = pg.Mesh(2)

    for p in pols:
        nodes = []
        for n in p.nodes():
            nn = poly.createNodeWithCheck(n.pos())
            if n.marker() != 0:
                nn.setMarker(n.marker())
            nodes.append(nn)

        for e in p.boundaries():
            poly.createEdge(nodes[e.node(0).id()], nodes[e.node(1).id()],
                            e.marker())

        if len(p.regionMarker()) > 0:
            for rm in p.regionMarker():
                poly.addRegionMarker(rm)

        if len(p.holeMarker()) > 0:
            for hm in p.holeMarker():
                poly.addHoleMarker(hm)

    return poly


def createParaDomain2D(*args, **kwargs):
    """API change here .. use createParaMeshPLC instead."""
    print("createParaDomain2D: API change: use createParaMeshPLC instead")
    return createParaMeshPLC(*args, **kwargs)


def createParaMeshPLC(sensors, paraDX=1, paraDepth=0, paraBoundary=2,
                      paraMaxCellSize=0.0, boundary=-1, boundaryMaxCellSize=0,
                      **kwargs):
    """Create a PLC mesh for an inversion parameter mesh.

    Create a PLC mesh for an inversion parameter mesh with for a given list of
    sensor positions.
    Sensor position assumed on the surface and must be sorted and unique.

    You can create a parameter mesh without sensors if you just set
    [xmin, xmax] as sensors.

    The PLC is a :gimliapi:`GIMLI::Mesh` and contain nodes, edges and
    two region markers, one for the parameters domain (marker=2) and
    a larger boundary around the outside (marker=1)

    TODO:

        * closed domains (boundary == 0)
        * additional topopoints
        * spline interpolations between sensorpoints or addpoints
        * subsurface sensors (partly .. see example)

    Parameters
    ----------
    sensors : [RVector3] | DataContainer with sensorPositions() | [xmin, xmax]
        Sensor positions. Must be sorted and unique in positive x direction.
        Depth need to be y-coordinate.
    paraDX : float [1]
        Relativ distance for refinement nodes between two electrodes (1=none),
        e.g., 0.5 means 1 additional node between two neighboring electrodes
        e.g., 0.33 means 2 additional equidistant nodes between two electrodes
    paraDepth : float, optional
        Maximum depth for parametric domain, 0 (default) means 0.4 * maximum
        sensor range.
    paraBoundary : float, optional
        Margin for parameter domain in absolute sensor distances. 2 (default).
    paraMaxCellSize: double, optional
        Maximum size for parametric size in m*m
    boundaryMaxCellSize: double, optional
        Maximum cells size in the boundary region in m*m
    boundary : float, optional
        Boundary width to be appended for domain prolongation in absolute
        para domain width.
        Values lover 0 force the boundary to be 4 times para domain width.

    Returns
    -------
    poly: :gimliapi:`GIMLI::Mesh`
        piecewise linear complex (PLC) containing nodes and edges

    Examples
    --------
    >>>  # no need to import matplotlib. pygimli's show does
    >>> import pygimli as pg
    >>> import pygimli.meshtools as plc
    >>> # Create the simplest paramesh PLC with a para box of 10 m without
    >>> # electrodes
    >>> p = plc.createParaMeshPLC([0,10])
    >>> # you can add subsurface electrodes now with
    >>> for z in range(1,4):
    ...     n = p.createNode((5,-z), -99)
    >>> ax,_ = pg.show(p)
    """
    noSensors = False
    if hasattr(sensors, 'sensorPositions'):  # obviously a DataContainer type
        sensors = sensors.sensorPositions()
    elif isinstance(sensors, np.ndarray):
        sensors = [pg.RVector3(s) for s in sensors]
    elif isinstance(sensors, list):
        if len(sensors) == 2:
            # guess we have just a desired Pbox with
            sensors = [pg.RVector3(sensors[0], 0.0),
                       pg.RVector3(sensors[1], 0.0)]
            noSensors = True
            paraBoundary = 0

    eSpacing = kwargs.pop('eSpacing', sensors[0].distance(sensors[1]))

    iz = 1
    xmin, ymin, zmin = sensors[0][0], sensors[0][1], sensors[0][2]
    xmax, ymax, zmax = xmin, ymin, zmin
    for e in sensors:
        xmin = min(xmin, e[0])
        xmax = max(xmax, e[0])
        ymin = min(ymin, e[1])
        ymax = max(ymax, e[1])
        zmin = min(zmin, e[2])
        zmax = max(zmax, e[2])

    if abs(ymin) < 1e-8 and abs(ymax) < 1e-8:
        iz = 2

    paraBound = eSpacing * paraBoundary

    if paraDepth == 0:
        paraDepth = 0.4 * (xmax - xmin)

    poly = pg.Mesh(2)
    # define para domain without surface
    n1 = poly.createNode([xmin - paraBound, sensors[0][iz]])
    n2 = poly.createNode([xmin - paraBound, sensors[0][iz] - paraDepth])
    n3 = poly.createNode([xmax + paraBound, sensors[-1][iz] - paraDepth])
    n4 = poly.createNode([xmax + paraBound, sensors[-1][iz]])

    if boundary < 0:
        boundary = 4

    bound = abs(xmax - xmin) * boundary
    if bound > paraBound:
        # define world without surface
        n11 = poly.createNode(n1.pos() - [bound, 0.])
        n12 = poly.createNode(n11.pos() - [0., bound + paraDepth])
        n14 = poly.createNode(n4.pos() + [bound, 0.])
        n13 = poly.createNode(n14.pos() - [0., bound + paraDepth])

        poly.createEdge(n1, n11, pg.MARKER_BOUND_HOMOGEN_NEUMANN)
        poly.createEdge(n11, n12, pg.MARKER_BOUND_MIXED)
        poly.createEdge(n12, n13, pg.MARKER_BOUND_MIXED)
        poly.createEdge(n13, n14, pg.MARKER_BOUND_MIXED)
        poly.createEdge(n14, n4, pg.MARKER_BOUND_HOMOGEN_NEUMANN)
        poly.addRegionMarker(n12.pos() + [1e-3, 1e-3], 1, boundaryMaxCellSize)

    poly.createEdge(n1, n2, 1)
    poly.createEdge(n2, n3, 1)
    poly.createEdge(n3, n4, 1)
    poly.addRegionMarker(n2.pos() + [1e-3, 1e-3], 2, paraMaxCellSize)

    # define surface
    nSurface = []
    nSurface.append(n1)
    if paraDX == 0.0:
        paraDX = 1.0

    if not noSensors:
        for i, e in enumerate(sensors):
            if iz == 2:
                e.rotateX(-math.pi / 2)
            if paraDX >= 0.5:
                nSurface.append(poly.createNode(e, pg.MARKER_NODE_SENSOR))
                if i < len(sensors) - 1:
                    e1 = sensors[i + 1]
                    if iz == 2:
                        e1.rotateX(-math.pi / 2)
                    nSurface.append(poly.createNode((e + e1) * 0.5))
                # print("Surface add ", e, el, nSurface[-2].pos(),
                #        nSurface[-1].pos())
            elif paraDX < 0.5:
                if i > 0:
                    e1 = sensors[i - 1]
                    if iz == 2:
                        e1.rotateX(-math.pi / 2)
                    nSurface.append(poly.createNode(e - (e - e1) * paraDX))
                nSurface.append(poly.createNode(e, pg.MARKER_NODE_SENSOR))
                if i < len(sensors) - 1:
                    e1 = sensors[i + 1]
                    if iz == 2:
                        e1.rotateX(-math.pi / 2)
                    nSurface.append(poly.createNode(e + (e1 - e) * paraDX))
                # print("Surface add ", nSurface[-3].pos(), nSurface[-2].pos(),
                #        nSurface[-1].pos())
    nSurface.append(n4)

    for i in range(len(nSurface) - 1, 0, -1):
        poly.createEdge(nSurface[i], nSurface[i - 1],
                        pg.MARKER_BOUND_HOMOGEN_NEUMANN)

    # print(poly)
    # pg.meshtools.writePLC(poly, "test.poly")
    # pg.show(poly)
    # pg.wait()
    return poly


def readPLC(filename, comment='#'):
    r"""Generic PLC reader.

    Read 2D :term:`Triangle` or 3D :term:`Tetgen` PLC files.

    Parameters
    ----------
    filename: string
        Filename *.poly

    comment: string ('#')
        String containing all characters that define a comment line. Identified
        lines will be ignored during import.

    Returns
    -------
    poly :
        :gimliapi:`GIMLI::Mesh`
    """
    with open(filename, 'r') as fi:
        content = fi.readlines()

    # Filter comment lines
    comment_lines = []
    for i, line in enumerate(content):
        if line[0] in comment:
            comment_lines.append(i)
    for j in comment_lines[::-1]:
        del(content[j])

    # Read header
    headerLine = content[0].split()

    if len(headerLine) != 4:
        raise Exception("Format unknown! header size != 4", headerLine)

    fromOne = 0
    nVerts = int(headerLine[0])
    dimension = int(headerLine[1])
    nPointsAttributes = int(headerLine[2])
    haveNodeMarker = int(headerLine[3])

    poly = pg.Mesh(dimension)

    # Nodes section
    for i in range(nVerts):
        row = content[1 + i].split()

        if len(row) == (1 + dimension + nPointsAttributes + haveNodeMarker):
            if i == 0:
                fromOne = int(row[0])
            if dimension == 2:
                n = poly.createNode((float(row[1]), float(row[2])))
            elif dimension == 3:
                n = poly.createNode((float(row[1]), float(row[2]),
                                     float(row[3])))
            if haveNodeMarker:
                n.setMarker(int(row[-1]))

        else:
            raise Exception("Poly file seams corrupt: node section line: " +
                            str(i) + " " + row)

    # Segment section
    row = content[1 + nVerts].split()

    if len(row) != 2:
        raise Exception("Format unknown for segment section " + row)

    nSegments = int(row[0])
    haveBoundaryMarker = int(row[1])

    if dimension == 2:
        for i in range(nSegments):
            row = content[2 + nVerts + i].split()

            if len(row) == (3 + haveBoundaryMarker):
                marker = 0
                if haveBoundaryMarker:
                    marker = int(row[3])

                poly.createEdge(
                    poly.node(int(row[1]) - fromOne),
                    poly.node(int(row[2]) - fromOne), marker)
    else:
        segment_offset = 0
        for i in range(nSegments):
            row = content[2 + nVerts + i + segment_offset].split()
            numBounds = int(row[0])
            numHoles = row[1]
            assert numHoles == '0', 'Can\'t handle Boundaries with holes yet'
            marker = 0
            if haveBoundaryMarker:
                marker = int(row[2])

            for k in range(numBounds):
                bound_row = content[2 + nVerts + i + segment_offset + 1]\
                    .split()
                ivec = [int(bound_row[1]), int(bound_row[2]),
                        int(bound_row[3]), int(bound_row[4])]
                poly.createBoundary(ivec, marker=marker)
                segment_offset += 1
        nSegments += segment_offset

    # Hole section
    row = content[2 + nVerts + nSegments].split()

    if len(row) != 1:
        raise Exception("Format unknown for hole section " + row)

    nHoles = int(row[0])
    for i in range(nHoles):
        row = content[3 + nVerts + nSegments + i].split()
        if len(row) == 3:
            poly.addHoleMarker([float(row[1]), float(row[2])])
        elif len(row) == 4 and dimension == 3:
            poly.addHoleMarker([float(row[1]), float(row[2]), float(row[3])])
        else:
            raise Exception("Poly file seams corrupt: hole section line (3):" +
                            row + " : " + str(i) + " " + str(len(row)))

    if (3 + nVerts + nSegments + nHoles) < len(content):
        # Region section
        row = content[3 + nVerts + nSegments + nHoles].split()

        if len(row) != 1:
            raise Exception("Format unknown for region section " + row)

        nRegions = int(row[0])

        for i in range(nRegions):
            row = content[4 + nVerts + nSegments + nHoles + i].split()
            if len(row) == 5:
                poly.addRegionMarker([float(row[1]), float(row[2])],
                                     marker=int(float(row[3])),
                                     area=float(row[4]))
            elif len(row) == 6 and dimension == 3:
                poly.addRegionMarker([float(row[1]), float(row[2]),
                                      float(row[3])],
                                     marker=int(float(row[4])),
                                     area=float(row[5]))
            else:
                raise Exception("Poly file seams corrupt: region section " +
                                "line (5): " + str(i) + " " + str(len(row)))

    return poly


def exportPLC(poly, fname, **kwargs):
    r"""General writer to save piece-wise linear complex (PLC) as poly file.

    Choose from poly.dimension() and forward appropriate to
    :gimliapi:`GIMLI::Mesh::exportAsTetgenPolyFile`
    and :py:mod:`pygimli.meshtools.writeTrianglePoly`

    Parameters
    ----------
    poly : :gimliapi:`GIMLI::Mesh`
        The polygon to be written.
    fname : string
        Filename of the file to write (\\*.n, \\*.e).

    Examples
    --------
    >>> import pygimli as pg
    >>> import tempfile, os
    >>> fname = tempfile.mktemp() # Create temporary string for filename.
    >>> world2d = pg.meshtools.createWorld(start=[-10, 0], end=[20, 0])
    >>> pg.meshtools.writePLC(world2d, fname)
    >>> read2d = pg.meshtools.readPLC(fname)
    >>> print(read2d)
    Mesh: Nodes: 4 Cells: 0 Boundaries: 4
    >>> world3d = pg.createGrid([0, 1], [0, 1], [-1, 0])
    >>> pg.meshtools.writePLC(world3d, fname)
    >>> os.remove(fname)
    """
    if poly.dimension() == 2:
        exportTrianglePoly(poly, fname, **kwargs)
    else:
        exportTetgenPoly(poly, fname, **kwargs)


def writePLC(*args, **kwargs):
    """
    Backward compatibility.
    Please use :py:mod:`pygimli.meshtools.exportPLC`.
    """
    return exportPLC(*args, **kwargs)


def exportTrianglePoly(poly, fname, float_format='.15e', **kwargs):
    r"""Write :term:`Triangle` poly.

    Write :term:`Triangle` :cite:`Shewchuk96b` ASCII file.
    See: ://www.cs.cmu.edu/~quake/triangle.html

    Parameters
    ----------
    poly : :gimliapi:`GIMLI::Mesh`
        mesh PLC holding nodes, edges, holes & regions
    fname : string
        Filename of the file to read (\\*.n, \\*.e)
    float_format : string
        format string for floats according to str.format()

    verbose : boolean [False]
        Be verbose during import.
    """
    if float_format[0] != '{':
        pfmt = '{:' + float_format + '}'
    else:
        pfmt = float_format
    with open(fname, 'w') as fid:
        fid.write('{:d}\t2\t0\t1\n'.format(poly.nodeCount()))
        nm = poly.nodeMarkers()
        bm = poly.boundaryMarkers()

        fmt = '{:d}' + ('\t' + pfmt) * 2 + '\t{:d}\n'
        for i, p in enumerate(poly.positions()):
            fid.write(fmt.format(i, p.x(), p.y(), nm[i]))
        fid.write('{:d}\t1\n'.format(poly.boundaryCount()))

        for i, b in enumerate(poly.boundaries()):
            fid.write('{:d}\t{:d}\t{:d}\t{:d}\n'.format(i, b.node(0).id(),
                                                        b.node(1).id(), bm[i]))
        fid.write('{:d}\n'.format(len(poly.holeMarker())))

        fmt = '{:d}' + ('\t' + pfmt) * 2 + '\n'
        for i, h in enumerate(poly.holeMarker()):
            fid.write(fmt.format(i, h.x(), h.y()))
        fid.write('{:d}\n'.format(len(poly.regionMarker())))

        fmt = '{:d}' + ('\t' + pfmt) * 3 + '\t{:.15e}\n'
        for i, r in enumerate(poly.regionMarker()):
            fid.write(fmt.format(i, r.x(), r.y(), r.marker(), r.area()))

    return


def writeTrianglePoly(*args, **kwargs):
    """ Backward compatibility.
    Please use :py:mod:`pygimli.meshtools.exportTrianglePoly`.
    """
    return exportTrianglePoly(*args, **kwargs)


def exportTetgenPoly(poly, filename, float_format='.12e', **kwargs):
    """
    Writes a given piecewise linear complex (mesh/poly ) into a Ascii file in
    :term:`Tetgen` .poly format.

    Parameters
    ----------

    filename: string
        Name in which the result will be written. The recommended file
        ending is '.poly'.

    poly: :gimliapi:`GIMLI::Mesh`
        Piecewise linear complex as :gimliapi:`GIMLI::Mesh` to be exported.

    float_format: format string ('.12e')
        Format that will be used to write float values in the Ascii file.
        Default is the exponential float form with a precision of 12 digits.

    """
    if filename[-5:] != '.poly':
        filename = filename + '.poly'
    polytxt = ''
    sep = '\t'  # standard tab seperated file
    assert poly.dim() == 3, 'Exit, only for 3D meshes.'
    boundary_marker = 1
    attribute_count = 0

    # Part 1/4: node list
    # intro line
    # <nodecount> <dimension (3)> <# of attributes> <boundary markers (0 or 1)>
    polytxt += '{0}{5}{1}{5}{2}{5}{3}{4}'.format(poly.nodeCount(), 3,
                                                 attribute_count,
                                                 boundary_marker,
                                                 os.linesep, sep)
    # loop over positions, attributes and marker(node)
    # <point idx> <x> <y> <z> [attributes] [boundary marker]
    point_str = '{:d}'  # index of the point
    for i in range(3):
        # coords as float with given precision
        point_str += sep + '{:%s}' % (float_format)
    point_str += sep + '{:d}' + os.linesep  # node marker
    for j, node in enumerate(poly.nodes()):
        fill = [node.id()]
        fill.extend([pos for pos in node.pos()])
        fill.append(node.marker())
        polytxt += point_str.format(*fill)

    # Part 2/4: boundary list
    # intro line
    # <# of facets> <boundary markers (0 or 1)>
    nBoundaries = poly.boundaryCount()
    # look for extra boundaries present in either the PLC or in kwargs
    extraBoundaries = []
    if 'extraBoundaries' in kwargs:
        extraBoundaries += kwargs.pop('extraBoundaries', [])

    if hasattr(poly, 'extraBoundaries'):
        extraBoundaries += poly.extraBoundaries

    if len(extraBoundaries) > 0:
        print("Detected ", len(extraBoundaries), " extra boundaries!")

    nBoundaries += len(extraBoundaries)
    polytxt += '{0:d}{2}1{1}'.format(nBoundaries, os.linesep, sep)
    # loop over facets, each facet can contain an arbitrary number of holes
    # and polygons, in our case, there is always one polygon per facet.
    for bound in poly.boundaries():
        # one line per facet
        # <# of polygons> [# of holes] [boundary marker]
        npolys = 1
        polytxt += '1{2}0{2}{0:d}{1}'.format(bound.marker(), os.linesep, sep)
        # inner loop over polygons
        # <# of corners> <corner 1> <corner 2> ... <corner #>
        for l in range(npolys):
            poly_str = '{:d}'.format(bound.nodeCount())
            for ind in bound.ids():
                poly_str += sep + '{:d}'.format(ind)

            polytxt += '{0}{1}'.format(poly_str, os.linesep)
        # inner loop over holes
        # not necessary yet ?! why is there an extra hole section?
        # because this is for 2D holes in facets only

    # part 2b: extra boundaries that cannot be part of mesh class
    for nodes in extraBoundaries:
        # <# of polygons> [# of holes] [boundary marker]
        npolys = 1
        polytxt += '1{2}0{2}{0:d}{1}'.format(111, os.linesep, sep)
        # <# of corners> <corner 1> <corner 2> ... <corner #>
        poly_str = '{:d}'.format(len(nodes))
        for ind in nodes:
            poly_str += sep + '{:d}'.format(ind)

        polytxt += '{0}{1}'.format(poly_str, os.linesep)

    # part 3/4: hole list
    # intro line
    # <# of holes>
    holes = poly.holeMarker()
    polytxt += '{:d}{}'.format(len(holes), os.linesep)
    # loop over hole markers
    # <hole #> <x> <y> <z>
    hole_str = '{:d}'
    for m in range(3):
        hole_str += sep + '{:%s}' % float_format

    hole_str += os.linesep
    for n, hole in enumerate(holes):
        polytxt += hole_str.format(n, *hole)

    # part 4/4: region attributes and volume constraints (optional)
    # intro line
    # <# of regions>
    regions = poly.regionMarker()
    polytxt += '{:d}{}'.format(len(regions), os.linesep)
    # loop over region markers
    # <region #> <x> <y> <z> <region number> <region attribute>
    region_str = '{:d}'
    for o in range(3):
        region_str += sep + '{:%s}' % (float_format)

    region_str += sep + '{:d}%s{:%s}' % (sep, float_format) + os.linesep
    for p, region in enumerate(regions):
        polytxt += region_str.format(p, region.x(), region.y(), region.z(),
                                     region.marker(),
                                     region.area())

    # writing file
    with open(filename, 'w') as out:
        out.write(polytxt)


def tetgen(filename, quality=1.2, preserveBoundary=False, verbose=False):
    """Create a mesh with :term:`Tetgen` from file.

    Create a :term:`Tetgen` :cite:`Si2004` mesh from a PLC.

    Forwards to system call tetgen, which must be known to your system.

    Parameters
    ----------
    filename: str

    quality: float [1.2]
        Refines mesh (to improve mesh quality). [1.1 ... ]

    preserveBoundary: bool [False]
        Preserve PLC boundary mesh

    verbose: bool [False]
        be verbose

    Returns
    -------
    mesh : :gimliapi:`GIMLI::Mesh`
    """
    filebody = filename.replace('.poly', '')
    syscal = 'tetgen -pazAC'
    syscal += 'q' + str(quality)

    if not verbose:
        syscal += 'Q'
    else:
        syscal += 'V'

    if preserveBoundary:
        syscal += 'Y'

    syscal += ' ' + filebody + '.poly'

    if verbose:
        print(syscal)

    system(syscal)
    system('meshconvert -it -BD -o ' + filebody + ' ' + filebody + '.1')
    try:
        os.remove(filebody + '.1.node')
        os.remove(filebody + '.1.ele')
        os.remove(filebody + '.1.face')
    except BaseException as e:
        print(e)
    mesh = pg.Mesh(filebody)
    return mesh


def polyAddVIP(filename, pos, marker=0, isRegionMarker=False,
               isHoleMarker=False, maxCellSize=0, verbose=False):
    """Add very important point (VIP) to a PLC file.

    Add very important point (VIP) to a PLC.
    Out of core wrapper for dcfemlib::polytools::polyAddVIP.

    If you wan add these points to a plc directly use
    :gimliapi:`GIMLI::Mesh::createNode`,
    :gimliapi:`GIMLI::Mesh::addRegionMarker` or
    :gimliapi:`GIMLI::Mesh::addHoleMarker`.

    Parameters
    ----------

    Returns
    -------
    """
    raise BaseException('obsolete use mesh methods directly')
#     syscal = "polyAddVIP -x " + str(pos[0]) + \
#         " -y " + str(pos[1]) + \
#         " -z " + str(pos[2])
#
#     if isHoleMarker:
#         syscal += " -H "
#     else:
#         syscal += " -m " + str(marker)
#
#     if isRegionMarker:
#         syscal += " -R "
#
#     if maxCellSize > 0:
#         syscal += " -a " + str(maxCellSize)
#
#     syscal += " " + filename
#
#     if verbose:
#         print(syscal)
#     system(syscal)
# # def polyAddVIP


def polyAddRectangle(filename, rect, marker=0, depth=0, clean=True):
    """Add horizontal plane to a PLC.

    Add horizontal plane to a PLC.
    Out of core wrapper for dcfemlib::polytools::polytools.
    Merge a meshed horizontal Rectangle with given marker[0] to
    a 3D PLC at a given depth [0] clean removes all out of core files

    Parameters
    ----------

    Returns
    -------
    """
    raise BaseException('obsolete use mesh methods directly')
#     rect.writeXY("__pad.xy", close=True)
#     system("polyCreateWorld -d2 -t __pad.xy -C __pad")
#     a = rect.area() / 29.0
#
#     system("dctriangle -a " + str(a) + " -q34.0 -S __pad")
#     system("polyCreateFacet -o __pad3d -m " + str(marker) + " __pad.bms")
#     system("polyTranslate -z " + str(depth) + " __pad3d")
#
#     # add node to the center of the rectangle
# #    system("polyAddVIP -x " + str(rect.start[0])
# #                    + " -y " + str(rect.start[1])
# #                    + " -m -1 __pad3d")
#     system("polyMerge " + filename + " __pad3d " + filename)
#
# #    system("polyCreateFacet -o __pad3d -m 1 __pad.bms")
# #    system("polyTranslate -z -0.1 __pad3d")
# #    system("polyMerge " + filename + " __pad3d " + filename)
#
#     if clean:
#         os.remove('__pad.xy')
#         os.remove('__pad.poly')
#         os.remove('__pad.bms')
#         os.remove('__pad3d.poly')
# def polyAddRectangle


def polyCreateWorld(filename, x=None, depth=None, y=None, marker=0,
                    maxCellSize=0, verbose=True):
    """Create the PLC of a default world.

    Create the PLC of a default world.
    Out of core wrapper for dcfemlib::polytools::polyCreateWorld

    Parameters
    ----------

    Returns
    -------
    """
    raise BaseException('obsolete use mesh methods directly')
    if depth is None:
        print("Please specify worlds depth.")
        return

    if x is None:
        print("Please specify worlds x dimension.")
        return

    dimension = 3
    z = depth

    if y is None:
        dimension = 2

    syscal = 'polyCreateWorld -d ' + str(dimension) \
        + ' -x ' + str(x) \
        + ' -y ' + str(y) \
        + ' -z ' + str(z) \
        + ' -m ' + str(marker) \

    if maxCellSize > 0:
        syscal += " -a " + str(maxCellSize)

    syscal = syscal + ' ' + filename

    if verbose:
        print(syscal)

    os.system(syscal)


def polyTranslate(filename, x=0.0, y=0.0, z=0.0):
    """Translate (move) a PLC.

    Out of core wrapper for dcfemlib::polytools.
    Spatial translate (move) the PLC (filename) by x, y and z

    Parameters
    ----------

    Returns
    -------
    """
    raise BaseException('obsolete use mesh methods directly')
    # system("polyTranslate " +
    #        " -x " + str(x) +
    #        " -y " + str(y) +
    #        " -z " + str(z) + " " + filename)


if __name__ == "__main__":
    pass
