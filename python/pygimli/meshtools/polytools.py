# -*- coding: utf-8 -*-
""" Tools to create or manage PLC """

import os
from os import system

import math
import numpy as np
import pygimli as pg

def polyCreateDefaultEdges_(poly, boundaryMarker=1, isClosed=True, **kwargs):
    """ INTERNAL """

    nEdges = poly.nodeCount()-1 + isClosed
    bm = None
    if hasattr(boundaryMarker, '__len__'):
        if len(boundaryMarker) == nEdges:
            bm = boundaryMarker
        else:
            raise Exception("marker length != nEdges",
                            len(boundaryMarker), nEdges)
    else:
        bm = [boundaryMarker] * nEdges

    for i in range(poly.nodeCount() - 1):
        poly.createEdge(poly.node(i), poly.node(i+1), bm[i])
    if isClosed:
        poly.createEdge(poly.node(poly.nodeCount()-1),
                        poly.node(0), bm[-1])

def createRectangle(start=None, end=None, pos=None, size=None, **kwargs):
    """
    Create rectangle polygon.

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
            Maximum cell size for the resulting triangle cells after mesh generation
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
    poly : gimliapi:`GIMLI::Mesh`
        The resulting polygon is a gimliapi:`GIMLI::Mesh`.

    Examples
    --------
    >>> from pygimli.meshtools import createRectangle
    >>> from pygimli.mplviewer import drawMesh
    >>> import matplotlib.pyplot as plt
    >>> rectangle = createRectangle(start=[4, -4], end=[6, -6], marker=4, area=0.1)
    >>>
    >>> fig, ax = plt.subplots()
    >>> drawMesh(ax, rectangle)
    >>> plt.show()
    """

    if not ((start and end) or (pos and size)):
        raise BaseException("createRectangle pls. give either start and end"
                            "OR pos and size." )
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
        poly.addHoleMarker(sPos + (ePos-sPos)*0.2)
    else:
        poly.addRegionMarker(sPos + (ePos-sPos)*0.2,
                             marker=kwargs.pop('marker', 1),
                             area=kwargs.pop('area', 0))

    if size is not None:
        poly.scale(size)
    if pos is not None:
        poly.translate(pos)



    polyCreateDefaultEdges_(poly, **kwargs)

    return poly

def createWorld(start, end, marker=1, area=0, layers=None):
    """
    Create simple rectangular world.

    Create simple rectangular world with appropriate boundary conditions.
    Surface boundary is set do pg.MARKER_BOUND_HOMOGEN_NEUMANN, i.e, -1 and
    inner surface is set to pg.MARKER_BOUND_MIXED, i.e., -2

    Parameters
    ----------
    start : [x, y]
        Upper/Left Corner
    end : [x, y]
        Lower/Right Corner
    marker : int
        Marker for the resulting triangle cells after mesh generation
    area : float
        Maximum cell size for the resulting triangle cells after mesh generation
    layers : [float]
        Add some layers to the world.

    Returns
    -------
    poly : gimliapi:`GIMLI::Mesh`
        The resulting polygon is a gimliapi:`GIMLI::Mesh`.

    Examples
    --------
    >>> from pygimli.meshtools import createWorld
    >>> from pygimli.mplviewer import drawMesh
    >>> import matplotlib.pyplot as plt
    >>> world = createWorld(start=[-5, 0], end=[5, -5], layers=[1,3,2,5])
    >>>
    >>> fig, ax = plt.subplots()
    >>> drawMesh(ax, world)
    >>> plt.show()
    """

    z = [start[1]]

    if layers is not None:
        z = z + layers

    z.append(end[1])
    rs = []
    poly = pg.Mesh(2)

    for i, depth in enumerate(z):
        n = poly.createNode([start[0], depth])
        if i > 0:
            poly.addRegionMarker(n.pos() + [0.2, 0.2], marker=i, area=area)

    for i, depth in enumerate(z[::-1]):
        poly.createNode([end[0], depth])

    polyCreateDefaultEdges_(poly,
                            boundaryMarker=[1]*(len(z)-1) + [3] + [2]*(len(z)-1) + [4])

    if layers is not None:
        for i in range(len(layers)):
            poly.createEdge(poly.node(i + 1),
                            poly.node(poly.nodeCount() - i - 2), 4 + i)

    return poly

def createCircle(pos=None, radius=1, segments=12, start=0, end=2.*math.pi,
                 **kwargs):

    """
    Create simple circle polygon.

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
            Maximum cell size for the resulting triangle cells after mesh generation
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
    poly : gimliapi:`GIMLI::Mesh`
        The resulting polygon is a gimliapi:`GIMLI::Mesh`.

    Examples
    --------
    >>> import matplotlib.pyplot as plt
    >>> from pygimli.mplviewer import drawMesh
    >>> import pygimli as pg
    >>> import math
    >>> from pygimli.meshtools import polytools as plc
    >>> c0 = plc.createCircle(pos=(-5.0, 0.0), radius=2, segments=6)
    >>> c1 = plc.createCircle(pos=(0.0, 0.0), segments=5, start=0, end=math.pi)
    >>> c2 = plc.createCircle(pos=(5.0, 0.0), segments=3, start=math.pi,
    ...                       end=1.5*math.pi, isClosed=False)
    >>> plc = plc.mergePLC([c0, c1, c2])
    >>>
    >>> fig, ax = plt.subplots()
    >>> drawMesh(ax, plc)
    >>> plt.show()
    """

    if pos == None:
        pos = [0.0, 0.0]

    poly = pg.Mesh(2)

    dPhi = (end - start) / (segments)
    nPhi = segments +1

    if abs((end%(2.*math.pi) - start)) < 1e-6:
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
        poly.addRegionMarker([0.0, 0.0],
                             marker=kwargs.pop('marker', 1),
                             area=kwargs.pop('area', 0))

    if hasattr(radius, '__len__'):
        poly.scale(radius)
    else:
        poly.scale([radius, radius])

    poly.translate(pos)

    polyCreateDefaultEdges_(poly, **kwargs)
    return poly

def createLine(start, end, segments, **kwargs):
    """
    Create simple line polygon.

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
    poly : gimliapi:`GIMLI::Mesh`
        The resulting polygon is a gimliapi:`GIMLI::Mesh`.

    Examples
    --------
    """
    poly = pg.Mesh(2)
    startPos = pg.RVector3(start)
    endPos = pg.RVector3(end)
    a = endPos - startPos

    dt = 1 / segments
    for i in range(0, segments + 1):
        if kwargs.pop('leftDirection', True):
            p = startPos + a * (dt * i)
        else:
            p = endPos -  a * (dt * i)

        poly.createNode(p)

    polyCreateDefaultEdges_(poly, isClosed=False, **kwargs)
    return poly

def mergePLC(pols):
    """
    Merge multiply polygons into a single polygon

    Common nodes and common edges will be checked and removed.

    Crossing or touching edges or Node/Edge intersections will NOT be
    recognized yet. -> TODO

    Parameters
    ----------
    pols: [gimliapi:`GIMLI::Mesh`]
        List of polygons that need to be merged

    Returns
    -------
    poly : gimliapi:`GIMLI::Mesh`
        The resulting polygon is a gimliapi:`GIMLI::Mesh`.

    Examples
    --------
    >>> from pygimli.meshtools import polytools as plc
    >>> from pygimli.meshtools import createMesh
    >>> from pygimli.mplviewer import drawMesh
    >>> import matplotlib.pyplot as plt
    >>> world = plc.createWorld(start=[-10, 0], end=[10, -10], marker=1)
    >>> c1 = plc.createCircle([-1, -4], radius=1.5, area=0.1, marker=2, segments=4)
    >>> c2 = plc.createCircle([-6, -5], radius=[1.5, 3.5], isHole=1)
    >>> r1 = plc.createRectangle(pos=[3, -5], size=[2, 2], marker=3)
    >>> r2 = plc.createRectangle(start=[4, -4], end=[6, -6], marker=4, area=0.1)
    >>> plc = plc.mergePLC([world, c1, c2, r1, r2])
    >>>
    >>> fig, ax = plt.subplots()
    >>> drawMesh(ax, plc)
    >>> drawMesh(ax, createMesh(plc))
    >>> plt.show()
    """
    poly = pg.Mesh(2)

    for p in pols:
        nodes = [poly.createNodeWithCheck(n.pos()) for n in p.nodes()]

        for e in p.boundaries():
            poly.createEdge(nodes[e.node(0).id()],
                            nodes[e.node(1).id()],
                            e.marker())

        if len(p.regionMarker()) > 0:
            for rm in p.regionMarker():
                poly.addRegionMarker(rm)

        if len(p.holeMarker()) > 0:
            for hm in p.holeMarker():
                poly.addHoleMarker(hm)

    return poly

def readPLC(filename):
    """
    Read 2D triangle POLY or 3D Tetgen PLC files.
    
    TODO 3D Tetgen PLC
     
    Parameters
    ----------
    filename: str
        Filename *.poly

    Returns
    -------
    poly : gimliapi:`GIMLI::Mesh`
        The resulting polygon is a gimliapi:`GIMLI::Mesh`.
        
    """
    with open(filename, 'r') as fi:
        content = fi.readlines()
    fi.close()
    
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
                n = poly.createNode((float(row[1]), float(row[2]), float(row[3])))
            if haveNodeMarker:
                n.setMarker(int(row[-1]))

        else:
            raise Exception("Poly file seams corrupt: node section line: "
                            + str(i) + " " + row)
      
    # Segment section
    row = content[1 + nVerts].split()

    if len(row) != 2:
        raise Exception("Format unknown for segment section " + row)
    
    nSegments = int(row[0])
    haveBoundaryMarker = int(row[1]);

    
    if dimension == 2:
        for i in range(nSegments):
            row = content[2 + nVerts + i].split()
    
            if len(row) == (3 + haveBoundaryMarker):
                marker = 0
                if haveBoundaryMarker:
                    marker = int(row[3])
                
                poly.createEdge(poly.node(int(row[1]) - fromOne),
                                poly.node(int(row[2]) - fromOne),
                                marker)
    else:
        raise Exception("Read segments for 3D tetgen format not yet supported")
    
    # Hole section
    row = content[2 + nVerts + nSegments].split()
    
    if len(row) != 1:
        raise Exception("Format unknown for hole section " + row)
    
    nHoles = int(row[0])
    for i in range(nHoles):
        row = content[3 + nVerts + nVerts].split()

        if len(row) == 3:
            poly.addHoleMarker([float(row[1]), float(row[2])])
        else:
            raise Exception("Poly file seams corrupt: hole section line (3): "
                            + str(i) + " " + str(len(row)))
      
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
                                     marker=int(float(row[3])), area=float(row[4]))
            else:
                raise Exception("Poly file seams corrupt: region section line (5): "
                                + str(i) + " " + str(len(row)))
    
    return poly
    

def tetgen(filename, quality=1.2, preserveBoundary=False, verbose=False):
    """
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
    mesh : gimliapi:`GIMLI::Mesh`
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
    except:
        None
    mesh = pg.Mesh(filebody)
    return mesh


def polyAddVIP(filename, pos, marker=0, isRegionMarker=False,
               isHoleMarker=False, maxCellSize=0, verbose=False):
    """
    Add very important point (VIP) to a PLC.

    Out of core wrapper for dcfemlib::polytools::polyAddVIP.

    Parameters
    ----------

    Returns
    -------
    """

    syscal = "polyAddVIP -x " + str(pos[0]) + \
        " -y " + str(pos[1]) + \
        " -z " + str(pos[2])

    if isHoleMarker:
        syscal += " -H "
    else:
        syscal += " -m " + str(marker)

    if isRegionMarker:
        syscal += " -R "

    if maxCellSize > 0:
        syscal += " -a " + str(maxCellSize)

    syscal += " " + filename

    if verbose:
        print(syscal)
    system(syscal)
# def polyAddVIP


def polyAddRectangle(filename, rect, marker=0, depth=0, clean=True):
    """
    Add horizontal plane to a PLC

    Out of core wrapper for dcfemlib::polytools::polytools.
    Merge a meshed horizontal Rectangle with given marker[0] to
    a 3D PLC at a given depth [0] clean removes all out of core files

    Parameters
    ----------

    Returns
    -------
    """

    rect.writeXY("__pad.xy", close=True)
    system("polyCreateWorld -d2 -t __pad.xy -C __pad")
    a = rect.area() / 29.0

    system("dctriangle -a " + str(a) + " -q34.0 -S __pad")
    system("polyCreateFacet -o __pad3d -m " + str(marker) + " __pad.bms")
    system("polyTranslate -z " + str(depth) + " __pad3d")

    # add node to the center of the rectangle
#    system("polyAddVIP -x " + str(rect.start[0])
#                    + " -y " + str(rect.start[1])
#                    + " -m -1 __pad3d")
    system("polyMerge " + filename + " __pad3d " + filename)

#    system("polyCreateFacet -o __pad3d -m 1 __pad.bms")
#    system("polyTranslate -z -0.1 __pad3d")
#    system("polyMerge " + filename + " __pad3d " + filename)

    if clean:
        os.remove('__pad.xy')
        os.remove('__pad.poly')
        os.remove('__pad.bms')
        os.remove('__pad3d.poly')
# def polyAddRectangle


def polyCreateWorld(filename, x=None, depth=None, y=None, marker=0,
                    maxCellSize=0, verbose=True):
    """
    Create the PLC of a default world.

    Out of core wrapper for dcfemlib::polytools::polyCreateWorld

    Parameters
    ----------

    Returns
    -------
    """

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


def polyTranslate(filename, x=0.0, y=0.0, z=0.0, verbose=True):
    """
    Translate (move) a PLC

    Out of core wrapper for dcfemlib::polytools.
    Spatial translate (move) the PLC (filename) by x, y and z

    Parameters
    ----------

    Returns
    -------
    """
    system("polyTranslate " +
           " -x " + str(x) +
           " -y " + str(y) +
           " -z " + str(z) + " " + filename)

if __name__ == "__main__":
    #from pygimli.meshtools import polytools as plc
    #from pygimli.meshtools import createMesh
    #from pygimli.mplviewer import drawMesh
    #import matplotlib.pyplot as plt
    #world = plc.createWorld(start=[-10, 0], end=[10, -10], marker=1)
    #c1 = plc.createCircle([-1, -4], radius=1.5, area=0.1, marker=2, segments=4)
    #c2 = plc.createCircle([-6, -5], radius=[1.5, 3.5], isHole=1)
    #r1 = plc.createRectangle(pos=[3, -5], size=[2, 2], marker=3)
    #r2 = plc.createRectangle(start=[4, -4], end=[6, -6], marker=4, area=0.1)
    #plc = plc.mergePLC([world, c1, c2, r1, r2])

    #fig, ax = plt.subplots()
    #drawMesh(ax, plc)
    #drawMesh(ax, createMesh(plc))
    #plt.show()

    #import matplotlib.pyplot as plt
    #import pygimli as pg
    #import math
    #from pygimli.meshtools import polytools as plc
    #c0 = plc.createCircle(pos=(-5.0, 0.0), radius=2, segments=6)
    #c1 = plc.createCircle(pos=(0.0, 0.0), segments=5, start=0, end=math.pi)
    #c2 = plc.createCircle(pos=(5.0, 0.0), segments=3, start=math.pi,
    #                      end=1.5*math.pi, isClosed=False)
    #fig, ax = plt.subplots()
    #pg.show([c0, c1, c2], axes=ax)
    #plt.show()
    pass
