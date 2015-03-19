# -*- coding: utf-8 -*-
""" Tools to create or manage PLC """

import os
from os import system

import numpy as np
import pygimli as pg

def polyCreateDefaultEdges_(poly, marker=1, closed=True):
    """ INTERNAL """
    
    nEdges = poly.nodeCount()-1 + closed
    bm = [marker]*nEdges
    if hasattr(marker, '__len__'):
        if len(marker) == nEdges:
            bm = marker
            
    for i in range(poly.nodeCount() - 1):
        poly.createEdge(poly.node(i), poly.node(i+1), bm[i])
    if closed:
        poly.createEdge(poly.node(poly.nodeCount()-1), 
                        poly.node(0), bm[-1])
        
def createRectangle(start=None, end=None, pos=None, size=None, marker=1, area=0, 
                    boundaryMarker=1, isHole=False):
    """
    Create rectangle polygon.
        
    Create rectangle with start position and a given size.
        
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
    marker : int
        Marker for the resulting triangle cells after mesh generation
    area : float
        Maximum cell size for the resulting triangle cells after mesh generation
    boundaryMarker : int
        Marker for the resulting boundary edges
    isHole :
        The Polygone will become a hole instead of a triangulation

    Returns
    -------
    poly : gimliapi:`GIMLI::Mesh`
        The resulting polygon is a gimliapi:`GIMLI::Mesh`.
        
    Examples
    --------
    TODO
    """
    
    if start is None:
        start = [-0.5, 0.5]
    if end is None:
        end = [0.5, -0.5]
    poly = pg.Mesh(2)
    
    poly.createNode(start)
    poly.createNode([start[0], end[1]])
    poly.createNode(end)
    poly.createNode([end[0], start[1]])
    
    if size is not None:
        poly.scale(size)
    if pos is not None:
        poly.translate(pos)
    
    if isHole:
        poly.addHoleMarker(poly.nodes()[0].pos() + [0.001, -0.001])
    else:
        poly.addRegionMarker(poly.nodes()[0].pos() + [0.001, -0.001], 
                             marker=marker, area=area)
    polyCreateDefaultEdges_(poly, marker=boundaryMarker, closed=True)
    return poly

def createWorld(start, end, marker=1, area=0):
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
        
    Returns
    -------
    poly : gimliapi:`GIMLI::Mesh`
        The resulting polygon is a gimliapi:`GIMLI::Mesh`.
        
    Examples
    --------
    TODO
    """
    return createRectangle(start, end, marker=marker, area=area,
                           boundaryMarker = [pg.MARKER_BOUND_MIXED,
                                             pg.MARKER_BOUND_MIXED,
                                             pg.MARKER_BOUND_MIXED,
                                             pg.MARKER_BOUND_HOMOGEN_NEUMANN])
    #poly = pg.Mesh(2)    
    
    #poly.createNode(start)
    #poly.createNode([start[0], end[1]])
    #poly.createNode(end)
    #poly.createNode([end[0], start[1]])

    #poly.addRegionMarker(poly.nodes()[0].pos() + [0.001, -0.001], 
                         #marker=marker, area=area)

    #for i in range(poly.nodeCount() - 1):
        #poly.createEdge(poly.node(i), poly.node(i+1), 
                        #pg.MARKER_BOUND_HOMOGEN_NEUMANN)

    #poly.createEdge(poly.node(poly.nodeCount()-1), 
                    #poly.node(0), pg.MARKER_BOUND_MIXED)
    #return poly

def createCircle(pos, radius, segments=12, marker=1, area=0,
                 boundaryMarker=1, leftDirection=True, isHole=False):
    """
    Create simple circle polygon.
        
    Parameters
    ----------
    pos : [x, y]
        Center position
    radius : float | [a,b]
        radius or halfaxes of the circle
    segments : int
        Discrete amount of segmens for the circle
    marker : int
        Marker for the resulting triangle cells after mesh generation
    area : float
        Maximum cell size for the resulting triangle cells after mesh generation
    boundaryMarker : int
        Marker for the resulting boundary edges
    leftDirection : bool
        Rotational direction
    isHole :
        The Polygone will become a hole instead of a triangulation
        
    Returns
    -------
    poly : gimliapi:`GIMLI::Mesh`
        The resulting polygon is a gimliapi:`GIMLI::Mesh`.
    
    Examples
    --------
    TODO
    """
    poly = pg.Mesh(2)    
    poly.createNode([0, 1.0])
    for i in range(1, segments):
        if leftDirection:
            xp = np.sin(-i * (2. * np.pi) / segments)
        else:
            xp = np.sin(i * (2. * np.pi) / segments)

        yp = np.cos(i * (2. * np.pi) / segments)
        poly.createNode([xp, yp])
    
    if hasattr(radius, '__len__'):
        poly.scale(radius)
    else:
        poly.scale([radius, radius])
    poly.translate(pos)
    if isHole:
        poly.addHoleMarker(poly.nodes()[0].pos() + [0.0, -0.001])
    else:
        poly.addRegionMarker(poly.nodes()[0].pos() + [0.0, -0.001], 
                             marker=marker, area=area)
    
    polyCreateDefaultEdges_(poly, boundaryMarker, closed=True)
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
    >>> import pygimli.polytools as plc
    >>> from pygimli.meshtools import createMesh
    >>> from pygimli.mplviewer import drawMesh
    >>> import matplotlib.pyplot as plt
    >>> world = plc.createWorld(start=[-10, 0], end=[10, -10], marker=1)
    >>> c1 = plc.createCircle([-1, -4], radius=1.5, area=0.1, marker=2)
    >>> c2 = plc.createCircle([-6, -5], radius=[1.5, 3.5], isHole=1)
    >>> r1 = plc.createRectangle([3, -5], size=[2, 2], marker=3)
    >>> r2 = plc.createRectangle([5, -5], size=[2, 2], marker=4, area=0.1)
    >>> plc = plc.mergePLC([world, c1, c2, r1, r2])
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
            poly.addRegionMarker(p.regionMarker()[0])
    
        if len(p.holeMarker()) > 0:
            poly.addHoleMarker(p.holeMarker()[0])
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
    from pygimli.polytools import *
    from pygimli.meshtools import createMesh
    from pygimli.mplviewer import drawMesh
    import matplotlib.pyplot as plt
    world = createWorldPolygon(start=[-10, 0], end=[10, -10], marker=1)
    c1 = createCirclePolygon([-1, -3], radius=1.5, area=0.1, marker=2)
    c2 = createCirclePolygon([-6, -5], radius=[1.5, 3.5], isHole=1)
    r1 = createRectanglePolygon([3, -5], size=[2, 2], marker=3)
    r2 = createRectanglePolygon([5, -5], size=[2, 2], marker=4, area=0.1)
    plc = mergePolygons([world, c1, c2, r1, r2])
    fig, ax = plt.subplots()
    drawMesh(ax, plc)
    drawMesh(ax, createMesh(plc))
    plt.show()