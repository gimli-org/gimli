# -*- coding: utf-8 -*-

import os
import pygimli as pg
import numpy as np

from pygimli.polytools import *

def createMesh(poly, quality=30, area=0.0,
               smooth=None, switches=None,
               regions=None, holes=None,
               verbose=False):
    """
    Create a mesh for a given PLC mesh.

    Create a mesh for a given PLC using
    :term:`triangle` or :term:`tetgen` if the gimli support for the
    meshgenerator is installed.
    The mesh poly need to be a valid PLC.

    If poly is a list coordinates a simple Delaunay mesh of the convex hull
    will be created.

    Tetgen support need to be implemented

    Parameters
    ----------
    poly: :gimliapi:`GIMLI::Mesh` or list
        * 2D or 3D gimli mesh that contains the PLC.
        * 2D mesh needs edges
        * 3D mesh needs ... to be implemented
        * List of x y pairs [[x0, y0], ... ,[xN, yN]]
        * PLC or list of PLC
    quality: float
        2D triangle quality sets a minimum angle constraint.
        Be careful with values above 34 degrees
    area: float
        2D maximum triangle size in m*m
    smooth: tuple
        [smoothing algorithm, number of iterations]
        0, no smoothing
        1, node center
        2, weighted node center
    switches: str
        force triangle to use the gives command switches

    Returns
    -------
    mesh: :gimliapi:`GIMLI::Mesh`
    """
    
    #poly == [pg.Mesh, ]
    if isinstance(poly, list):
        if isinstance(poly[0], pg.Mesh):
            return createMesh(mergePolygons(poly),
                              quality, area, smooth, switches, verbose)
    #poly == [pos, pos, ]    
    if isinstance(poly, list) or isinstance(poly, type(zip)):
        delPLC = pg.Mesh(2)
        for p in poly:
            delPLC.createNode(p[0], p[1], 0.0)
        return createMesh(delPLC, switches='-zeY')

    #poly == Mesh
    if poly.dim() == 2:
        if poly.nodeCount() == 0:
            raise Exception("No nodes in poly to create a valid mesh")

        tri = pg.TriangleWrapper(poly)

        if switches is None:
            # -D Conforming delaunay
            # -F Uses Steven Fortune's sweepline algorithm
            # no -a here ignores per region area
            switches = 'pzaeA'

            if area > 0:
                switches += 'a' + str(area)

            switches += 'q' + str(quality)

        if not verbose:
            switches += 'Q'

        # print(switches)
        tri.setSwitches(switches)
        mesh = tri.generate()

        if smooth is not None:
            mesh.smooth(nodeMoving=True,
                        edgeSwapping=False,
                        smoothFunction=smooth[0],
                        smoothIteration=smooth[1])
        return mesh

    else:
        raise('not yet implemented')

def readGmsh(fname, verbose=False):
    """
    Read :term:`Gmsh` ASCII file and return instance of GIMLI::Mesh class.

    Parameters
    ----------
    fname : string
        Filename of the file to read (\\*.msh). The file must conform
        to the `MSH ASCII file version 2 <http://geuz.org/gmsh/doc/
        texinfo/gmsh.html#MSH-ASCII-file-format>`_ format.
    verbose : boolean, optional
        Be verbose during import.

    Notes
    -----
    Physical groups specified in Gmsh are interpreted as follows:

    - Points with the physical number 99 are interpreted as sensors.
    - Physical Lines and Surfaces define boundaries in 2D and 3D, respectively.
        - Physical Number 1: homogeneous Neumann condition
        - Physical Number 2: mixed boundary condition
        - Physical Number 3: homogeneous Dirichlet condition
        - Physical Number 4: Dirichlet condition
    - Physical Surfaces and Volumes define regions in 2D and 3D, respectively.
        - Physical Number 1: No inversion region
        - Physical Number >= 2: Inversion region

    """
    inNodes, inElements, ncount, ecount = 0, 0, 0, 0
    fid = open(fname)
    if verbose:
        print('Reading %s... \n' % fname)

    for line in fid:

        if line[0] == '$':
            if line.find('Nodes') > 0:
                inNodes = 1
            if line.find('EndNodes') > 0:
                inNodes = 0
            if line.find('Elements') > 0:
                inElements = 1
            if line.find('EndElements') > 0:
                inElements = 0

        else:
            if inNodes == 1:
                if len(line.split()) == 1:
                    nodes = np.zeros((int(line), 3))
                    if verbose:
                        print('  Nodes: %s' % int(line))
                else:
                    nodes[ncount, :] = np.array(line.split(), 'float')[1:]
                    ncount += 1

            elif inElements == 1:
                if len(line.split()) == 1:
                    if verbose:
                        print('  Entries: %s' % int(line))
                    points, lines, triangles, tets = [], [], [], []

                else:
                    entry = list(map(int, line.split()))[1:]

                    if entry[0] == 15:
                        points.append((entry[-2], entry[-3]))
                    elif entry[0] == 1:
                        lines.append((entry[-2], entry[-1], entry[2]))
                    elif entry[0] == 2:
                        triangles.append((entry[-3], entry[-2],
                                          entry[-1], entry[2]))
                    elif entry[0] == 4:
                        tets.append((entry[-4], entry[-3], entry[-2],
                                     entry[-1], entry[2]))
    fid.close()

    lines = np.asarray(lines)
    triangles = np.asarray(triangles)
    tets = np.asarray(tets)

    if verbose:
        print('    Points: %s' % len(points))
        print('    Lines: %s' % len(lines))
        print('    Triangles: %s' % len(triangles))
        print('    Tetrahedra: %s \n' % len(tets))
        print('Creating mesh object... \n')

    # check dimension
    if len(tets) == 0:
        dim, bounds, cells = 2, lines, triangles
        zero_dim = np.abs(nodes.sum(0)).argmin()  # identify zero dimension
    else:
        dim, bounds, cells = 3, triangles, tets
    if verbose:
        print('  Dimension: %s-D' % dim)

    # creating instance of GIMLI::Mesh class
    mesh = pg.Mesh(dim)

    # replacing boundary markers (gmsh does not allow negative phys. regions)
    bound_marker = (pg.MARKER_BOUND_HOMOGEN_NEUMANN, pg.MARKER_BOUND_MIXED,
                    pg.MARKER_BOUND_HOMOGEN_DIRICHLET,
                    pg.MARKER_BOUND_DIRICHLET)
    for i in range(4):
        bounds[:, dim][bounds[:, dim] == i + 1] = bound_marker[i]

    # account for CEM markers
    bounds[:, dim][bounds[:, dim] >= 10000] *= -1

    if verbose:
        bound_types = np.unique(bounds[:, dim])
        regions = np.unique(cells[:, dim + 1])
        print('  Regions: %s ' % len(regions) + str(tuple(regions)))
        print('  Boundary types: %s ' % len(bound_types) +
              str(tuple(bound_types)))

    for node in nodes:
        if dim == 2:
            mesh.createNode(node[0], node[3 - zero_dim], 0)
        else:
            mesh.createNode(node)

    for cell in cells:
        if dim == 2:
            mesh.createTriangle(
                mesh.node(int(cell[0] - 1)), mesh.node(int(cell[1] - 1)),
                mesh.node(int(cell[2] - 1)), marker=int(cell[3]))
        else:
            mesh.createTetrahedron(
                mesh.node(int(cell[0] - 1)), mesh.node(int(cell[1] - 1)),
                mesh.node(int(cell[2] - 1)), mesh.node(
                    int(cell[3] - 1)),
                marker=int(cell[4]))

    mesh.createNeighbourInfos()

    for bound in bounds:
        if dim == 2:
            mesh.createEdge(
                mesh.node(int(bound[0] - 1)),
                mesh.node(int(bound[1] - 1)),
                marker=int(bound[2]))
        else:
            mesh.createTriangleFace(
                mesh.node(int(bound[0] - 1)), mesh.node(int(bound[1] - 1)),
                mesh.node(int(bound[2] - 1)), marker=int(bound[3]))

    # assign marker to corresponding nodes (sensors, reference nodes, etc.)
    if len(points) > 0:
        for point in points:
            mesh.node(point[0] - 1).setMarker(-point[1])

    if verbose:
        if len(points) > 0:
            points = np.asarray(points)
            node_types = np.unique(points[:, 1])
            print('  Marked nodes: %s ' % len(points) + str(tuple(node_types)))
        print('\nDone. \n')
        print('  ' + str(mesh))

    return mesh


def readTriangle(fname, verbose=False):
    """
    Read :term:`Triangle` :cite:`Shewchuk96b` mesh.

    Read :term:`Triangle` :cite:`Shewchuk96b` ASCII files and return a instance
    of GIMLI::Mesh class.
    See: ://www.cs.cmu.edu/~quake/triangle.html

    Parameters
    ----------
    fname : string
        Filename of the file to read (\\*.n, \\*.e)

    verbose : boolean, optional
        Be verbose during import.

    """

    raise("implement me!")
    os.system('meshconvert -d2 ' + fname)
    return pg.Mesh(2)


def readTetgen(fname, verbose=False):
    """
    Read :term:`Tetgen` :cite:`Si2004` mesh.

    Read :term:`Tetgen` :cite:`Si2004` ASCII files and return instance
    of GIMLI::Mesh class.
    See: http://tetgen.org/

    Parameters
    ----------
    fname : string
        Filename of the file to read (\\*.n, \\*.e \\*.f)

    verbose : boolean, optional
        Be verbose during import.

    """
    raise("implement me!")
    os.system('meshconvert -d3 -D ..' + fname)
    return pg.Mesh(3)


def readHydrus2dMesh(fname='MESHTRIA.TXT'):
    """
    Import mesh from Hydrus2D.

    Parameters
    ----------
    fname : str, optional
        Filename of Hydrus output file.

    See Also
    --------
    readHydrus3dMesh : Similar routine for three-dimensional meshes.

    References
    ----------
    .. [1] http://www.pc-progress.com/en/Default.aspx?h3d-description
    """
    fid = open(fname)
    line = fid.readline().split()
    nnodes = int(line[1])
    ncells = int(line[3])
    mesh = pg.Mesh()
    for i in range(nnodes):
        line = fid.readline().split()
        mesh.createNode(
            pg.RVector3(float(line[1]) / 100., float(line[2]) / 100., 0.))

    for i in range(3):
        line = fid.readline()

    for i in range(ncells):
        line = fid.readline().split()
        if len(line) == 4:
            mesh.createTriangle(
                mesh.node(int(line[1]) - 1),
                mesh.node(int(line[2]) - 1),
                mesh.node(int(line[3]) - 1),
                1)
        elif len(line) == 5:
            mesh.createTetrahedron(
                mesh.node(int(line[1]) - 1), mesh.node(int(line[2]) - 1),
                mesh.node(int(line[3]) - 1), mesh.node(int(line[4]) - 1), 1)

    fid.close()
    return mesh


def readHydrus3dMesh(filename='MESHTRIA.TXT'):
    """
    Import mesh from Hydrus3D.

    Parameters
    ----------
    fname : str, optional
        Filename of Hydrus output file.

    See Also
    --------
    readHydrus2dMesh : Similar routine for two-dimensional meshes.

    References
    ----------
    .. [1] http://www.pc-progress.com/en/Default.aspx?h3d-description
    """
    f = open(filename, 'r')
    for i in range(6):
        line1 = f.readline()

    nnodes = int(line1.split()[0])
    ncells = int(line1.split()[1])
    print(nnodes, ncells)
    line1 = f.readline()
    nodes = []
    dx = 0.01
    mesh = pg.Mesh()
    for ni in range(nnodes):
        pos = f.readline().split()
        p = pg.RVector3(
            float(pos[1]) * dx,
            float(pos[2]) * dx,
            float(pos[3]) * dx * (-1.))
        n = mesh.createNode(p)
        nodes.append(n)

    line1 = f.readline()
    line1 = f.readline()
    cells = []
    for ci in range(ncells):
        pos = f.readline().split()
        i, j, k, l = int(pos[1]), int(pos[2]), int(pos[3]), int(pos[4]),
        c = mesh.createTetrahedron(
            nodes[i - 1],
            nodes[j - 1],
            nodes[k - 1],
            nodes[l - 1])
        cells.append(c)

    f.close()
    return mesh


def transform2DMeshTo3D(mesh, x, y, z=None):
    """
    Transform a 2D mesh into 3D coordinates using a point list (e.g. from GPS)

    Parameters
    ----------
    mesh: GIMLi::Mesh
    x,y: array of x/y positions along 2d profile
    z: optional height to add (topographical correction if computed flat earth)

    See Also
    --------

    References
    ----------
    """

    # get mesh node positions
    mt, mz = pg.x(mesh.positions()), pg.y(mesh.positions())  # mesh tape and z
    # compute length of reference points along tape
    pt = np.hstack((0., np.cumsum(np.sqrt(np.diff(x)**2 + np.diff(y)**2))))
    #  interpolate node positions from tape to x/y using tape positions
    mx = np.interp(mt, pt, x)
    my = np.interp(mt, pt, y)
    # compute z offset by interpolating z
    if z is None:
        oz = np.zeros(len(mt))
    else:
        oz = np.interp(mt, pt, z)

    # set the positions in the mesh
    for i, node in enumerate(mesh.nodes()):
        node.setPos(pg.RVector3(mx[i], my[i], mz[i] + oz[i]))


def rot2DGridToWorld(mesh, start, end):
    """
    ..

    todo:: Complete Documentation. ...rotate a given 2D grid in...
    """
    mesh.rotate(pg.degToRad(pg.RVector3(-90.0, 0.0, 0.0)))

    src = pg.RVector3(0.0, 0.0, 0.0).norm(pg.RVector3(0.0, 0.0, -10.0),
                                          pg.RVector3(10.0, 0.0, -10.0))
    dest = start.norm(start - pg.RVector3(0.0, 0.0, 10.0), end)

    q = pg.getRotation(src, dest)
    rot = pg.RMatrix(4, 4)
    q.rotMatrix(rot)
    mesh.transform(rot)
    mesh.translate(start)


def merge2Meshes(m1, m2):
    """Merge two meshes into one new mesh and return combined mesh."""

    mesh = pg.Mesh(m1)

    for c in m2.cells():
        mesh.copyCell(c)

    for b in m2.boundaries():
        mesh.copyBoundary(b)

    for key in list(mesh.exportDataMap().keys()):
        d = mesh.exportDataMap()[key]
        print(d)
        d.resize(mesh.cellCount())
        d.setVal(m1.exportDataMap()[key], 0, m1.cellCount())
        d.setVal(
            m2.exportDataMap()[key],
            m1.cellCount(),
            m1.cellCount() + m2.cellCount())
        mesh.addExportData(key, d)

    return mesh


def mergeMeshes(meshlist):
    """
    Merge several meshes into one new mesh and return the new mesh.

    Parameters
    ----------
    meshlist : list
        List of at least two meshes to be merged.

    See Also
    --------
    merge2Meshes
    """

    if not isinstance(meshlist, list):
        raise "argument meshlist is no list"

    if len(meshlist) < 2:
        raise "to few meshes in meshlist"

    mesh = meshlist[0]

    for m in range(1, len(meshlist)):
        mesh = merge2Meshes(mesh, meshlist[m])

    return mesh


def createParaDomain2D(sensors, paraDX=1, paraDepth=0,
                       paraBoundary=2, paraMaxCellSize=0, boundary=-1,
                       verbose=False, *args, **kwargs):
    """
    Create a PLC mesh for an inversion parameter mesh.

    Create a PLC mesh for an inversion parameter mesh for a given list of
    sensor positions.
    Sensor position assumed on the surface and must be sorted and unique.

    The PLC is a :gimliapi:`GIMLI::Mesh` and contain nodes, edges and
    two region markers, one for the parameters domain (marker=2) and
    a larger boundary around the outside (marker=1)

    TODO
    ----
        * closed domains (boundary == 0)
        * additional topopoints
        * spline interpolations between sensorpoints or addpoints
        * subsurface sensors

    Parameters
    ----------
    sensors : list of RVector3 objects
        Sensor positions. Must be sorted and unique in positive x direction.
        Depth need to be y-coordinate.
    paraDX : float
        Relativ distance for refinement nodes between two electrodes,
        e.g., 0.5 means 1 additional node in the middle between two electrodes
        e.g., 0.33 means 2 additional node evenly spaced between two electrodes
    paraDepth : float, optional
        Maximum depth for parametric domain, 0 (default) means 0.4 * maximum
        sensor range.
    paraBoundary : float, optional
        Margin for parameter domain in absolute sensor distances. 2 (default).
    paraMaxCellSize: double, optional
        Maximum size for parametric size in m*m
    boundary : float, optional
        Boundary width to be appended for domain prolongation in absolute
        para domain width.
        Values lover 0 force the boundary to be 4 times para domain width.

    Returns
    -------
    poly: :gimliapi:`GIMLI::Mesh`

    """

    eSpacing = sensors[0].distance(sensors[1])

    xmin, ymin = sensors[0][0], sensors[0][1]
    xmax = xmin
    ymax = ymin
    for e in sensors:
        xmin = min(xmin, e[0])
        xmax = max(xmax, e[0])
        ymin = min(ymin, e[1])
        ymax = max(ymax, e[1])

    paraBound = eSpacing * paraBoundary

    if paraDepth == 0:
        paraDepth = 0.4 * (xmax - xmin)

    poly = pg.Mesh(2)
    # define para domain without surface
    n1 = poly.createNode([xmin - paraBoundary, sensors[0][1]])
    n2 = poly.createNode([xmin - paraBoundary, sensors[0][1] - paraDepth])
    n3 = poly.createNode([xmax + paraBoundary, sensors[-1][1] - paraDepth])
    n4 = poly.createNode([xmax + paraBoundary, sensors[-1][1]])

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
        poly.addRegionMarker(n12.pos() + [1e-3, 1e-3], 1)

    poly.createEdge(n1, n2, 1)
    poly.createEdge(n2, n3, 1)
    poly.createEdge(n3, n4, 1)
    poly.addRegionMarker(n2.pos() + [1e-3, 1e-3], 2, paraMaxCellSize)

    # define surface
    nSurface = []
    nSurface.append(n1)
    for i, e in enumerate(sensors):
        if paraDX >= 0.5:
            nSurface.append(poly.createNode(e, pg.MARKER_NODE_SENSOR))
            if (i < len(sensors) - 1):
                nSurface.append(poly.createNode((sensors[i + 1] + e) * 0.5))
            # print("Surface add ", e, el, nSurface[-2].pos(),
            #        nSurface[-1].pos())
        elif paraDX < 0.5:
            if (i > 0):
                nSurface.append(
                    poly.createNode(e - (e - sensors[i - 1]) * paraDX))
            nSurface.append(poly.createNode(e, pg.MARKER_NODE_SENSOR))
            if (i < len(sensors) - 1):
                nSurface.append(
                    poly.createNode(e + (sensors[i + 1] - e) * paraDX))
            # print("Surface add ", nSurface[-3].pos(), nSurface[-2].pos(),
            #        nSurface[-1].pos())
    nSurface.append(n4)

    for i in range(len(nSurface) - 1, 0, -1):
        poly.createEdge(nSurface[i], nSurface[i - 1],
                        pg.MARKER_BOUND_HOMOGEN_NEUMANN)

#     for n in poly.nodes():
#        print(n.pos())
#    sys.exit()
    return poly


def createParaMesh2dGrid(sensors, paraDX=1, paraDZ=1, paraDepth=0, nLayers=11,
                         boundary=-1, paraBoundary=2, verbose=False, *args,
                         **kwargs):
    """
    Create a grid style mesh for an inversion parameter mesh.

    Return parameter grid for a given list of sensor positions.

    Parameters
    ----------
    sensors : list of RVector3 objects
        Sensor positions. Must be sorted in positve x direction
    paraDX : float, optional
        Horizontal distance between sensors, relative regarding sensor
        distance. Value must be greater than 0 otherwise 1 is assumed.
    paraDZ : float, optional
        Vertical distance to the first depth layer, relative regarding sensor
        distance. Value must be greater than 0 otherwise 1 is assumed.
    paraDepth : float, optional
        Maximum depth for parametric domain, 0 (default) means 0.4 * maxmimum
        sensor range.
    nLayers : int, optional
        Number of depth layers.
    boundary : int, optional
        Boundary width to be appended for domain prolongation in absolute
        para domain width.
        Values lover 0 force the boundary to be 4 times para domain width.
    paraBoundary : int, optional
        Offset for parameter domain boundary in absolute sensor distance.
        2 (default).
    verbose : boolean, optional
        Be verbose.

    Returns
    -------
    mesh: :gimliapi:`GIMLI::Mesh`

    """

    mesh = pg.Mesh(2)

    # maybe separate x y z and sort
    if isinstance(sensors, np.ndarray):
        sensors = [pg.RVector3(s, 0) for s in sensors]

    sensorX = pg.x(sensors)

    eSpacing = abs(sensorX[1] - sensorX[0])

    xmin = min(sensorX) - paraBoundary * eSpacing
    xmax = max(sensorX) + paraBoundary * eSpacing

    if paraDX == 0:
        paraDX = 1.
    if paraDZ == 0:
        paraDZ = 1.

    dx = eSpacing * paraDX
    dz = eSpacing * paraDZ

    if paraDepth == 0:
        paraDepth = 0.4 * (xmax - xmin)

    x = pg.utils.grange(xmin, xmax, dx=dx)

    y = -pg.increasingRange(dz, paraDepth, nLayers)

    mesh.createGrid(x, y)

    list(map(lambda cell: cell.setMarker(2), mesh.cells()))

    paraXLimits = [xmin, xmax]
#    paraYLimits = [min(y), max(y)]  # not used

    if boundary < 0:
        boundary = abs((paraXLimits[1] - paraXLimits[0]) * 4.0)

    mesh = pg.meshtools.appendTriangleBoundary(mesh,
                                               xbound=boundary,
                                               ybound=boundary,
                                               marker=1, *args, **kwargs)

    return mesh
