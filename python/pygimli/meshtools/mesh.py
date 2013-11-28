# -*- coding: utf-8 -*-

import pygimli as g
import numpy as np


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
        - Physical Number 1: homogenous Neumann condition
        - Physical Number 2: mixed boundary condition
        - Physical Number 3: homogeneous Dirichlet condition
        - Physical Number 4: Dirichlet condition
    - Physical Surfaces and Volumes define regions in 2D and 3D, respectively.
        - Physical Number 1: No inversion region
        - Physical Number >= 2: Inversion region

    References
    ----------
    .. [1] C. Geuzaine and J.-F. Remacle. Gmsh: a three-dimensional finite
           element mesh generator with built-in pre- and post-processing
           facilities. International Journal for Numerical Methods in
           Engineering 79(11), pp. 1309-1331, 2009.
    """
    inNodes, inElements, ncount, ecount = 0, 0, 0, 0
    fid = open(fname)
    if verbose:
        print 'Reading %s... \n' % fname

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
                        print '  Nodes: %s' % int(line)
                else:
                    nodes[ncount, :] = np.array(line.split(), 'float')[1:]
                    ncount += 1

            elif inElements == 1:
                if len(line.split()) == 1:
                    if verbose:
                        print '  Entries: %s' % int(line)
                    points, lines, triangles, tets = [], [], [], []

                else:
                    entry = map(int, line.split())[1:]

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
        print '    Points: %s' % len(points)
        print '    Lines: %s' % len(lines)
        print '    Triangles: %s' % len(triangles)
        print '    Tetrahedra: %s \n' % len(tets)
        print 'Creating mesh object... \n'

    # check dimension
    if len(tets) == 0:
        dim, bounds, cells = 2, lines, triangles
        zero_dim = np.abs(nodes.sum(0)).argmin()  # identify zero dimension
    else:
        dim, bounds, cells = 3, triangles, tets
    if verbose:
        print '  Dimension: %s-D' % dim

    # creating instance of GIMLI::Mesh class
    mesh = g.Mesh(dim)

    # replacing boundary markers (gmsh does not allow negative physical
    # regions)
    bound_marker = (g.MARKER_BOUND_HOMOGEN_NEUMANN, g.MARKER_BOUND_MIXED,
                    g.MARKER_BOUND_HOMOGEN_DIRICHLET, g.MARKER_BOUND_DIRICHLET)
    for i in range(4):
        bounds[:, dim][bounds[:, dim] == i + 1] = bound_marker[i]

    if verbose:
        bound_types = np.unique(bounds[:, dim])
        regions = np.unique(cells[:, dim + 1])
        print '  Regions: %s ' % len(regions) + str(tuple(regions))
        print '  Boundary types: %s ' % len(bound_types) + str(tuple(bound_types))

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
        points = np.asarray(points)
        node_types = np.unique(points[:, 1])
        print '  Marked nodes: %s ' % len(points) + str(tuple(node_types))
        print '\nDone. \n'
        print '  ' + str(mesh)

    return mesh


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
    mesh = g.Mesh()
    for i in range(nnodes):
        line = fid.readline().split()
        mesh.createNode(
            g.RVector3(float(line[1]) / 100.,
                       float(line[2]) / 100.,
                       0.))

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
    print nnodes, ncells
    line1 = f.readline()
    nodes = []
    dx = 0.01
    mesh = g.Mesh()
    for ni in range(nnodes):
        pos = f.readline().split()
        p = g.RVector3(
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

def transform2DMeshTo3D( mesh, x, y, z=None ):
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
    mt, mz = g.x( mesh.positions() ), g.y( mesh.positions() ) # mesh tape and z
    # compute length of reference points along tape
    pt = np.hstack( (0., np.cumsum( np.sqrt( np.diff( x )**2 + np.diff( y )**2 ) ) ) )
    #  interpolate node positions from tape to x/y using tape positions
    mx = np.interp( mt, pt, x )
    my = np.interp( mt, pt, y )
    # compute z offset by interpolating z
    if z is None:
        oz = np.zeros( len(mt) )
    else:
        oz = np.interp( mt, pt, z )
    
    # set the positions in the mesh
    for i, node in enumerate( mesh.nodes() ):
        node.setPos( g.RVector3( mx[i], my[i], mz[i]+oz[i] ) )

def rot2DGridToWorld(mesh, start, end):
    """
    ..

    todo:: Complete Documentation. ...rotate a given 2D grid in...
    """
    mesh.rotate(g.degToRad(g.RVector3(-90.0, 0.0, 0.0)))

    src = g.RVector3(0.0, 0.0, 0.0).norm(g.RVector3(0.0, 0.0, -10.0),
                                         g.RVector3(10.0, 0.0, -10.0))
    dest = start.norm(start - g.RVector3(0.0, 0.0, 10.0), end)

    q = g.getRotation(src, dest)
    rot = g.RMatrix(4, 4)
    q.rotMatrix(rot)
    mesh.transform(rot)
    mesh.translate(start)


def merge2Meshes(m1, m2):
    """Merge two meshes into one new mesh and return combined mesh."""

    mesh = g.Mesh(m1)

    for c in m2.cells():
        mesh.copyCell(c)

    for b in m2.boundaries():
        mesh.copyBoundary(b)

    for key in mesh.exportDataMap().keys():
        d = mesh.exportDataMap()[key]
        print d
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


def createParaMesh2dGrid(sensors, paraDX=1, paraDZ=1, paraDepth=0, nLayers=11,
                         boundary=-1, paraBoundary=2, verbose=False, *args,
                         **kwargs):
    """
    Return parameter mesh for a given list of sensor positions.

    Parameters
    ----------
    sensors : list of RVector3 objects
        Sensor positions.
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
        Triangle boundary to be appended for domain prolongation values lover 0
        force boundary to be 4 times para domain width.
    paraBoundary : int, optional
        Offset for parameter domain in absolute sensor distance 2 (default).
    verbose : boolean, optional
        Be verbose.
    """
    mesh = g.Mesh(2)

    # maybe separate x y z and sort
    sensorX = g.x(sensors)
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

    x = g.utils.grange(xmin, xmax, dx=dx)

    y = -g.increasingRange(dz, paraDepth, nLayers)

    mesh.createGrid(x, y)

    map(lambda cell: cell.setMarker(2), mesh.cells())

    paraXLimits = [xmin, xmax]
    paraYLimits = [min(y), max(y)]

    if boundary < 0:
        boundary = abs((paraXLimits[1] - paraXLimits[0]) * 4.0)

    mesh = g.meshtools.appendTriangleBoundary(
        mesh, xbound=boundary, ybound=boundary, marker=1,
        *args, **kwargs)

    return mesh
