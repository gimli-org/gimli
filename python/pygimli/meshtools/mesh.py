# -*- coding: utf-8 -*-
"""Generally mesh generation and maintenance."""

import numpy as np
import os

import pygimli as pg


def createMesh(poly, quality=30, area=0.0, smooth=None, switches=None,
               verbose=False, **kwargs):
    """Create a mesh for a given geometry polygon.

    The mesh is created by :term:`triangle` or :term:`tetgen` if the
    gimli support for these mesh generators are installed.
    The geometry needs to contain nodes and boundaries and should be valid
    in the sense that the boundaries are non intersecting.

    If poly is a list of coordinates a simple Delaunay mesh of the convex hull
    will be created.
    TODO: Tetgen support need to be implemented

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
        Be careful with values above 34 degrees.
    area: float
        2D maximum triangle size in m*m
    smooth: tuple
        [smoothing algorithm, number of iterations]
        0, no smoothing
        1, node center
        2, weighted node center
    switches: str
        Force triangle to use the gives command switches.

    Returns
    -------
    mesh: :gimliapi:`GIMLI::Mesh`

    Examples
    --------
    >>> # no need to import matplotlib. pygimli's show does
    >>> import pygimli as pg
    >>> import pygimli.meshtools as mt
    >>> rect = mt.createRectangle(start=[0, 0], end=[4, 1])
    >>> ax, _ = pg.show(mt.createMesh(rect, quality=10))
    >>> ax, _ = pg.show(mt.createMesh(rect, quality=33))
    >>> ax, _ = pg.show(mt.createMesh(rect, quality=33, area=0.01))
    >>> pg.wait()
    """
    #  poly == [pg.Mesh, ]
    if isinstance(poly, list):
        if isinstance(poly[0], pg.Mesh):
            return createMesh(
                pg.meshtools.mergePLC(poly), quality, area, smooth, switches,
                verbose)
    # poly == [pos, pos, ]
    if isinstance(poly, list) or isinstance(poly, type(zip)):
        delPLC = pg.Mesh(2)
        for p in poly:
            delPLC.createNode(p[0], p[1], 0.0)
        return createMesh(delPLC, switches='-zeY')

    # poly == Mesh
    if poly.dim() == 2:
        if poly.nodeCount() == 0:
            raise Exception("No nodes in poly to create a valid mesh")

        tri = pg.TriangleWrapper(poly)

        if switches is None:
            # -D Conforming delaunay
            # -F Uses Steven Fortune's sweepline algorithm
            # no -a here ignores per region area
            switches = 'pazeA'

            if area > 0:
                switches += 'a' + str(area)
                pass
            else:
                switches += 'a'

            # switches = switches.replace('.', ',')
            switches += 'q' + str(quality)

        if not verbose:
            switches += 'Q'

        if verbose:
            print(switches)

        tri.setSwitches(switches)
        mesh = tri.generate()

        if smooth is not None:
            mesh.smooth(nodeMoving=kwargs.pop('node_move', False),
                        edgeSwapping=False, smoothFunction=smooth[0],
                        smoothIteration=smooth[1])
        return mesh

    else:
        raise Exception('not yet implemented')


def refineQuad2Tri(mesh, style=1):
    """Refine mesh of quadrangles into a mesh of triangle cells.

        TODO mixed meshes

    Parameters
    ----------
    mesh : :gimliapi:`GIMLI::Mesh`
        Mesh containing quadrangle cells.

    style: int [1]
        * 1 bisect each quadrangle into 2 triangles
        * 2 bisect each quadrangle into 4 triangles

    Returns
    -------
    ret : :gimliapi:`GIMLI::Mesh`
        Mesh containing triangle cells.

    Examples
    --------
    >>> # no need to import matplotlib. pygimli's show does
    >>> import pygimli as pg
    >>> import pygimli.meshtools as mt
    >>> quads = pg.createGrid(range(10), range(10))
    >>> ax, _ = pg.show(quads)
    >>> ax, _ = pg.show(mt.refineQuad2Tri(quads, style=1))
    >>> ax, _ = pg.show(mt.refineQuad2Tri(quads, style=2))
    >>> pg.wait()

    """
    out = pg.Mesh(2)
    newNode = None

    for n in mesh.nodes():
        out.createNode(n.pos())

    for c in mesh.cells():

        if style == 1:
            out.createCell([c.node(0).id(), c.node(1).id(), c.node(2).id()])
            out.createCell([c.node(0).id(), c.node(2).id(), c.node(3).id()])

        elif style == 2:
            newNode = out.createNodeWithCheck(c.center())

            for i in range(4):
                out.createCell([c.node(i).id(), c.node((i + 1) % 4).id(),
                                newNode.id()])

        for i in range(c.boundaryCount()):
            b = c.boundary(i)
            if b.marker() != 0:
                out.createBoundary([b.node(0).id(), b.node(1).id()],
                                   b.marker())

    out.createNeighbourInfos()

    return out


def readGmsh(fname, verbose=False):
    r"""Read :term:`Gmsh` ASCII file and return instance of GIMLI::Mesh class.

    Parameters
    ----------
    fname : string
        Filename of the file to read (\\*.msh). The file must conform
        to the `MSH ASCII file version 2
        <http://gmsh.info/doc/texinfo/gmsh.html#MSH-ASCII-file-format>`_ format
    verbose : boolean, optional
        Be verbose during import.

    Notes
    -----
    Physical groups specified in Gmsh are interpreted as follows:

    - Points with the physical number 99 are interpreted as sensors.
    - Physical Lines and Surfaces define boundaries in 2D and 3D, respectively.
        - Physical Number 1: Homogeneous Neumann condition
        - Physical Number 2: Mixed boundary condition
        - Physical Number 3: Homogeneous Dirichlet condition
        - Physical Number 4: Dirichlet condition
    - Physical Surfaces and Volumes define regions in 2D and 3D, respectively.
        - Physical Number 1: No inversion region
        - Physical Number >= 2: Inversion region

    Examples
    --------
    >>> import tempfile, os
    >>> from pygimli.meshtools import readGmsh
    >>> gmsh = '''
    ... $MeshFormat
    ... 2.2 0 8
    ... $EndMeshFormat
    ... $Nodes
    ... 3
    ... 1 0 0 0
    ... 2 0 1 0
    ... 3 1 1 0
    ... $EndNodes
    ... $Elements
    ... 7
    ... 1 15 2 0 1 1
    ... 2 15 2 0 2 2
    ... 3 15 2 0 3 3
    ... 4 1 2 0 1 2 3
    ... 5 1 2 0 2 3 1
    ... 6 1 2 0 3 1 2
    ... 7 2 2 0 5 1 2 3
    ... $EndElements
    ... '''
    >>> fname = tempfile.mktemp()
    >>> with open(fname, "w") as f:
    ...     f.writelines(gmsh)
    >>> mesh = readGmsh(fname)
    >>> print(mesh)
    Mesh: Nodes: 3 Cells: 1 Boundaries: 3
    >>> os.remove(fname)
    """
    inNodes, inElements, ncount = 0, 0, 0
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
                    entry = [int(e_) for e_ in line.split()][1:]

                    if entry[0] == 15:
                        points.append((entry[-2], entry[-3]))
                    elif entry[0] == 1:
                        lines.append((entry[-2], entry[-1], entry[2]))
                    elif entry[0] == 2:
                        triangles.append((entry[-3], entry[-2], entry[-1],
                                          entry[2]))
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

    if bounds.any():
        for i in range(4):
            bounds[:, dim][bounds[:, dim] == i + 1] = bound_marker[i]

        # account for CEM markers
        bounds[:, dim][bounds[:, dim] >= 10000] *= -1

        if verbose:
            bound_types = np.unique(bounds[:, dim])
            print('  Boundary types: %s ' % len(bound_types) + str(tuple(
                bound_types)))
    else:
        print("WARNING: No boundary conditions found.",
              "Setting Neumann on the outer edges by default.")

    if verbose:
        regions = np.unique(cells[:, dim + 1])
        print('  Regions: %s ' % len(regions) + str(tuple(regions)))

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
                mesh.node(int(cell[2] - 1)), mesh.node(int(cell[3] - 1)),
                marker=int(cell[4]))

    mesh.createNeighbourInfos()

    # Set Neumann on outer edges by default (can be overriden by Gmsh info)
    for b in mesh.boundaries():
        if not b.leftCell() or not b.rightCell():
            b.setMarker(pg.MARKER_BOUND_HOMOGEN_NEUMANN)

    for bound in bounds:
        if dim == 2:
            mesh.createEdge(
                mesh.node(int(bound[0] - 1)), mesh.node(int(bound[1] - 1)),
                marker=int(bound[2]))
        else:
            mesh.createTriangleFace(
                mesh.node(int(bound[0] - 1)), mesh.node(int(bound[1] - 1)),
                mesh.node(int(bound[2] - 1)), marker=int(bound[3]))

    # assign marker to corresponding nodes (sensors, reference nodes, etc.)
    if points:
        for point in points:
            mesh.node(point[0] - 1).setMarker(-point[1])

    if verbose:
        if points:
            points = np.asarray(points)
            node_types = np.unique(points[:, 1])
            print('  Marked nodes: %s ' % len(points) + str(tuple(node_types)))
        print('\nDone. \n')
        print('  ' + str(mesh))

    return mesh


def readTriangle(fname, verbose=False):
    r"""Read :term:`Triangle` :cite:`Shewchuk96b` mesh.

    Read :term:`Triangle` :cite:`Shewchuk96b` ASCII mesh files and return an
    instance of GIMLI::Mesh class.
    See: ://www.cs.cmu.edu/~quake/triangle.html

    Parameters
    ----------
    fname : string
        Filename of the file to read (\\*.n, \\*.e)

    verbose : boolean, optional
        Be verbose during import.

    """
    raise Exception("implement me!" + fname + str(verbose))
    # os.system('meshconvert -d2 ' + fname)
    # return pg.Mesh(2)


def readTetgen(fname, comment='#', verbose=True,
               default_cell_marker=0, load_faces=True,
               quadratic=False):
    """
    Reads and converts a mesh from the basic :term:`Tetgen` output.

    Read :term:`Tetgen` :cite:`Si2004` ASCII files and return instance
    of :gimliapi:`GIMLI::Mesh` class.
    See: http://tetgen.org/

    Parameters
    ----------

    fname: str
        Base name of the tetgen output, without ending. All additional files
        (.ele and .face respectively) has to be the same basename.

    comment: str ('#')
        String consisting of all symbols that indicates a comment in the input
        files. Standard for tetgen files is the '#'.

    verbose: boolean (True)
        Enables console output during the import process.

    default_cell_marker: int (0)
        :term:`Tetgen` files can contain cell markers, but doesn't have to.
        If no marker are found, the given interger is used.

    load_faces:
        Optional decision weather the faces of the :term:`Tetgen` output are
        loaded or not. Note that without the -f in during the tetgen call,
        the faces in the .face file will only contain the faces of the original
        input poly file and not all faces. If only part of the faces are
        imported a createNeighbourInfos call of the mesh will fail.

    quadratic: boolean (False)
        Returns a P2 refined mesh when True (to be removed, as soon as direct
        import of quadratic meshs is possible).

    Returns
    -------
    mesh: :gimliapi:`GIMLI::Mesh`
    """
    mesh = pg.Mesh(3)

    # Part 1/3: Nodes, essential
    with open(fname + '.node', 'r') as node_in:
        node_lines = pg.utils.filterLinesByCommentStr(node_in.readlines(),
                                                      comment)
    node_1 = node_lines[0].split()
    node_count = int(node_1[0])
    assert int(node_1[1]) == 3, 'Wrong dim: {}, should be 3.'.format(node_1[1])
    number_node_attr = int(node_1[2])
    node_markers = int(node_1[3])
    node_attributes = []
    for i in range(number_node_attr):
        node_attributes.append([])

    for n in range(node_count):
        node_n = node_lines[n + 1].split()
        if node_markers:
            node_marker_n = int(node_n[-1])
        else:
            node_marker_n = n
        mesh.createNode([float(node_n[1]),
                         float(node_n[2]),
                         float(node_n[3])],
                        marker=node_marker_n)
        for m in range(number_node_attr):
            node_attributes[m].append(float(node_n[4 + m]))

    for k in range(number_node_attr):
        if verbose:
            print('Add node data to mesh.')
        mesh.addData('node_data_{}'.format(k + 1), node_attributes[k])

    # Part 2/3: Tetrahedrons, optional
    if os.path.exists(fname + '.ele'):
        if verbose:
            print('Found .ele file. Adding cells and cell marker.')
        with open(fname + '.ele', 'r') as cell_in:
            cell_lines = pg.utils.filterLinesByCommentStr(cell_in.readlines(),
                                                          comment)
        cell_1 = cell_lines[0].split()
        cell_count = int(cell_1[0])

        nodes_per_cell = int(cell_1[1])  # 4 or 10
        if nodes_per_cell == 10:
            quadratic = True
            raise Exception('Cannot import quadratic meshes directly yet.')

        cell_markers = int(cell_1[2])
        for c in range(cell_count):
            cell_n = cell_lines[c + 1].split()
            if cell_markers:
                cell_marker_n = int(cell_n[-1])
            else:
                cell_marker_n = default_cell_marker
            mesh.createCell([int(ind) for ind in cell_n[1:5]],
                            marker=cell_marker_n)
            # in order to import quadratic meshes directly, i ned the sorting
            # of the node indices
#           mesh.createCell([int(ind) for ind in cell_n[1:nodes_per_cell + 1]],
#                           marker=cell_marker_n)

    # Part 3/3: Boundaries and Marker, optional
    if os.path.exists(fname + '.face') and load_faces:
        if verbose:
            print('Found .face file. Adding boundaries and boundary marker.')
        with open(fname + '.face', 'r') as face_in:
            face_lines = pg.utils.filterLinesByCommentStr(face_in.readlines(),
                                                          comment)

        face_1 = face_lines[0].split()
        face_count = int(face_1[0])
        face_markers = int(face_1[1])

        for f in range(face_count):
            face_n = face_lines[f + 1].split()
            if face_markers:
                face_marker_n = int(face_n[-1])
            else:
                face_marker_n = 0
            mesh.createBoundary(
                [int(ind) for ind in face_n[1:4]],
                marker=face_marker_n)
            # quadratic
#            mesh.createBoundary(
#                [int(ind) for ind in face_n[1:len(face_n) - face_markers]],
#                marker=face_marker_n)

    if quadratic:
        mesh = mesh.createP2()

    if verbose:
        print(mesh)

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
    .. http://www.pc-progress.com/en/Default.aspx?h3d-description
    """
    fid = open(fname)
    line = fid.readline().split()
    nnodes = int(line[1])
    ncells = int(line[3])
    mesh = pg.Mesh()
    for _ in range(nnodes):
        line = fid.readline().split()
        mesh.createNode(pg.RVector3(
            float(line[1]) / 100., float(line[2]) / 100., 0.))

    for _ in range(3):
        line = fid.readline()

    for _ in range(ncells):
        line = fid.readline().split()
        if len(line) == 4:
            mesh.createTriangle(
                mesh.node(int(line[1]) - 1), mesh.node(int(line[2]) - 1),
                mesh.node(int(line[3]) - 1), 1)
        elif len(line) == 5:
            mesh.createTetrahedron(
                mesh.node(int(line[1]) - 1), mesh.node(int(line[2]) - 1),
                mesh.node(int(line[3]) - 1), mesh.node(int(line[4]) - 1), 1)

    fid.close()
    mesh.createNeighbourInfos()
    return mesh


def readHydrus3dMesh(filename='MESHTRIA.TXT'):
    """Import mesh from Hydrus3D.

    Parameters
    ----------
    fname : str, optional
        Filename of Hydrus output file.

    See Also
    --------
    readHydrus2dMesh : Similar routine for two-dimensional meshes.

    References
    ----------
    .. http://www.pc-progress.com/en/Default.aspx?h3d-description
    """
    f = open(filename, 'r')
    for i in range(6):
        line1 = f.readline()

    nnodes = int(line1.split()[0])
    ncells = int(line1.split()[1])
    # print(nnodes, ncells)

    line1 = f.readline()
    nodes = []
    dx = 0.01
    mesh = pg.Mesh()
    for _ in range(nnodes):
        pos = f.readline().split()
        p = pg.RVector3(
            float(pos[1]) * dx, float(pos[2]) * dx, float(pos[3]) * dx * (-1.))
        n = mesh.createNode(p)
        nodes.append(n)

    line1 = f.readline()
    line1 = f.readline()
    cells = []
    for _ in range(ncells):
        pos = f.readline().split()
        i, j, k, l = int(pos[1]), int(pos[2]), int(pos[3]), int(pos[4]),
        c = mesh.createTetrahedron(nodes[i - 1], nodes[j - 1], nodes[k - 1],
                                   nodes[l - 1])
        cells.append(c)

    f.close()
    mesh.createNeighbourInfos()
    return mesh


def readGambitNeutral(filename, verbose=False):
    r"""Import Gambit Neutral meshes *.neu.

    See. https://www.sharcnet.ca/Software/Gambit/html/users_guide/ug01.htm

    Not fully implemented. If needed, we can improve this importer just send us
    an example file.

    Parameters
    ----------
    fname : string
        Filename of the file to read (\\*.n, \\*.e \\*.f)

    verbose : boolean, optional
        Be verbose during import.
    """
    with open(filename, 'r') as fi:
        content = fi.readlines()
    fi.close()

    mesh = pg.Mesh(2)

    for i, line in enumerate(content):
        if 'ENDOFSECTION' in line:
            break

        try:
            nVerts = int(content[i - 1].split()[0])
            nElements = int(content[i - 1].split()[1])
        except:
            raise Exception("Cannot interpret GAMBIT Neutral header: " +
                            content[0:i])

    for i, line in enumerate(content):
        if 'NODAL COORDINATES' in line:
            break
    for j in range(nVerts):
        vertx = float(content[i + j + 1].split()[1])
        verty = float(content[i + j + 1].split()[2])
        mesh.createNode((vertx, verty))

    for i, line in enumerate(content):
        if 'ELEMENTS/CELLS' in line:
            break
    for j in range(nElements):
        nNodes = int(content[i + j + 1].split()[1])
        nodes = []
        for k in range(nNodes):
            nodes.append(int(content[i + 1 + j].split()[3 + k]) - 1)
        mesh.createCell(nodes)

    if verbose:
        print("Gambit neutral file imported: ", mesh)

    mesh.createNeighbourInfos()
    return mesh


def convertHDF5Mesh(h5Mesh, group='mesh', indices='cell_indices',
                    pos='coordinates', cells='topology', marker='values',
                    marker_default=0, dimension=3, verbose=True,
                    useFenicsIndices=False):
    """
    Converts instance of a hdf5 mesh to a :gimliapi:`GIMLI::Mesh`.
    For full documentation please see :py:mod:`pygimli:meshtools:readHDF5Mesh`.
    """
    # open mesh containing group inside hdf file
    inmesh = h5Mesh.get(group)

    # get node positions
    mesh_pos = inmesh[pos][:]

    # get cells
    if useFenicsIndices:
        # case 1/2: using fenics specific indices
        # TODO why extra indices? cause this is really slow
        # figure out the nature and purpose of the extra fenics indices
        mesh_indices = inmesh[indices]
        try:
            mesh_cells = inmesh[cells][mesh_indices]
        except IndexError as IE:
            print('Fenics Indices arren\'t just in arbitrary order in range\
(0, cellCount) as expected. Need Fix.')
            raise IE
    else:
        # case 2/2: indices implicit: [0, ... cellCount)
        mesh_cells = inmesh[cells][:]

    # get marker
    try:
        # case 1/3: marker found in hdf file
        mesh_marker = inmesh[marker][:]
    except KeyError:
        if isinstance(marker_default, int):
            # case 2/3: marker not found, given as int
            mesh_marker = np.ones(len(mesh_cells), dtype=int) * marker_default
        else:
            # case 3/3: marker not found, given as array_like
            mesh_marker = marker_default

    # start building pygimli mesh
    mesh = pg.Mesh(dimension)
    for node in mesh_pos:
        mesh.createNode(node)
    for i, cell in enumerate(mesh_cells):
        mesh.createCell(pg.IndexArray(cell), marker=int(mesh_marker[i]))
    mesh.createNeighbourInfos()
    if verbose:
        print('converted mesh:', mesh)
    return mesh


def readHDF5Mesh(filename, group='mesh', indices='cell_indices',
                 pos='coordinates', cells='topology', marker='values',
                 marker_default=0, dimension=3, verbose=True,
                 useFenicsIndices=False):
    '''
    Function for loading a mesh from HDF5 file format.

    Returns an instance of :gimliapi:`GIMLI::Mesh` class.
    Default values for keywords are suited for :term:`FEniCS` syntax
    .h5 meshes.

    Requirements: h5py module

    TODO:
        * Fenics hdf5 meshs doesn't have boundary markers.

    Parameters
    ----------

    filename: string
        Name of the mesh that has to be transformed into :term:`pyGIMLi`
        format.

    group: string ['domains']
        hdf group that contains the mesh informations (see other keyword
        arguments). Default is 'domains' for :term:`FEniCS` compatibility.

    indices: string ['cell_indices']
        Key for the part of the hdf file containing the indices of the cells.

    pos: string ['coordinates']
        Key for the part of the hdf file containing the nodepositions.

    cells: string ['topology']
        Key for the part of the hdf file containing the array which defies the
        cells. Usually of shape (cellCount, 3) for 2D meshs or (cellCount, 4)
        for 3D tetrahedra meshes. For each cell the indices of the
        corresponding node indices is given.

    marker: string ['values']
        If marker is part of the hdf data container, the corresponding array
        is used as identifier for the cell markers. If not found, the cell
        markers will be set to marker_default.

    marker_default: int or array [0]
        Default marker if no markers are found in the hdf file. If array, size
        has to match the cellCount of the mesh.

    dimension: int [3]
        Dimension of the in/outpu mesh, no own check for dimensions yet.
        Fixed on 3 for now.

    Yields
    ------

    mesh:
        :gimliapi:`GIMLI::Mesh`

    '''
    h5py = pg.io.opt_import('h5py',
                            requiredFor='import mesh in .h5 data format')
    h5 = h5py.File(filename, 'r')
    if verbose:
        print('loaded hdf5 mesh:', h5)

    mesh = convertHDF5Mesh(h5, group=group, indices=indices,
                           pos=pos, cells=cells, marker=marker,
                           marker_default=marker_default, dimension=dimension,
                           verbose=verbose, useFenicsIndices=useFenicsIndices)

    h5.close()
    return mesh


def readFenicsHDF5Mesh(filename, group='mesh', verbose=True):
    """
    Reads :term:`FEniCS` mesh from file format .h5 and returns a
    :gimliapi:`GIMLI::Mesh`.
    """
    mesh = readHDF5Mesh(filename, group=group, indices='cell_indices',
                        pos='coordinates', cells='topology', marker='values',
                        marker_default=0, dimension=3, verbose=verbose,
                        useFenicsIndices=False)
    return mesh


def exportHDF5Mesh(mesh, exportname, group='mesh', indices='cell_indices',
                   pos='coordinates', cells='topology', marker='values'):
    '''
    Writes given in a hdf5 format file.

    3D tetrahedron meshes only! Boundary markers are ignored.

    Keywords are explained in :py:mod:`pygimli.meshtools.readHDFS`
    '''
    h5py = pg.io.opt_import('h5py',
                            requiredFor='export mesh in .h5 data format')
    if not isinstance(mesh, pg.Mesh):
        mesh = pg.Mesh(mesh)

    # prepare output for writing in hdf data container
    pg_pos = mesh.positions()
    mesh_pos = np.array((np.array(pg.x(pg_pos)),
                         np.array(pg.y(pg_pos)),
                         np.array(pg.z(pg_pos)))).T

    mesh_cells = np.zeros((mesh.cellCount(), 4))  # hard coded for tetrahedrons
    for i, cell in enumerate(mesh.cells()):
        mesh_cells[i] = cell.ids()

    mesh_indices = np.arange(0, mesh.cellCount() + 1, 1, dtype=np.int64)
    mesh_markers = np.array(mesh.cellMarkers())

    with h5py.File(exportname, 'w') as out:
        # writing indices
        idx_name = '{}/{}'.format(group, indices)
        out.create_dataset(idx_name, data=mesh_indices, dtype=np.int64)
        # writing node positions
        pos_name = '{}/{}'.format(group, pos)
        out.create_dataset(pos_name, data=mesh_pos, dtype=float)
        # writing cells via indices
        cells_name = '{}/{}'.format(group, cells)
        out.create_dataset(cells_name, data=mesh_cells, dtype=np.int64)
        # writing marker
        marker_name = '{}/{}'.format(group, marker)
        out.create_dataset(marker_name, data=mesh_markers, dtype=np.uint64)
        out[group][cells].attrs['celltype'] = np.array(('tetrahedron'))
        out[group][cells].attrs['partition'] = np.array([0], dtype=np.uint64)
    return True


def exportFenicsHDF5Mesh(mesh, exportname, group='mesh'):
    """
    Exports Gimli mesh in HDF5 format suitable for Fenics.

    Equivalent to calling the function
    :py:mod:`pygimli.meshtools.exportHDF5Mesh(mesh, exportname, group=group,
    indices='cell_indices', pos='coordinates', cells='topology',
    marker='values')`.

    Parameters
    ----------

    mesh: :gimliapi:GIMLI::Mesh`
        Mesh to be saved.

    exportname: string
        Name under which the mesh is saved.

    group: string ('mesh')
        Identification string under which the mesh is saved. Important for
        :term:`FEniCS` to know.

    """
    return exportHDF5Mesh(mesh, exportname, group=group,
                          indices='cell_indices',
                          pos='coordinates', cells='topology', marker='values')


def transform2DMeshTo3D(mesh, x, y, z=None):
    """2D mesh into 3D coordinates.

    Transform a 2D mesh into 3D coordinates using a point list (e.g. from GPS)

    Parameters
    ----------
    mesh: :gimliapi:`GIMLI::Mesh`

    x,y: [float]
        array of x/y positions along 2d profile

    z: [float]
        optional height to add (topographical correction on top of flat earth)

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
    """Rotate a 2D mesh into 3D world coordinates.

    todo:: Complete Documentation. ...rotate a given 2D grid in...
    """
    mesh.rotate(pg.degToRad(pg.RVector3(-90.0, 0.0, 0.0)))

    src = pg.RVector3(0.0, 0.0, 0.0).norm(
        pg.RVector3(0.0, 0.0, -10.0), pg.RVector3(10.0, 0.0, -10.0))
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
        d.resize(mesh.cellCount())
        d.setVal(m1.exportDataMap()[key], 0, m1.cellCount())
        d.setVal(m2.exportDataMap()[key], m1.cellCount(),
                 m1.cellCount() + m2.cellCount())
        mesh.addExportData(key, d)

    return mesh


def mergeMeshes(meshlist):
    """Merge several meshes into one new mesh and return the new mesh.

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
        raise Exception("argument meshlist is no list")

    if len(meshlist) < 2:
        raise Exception("to few meshes in meshlist")

    mesh = meshlist[0]

    for m in range(1, len(meshlist)):
        mesh = merge2Meshes(mesh, meshlist[m])

    return mesh


def createParaMesh(*args, **kwargs):
    """Create parameter mesh from list of sensor positions.

    Create parameter mesh from list of sensor positions.

    Parameters
    ----------
    sensors : list of RVector3 objects
        Sensor positions. Must be sorted and unique in positive x direction.
        Depth need to be y-coordinate.
    paraDX : float
        Relativ distance for refinement nodes between two electrodes (1=none),
        e.g., 0.5 means 1 additional node between two neighboring electrodes
        e.g., 0.33 means 2 additional equidistant nodes between two electrodes
    paraDepth : float, optional
        Maximum depth for parametric domain, 0 (default) means 0.4 * maximum
        sensor range.
    paraBoundary : float, optional
        Margin for parameter domain in absolute sensor distances. 2 (default).
    paraMaxCellSize: double [0], optional
        Maximum size for parametric size in m*m (0-no constraint)
    boundaryMaxCellSize: double [0], optional
        Maximum cells size in the boundary region in m*m (0-no constraint)
    boundary : float, optional
        Boundary width to be appended for domain prolongation in absolute
        para domain width.
        Values lover 0 force the boundary to be 4 times para domain width.

    Returns
    -------
    poly: :gimliapi:`GIMLI::Mesh`

    """
    plc = pg.meshtools.createParaMeshPLC(*args, **kwargs)
    kwargs.pop('paraMaxCellSize', 0)
    kwargs.pop('boundaryMaxCellSize', 0)
    mesh = createMesh(plc, **kwargs)
    return mesh


def createParaMesh2dGrid(*args, **kwargs):
    """API change here .. use createParaMesh2DGrid instead."""
    print("createParaMesh2dGrid: API change: pls use createParaMesh2DGrid")
    return createParaMesh2DGrid(*args, **kwargs)


def createParaMesh2DGrid(sensors, paraDX=1, paraDZ=1, paraDepth=0, nLayers=11,
                         boundary=-1, paraBoundary=2, **kwargs):
    """Create a grid style mesh for an inversion parameter mesh.

    Create a grid style mesh for an inversion parameter mesh.
    Return parameter grid for a given list of sensor positions.
    Uses and forwards arguments to
    :py:mod:`pygimli.meshtools.appendTriangleBoundary`.

    Parameters
    ----------
    sensors : list of RVector3 objects or data container with sensorPositions
        Sensor positions. Must be sorted in positive x direction
    paraDX : float, optional
        Horizontal distance between sensors, relative regarding sensor
        distance. Value must be greater than 0 otherwise 1 is assumed.
    paraDZ : float, optional
        Vertical distance to the first depth layer, relative regarding sensor
        distance. Value must be greater than 0 otherwise 1 is assumed.
    paraDepth : float, optional
        Maximum depth for parametric domain, 0 (default) means 0.4 * maximum
        sensor range.
    nLayers : int, optional [11]
        Number of depth layers.
    boundary : int, optional [-1]
        Boundary width to be appended for domain prolongation in absolute
        para domain width.
        Values lower than 0 force the boundary to be 4 times para domain width.
    paraBoundary : int, optional [2]
        Offset to the parameter domain boundary in absolute sensor spacing.

    Returns
    -------
    mesh: :gimliapi:`GIMLI::Mesh`

    Examples
    --------
    >>> import pygimli as pg
    >>> import matplotlib.pyplot as plt
    >>>
    >>> from pygimli.meshtools import createParaMesh2DGrid
    >>> mesh = createParaMesh2DGrid(sensors=pg.RVector(range(10)),
    ...                             boundary=1, paraDX=1,
    ...                             paraDZ=1, paraDepth=5)
    >>> ax, _ = pg.show(mesh, mesh.cellMarkers(), alpha=0.3, cmap="summer",
    ...                 hold=True)
    >>> ax, _ = pg.show(mesh, ax=ax)
    """
    mesh = pg.Mesh(2)

    # maybe separate x y z and sort
    if isinstance(sensors, np.ndarray) or isinstance(sensors, pg.RVector):
        sensors = [pg.RVector3(s, 0) for s in sensors]

    if isinstance(sensors, pg.DataContainer):
        sensors = sensors.sensorPositions()

    sensorX = pg.x(sensors)

    eSpacing = abs(sensorX[1] - sensorX[0])

    xmin = min(sensorX) - paraBoundary * eSpacing
    xmax = max(sensorX) + paraBoundary * eSpacing

    if paraDX == 0:
        paraDX = 1.
    if paraDZ == 0:
        paraDZ = 1.

    dx = paraDX
    dz = paraDZ
    if eSpacing > 0:
        dx = eSpacing * paraDX
        # dz = eSpacing * paraDZ  # not really making sense

    if paraDepth == 0:
        paraDepth = 0.4 * (xmax - xmin)

    # print(xmin, xmax, dx)
    x = pg.utils.grange(xmin, xmax, dx=dx)

    y = -pg.increasingRange(dz, paraDepth, nLayers)

    mesh.createGrid(x, y)
    mesh.setCellMarkers([2] * mesh.cellCount())

    paraXLimits = [xmin, xmax]
    #    paraYLimits = [min(y), max(y)]  # not used

    if boundary < 0:
        boundary = abs((paraXLimits[1] - paraXLimits[0]) * 4.0)

    mesh = pg.meshtools.appendTriangleBoundary(
        mesh, xbound=boundary, ybound=boundary, marker=1, **kwargs)

    return mesh

if __name__ == "__main__":
    pass
