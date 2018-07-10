#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""TODO DOCUMENT ME"""

from copy import deepcopy

import numpy as np

import pygimli as pg

from pygimli.utils import unique


def parseArgToArray(arg, nDof, mesh=None, userData=None):
    """
    Parse array related arguments to create a valid value array.

    Parameters
    ----------
    arg : float | int | iterable | callable
        The target array value that will be converted to an array.

        If arg is a callable with it must fulfill:

        :: arg(cell|node|boundary, userData=None)

        Where MeshEntity is one of
        :gimliapi:`GIMLI::Cell` ,
        :gimliapi:`GIMLI::Node` or
        :gimliapi:`GIMLI::Boundary`
        depeding on nDof, where nDof is mesh.cellCount(),
        mesh.nodeCount() or mesh.boundaryCount(),
        respectively.

    nDof : int | [int]
        Desired array size.

    mesh : :gimliapi:`GIMLI::Mesh`
        Used if arg is callable

    userData : class
        Used if arg is callable

    Returns
    -------

    ret : :gimliapi:`GIMLI::RVector`
        Array of desired length filled with the appropriate values.
    """
    nDofs = nDof

    if not hasattr(nDof, '__len__'):
        nDofs = [nDof]

    try:
        return pg.RVector(nDofs[0], float(arg))
    except BaseException as _:
        pass

    if hasattr(arg, '__len__'):
        if isinstance(arg, np.ndarray):
            if len(arg) == nDofs[0]:
                return arg
            else:
                raise Exception('Given array does not have requested (' +
                                str(nDof) + ') size (' +
                                str(len(arg)) + ')')

        for n in nDofs:
            if len(arg) == n:
                return arg

        try:
            # [marker, val] || [[marker, val]]
            return parseMapToCellArray(arg, mesh)
        except:
            raise Exception("Array 'arg' has the wrong size: " +
                            str(len(arg)) + " != " + str(nDof))
    elif hasattr(arg, '__call__'):
        ret = pg.RVector(nDofs[0], 0.0)

        if not mesh:
            raise Exception("Please provide a mesh for the callable"
                            "argument to parse ")

        if nDofs[0] == mesh.nodeCount():
            for n in mesh.nodes():
                if userData:
                    ret[n.id()] = arg(node=n, userData=userData)
                else:
                    ret[n.id()] = arg(node=n)
        elif nDofs[0] == mesh.cellCount():
            for c in mesh.cells():
                if userData:
                    ret[c.id()] = arg(cell=c, userData=userData)
                else:
                    ret[c.id()] = arg(cell=c)
        elif nDofs[0] == mesh.boundaryCount():
            for b in mesh.boundaries():
                if userData:
                    ret[b.id()] = arg(boundary=b, userData=userData)
                else:
                    ret[b.id()] = arg(boundary=b)
        else:
            raise Exception("Cannot parse callable argument " + str(ndof) +
                            " nodes: " + str(mesh.nodeCount()) +
                            " cells: " + str(mesh.cellCount()))

        return ret
    raise Exception("Cannot parse argument type " + str(type(arg)))


def generateBoundaryValue(boundary, arg, time=0.0, userData=None):
    """
    Generate a value for the given Boundary.

    Parameters
    ----------
    boundary : :gimliapi:`GIMLI::Boundary` or list of ..
        The related boundary.

    arg : convertible | iterable | callable or list of ..

        - convertible into float
        - iterable of minimum length = boundary.id()
        - callable generator function

        If arg is a callable it must fulfill:

        :: arg(boundary=:gimliapi:`GIMLI::Boundary`, time=0.0, userData=None)

        The callable function arg have to return appropriate values
        for all nodes of the boundary or one scalar value for all nodes.

    Returns
    -------
    val : float or list of ..

    Examples
    --------
    >>> import pygimli as pg
    >>>
    >>> # def uBoundary(boundary):
    >>> #    u = []
    >>> #    for n in boundary.nodes():
    >>> #        u.append(n[0] + n[1])
    >>> #    return u
    """
    if hasattr(boundary, '__len__'):
        vals = np.zeros(len(boundary))
        for i, b in enumerate(boundary):
            vals[i] = generateBoundaryValue(b, arg[i], time, userData)
        return vals
    val = 0.

    #print(arg, callable(arg))
    if hasattr(arg, '__call__'):
        kwargs = dict()
        kwargs['boundary'] = boundary

        if time != 0.0 and time is not None:
            kwargs['time'] = time
        if userData:
            kwargs['userData'] = userData
        try:
            val = arg(**kwargs)
        except Exception as e:
            print(arg, "(", kwargs, ")")
            pg.logger.critical(e)
            raise Exception("Wrong arguments for callback function.")

    elif hasattr(arg, '__len__'):
        if callable(arg[0]):
            kwargs = arg[1]
            if time != 0.0 and time is not None:
                kwargs['time'] = time
            val = arg[0](boundary=boundary, **kwargs)
        else:
            val = generateBoundaryValue(boundary, arg[boundary.id()], userData)
    else:
        try:
            val = float(arg)
        except ValueError:
            raise arg
    return val


def parseArgPairToBoundaryArray(pair, mesh):
    """
    Parse boundary related pair argument to
    [ :gimliapi:`GIMLI::Boundary`, value|callable ] list.

    Parameters
    ----------

    pair : tuple

        - [marker, arg]
        - [[marker, ...], arg]
        - [boundary, arg]
        - [[boundary,...], arg]
        - [node, arg]
        - [marker, callable, *kwargs]
        - [[marker, ...], callable, *kwargs]

        arg will be parsed by
        :py:mod:`pygimli.solver.solver.generateBoundaryValue`
        and distributed to each boundary.
        Callable functions will be executed at run time.

    mesh : :gimliapi:`GIMLI::Mesh`
        Used to find boundaries by marker.

    Returns
    -------

    boundaries : list()
        [ :gimliapi:`GIMLI::Boundary`, value|callable ]
    """
    boundaries = []
    bounds = []

    #print('*'*30, pair)
    if isinstance(pair[0], int):
        bounds = mesh.findBoundaryByMarker(pair[0])
    elif isinstance(pair[0], list):
        # [[,,..], ]
        for b in pair[0]:
            for bi in mesh.boundaries(pg.find(mesh.boundaryMarkers() == b)):
                bounds.append(bi)

    elif isinstance(pair[0], pg.stdVectorBounds):
        bounds = pair[0]
    elif isinstance(pair[0], pg.Boundary):
        boundaries.append(pair)
        return boundaries
    elif isinstance(pair[0], pg.Node):
        boundaries.append(pair)
        return boundaries

    for b in bounds:
        val = None
        #print('-'*50)
        #print(b, pair[1], callable(pair[1]))
        #print('+'*50)
        if callable(pair[1]):
            # don't execute the callable here
            # we want to call them at runtime
            if len(pair) > 2:
                val = pair[1:]
            else:
                val = pair[1]
        else:
            val = generateBoundaryValue(b, pair[1])
        boundaries.append([b, val])

    #print('#'*30)
    return boundaries


def parseArgToBoundaries(args, mesh):
    """
    Parse boundary related arguments to create a valid boundary value list:
    [ :gimliapi:`GIMLI::Boundary`, value|callable ]

    Parameters
    ----------
    args : callable | pair | [pair, ...]
        If args is just a callable than every boundary will be evaluated
        at runtime with this function as args(boundary).
        Else see :py:mod:`pygimli.solver.solver.parseArgPairToBoundaryArray`

    mesh : :gimliapi:`GIMLI::Mesh`
        Used to find boundaries by marker

    Returns
    -------
    boundaries : list()
        [ :gimliapi:`GIMLI::Boundary`, value|callable ]

    Examples
    --------
    >>> # no need to import matplotlib. pygimli's show does
    >>> import pygimli as pg
    >>> mesh = pg.meshtools.createWorld([0, 0], [1, -1], worldMarker=0)
    >>> ax, _ = pg.show(mesh, boundaryMarker=True)
    >>> # all edges with marker 1 get value 1.0
    >>> b = pg.solver.parseArgToBoundaries([1, 1.0], mesh)
    >>> print(len(b))
    1
    >>> # same as above with marker 2 get value 2
    >>> b = pg.solver.parseArgToBoundaries([[1, 1.0], [2, 2.0]], mesh)
    >>> print(len(b))
    2
    >>> # same as above with marker 3 get value 3
    >>> b = pg.solver.parseArgToBoundaries([[1, 1.], [2, 2.], [3, 3.]], mesh)
    >>> print(len(b))
    3
    >>> # edges with marker 1 and 2 get value 1
    >>> b = pg.solver.parseArgToBoundaries([[1, 2], 1.0], mesh)
    >>> print(len(b))
    2
    >>> b = pg.solver.parseArgToBoundaries([[1, 2, 3], 1.0], mesh)
    >>> print(len(b))
    3
    >>> b = pg.solver.parseArgToBoundaries([[[1, 2, 3], 1.0], [4, 4.0]], mesh)
    >>> print(len(b))
    4
    >>> b = pg.solver.parseArgToBoundaries([mesh.node(0), 0.0], mesh)
    >>> print(len(b))
    1
    >>> def bCall(boundary):
    ...     u = []
    ...     for i, n in enumerate(boundary.nodes()):
    ...         u.append(i)
    ...     return u
    >>> b = pg.solver.parseArgToBoundaries([1, bCall], mesh)
    >>> print(len(b),b[0][1](b[0][0]))
    1 [0, 1]
    >>> def bCall(boundary, a1, a2):
    ...     return a1 + a2
    >>> b = pg.solver.parseArgToBoundaries([1, bCall, {'a1':2, 'a2':3}], mesh)
    >>> print(len(b), b[0][1][0](b[0][0], **b[0][1][1]))
    1 5
    >>> b = pg.solver.parseArgToBoundaries([[1, bCall, {'a1':1, 'a2':2}],
    ...                                     [2, bCall, {'a1':3, 'a2':4}]], mesh)
    >>> print(len(b), b[0][1][0](b[0][0], **b[0][1][1]))
    2 3
    >>> b = pg.solver.parseArgToBoundaries([[1,2], bCall, {'a1':4, 'a2':5}], mesh)
    >>> print(len(b), b[1][1][0](b[1][0], **b[1][1][1]))
    2 9
    >>> pg.wait()
    """
    boundaries = []
    if isinstance(args, list):
        #print('!'*100)
        #print(args)
        if len(args) == 2:
            # if isinstance(args[0], list):
                # print('~'*10, '[[,],]')
                # boundaries += parseArgPairToBoundaryArray(args, mesh)
            # else:
            try:
                    # [[,], [,]]
                if len(args[0]) == 2 and len(args[1]) == 2:
                    #print('~'*10, '[[,],[,]]')
                    boundaries += parseArgPairToBoundaryArray(args[0], mesh)
                    boundaries += parseArgPairToBoundaryArray(args[1], mesh)
                elif len(args[0]) == 3 and callable(args[0][1]):
                    #print('~'*10, '[[marker,callable,*kwrags],[,]]')
                    boundaries += parseArgPairToBoundaryArray(args[0], mesh)
                    boundaries += parseArgPairToBoundaryArray(args[1], mesh)
                else:
                    #print('~'*10, '??[[,]]')
                    boundaries += parseArgPairToBoundaryArray(args, mesh)
            except BaseException as _:
                # [,]
                # print('~'*10, '[,]')
                boundaries += parseArgPairToBoundaryArray(args, mesh)
        elif len(args) == 3 and callable(args[1]):
            # print('~'*10, '[[,], callable, kwargs)
            boundaries += parseArgPairToBoundaryArray(args, mesh)
        else:
            # print('~'*10, '[[,], [,], ...]')
            # [[,], [,], ...]
            for a in args:
                boundaries += parseArgPairToBoundaryArray(a, mesh)

    elif hasattr(args, '__call__'):
        for b in mesh.boundaries():
            if not b.leftCell() or not b.rightCell():
                # if args(b) is not None:
                boundaries.append([b, args])

    elif isinstance(args, float) or isinstance(args, int):
        for b in mesh.boundaries():
            if not b.leftCell() or not b.rightCell():
                # if args(b) is not None:
                boundaries.append([b, args])

    else:
        raise Exception('cannot interprete boundary token', args)

    return boundaries


def parseMapToCellArray(attributeMap, mesh, default=0.0):
    """
    Parse a value map to cell attributes.

    A map should consist of pairs of marker and value.
    A marker is an integer and corresponds to the cell.marker().

    Parameters
    ----------
    mesh : :gimliapi:`GIMLI::Mesh`
        For each cell of mesh a value will be returned.

    attributeMap : list | dict
        List of pairs [marker, value] ] || [[marker, value]],
        or dictionary with marker keys

    default : float [0.0]
        Fill all unmapped atts to this default.

    Returns
    -------
    atts : array
        Array of length mesh.cellCount()
    """

    atts = pg.RVector(mesh.cellCount(), default)

    if isinstance(attributeMap, dict):
        raise Exception("Please implement me!")
    elif hasattr(attributeMap, '__len__'):
        if not hasattr(attributeMap[0], '__len__'):
            # assuming [marker, value]
            attributeMap = [attributeMap]

        for pair in attributeMap:
            if hasattr(pair, '__len__'):
                idx = pg.find(mesh.cellMarkers() == pair[0])
                if len(idx) == 0:
                    print("Warning! parseMapToCellArray: cannot find marker " +
                          str(pair[0]) + " within mesh.")
                else:
                    #print('---------------------')
                    #print(atts, idx, pair[1], type(pair[1]), float(pair[1]))
                    if isinstance(pair[1], np.complex):
                        #print('+++++++++++++++++')
                        if not isinstance(atts, pg.CVector):
                            atts = pg.toComplex(atts)
                        atts[idx] = pair[1]
                    else:
                        atts[idx] = float(pair[1])
            else:
                raise Exception("Please provide a list of [int, value] pairs" +
                                str(pair))
    else:
        print("attributeMap: ", attributeMap)
        raise Exception("Cannot interpret attributeMap!")

    return atts


def grad(mesh, u, r=None):
    r"""
    Return the discrete interpolated gradient :math:`\mathbf{v}`
    for a given scalar field :math:`\mathbf{u}`.

    .. math::
        \mathbf{v}(\mathbf{r}_{\mathcal{C}})
        &=
        \nabla u(\mathbf{r}_{\mathcal{N}})
        \\
        (\mathbf{v_x}(\mathbf{r}_{\mathcal{C}}),
         \mathbf{v_y}(\mathbf{r}_{\mathcal{C}}),
         \mathbf{v_z}(\mathbf{r}_{\mathcal{C}}))^{\text{T}}
        &=
        \left(\frac{\partial u}{\partial x},
         \frac{\partial u}{\partial y},
         \frac{\partial u}{\partial z}\right)^{\text{T}}

    With :math:`\mathcal{N}=\cup_{i=0}^{N} \text{Node}_i`,
    :math:`\mathcal{C}=\cup_{j=0}^{M} \text{Cell}_j`,
    :math:`\mathbf{u}=\{u(\mathbf{r}_{i})\}\in I\!R` and
    :math:`\mathbf{r}_{i} = (x_i, y_i, z_i)^{\text{T}}`

    The discrete scalar field
    :math:`\mathbf{u}(\mathbf{r}_{\mathcal{N}})`
    need to be defined for each node position :math:`\mathbf{r}_{\mathcal{N}}`.
    The resulting vector field :math:`\mathbf{v}(\mathbf{r}_{\mathcal{C}})`
    is defined for each cell center position :math:`\mathbf{r}_{\mathcal{C}}`.
    If you need other positions than the cell center,
    provide an appropriate array of coordinates :math:`\mathbf{r}`.

    See also
    --------
    GIMLI::Mesh::cellDataToBoundaryGradient
    GIMLI::Mesh::boundaryDataToCellGradient

    Parameters
    ----------
    mesh : :gimliapi:`GIMLI::Mesh`
        Discretization base, interpolation will be performed via finite element
        base shape functions.

    u : array | callable
        Scalar field per mesh node position or an appropriate
        callable([[x,y,z]])

    r : ndarray((M, 3)) [mesh.cellCenter()]
        Alternative target coordinates :math:`\mathbf{r} for the resulting
        gradient field. i.e., the positions where the vector field is defined.
        Default are all cell centers.

    Returns
    -------
    v : ndarray((M, 3))
        Resulting vector field defined on
        :math:`\mathbf{v}(\mathbf{r}_{\mathcal{C}})`.
        M is number of cells or length of given alternative coordinates r.

    Examples
    --------
    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> import pygimli as pg
    >>> fig, ax = plt.subplots()
    >>> mesh = pg.createGrid(x=np.linspace(0, 1, 20), y=np.linspace(0, 1, 20))
    >>> u = lambda p: pg.x(p)**2 * pg.y(p)
    >>> _ = pg.show(mesh, u(mesh.nodeCenters()), ax=ax)
    >>> _ = pg.show(mesh, [2.*pg.y(mesh.cellCenters())*pg.x(mesh.cellCenters()),
    ...             pg.x(mesh.cellCenters())**2], ax=ax)
    >>> _ = pg.show(mesh, pg.solver.grad(mesh, u), ax=ax, color='w',
    ...             linewidth=0.4)
    >>> plt.show()
    """

    if r is None:
        r = mesh.cellCenters()

    uv = u
    if callable(u) and not isinstance(u, pg.RVector):
        uv = u(mesh.nodeCenters())

    if len(uv) == mesh.cellCount():
        uv = pg.meshtools.cellDataToNodeData(mesh, uv)

    v = np.ndarray((len(r), 3))

    for i, _ in enumerate(v):
        c = mesh.findCell(r[i])
        if c:
            v[i] = c.grad(r[i], uv)

    return v


def div(mesh, v):
    r"""Return the discrete interpolated divergence field.

        Return the discrete interpolated divergence field. :math:`\mathbf{u}`
        for each cell for a given vector field :math:`\mathbf{v}`.
        First order integration via boundary center.

    .. math::
        d(cells) & = \nabla\cdot\vec{v} \\
        d(c_i) & = \sum_{j=0}^{N_B}\vec{v}_{B_j} \cdot \vec{n}_{B_j}

    Parameters
    ----------
    mesh : :gimliapi:`GIMLI::Mesh`
        Discretization base, interpolation will be performed via finite element
        base shape functions.

    V : array(N,3) | R3Vector
        Vector field at cell centers or boundary centers

    Returns
    -------
    d : array(M)
        Array of divergence values for each cell in the given mesh.

    Examples
    --------
    >>> import pygimli as pg
    >>> import numpy as np
    >>> v = lambda p: p
    >>> mesh = pg.createGrid(x=np.linspace(0, 1, 4))
    >>> print(pg.round(pg.solver.div(mesh, v(mesh.boundaryCenters())), 1e-5))
    <class 'pygimli.core._pygimli_.RVector'> 3 [1.0, 1.0, 1.0]
    >>> print(pg.round(pg.solver.div(mesh, v(mesh.cellCenters())), 1e-5))
    <class 'pygimli.core._pygimli_.RVector'> 3 [0.5, 1.0, 0.5]
    >>> mesh = pg.createGrid(x=np.linspace(0, 1, 4),
    ...                      y=np.linspace(0, 1, 4))
    >>> print(pg.round(pg.solver.div(mesh, v(mesh.boundaryCenters())), 1e-5))
    <class 'pygimli.core._pygimli_.RVector'> 9 [2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0]
    >>> divCells = pg.solver.div(mesh, v(mesh.cellCenters()))
    >>> # divergence from boundary values are exact where the divergence from
    >>> # interpolated cell center values are wrong due to boundary interp
    >>> print(sum(divCells))
    12.0
    >>> mesh = pg.createGrid(x=np.linspace(0, 1, 4),
    ...                      y=np.linspace(0, 1, 4),
    ...                      z=np.linspace(0, 1, 4))
    >>> print(sum(pg.solver.div(mesh, v(mesh.boundaryCenters()))))
    81.0
    >>> divCells = pg.solver.div(mesh, v(mesh.cellCenters()))
    >>> print(sum(divCells))
    54.0
    """
    mesh.createNeighbourInfos()
    d = None
    if hasattr(v, '__len__'):
        if len(v) == mesh.boundaryCount():
            d = mesh.divergence(v)
        elif len(v) == mesh.nodeCount():
            d = mesh.divergence(pg.meshtools.nodeDataToBoundaryData(mesh, v))
        elif len(v) == mesh.cellCount():
            CtB = mesh.cellToBoundaryInterpolation()
            d = mesh.divergence(np.array([CtB*pg.x(v),
                                          CtB*pg.y(v),
                                          CtB*pg.z(v)]).T)
        else:
            print(len(v), mesh)
            raise BaseException("implement me")
    elif callable(v):
        raise BaseException("implement me")

    return d


def divergence(mesh, F=None, normMap=None, order=1):
    """Divergence for callable function F((x,y,z)).

    MOVE THIS to a better place

    Divergence for callable function F((x,y,z)). Return sum div over boundary.

    Parameters
    ----------

    Returns
    -------
    """
    if F is None:
        F = lambda r: r

    diverg = 0
    directionCheck = False

    if mesh.cellCount() > 0:
        directionCheck = True

    bNorms = None
    if normMap is not None:
        bNorms = np.zeros((mesh.boundaryCount(), 2))
        for pair in normMap:
            bounds = mesh.findBoundaryByMarker(pair[0])
            for b in bounds:
                bN = [0.0, 0.0]
                if not isinstance(pair[1][0], str):
                    bN[0] = pair[1][0]
                if not isinstance(pair[1][1], str):
                    bN[1] = pair[1][1]

                bNorms[b.id()] = bN

    for b in mesh.boundaries():

        if directionCheck:
            if b.leftCell() is None and b.rightCell() is None:
                # print(b.id(), b.leftCell(), b.rightCell())
                sw = pg.Stopwatch(True)
                mesh.createNeighbourInfos()
                print("NeighbourInfos()", sw.duration(True))
                # return gauss(grid, F)

            # don't calc for inner boundaries
            if b.leftCell() is not None and b.rightCell() is not None:
                continue

        tmpdiv = 0
        shape = b.shape()

        if order == 1:
            if bNorms is not None:
                tmpdiv = shape.norm().dot(bNorms[b.id()]) * shape.domainSize()
            else:
                tmpdiv = shape.norm().dot(
                    F(shape.center())) * shape.domainSize()
        else:
            weights = pg.IntegrationRules.instance().weights(shape, order)
            abscissa = pg.IntegrationRules.instance().abscissa(shape, order)

            for i, p in enumerate(abscissa):
                rPos = shape.xyz(p)
                tmpdiv += shape.norm().dot(F(rPos)) * \
                    weights[i] * shape.domainSize()

        if directionCheck and b.leftCell() is None:
            tmpdiv *= -1
            # raise Exception("invalid mesh: left is None .. every
            # boundary need leftCell")

        diverg += tmpdiv

    return diverg


def triDiagToeplitz(dom, a, l, r, start=0, end=-1):
    """Create tri-diagonal Toeplitz matrix."""
    A = pg.RSparseMapMatrix(dom, dom)

    if end == -1:
        end = dom

    for i in range(start, end):
        A.addVal(i, i, a)
        if i > start:
            A.addVal(i, i - 1, l)

        if i < end - 1:
            A.addVal(i, i + 1, r)
    return A


def identity(dom, start=0, end=-1, scale=1):
    """Create identity matrix."""
    A = pg.RSparseMapMatrix(dom, dom)

    if end == -1:
        end = dom

    for i in range(start, end):
        if hasattr(scale, '__len__'):
            A.addVal(i, i, scale[i])
        else:
            A.addVal(i, i, scale)
    return A


def showSparseMatrix(A, full=False):
    """Show the content of a sparse matrix.

    Parameters
    ----------
    A : :gimliapi:`GIMLI::SparseMatrix` | :gimliapi:`GIMLI::SparseMapMatrix`
        Matrix to be shown.
    full : bool [False]
        Show as dense matrix.
    """
    S = A
    if isinstance(A, pg.RSparseMapMatrix):
        return showSparseMatrix(pg.SparseMatrix(A), full)
    else:
        rows = S.vecRowIdx()
        cols = S.vecColPtr()
        vals = S.vecVals()

        if full:
            Sd = pg.RMatrix(S.rows(), S.cols())

        for i in range(S.rows()):
            for j in range(cols[i], cols[i + 1]):
                if full:
                    Sd[i].setVal(vals[j], rows[j])
                else:
                    print(i, rows[j], vals[j])

        if full:
            print(np.array(Sd))


def linsolve(A, b, verbose=False):
    """
    DEPRECATED wrong name style
    """
    pg.deprecated('linsolve', 'linSolve')
    return linSolve(A, b, verbose)

def linSolve(A, b, verbose=False):
    r"""Direct solution after :math:`\textbf{x}` using core LinSolver.

    .. math::
        \textbf{A}\textbf{x} = \textbf{b}

    If :math:`\textbf{A}` is symmetric, sparse and positive definite.

    Parameters
    ----------
    A : :gimliapi:`GIMLI::RSparseMatrix` | :gimliapi:`GIMLI::RSparseMapMatrix`|
        numpy.array
        System matrix. Need to be symmetric, sparse and positive definite.

    b : iterable array
        Right hand side of the equation.

    verbose : bool [False]
        Be verbose.

    Returns
    -------

    x : :gimliapi:`GIMLI::RVector`
        Solution vector
    """
    x = pg.RVector(len(b), .0)

    if isinstance(A, pg.RSparseMatrix):
        solver = pg.LinSolver(A, verbose=verbose)
        solver.solve(b, x)
        return x
    elif isinstance(A, pg.RSparseMapMatrix):
        S = pg.RSparseMatrix(A)
        solver = pg.LinSolver(S, verbose=verbose)
        solver.solve(b, x)
        return x
    elif isinstance(A, np.ndarray):
        return np.linalg.solve(A, b)

    import scipy.sparse
    #pg.optImport('scipy.sparse')
    #print(type(A))
    if type(A) == scipy.sparse.csr.csr_matrix:
        return scipy.sparse.linalg.spsolve(A, b)
    else:
        raise StandardException("Don't know how to lineare solve a system"
                                " with matrixtype:" + type(A))

    return x

def assembleForceVector(mesh, f, userData=None):
    """
    DEPRECATED use assembleLoadVector instead
    """
    return assembleLoadVector(mesh, f, userData)

def assembleLoadVector(mesh, f, userData=None):
    """Create right hand side vector based on the given mesh and load
    or force values.

    Create right hand side vector based on the given mesh and load or force
    values.

    Parameters
    ----------
    f: float, array, callable(cell, [userData]), [...]

        Force Values
        float -> ones(mesh.nodeCount()) * vals,
        for each node [0 .. mesh.nodeCount()]
        for each cell [0 .. mesh.cellCount()]
    """
    if isinstance(f, list) or hasattr(f, 'ndim'):
        if isinstance(f, list):

            rhs = np.zeros((len(f), mesh.nodeCount()))
            for i, fi in enumerate(f):
                userData['i'] = i
                rhs[i] = assembleForceVector(mesh, fi, userData)

            return rhs

        elif f.ndim == 2:
            # assume rhs [n, nNodes] array is already a valid
            return f

    rhs = pg.RVector(mesh.nodeCount(), 0)

    if hasattr(f, '__call__') and not isinstance(f, pg.RVector):
        for c in mesh.cells():
            if userData is not None:
                f(c, rhs, userData)
            else:
                f(c, rhs)
    else:
        fArray = None
        if hasattr(f, '__len__'):
            if len(f) == mesh.cellCount() or len(f) == mesh.nodeCount():
                fArray = f

        if fArray is None:
            fArray = parseArgToArray(f, mesh.cellCount(), mesh, userData)

        if len(fArray) == mesh.cellCount():
            b_l = pg.ElementMatrix()

            for c in mesh.cells():
                if fArray[c.id()] != 0.0:
                    b_l.u(c)
                    rhs.add(b_l, fArray[c.id()])

#            print("test reference solution:")
#            rhsRef = pg.RVector(mesh.nodeCount(), 0)
#            for c in mesh.cells():
#                b_l.u(c)
#                for i, idx in enumerate(b_l.idx()):
#                    rhsRef[idx] += b_l.row(0)[i] * fArray[c.id()]
#            np.testing.assert_allclose(rhs, rhsRef)
#            print("Remove revtest in assembleForceVector after check")

        elif len(fArray) == mesh.nodeCount():
            fA = pg.RVector(fArray)
            b_l = pg.ElementMatrix()
            for c in mesh.cells():
                b_l.u(c)
                # rhs.addVal(b_l.row(0) * fArray[b_l.idx()], b_l.idx())
                rhs.add(b_l, fA)
            # print("test reference solution:")
            # rhsRef = pg.RVector(mesh.nodeCount(), 0)
            # for c in mesh.cells():
            #     b_l.u(c)
            #     for i, idx in enumerate(b_l.idx()):
            #         rhsRef[idx] += b_l.row(0)[i] * fArray[idx]
            # np.testing.assert_allclose(rhs, rhsRef)
            # print("Remove revtest in assembleForceVector after check")

            # rhs = pg.RVector(fArray)
        else:
            raise Exception("Forcevector have the wrong size: " +
                            str(len(fArray)))

    return rhs


def assembleNeumannBC(S, boundaryPairs, rhs=None, time=0.0, userData=None):
    r"""Apply Neumann condition to the system matrix S.

    Apply Neumann condition to the system matrix S.
    .. math::
        \frac{\partial u(\arr{r}, t)}{\partial\textbf{n}}
        = \textbf{n}\grad u(\arr{r}, t) = g \quad\text{with}\quad\arr{r}
        \quad\text{on}\quad \partial\Omega

    Parameters
    ----------

    S : :gimliapi:`GIMLI::RSparseMatrix`
        System matrix of the system equation.

    boundaryPair : list()
        List of pairs [ :gimliapi:`GIMLI::Boundary`, g ].
        The value g will assigned to the nodes of the boundaries.
        Later assignment overwrites prior.

        :math:`g` need to be a scalar value (float or int) or
        a value generator function that will be executed at run time.
        See :py:mod:`pygimli.solver.solver.parseArgToBoundaries`

        See tutorial section for an example,
        e.g., Modeling with Boundary Conditions

    rhs :
        TODO DOCU
    time : float
        Will be forwarded to value generator.

    userData : class
        Will be forwarded to value generator.
    """
    Se = pg.ElementMatrix()

    if not hasattr(boundaryPairs, '__getitem__'):
        raise BaseException("Boundary pairs need to be a list of "
                            "[boundary, value]")

    for pair in boundaryPairs:
        boundary = pair[0]

        val = pair[1]

        g = generateBoundaryValue(boundary, val, time, userData)

        if g is not 0.0 and g is not None:

            #Se.u(boundary)
            #if rhs is not None:
                #rhs.add(Se, g)

            # check for robyn condition
            Se.u2(boundary)
            Se *= g

def _assembleUDirichlet(S, rhs, uDirIndex, uDirchlet):
    """This should be moved directly into the core"""

    if rhs is not None:
        udirTmp = pg.RVector(S.rows(), 0.0)
        udirTmp.setVal(uDirchlet, uDirIndex)
        rhs -= S * udirTmp

    for i in uDirIndex:
        S.cleanRow(i)
        S.cleanCol(i)
        S.setVal(i, i, 1.0)

    if rhs is not None:
        rhs[uDirIndex] = uDirchlet
        #rhs.setVal(uDirchlet, uDirIndex)


def assembleDirichletBC(S, boundaryPairs, rhs=None, time=0.0, userData=None,
                        nodePairs=None):
    r"""Apply Dirichlet boundary condition.

    Apply Dirichlet boundary condition to the system matrix S and rhs vector.
    The right hand side values for h can be given for each boundary
    element individually by setting proper boundary pair arguments.

    .. math::
        u(\textbf{r}, t) = h
        \quad\text{for}\quad\textbf{r}\quad\text{on}\quad\delta\Omega=\Gamma_{\text{Dirichlet}}

    Parameters
    ----------
    S : :gimliapi:`GIMLI::RSparseMatrix`
        System matrix of the system equation.

    boundaryPair : list()
        List of pairs [ :gimliapi:`GIMLI::Boundary`, h ].
        The value :math:`h` will assigned to the nodes of the boundaries.
        Later assignment overwrites prior.

        :math:`h` need to be a scalar value (float or int) or
        a value generator function that will be executed at runtime.
        See :py:mod:`pygimli.solver.solver.parseArgToBoundaries`
        and :ref:`tut:modelling_bc` for example syntax,


    nodePairs : list()
        List of pairs [ nodeID, uD ].
        The value uD will assigned to the nodes given there ids.
        This node value settings will overwrite any prior settings due to
        boundaryPair.

    rhs : :gimliapi:`GIMLI::RVector`
        Right hand side vector of the system equation will bet set to
        :math:`u_{\text{D}}`

    time : float
        Will be forwarded to value generator.

    userData : class
        Will be forwarded to value generator.

    """

    if not hasattr(boundaryPairs, '__getitem__'):
        raise BaseException("Boundary pairs need to be a list of "
                            "[boundary, value]")
    uDirNodes = []
    uDirVal = dict()

    for pair in boundaryPairs:
        boundary = pair[0]
        val = pair[1]
        uD = generateBoundaryValue(boundary, val, time, userData)

        if uD is not None:
            if isinstance(boundary, pg.Node):
                n = boundary
                uDirNodes.append(n)
                uDirVal[n.id()] = uD
            else:
                for i, n in enumerate(boundary.nodes()):
                    uDirNodes.append(n)
                    if hasattr(uD, '__iter__'):
                        uDirVal[n.id()] = uD[i]
                    else:
                        uDirVal[n.id()] = uD

    if len(uDirNodes) == 0 and nodePairs is None:
        return

    uniqueNodes = unique(uDirNodes)

    uDirchlet = []
    uDirIndex = []

    for i, n in enumerate(uniqueNodes):
        uDirIndex.append(n.id())
        uDirchlet.append(uDirVal[n.id()])

    if nodePairs is not None:
        #print("nodePairs", nodePairs)

        if len(nodePairs) == 2 and type(nodePairs[0]) == int :
            # assume a single Node [NodeId, val]
            nodePairs = [nodePairs]

        for i, [n, val] in enumerate(nodePairs):
            uDirIndex.append(n)
            if hasattr(val, '__call__'):
                raise("callabe node pairs need to be implement.")
            uDirchlet.append(val)

    _assembleUDirichlet(S, rhs, uDirIndex, uDirchlet)


def assembleNeumannBC(S, boundaryPairs, rhs=None, time=0.0, userData=None):
    r"""Apply Neumann condition to the system matrix S.

    Apply Neumann condition to the system matrix S.
    The right hand side values for g can be given for each boundary
    element individually by setting proper boundary pair arguments.

    .. math::
        \frac{\partial u(\textbf{r}, t)}{\partial\textbf{n}}
        = \textbf{n}\nabla u(\textbf{r}, t) = g
        \quad\text{for}\quad\textbf{r}\quad\text{on}\quad\delta\Omega=\Gamma_{\text{Neumann}}

    Parameters
    ----------

    S : :gimliapi:`GIMLI::SparseMatrix`
        System matrix of the system equation.

    boundaryPair : list()
        List of pairs [ :gimliapi:`GIMLI::Boundary`, g ].
        The value :math:`g` will assigned to the nodes of the boundaries.
        Later assignment overwrites prior.

        :math:`g` need to be a scalar value (float or int) or
        a value generator function that will be executed at run time.

        See :py:mod:`pygimli.solver.solver.parseArgToBoundaries`
        and :ref:`tut:modelling_bc` for example syntax,

    rhs :
        TODO
    time : float
        Will be forwarded to value generator.

    userData : class
        Will be forwarded to value generator.
    """

    if rhs is None:
        raise BaseException("Neumann Boundary condition needs rhs vector.")

    if not hasattr(boundaryPairs, '__getitem__'):
        raise BaseException("Boundary pairs need to be a list of "
                            "[boundary, value]")

    Se = pg.ElementMatrix()

    for pair in boundaryPairs:
        boundary = pair[0]
        val = pair[1]
        g = generateBoundaryValue(boundary, val, time, userData)

        if g is not 0.0 and g is not None:
            Se.u(boundary)
            rhs.add(Se, g)

def assembleRobinBC(S, boundaryPairs, rhs=None, time=0.0, userData=None):
    r"""Apply Robin boundary condition.

    Apply Robin boundary condition to the system matrix S and rhs vector
    (if needed for b != 1 and g != 0).

    .. math::
        \alpha u(\textbf{r}, t) + \beta\frac{\partial u(\textbf{r}, t)}{\partial\textbf{n}}
        = \gamma
        \quad\text{for}\quad\textbf{r}\quad\text{on}\quad\delta\Omega=\Gamma_{\text{Robin}}\\
        \quad\text{currently only with}\quad \beta = 1 \quad\text{and}\quad \gamma = 0

    TODO
        * b!=1 and g!=0 variable
        * check for b = 0 and move to dirichlet

    Parameters
    ----------

    S : :gimliapi:`GIMLI::SparseMatrix`
        System matrix of the system equation.

    boundaryPair : list()
        List of pairs [ :gimliapi:`GIMLI::Boundary`, :math:`\alpha` ].
        The value :math:`\alpha` will assigned to the nodes of the boundaries.
        Later assignment overwrites prior.

        :math:`\alpha` needs to be a scalar value (float or int) or
        a value generator function that will be executed at run time.
        See :py:mod:`pygimli.solver.solver.parseArgToBoundaries`
        and :ref:`tut:modelling_bc` for example syntax,

    time : float
        Will be forwarded to value generator.

    userData : class
        Will be forwarded to value generator.
    """

    if not hasattr(boundaryPairs, '__getitem__'):
        raise BaseException("Boundary pairs need to be a list of "
                            "[boundary, value]")

    Sp = pg.ElementMatrix()
    Sq = pg.ElementMatrix()

    #if isinstance(rhs, np.ndarray):
        #rhs = pg.RVector(rhs)

    for pair in boundaryPairs:
        boundary = pair[0]
        val = pair[1]
        ### p = alpha / alpha
        p = generateBoundaryValue(boundary, val, time, userData)
        #### p = gamma / alpha
        #p = 20.
        #q = -41.0
        q = None

        if p is not 0.0 and p is not None:
            Sp.u2(boundary)
            S.add(Sp, p)
            #Sp *= p
            #S += Sp

            if q != None:
                Sq.u(boundary)
                rhs.add(Sq, -p*q)

def assembleBC_(bc, mesh, S, rhs, time=None, userData=None):
    r"""Shortcut to apply all boundary conditions.

    This is a helper function for the solver call.
    Shortcut to apply all boundary conditions will only forward to
    appropriate assemble functions.
    """

    ## we can't iterate because we want the following fixed order
    bct = dict(bc)
    if 'Neumann' in bct:
        assembleNeumannBC(S, parseArgToBoundaries(bct.pop('Neumann'), mesh),
                          rhs=rhs, time=time, userData=userData)
    if 'Robin' in bct:
        assembleRobinBC(S, parseArgToBoundaries(bct.pop('Robin'), mesh),
                        rhs=rhs, time=time, userData=userData)
    if 'Dirichlet' in bct:
        assembleDirichletBC(S, parseArgToBoundaries(bct.pop('Dirichlet'), mesh),
                            rhs=rhs, time=time, userData=userData)
    if 'Node' in bct:
        assembleDirichletBC(S, [], nodePairs=bct.pop('Node'),
                            rhs=rhs, time=time, userData=userData)

    if len(bct.keys()) > 0:
        pg.logger.warn("Unknown boundary condition[s]" + \
                       str(bct.keys()) + " will be ignored")

def createLoadVector(mesh, f, userData=None):
    return assembleLoadVector(mesh, f, userData)

def createStiffnessMatrix(mesh, a=None):
    r"""Create the Stiffness matrix.

    Calculates the Stiffness matrix :math:`{\bf S}` for the given mesh scaled
    with the per cell values a.

    ..math::
            ...

    Parameters
    ----------
    mesh : :gimliapi:`GIMLI::Mesh`
        Arbitrary mesh to calculate the stiffness for.
        Type of base and shape functions depends on the cell types.

    a : array, either complex or real, callable
        Per cell values., e.g., physical parameter. If None given default is 1.

    Returns
    -------
    A : :gimliapi:`GIMLI::RSparseMatrix`
        Stiffness matrix
    """

    if a is None:
        a = pg.RVector(mesh.cellCount(), 1.0)

    A = None

    if isinstance(a[0], float) or isinstance(a[0], np.float64):
        A = pg.RSparseMatrix()
        A.fillStiffnessMatrix(mesh, a)
        return A
    else:
        A = pg.CSparseMatrix()

    # create matrix structure regarding the mesh
    A.buildSparsityPattern(mesh)

    # define a local element matrix
    A_l = pg.ElementMatrix()
    for c in mesh.cells():
        A_l.ux2uy2uz2(c)
        A.add(A_l, scale=a[c.id()])

#        if c.id() == 0:
#            print(c.id(), A_l)
    return A


def createMassMatrix(mesh, b=None):
    r"""Create the mass matrix.

    Calculates the Mass matrix (Finite element identity matrix)
    the given mesh.

    ..math::
            ...

    Parameters
    ----------
    mesh : :gimliapi:`GIMLI::Mesh`

        Arbitrary mesh to calculate the mass element matrix.
        Type of base and shape functions depends on the cell types.

    b : array
        Per cell values. If None given default is 1.

    Returns
    -------
    A : :gimliapi:`GIMLI::RSparseMatrix`
        Mass element matrix
    """

    # need callable here
    if b is None:
        b = pg.RVector(mesh.cellCount(), 1.0)

    B = pg.RSparseMatrix()
    B.fillMassMatrix(mesh, b)
    return B

    # create matrix structure regarding the mesh
    # B.buildSparsityPattern(mesh)
    # define a local element matrix
    # B_l = pg.ElementMatrix()
    # for c in mesh.cells():
    #    B_l.u2(c)
    # # check if b[i] == B*b
    #    B_l *= b[c.id()]
    #    B += B_l
    # return B


def _feNorm(u, A):
    """Create a norm within a Finite Element space.

    Create the Finite Element Norm with a preassembled system matrix.
    """
    return np.sqrt(pg.dot(u, A.mult(u)))
    

def L2Norm(u, M=None, mesh=None):
    r"""Create Lebesgue (L2) norm for the finite element space.

    Finde the L2 Norm for a solution in the finite element space.
    :math:`u` exact solution
    :math:`{\bf M}` Mass matrix, i.e., Finite element identity matrix.

    .. math::

        L2(f(x)) = || f(x) ||_{L^2} & = (\int |f(x)|^2 \d x)^{1/2} \\
                 & \approx h (\sum |f(x)|^2 )^{1/2}
        L2(u) = || u ||_{L^2} & = (\int |u|^2 \d x)^{1/2} \\
                 & \approx (\sum M (u)) ^{1/2}
        e_{L2_rel} = \frac{L2(u)}{L2(u)} & = 
                \frac{(\sum M(u))^{1/2}}{(\sum M u)^{1/2}}

    The error for any approximated solution :math:`u_h` correlates to the L2 
    norm of 'L2Norm(u - u_h, M)'. If you like relative values, you can also 
    normalize this error with 'L2Norm(u - u_h, M) / L2Norm(u, M)*100'.

    Parameters
    ----------
    u : iterable
        Node based value to compute the L2 norm for.

    M : Matrix
        Mass element matrix.
        
    mesh : :gimliapi:`GIMLI::Mesh`
        Mesh with the FE space to generate M if necessary.

    Returns
    -------
    ret : float
        :math:`L2(u)` norm.

    """
    if M is None and mesh is not None:
        M = solver.createMassMatrix(mesh)

    if M is None:
        # M is Identity matrix
        return np.sqrt(pg.dot(u, u))

    return _feNorm(u, M)


def H1Norm(u, S=None, mesh=None):
    r"""Create (H1) norm for the finite element space.

    Parameters
    ----------
    u : iterable
        Node based value to compute the H1 norm for.

    S : Matrix
        Stiffness matrix.
        
    mesh : :gimliapi:`GIMLI::Mesh`
        Mesh with the FE space to generate S if necessary.

    Returns
    -------
    ret : float
        :math:`H1(u)` norm.

    """
    if S is None and mesh is not None:
        S = solver.createStiffnessMatrix(mesh)

    if S is None:
        raise StandardException("Need S or mesh here to calculate H1Norm")
        
    return _feNorm(u, S)


def solve(mesh, **kwargs):
    r"""Solve partial differential equation.

    This is a syntactic sugar convenience function for solving partial
    differential equation on a given mesh.
    Using the best guess method for the given parameter.
    Currently only Finite Element calculation using
    :py:mod:`pygimli.solver.solveFiniteElements`

    TODO
    :py:mod:`pygimli.solver.solveFiniteVolume`

    """
    return solveFiniteElements(mesh, **kwargs)


def solveFiniteElements(mesh, a=1.0, b=0.0, f=0.0, bc=None,
                        times=None, userData=None,
                        verbose=False, **kwargs):
    r"""Solve partial differential equation with Finite Elements.

    This is a syntactic sugar convenience function for using the Finite Element
    functionality of the library core to solve partial differential equation (PDE)
    that match the following form:

    .. math::

        \frac{\partial u}{\partial t} & = \nabla\cdot(a \nabla u) + b u + f(\mathbf{r},t)~~|~~\Omega_{\text{Mesh}}\\
        u & = h~~|~~\Gamma_{\text{Dirichlet}}\\
        \frac{\partial u}{\partial \mathbf{n}} & = g~~|~~\Gamma_{\text{Neumann}}\\
        \alpha u + \beta\frac{\partial u}{\partial \mathbf{n}} & = \gamma~~|~~\Gamma_{\text{Robin}}


    for the scalar solution :math:`u(\mathbf{r}, t)` at each node of a given mesh.
    The Domain :math:`\Omega` and the Boundary :math:`\Gamma` are defined
    through the mesh with appropriate boundary marker.

    TODO:

        * unsteady ub and dub
        * 'Infinity' Boundary condition (u vanishes at infinity)
        * 'Cauchy' Boundary condition
        (guaranties u and du on same boundary, will never work here because the problem
        becomes ill posed and would need some inverse strategy to solve.)

    Parameters
    ----------
    mesh : :gimliapi:`GIMLI::Mesh`
        Mesh represents spatial discretization of the calculation domain

    a   : value | array | callable(cell, userData)
        Cell values

    b   : value | array | callable(cell, userData)
        Cell values

    f : value | array(cells) | array(nodes) | callable(args, kwargs)
        force values

    bc : dict()
        Dictionary of boundary conditions.
        Current supported boundary conditions are by dictionary keys:
        'Dirichlet', 'Neumann', 'Robin'.

        The dictionary can contain multiple "key: Arg"
        Arg will be parsed by :py:mod:`pygimli.solver.solver.parseArgPairToBoundaryArray`

        If the dictionary key is 'Node' then fixed values for single node
        indices can by be given. e.g., bc={'Node': [nodeID, value]}.
        Note this is only a shortcut for
        bc={'Dirichlet': [mesh.node(nodeID), value]}.

    times : array [None]
        Solve as time dependent problem for the given times.

    **kwargs:

        u0 : value | array | callable(pos, userData)
            Node values

        uB : value | array | callable(pos, userData)
            DEPRECATED use bc={'Dirichlet' | uB}

        uN : list([node, value])
            DEPRECATED use bc={'Node' | uN}

        duB : value | array | callable(pos, userData)
            DEPRECATED use bc={'Neumann' | duB}

        theta : float [1]
            - :math:`theta = 0` means explicit Euler, maybe stable for
            :math:`\Delta t \quad\text{near}\quad h`
            - :math:`theta = 0.5`, Crank-Nicolsen, maybe instable
            - :math:`theta = 1`, implicit Euler

            If unsure choose :math:`\theta = 0.5 + \epsilon`, which is probably stable.

        stats : bool
            Give some statistics.

        progress : bool
            Give some calculation progress.

        assembleOnly : bool
            if set stops after matrix asssemblation

        ws : dict
            The WorkSpace is a dictionary that will get
            some temporary data during the calculation.
            Any keyvalue 'u' in the dictionary will be used for the resulting array.

    Returns
    -------
    u : array
        Returns the solution u either 1,n array for stationary problems or
        for m,n array for m time steps

    See also
    --------
        :ref:`tut:modelling` and :py:mod:`pygimli.solver.solve`

    Examples
    --------
    >>> import pygimli as pg
    >>> from pygimli.meshtools import polytools as plc
    >>> from pygimli.mplviewer import drawField, drawMesh
    >>> import matplotlib.pyplot as plt
    >>> world = plc.createWorld(start=[-10, 0], end=[10, -10],
    ...                         marker=1, worldMarker=False)
    >>> c1 = plc.createCircle(pos=[0.0, -5.0], radius=3.0, area=.1, marker=2)
    >>> mesh = pg.meshtools.createMesh([world, c1], quality=34.3)
    >>> u = pg.solver.solveFiniteElements(mesh, a=[[1, 100], [2, 1]],
    ...                                   bc={'Dirichlet':[[4, 1.0], [2, 0.0]]})
    >>> fig, ax = plt.subplots()
    >>> pc = drawField(ax, mesh, u)
    >>> drawMesh(ax, mesh)
    >>> plt.show()
    """
    if bc is None:
        bc={}
    if 'uB' in kwargs:
        pg.deprecated('bc arg uB', "bc={'Dirichlet': uB}")
        bc['Dirichlet'] = kwargs.pop('uB')
    if 'duB' in kwargs:
        pg.deprecated('bc arg duB', "bc={'Neumann': duB}")
        bc['Neumann'] = kwargs.pop('duB')
    if 'uN' in kwargs:
        pg.deprecated('bc arg uN', "bc={'Node': duB}")
        bc['Node'] = kwargs.pop('uN')

    workSpace = kwargs.pop('ws', dict())
    debug = kwargs.pop('debug', False)
    stats = kwargs.pop('stats', False)

    mesh.createNeighbourInfos()
    if verbose:
        print("Mesh: ", str(mesh))

    dof = mesh.nodeCount()

    swatch = pg.Stopwatch(True)
    swatch2 = pg.Stopwatch(True)

    # check for material parameter
    a = parseArgToArray(a, nDof=mesh.cellCount(), mesh=mesh, userData=userData)

    S = createStiffnessMatrix(mesh, a)
    M = createMassMatrix(mesh)

    A = None

    if b != 0:
        b = parseArgToArray(b, nDof=mesh.cellCount(), mesh=mesh, userData=userData)
        A = S + M * b
    else:
        A = S
    
    if debug:
        print("5: ", swatch2.duration(True))

    if times is None:
        rhs = assembleForceVector(mesh, f, userData=userData)

        assembleBC_(bc, mesh, S, rhs, time=None, userData=userData)
        
        if debug:
            print("6c: ", swatch2.duration(True))

        # create result array
        u = None
        if 'u' in workSpace:
            u = workSpace['u']

        singleForce = True
        if hasattr(rhs, 'ndim'):
            if rhs.ndim == 2:
                singleForce = False
                if u is None:
                    u = np.zeros(rhs.shape)
        else:
            if isinstance(a[0], complex):
                if u is None:
                    u = pg.CVector(rhs.size(), 0.0)
                rhs = pg.toComplex(rhs)
            else:
                if u is None:
                    u = pg.RVector(rhs.size(), 0.0)

        assembleTime = swatch.duration(True)
        if stats:
            stats.assembleTime = assembleTime

        if verbose:
            print(("Asssemblation time: ", assembleTime))

        # showSparseMatrix(S)

        workSpace['S'] = S
        workSpace['M'] = M
        workSpace['A'] = A
        workSpace['rhs'] = rhs

        if 'assembleOnly' in kwargs:
            return

        solver = pg.LinSolver(False)
        solver.setMatrix(S, 0)

        if singleForce:
            u = solver.solve(rhs)
        else:
            for i, r in enumerate(rhs):
                if isinstance(u[i], pg.RVector):
                    solver.solve(r, u[i])
                else:
                    u[i] = solver.solve(r)

        solverTime = swatch.duration(True)
        if verbose:
            if stats:
                stats.solverTime = solverTime
            print(("Solving time: ", solverTime))


        if len(kwargs.keys()) > 0:
            print("Warning! Unused arguments", **kwargs)
        return u

    else: # times given

        pg.solver.checkCFL(times, mesh, max(a))

        if debug:
            print("start TL", swatch.duration())

        M = createMassMatrix(mesh)
        F = assembleForceVector(mesh, f)

        u0 = np.zeros(dof)
        if 'u0' in kwargs:
            u0 = parseArgToArray(kwargs['u0'], dof, mesh, userData)
        
        progress = None
        if 'progress' in kwargs:
            from pygimli.utils import ProgressBar
            progress = ProgressBar(its=len(times), width=40, sign='+')

        theta = kwargs.pop('theta', 1.0)
        dynamic = kwargs.pop('dynamic', False)

        if not dynamic:
            S = createStiffnessMatrix(mesh, a)
            assembleBC_(bc, mesh, S, F, time=0.0, userData=userData)
            return crankNicolson(times, theta, S, M, F, u0=u0, progress=progress)

        rhs = np.zeros((len(times), dof))
        # rhs kann zeitabhängig sein ..wird hier nicht berücksichtigt
        rhs[:] = F  # this is slow: optimize

        if debug:
            print("rhs", swatch.duration())
        U = np.zeros((len(times), dof))
        U[0, :] = u0

        # init state
        u = pg.RVector(dof, 0.0)

        if debug:
            print("u0", swatch.duration())

        measure = 0.
        for n in range(1, len(times)):
            swatch.reset()

            dt = times[n] - times[n - 1]

            # previous timestep
            # print "i: ", i, dt, U[i - 1]

            swatch.reset()
            # (A + a*B)u is fastest,
            # followed by A*u + (B*u)*a and finally A*u + a*B*u and
            b = (M + (dt * (theta - 1.)) * S) * U[n - 1] + \
                dt * ((1.0 - theta) * rhs[n - 1] + theta * rhs[n])

            # print ('a',swatch.duration(True))
            # b = M * U[n - 1] - (A * U[n - 1]) * (dt*(1.0 - theta)) + \
            # dt * ((1.0 - theta) * rhs[n - 1] + theta * rhs[n])

            # print ('b',swatch.duration(True))

            # b = M * U[n - 1] - (dt*(1.0 - theta)) * A * U[n - 1] + \
            # dt * ((1.0 - theta) * rhs[n - 1] + theta * rhs[n])
            # print ('c',swatch.duration(True))

            measure += swatch.duration()

            A = M + S * dt * theta

            assembleBC_(bc, mesh, S, b, time=times[n], userData=userData)
            
            # u = S/b
            t_prep = swatch.duration(True)
            solver = pg.LinSolver(S, verbose)
            solver.solve(b, u)

            if 'plotTimeStep' in kwargs:
                kwargs['plotTimeStep'](u, times[n])

            U[n, :] = np.asarray(u)

            if progress:
                progress.update(n,
                            ' t_prep: ' + str(round(t_prep, 5)) + 's' + \
                            ' t_step: ' + str(round(swatch.duration(),5)) + 's')

        if debug:
            print("Measure(" + str(len(times)) + "): ",
                  measure, measure / len(times))
        return U

def checkCFL(times, mesh, vMax):
    """Check Courant-Friedrichs-Lewy condition.

    For advection and flow problems. CFL Number should be lower then 1 to 
    ensure stability.

    Parameters
    ----------
    """
    if times is not None:
        dt = times[1] - times[0]
        dx = 0.0
        
        if mesh.dimension() == 1:
            dx = min(mesh.cellSizes())
        else:
            dx = min(mesh.boundarySizes())

        c = vMax * dt / dx
        if c > 1:
            print("Courant-Friedrichs-Lewy Number:", c,
                  "but sould be lower 1 to ensure movement inside a cell "
                  "per timestep. ("
                  "vmax =", vMax,
                  "dt =", dt,
                  "dx =", dx,
                  "dt <", dx/vMax,
                  " | N > ", int((times[-1]-times[0])/(dx/vMax))+1, ")")
    return c

def crankNicolson(times, theta, S, I, f, u0=None, progress=None, debug=None):
    """
        S = constant over time
        f = constant over time
    """
    if len(times) < 2:
        raise BaseException("We need at least 2 times for "
                            "Crank-Nicolsen time discretization." + str(len(times)))
    sw = pg.Stopwatch(True)

    if u0 is None:
        u0 = np.zeros(len(f))

    u = np.zeros((len(times), len(f)))
    u[0, :] = u0
    dt = times[1] - times[0]

    rhs = np.zeros((len(times), len(f)))
    rhs[:] = f

    timeAssemble = []
    timeSolve = []

    timeMeasure = False
    if progress:
        timeMeasure = True

    A = I + S * (dt * theta)
    solver = pg.LinSolver(A, verbose=False)
    St = I - S * dt # cache what is possible the theta=0
    for n in range(1, len(times)):

        if timeMeasure:
            pg.tic()

#        pg.tic()
        #bRef = (I + (dt * (theta - 1.)) * S) * u[n - 1] + \
           #dt * ((1.0 - theta) * rhs[n - 1] + theta * rhs[n])
#        pg.toc()
#
#        pg.tic()
        #b = I * u[n - 1] + ((dt * (theta - 1.)) * S) * u[n - 1] + \
           #dt * ((1.0 - theta) * rhs[n - 1] + theta * rhs[n])
#        pg.toc()
#
#        pg.tic()
    
        if theta == 0:
            b = St * u[n - 1] + dt * rhs[n - 1]
        else:
            b = I * u[n - 1] + S.mult(dt * (theta - 1.) * u[n - 1]) + \
                dt * ((1.0 - theta) * rhs[n - 1] + theta * rhs[n])

#        pg.toc()
#        print(np.linalg.norm(b-b1))
        #np.testing.assert_allclose(bRef, b)

        if timeMeasure:
            timeAssemble.append(pg.dur())

        if timeMeasure:
            pg.tic()

        u[n, :] = solver.solve(b)

        if timeMeasure:
            timeSolve.append(pg.dur())

        # A = (I + dt * theta * S)
        # u[n, : ] = linsolve(A, b)

        if progress:
            progress.update(n,
                            ' t_prep: ' + str(round(timeAssemble[-1], 5)) + 's' + \
                            ' t_step: ' + str(round(timeSolve[-1], 5)) + 's')

        #if verbose and (n % verbose == 0):
            ## print(min(u[n]), max(u[n]))
            #print("timesteps:", n, "/", len(times),
                  #'runtime:', sw.duration(), "s",
                  #'assemble:', np.mean(timeAssemble),
                  #'solve:', np.mean(timeSolve))
    return u


class RungeKutta(object):
    """TODO DOCUMENT ME"""
    rk4a = [0.0,
            -567301805773.0/1357537059087.0,
            -2404267990393.0/2016746695238.0,
            -3550918686646.0/2091501179385.0,
            -1275806237668.0/842570457699.0]
    rk4b = [1432997174477.0/9575080441755.0,
            5161836677717.0/13612068292357.0,
            1720146321549.0/2090206949498.0,
            3134564353537.0/4481467310338.0,
            2277821191437.0/14882151754819.0]
    rk4c = [0.0,
            1432997174477.0/9575080441755.0,
            2526269341429.0/6820363962896.0,
            2006345519317.0/3224310063776.0,
            2802321613138.0/2924317926251.0]

    def __init__(self, solver, verbose=False):
        """TODO DOCUMENT_ME"""
        self.solver = solver
        self.verbose = verbose
        self.order = 5
        self.dt = None
        self.time = 0
        self.tMax = None
        self.u = None
        self.resu = None
        self.nSteps = 0

    def run(self, u0, dt, tMax=1):
        """TODO DOCUMENT_ME"""
        self.start(u0, dt, tMax)

        for _ in range(self.nSteps):
            self.step()

        return self.u

    def start(self, u0, dt, tMax=1):
        """TODO DOCUMENT_ME"""
        self.nSteps = int(np.ceil(tMax/dt))
        self.dt = dt
        self.time = 0
        self.tMax = tMax
        self.u = deepcopy(u0)
        self.resu = deepcopy(u0)

        if isinstance(self.resu, list):
            for r in self.resu:
                r *= 0.0
        else:
            self.resu *= 0.0

    def step(self):
        """TODO DOCUMENT ME"""
        if self.time + self.dt > self.tMax:
            self.dt = self.tMax - self.time

        if self.order == 1:
            # explicit Euler
            k1 = self.solver.explicitRHS(self.u,
                                         self.time)
            self.u += self.dt * k1

        elif self.order == 3:
            k1 = self.solver.explicitRHS(self.u, self.time)
            k1 = self.u + self.dt * k1

            k2 = self.solver.explicitRHS(k1, self.time)
            k2 = (3*self.u + k1 + self.dt*k2)/4

            k3 = self.solver.explicitRHS(k2, self.time)

            self.u = (self.u + 2*k2 + 2*self.dt*k3)/3

        elif self.order == 4:
            # classical 4 step Runga-Kutta rk4
            k1 = self.solver.explicitRHS(self.u,
                                         self.time)
            k2 = self.solver.explicitRHS(self.u + self.dt/2 * k1,
                                         self.time + self.dt/2)
            k3 = self.solver.explicitRHS(self.u + self.dt/2 * k2,
                                         self.time + self.dt/2)
            k4 = self.solver.explicitRHS(self.u + self.dt * k3,
                                         self.time + self.dt)
            self.u += 1./6. * self.dt * (k1 + 2.0 * k2 + 2.0 * k3 + k4)

        elif self.order == 5:
            # low storage Version of rk4
            for jRK in range(5):
                tLocal = self.time + self.rk4c[jRK] * self.dt

                rhs = self.solver.explicitRHS(self.u, tLocal)

                if isinstance(self.resu, list):
                    for i in range(len(self.resu)):
                        self.resu[i] = self.rk4a[jRK] * self.resu[i] + \
                                       self.dt * rhs[i]

                        self.u[i] += self.rk4b[jRK] * self.resu[i]
                else:
                    self.resu = self.rk4a[jRK] * self.resu + self.dt * rhs
                    self.u += self.rk4b[jRK] * self.resu

        self.time += self.dt
        return self.u


if __name__ == "__main__":
    pass
