#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""TODO DOCUMENT ME"""
from copy import deepcopy

import numpy as np
import pygimli as pg
from pygimli.utils import unique


def parseDictKey_(key, markers):
    """ Parse dictionary key of type str to marker list.

    Utility function to parse a dictionary key string into a valid list of markers containing in a given markers list.

    Parameters
    ----------
    key: str
        Supported are
        - 'm1': Single marker
        - 'm1,m2': Comma separated list
        - ':': Slice wildcard
        - 'start:stop:step': Slice like syntax

    markers: [int]
        List of integers, e.g., cell or boundary markers

    Returns
    -------
    mas: [int]
        List of integers described by key
    """
    markers = pg.unique(markers)
    if ',' in key:
        mas = [int(k) for k in key.split(',')]
    elif ':' in key:
        sse = key.split(':')

        start = markers[0]
        stop = markers[-1] + 1
        step = 1

        if len(sse) > 0:
            try:
                start = int(sse[0])
            except:
                pass
        if len(sse) > 1:
            try:
                stop = int(sse[1])
            except:
                pass
        if len(sse) > 2:
            try:
                step = int(sse[2])
            except:
                pass

        mas = list(range(start, stop, step))
        return mas
    else:
        mas = [int(key)]

    return [m for m in mas if m in markers]


def cellValues(mesh, arg, **kwargs):
    """Get a value for each cell.

    Returns a array or vector of length mesh.cellCount() based on arg.
    The preferable arg is a dictionary for the cell marker and the appropriate cell value. The designated value can be calculated using a
    callable(cell, **kwargs), which is called on demand

    Parameters
    ----------
    arg : float | int | complex | ndarray | iterable | callable | dict
        Argument to be parsed as cell data.
        If arg is a dictionary, the dict key will be interpreted as cell marker:

        Dictionary is key: value. Value can be float, int, complex or ndarray. The last for anistropic or elastic tensors.

        Key can be integer for cell marker or str, which will be interpreted as splice or list. See examples or parseDictKey_.

    mesh : :gimliapi:`GIMLI::Mesh`
        Used if arg is callable

    userData : class
        Used if arg is callable

    Returns
    -------
    ret : :gimliapi:`GIMLI::RVector` | ndarray(mesh.cellCount(), xx )
        Array of desired length filled with the appropriate values.

    Examples
    --------
    >>> import pygimli as pg
    >>>
    >>> mesh = pg.createGrid(x=range(5))
    >>> mesh.setCellMarkers([1, 1, 2, 2])
    >>> print(mesh.cellCount())
    4
    >>> print(pg.solver.cellValues(mesh, [1, 2, 3, 4]))
    [1, 2, 3, 4]
    >>> print(pg.solver.cellValues(mesh, {1:1.0, 2:10}))
    [1.0, 1.0, 10, 10]
    >>> print(pg.solver.cellValues(mesh, {':':2.0}))
    [2.0, 2.0, 2.0, 2.0]
    >>> print(pg.solver.cellValues(mesh, {'0:2':3.0}))
    [3.0, 3.0, None, None]
    >>> print(np.array(pg.solver.cellValues(mesh, {'1:3' : np.diag([1.0, 2.0])})))
    [[[1. 0.]
      [0. 2.]]
    <BLANKLINE>
     [[1. 0.]
      [0. 2.]]
    <BLANKLINE>
     [[1. 0.]
      [0. 2.]]
    <BLANKLINE>
     [[1. 0.]
      [0. 2.]]]
    >>> print(np.array(pg.solver.cellValues(mesh, {':' : pg.core.CMatrix(2, 2)})))
    [[[0.+0.j 0.+0.j]
      [0.+0.j 0.+0.j]]
    <BLANKLINE>
     [[0.+0.j 0.+0.j]
      [0.+0.j 0.+0.j]]
    <BLANKLINE>
     [[0.+0.j 0.+0.j]
      [0.+0.j 0.+0.j]]
    <BLANKLINE>
     [[0.+0.j 0.+0.j]
      [0.+0.j 0.+0.j]]]
    >>> print(pg.solver.cellValues(mesh, {'1,2':1 + 1j*2.0}))
    [(1+2j), (1+2j), (1+2j), (1+2j)]
    >>> def cellVal(c, b=1):
    ...     return c.center()[0]*b
    >>> t = pg.solver.cellValues(mesh, {':' : cellVal})
    >>> print([t[c.id()](c) for c in mesh.cells()])
    [0.5, 1.5, 2.5, 3.5]
    """
    if isinstance(arg, dict):

        try:
            val = list(arg.values())[0]
        except:
            pg.error("Can't interpret empty dictionary:", arg)
            val = 1.0

        ret = [None] * mesh.cellCount()

        for key, val in arg.items():
            if isinstance(key, str):
                mas = parseDictKey_(key, mesh.cellMarkers())

                for m in mas:
                    for i in pg.find(mesh.cellMarkers() == m):
                        ret[i] = val
            else:
                for i in pg.find(mesh.cellMarkers() == key):
                    ret[i] = val

        return ret

    # if arg have allready correct size
    if hasattr(arg, '__len__'):
        if len(arg) == mesh.cellCount():
            return arg

    # if arg if scalar or globel data type, ndarray or Matrix but not the right size assume global tensor
    if isinstance(arg, np.ndarray) or \
        isinstance(arg, pg.core.RMatrix) or \
        isinstance(arg, pg.core.CMatrix) or \
        isinstance(arg, float) or \
        isinstance(arg, int) or \
        isinstance(arg, complex):
        return [arg]*mesh.cellCount()

    return parseArgToArray(arg,
                           nDof=mesh.cellCount(),
                           mesh=mesh, **kwargs)


def parseArgToArray(arg, nDof, mesh=None, userData={}):
    """
    Parse array related arguments to create a valid value array.

    Parameters
    ----------
    arg : float | int | iterable | callable
        The target array value that will be converted to an array.

        If arg is a callable with it must fulfill:

        :: arg(cell|node|boundary, userData={})

        Where MeshEntity is one of
        :gimliapi:`GIMLI::Cell` ,
        :gimliapi:`GIMLI::Node` or
        :gimliapi:`GIMLI::Boundary`
        depending on nDof, where nDof is mesh.cellCount(),
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
    pg.warn('check if obsolete: parseArgToArray')
    if not hasattr(nDof, '__len__'):
        nDof = [nDof]

    try:
        return pg.Vector(nDof[0], float(arg))
    except BaseException as _:
        pass

    if hasattr(arg, '__len__'):
        if isinstance(arg, np.ndarray):
            if len(arg) == nDof[0]:
                return arg
            else:
                raise Exception('Given array does not have requested (' +
                                str(nDof) + ') size (' +
                                str(len(arg)) + ')')

        for n in nDof:
            if len(arg) == n:
                return arg

        try:
            # [marker, val] || [[marker, val]]
            return parseMapToCellArray(arg, mesh)
        except:
            raise Exception("Array 'arg' has the wrong size: " +
                            str(len(arg)) + " != " + str(nDof))
    elif hasattr(arg, '__call__'):
        ret = pg.Vector(nDof[0], 0.0)

        if not mesh:
            raise Exception("Please provide a mesh for the callable"
                            "argument to parse ")

        if nDof[0] == mesh.nodeCount():
            for n in mesh.nodes():
                if userData:
                    ret[n.id()] = arg(node=n, userData=userData)
                else:
                    ret[n.id()] = arg(node=n)
        elif nDof[0] == mesh.cellCount():
            for c in mesh.cells():
                if userData:
                    ret[c.id()] = arg(cell=c, userData=userData)
                else:
                    ret[c.id()] = arg(cell=c)
        elif nDof[0] == mesh.boundaryCount():
            for b in mesh.boundaries():
                if userData:
                    ret[b.id()] = arg(boundary=b, userData=userData)
                else:
                    ret[b.id()] = arg(boundary=b)
        else:
            raise Exception("Cannot parse callable argument " + str(nDof) +
                            " nodes: " + str(mesh.nodeCount()) +
                            " cells: " + str(mesh.cellCount()))

        return ret
    raise Exception("Cannot parse argument type " + str(type(arg)))


def generateBoundaryValue(boundary, arg, time=0.0, userData={}):
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

        :: arg(boundary=:gimliapi:`GIMLI::Boundary`, time=0.0, userData={})

        The callable function arg have to return appropriate values
        for all nodes of the boundary or one scalar value for all nodes.

    Returns
    -------
    val : float or list of ..
    """
    if hasattr(boundary, '__len__'):
        pg.deprecated('bad design')
        values = np.zeros(len(boundary))
        for i, b in enumerate(boundary):
            values[i] = generateBoundaryValue(b, arg[i], time, userData)
        return values
    val = 0.

    if callable(arg):
        kwargs = dict()
        if time != 0.0 and time is not None:
            kwargs['time'] = time
        if userData is not None and userData.keys():
            kwargs['userData'] = userData
        try:
            # val(boundary, time=0, userData={})
            val = arg(boundary, **kwargs)
        except BaseException as e:
            print(arg, "(", kwargs, ")")
            pg.critical("Wrong arguments for callback function.", e)

    elif hasattr(arg, '__len__'):
        if callable(arg[0]):
            kwargs = arg[1]
            if time != 0.0 and time is not None:
                kwargs['time'] = time
            val = arg[0](boundary=boundary, **kwargs)
        else:
            if len(arg) < 4:
                # we assume vector valued arguments
                return arg
            else:
                pg.warn('necessary?? generateBoundaryValue')
                val = generateBoundaryValue(boundary, arg[boundary.id()],
                userData)
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
        - [marker, [callable, *kwargs]]
        - [marker, [arg_x, arg_y, arg_z]]
        - [boundary, arg]
        - ['*', arg]
        - [node, arg]
        - [[marker, ...], arg] (REMOVE ME because of bad design)
        - [[boundary,...], arg]  (REMOVE ME because of bad design)
        - [marker, callable, *kwargs] (REMOVE ME because of bad design)
        - [[marker, ...], callable, *kwargs]  (REMOVE ME because of bad design)

        arg will be parsed by
        :py:mod:`pygimli.solver.solver.generateBoundaryValue`
        and distributed to each boundary.
        Callable functions will be executed at run time.
        '*' will be interpreted as
        'all boundary elements with one neighboring cell'

    mesh : :gimliapi:`GIMLI::Mesh`
        Used to find boundaries by marker.

    Returns
    -------

    bc : list()
        [:gimliapi:`GIMLI::Boundary`, value|callable]
    """
    bc = []
    bounds = []
    if isinstance(pair[1], list):
        #  [marker, [callable, *kwargs]]
        if callable(pair[1][0]):
            pair = [pair[0]] + pair[1]

    if pair[0] == '*':
        mesh.createNeighbourInfos()
        for b in mesh.boundaries():
            if b.leftCell() is not None and b.rightCell() is None:
                bounds.append(b)
    elif isinstance(pair[0], int):
        bounds = mesh.findBoundaryByMarker(pair[0])
    elif isinstance(pair[0], pg.core.Node):
        bc.append(pair)
        return bc

    # bad Design .. need to remove
    elif isinstance(pair[0], list):
        pg.deprecated('bad design')
        # [[,,..], ]
        for b in pair[0]:
            for bi in mesh.boundaries(pg.find(mesh.boundaryMarkers() == b)):
                bounds.append(bi)

    elif isinstance(pair[0], pg.core.stdVectorBounds):
        pg.deprecated('bad design')
        pg.warn('inuse? pair[0], pg.core.stdVectorBounds)')#20200115
        bounds = pair[0]
    elif isinstance(pair[0], pg.core.Boundary):
        pg.warn('inuse? isinstance(pair[0], pg.core.Boundary)')#20200115
        bc.append(pair)
        return bc

    for b in bounds:
        val = None
        # print('-'*50)
        # print(b, pair[1], callable(pair[1]))
        # print('+'*50)
        if callable(pair[1]):
            # don't execute the callable here
            # we want to call them at runtime
            if len(pair) > 2:
                val = pair[1:]
            else:
                val = pair[1]
        else:
            # this will be executed
            val = generateBoundaryValue(b, pair[1])
        bc.append([b, val])

    #print('#'*30)
    return bc


def parseArgToBoundaries(args, mesh):
    """
    Parse boundary related arguments to create a valid boundary value list:
    [ :gimliapi:`GIMLI::Boundary`, value|callable ]

    TODO
    ----
    - callable dynamic at runtime

    Parameters
    ----------
    args : dict, float, callable
        Dictionary is prefered (key=value|callable).
        If args is just a callable or float every outer boundary will be processed with args.

        List pairs will be removed or not correct parsed for vector valued problems.
        Callable will be evaluated at runtime. See examples.
        Else see :py:mod:`pygimli.solver.solver.parseArgPairToBoundaryArray`

    mesh : :gimliapi:`GIMLI::Mesh`
        Used to find boundaries by marker

    Returns
    -------
    boundaries : list()
        [ :gimliapi:`GIMLI::Boundary`, value|callable ]

    Examples
    --------
    >>> # no need to import matplotlib. pygimli show does
    >>> import pygimli as pg
    >>> import pygimli.meshtools as mt
    >>> plc = mt.createWorld([0, 0], [1, -1], worldMarker=0)
    >>> ax, _ = pg.show(plc, boundaryMarker=True)
    >>> mesh = mt.createMesh(plc)
    >>> # all four outer boundaries get value = 1.0
    >>> b = pg.solver.parseArgToBoundaries(1.0, mesh)
    >>> print(len(b))
    4
    >>> # all edges with marker 1 get value = 1.0
    >>> b = pg.solver.parseArgToBoundaries({1: 1.0}, mesh)
    >>> print(len(b))
    1
    >>> # same as above with marker 2 get value 2
    >>> b = pg.solver.parseArgToBoundaries({'1': 1.0, 2 : 2.0}, mesh)
    >>> print(len(b))
    2
    >>> # same as above with marker 3 get value 3
    >>> b = pg.solver.parseArgToBoundaries({1:1., 2:2., 3:3.}, mesh)
    >>> print(len(b))
    3
    >>> # Boundary values for vector valued problem
    >>> b = pg.solver.parseArgToBoundaries({1:[1.0, 1.0]}, mesh)
    >>> print(len(b), b[0][1])
    1 [1.0, 1.0]
    >>> # edges with marker 1 and 2 get value 1
    >>> b = pg.solver.parseArgToBoundaries({'1,2':1.0}, mesh)
    >>> print(len(b))
    2
    >>> b = pg.solver.parseArgToBoundaries({'1, 2, 3': 1.0}, mesh)
    >>> print(len(b))
    3
    >>> b = pg.solver.parseArgToBoundaries({'1:4':1.0, 4:4.0}, mesh)
    >>> print(len(b))
    4
    >>> b = pg.solver.parseArgToBoundaries({mesh.node(0):0.0}, mesh)
    >>> print(len(b))
    1
    >>> def bCall(boundary):
    ...     u = []
    ...     for i, n in enumerate(boundary.nodes()):
    ...         u.append(i)
    ...     return u
    >>> b = pg.solver.parseArgToBoundaries({1:bCall}, mesh)
    >>> print(len(b),b[0][1](b[0][0]))
    1 [0, 1]
    >>> def bCall(boundary, a1, a2):
    ...     return a1 + a2
    >>> b = pg.solver.parseArgToBoundaries({1: [bCall, {'a1':2, 'a2':3}]}, mesh)
    >>> print(len(b), b[0][1][0](b[0][0], **b[0][1][1]))
    1 5
    >>> b = pg.solver.parseArgToBoundaries({1: [bCall, {'a1':1, 'a2':2}],
    ...                                     2: [bCall, {'a1':3, 'a2':4}]}, mesh)
    >>> print(len(b), b[0][1][0](b[0][0], **b[0][1][1]))
    2 3
    >>> b = pg.solver.parseArgToBoundaries({'1,2': [bCall, {'a1':4, 'a2':5}]}, mesh)
    >>> print(len(b), b[1][1][0](b[1][0], **b[1][1][1]))
    2 9
    >>> pg.wait()
    """
    boundaries = list()

    if isinstance(args, dict):

        try:
            val = list(args.values())[0]
        except:
            pg.error("Can't interprete empty dictionary:", args)

        for key, val in args.items():
            if isinstance(key, str) and key != '*':
                markers = parseDictKey_(key, mesh.boundaryMarkers())

                for m in markers:
                    boundaries += parseArgPairToBoundaryArray([m, val], mesh)
            else:
                boundaries += parseArgPairToBoundaryArray([key, val], mesh)

        return boundaries

    if hasattr(args, '__call__') or \
        isinstance(args, float) or isinstance(args, int):
        return parseArgToBoundaries({'*': args}, mesh)


    elif isinstance(args, list):
        pg.warn('DEPRECATED by bad design [parseArgToBoundaries(lists)] check if obsolete') # 20200115
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

    else:
        raise Exception('cannot interpret boundary token', args)

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
        Fill all unmapped attributes to this default.

    Returns
    -------
    att : array
        Array of length mesh.cellCount()
    """
    pg.warn('check if obsolete: parseMapToCellArray')
    att = pg.Vector(mesh.cellCount(), default)

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
                    pg.warn("parseMapToCellArray: cannot find marker " +
                          str(pair[0]) + " within mesh.")
                else:
                    #print('---------------------')
                    #print(att, idx, pair[1], type(pair[1]), float(pair[1]))
                    if isinstance(pair[1], np.complex):
                        if not isinstance(att, pg.CVector):
                            att = pg.math.toComplex(att)
                        att.setVal(val=pair[1], ids=idx)
                    else:
                        att.setVal(val=float(pair[1]), ids=idx)
            else:
                raise Exception("Please provide a list of [int, value] pairs" +
                                str(pair))
    else:
        print("attributeMap: ", attributeMap)
        raise Exception("Cannot interpret attributeMap!")

    return att


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
    >>> _ = pg.show(mesh, u(mesh.positions()), ax=ax)
    >>> _ = pg.show(mesh, [2.*pg.y(mesh.cellCenters())*pg.x(mesh.cellCenters()),
    ...             pg.x(mesh.cellCenters())**2], ax=ax)
    >>> _ = pg.show(mesh, pg.solver.grad(mesh, u), ax=ax, color='w',
    ...             linewidth=0.4)
    >>> plt.show()
    """
    if r is None:
        r = mesh.cellCenters()

    uv = u
    if callable(u) and not isinstance(u, pg.Vector):
        uv = u(mesh.positions())

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
    >>> print(pg.math.round(pg.solver.div(mesh, v(mesh.boundaryCenters())), 1e-5))
    3 [1.0, 1.0, 1.0]
    >>> print(pg.math.round(pg.solver.div(mesh, v(mesh.cellCenters())), 1e-5))
    3 [0.5, 1.0, 0.5]
    >>> mesh = pg.createGrid(x=np.linspace(0, 1, 4),
    ...                      y=np.linspace(0, 1, 4))
    >>> print(pg.math.round(pg.solver.div(mesh, v(mesh.boundaryCenters())), 1e-5))
    9 [2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0]
    >>> divCells = pg.solver.div(mesh, v(mesh.cellCenters()))
    >>> # divergence from boundary values are exact where the divergence from
    >>> # interpolated cell center values are wrong due to interpolation to the boundary
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
    mesh.createNeighborInfos()
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


def divergence(mesh, func=None, normMap=None, order=1):
    """Divergence for callable function func((x,y,z)).

    MOVE THIS to a better place

    Divergence for callable function func((x,y,z)).
    Return sum div over boundary.

    Parameters
    ----------

    Returns
    -------
    """
    if func is None:
        func = lambda r: r

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
                sw = pg.core.Stopwatch(True)
                mesh.createNeighborInfos()
                print("NeighborInfos()", sw.duration(True))
                # return gauss(grid, F)

            # don't calc for inner boundaries
            if b.leftCell() is not None and b.rightCell() is not None:
                continue

        divS = 0
        shape = b.shape()

        if order == 1:
            if bNorms is not None:
                divS = shape.norm().dot(bNorms[b.id()]) * shape.domainSize()
            else:
                divS = shape.norm().dot(
                    func(shape.center())) * shape.domainSize()
        else:
            weights = pg.core.IntegrationRules.instance().weights(shape, order)
            abscissa = pg.core.IntegrationRules.instance().abscissa(shape, order)

            for i, p in enumerate(abscissa):
                rPos = shape.xyz(p)
                divS += shape.norm().dot(func(rPos)) * \
                    weights[i] * shape.domainSize()

        if directionCheck and b.leftCell() is None:
            divS *= -1
            # raise Exception("invalid mesh: left is None .. every
            # boundary need leftCell")

        diverg += divS

    return diverg


def triDiagToeplitz(dom, a, l, r, start=0, end=-1):
    """Create tri-diagonal Toeplitz matrix."""
    A = pg.matrix.SparseMapMatrix(dom, dom)

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
    A = pg.matrix.SparseMapMatrix(dom, dom)

    if end == -1:
        end = dom

    for i in range(start, end):
        if hasattr(scale, '__len__'):
            A.addVal(i, i, scale[i])
        else:
            A.addVal(i, i, scale)
    return A


def showSparseMatrix(mat, full=False):
    """Show the content of a sparse matrix.

    Parameters
    ----------
    mat : :gimliapi:`GIMLI::SparseMatrix` | :gimliapi:`GIMLI::SparseMapMatrix`
        Matrix to be shown.
    full : bool [False]
        Show as dense matrix.
    """

    if isinstance(mat, pg.matrix.SparseMapMatrix):
        m_ = pg.matrix.SparseMatrix(mat)
        return showSparseMatrix(m_, full)
    else:
        rows = mat.vecRowIdx()
        cols = mat.vecColPtr()
        vals = mat.vecVals()

        matD = None
        if full:
            matD = pg.Matrix(mat.rows(), mat.cols())

        for i in range(mat.rows()):
            for j in range(cols[i], cols[i + 1]):
                if full:
                    matD[i, rows[j]] = vals[j]
                else:
                    print(i, rows[j], vals[j])

        if full:
            print(np.array(matD))


def linSolve(mat, b, solver=None, verbose=False):
    r"""Direct solution after :math:`\textbf{x}` using core LinSolver.

    .. math::
        \textbf{A}\textbf{x} = \textbf{b}

    If :math:`\textbf{A}` is symmetric, sparse and positive definite.

    Parameters
    ----------
    mat : :gimliapi:`GIMLI::RSparseMatrix` | :gimliapi:`GIMLI::RSparseMapMatrix` |
        numpy.array
        System matrix. Need to be symmetric, sparse and positive definite.

    b : iterable array
        Right hand side of the equation.

    solver : str [None]
        Try to choose a solver, 'pg' for pygimli core cholmod or umfpack.
        'np' for numpy linalg or scipy.sparse.linalg.
        Automatic choosing if solver is None depending on matrixtype.

    verbose : bool [False]
        Be verbose.

    Returns
    -------

    x : :gimliapi:`GIMLI::RVector`
        Solution vector
    """
    x = pg.Vector(len(b), .0)

    if solver is None:
        if isinstance(mat, pg.matrix.SparseMatrix) or \
           isinstance(mat, pg.matrix.SparseMapMatrix) or \
            isinstance(mat, pg.matrix.BlockMatrix) or \
            isinstance(mat, pg.matrix.CSparseMatrix):
            solver = 'pg'

    if solver == 'pg':
        _m = mat
        if isinstance(mat, pg.matrix.CSparseMatrix):
            x = pg.CVector(len(b), 0)
        elif isinstance(mat, pg.matrix.SparseMatrix):
            pass
        elif isinstance(mat, pg.matrix.SparseMapMatrix):
            _m = pg.matrix.SparseMatrix(mat)
        elif isinstance(mat, pg.matrix.BlockMatrix):
            _m = mat.sparseMapMatrix()
        else:
            pg.critical("Solver '" + solver + "' does not know how to "
                        "solve linear system with matrixtype:" + mat)

        ls = pg.core.LinSolver(_m, verbose=verbose)
        ls.solve(b, x)
    else:

        if isinstance(mat, np.ndarray):
            return np.linalg.solve(mat, b)

        scipy = pg.optImport('scipy')
        #import scipy.sparse
        from scipy.sparse.linalg import spsolve

        if isinstance(mat, scipy.sparse.csr.csr_matrix) or \
            isinstance(mat, scipy.sparse.coo.coo_matrix):
            if verbose:
                pg.info("linSolve use scipy.sparse")
            return spsolve(mat, b)

        return linSolve(pg.utils.sparseMatrix2csr(mat), b,
                        solver='np', verbose=verbose)

    return x


def assembleLoadVector(mesh, f, userData={}):
    r"""Assemble the load vector. See createLoadVector.
    Maybe we will remove this
    """
    pg.deprecate('createLoadVector') # 20200115
    return createLoadVector(mesh, f, userData)


def createForceVector(mesh, f, userData={}):
    """ Create a right hand side vector for vector valued solutions.

    Parameters
    ----------
    f: [ convertable ]
        List of rhs side options. Must be convertable to createLoadVector.
        See :py:mod:`createLoadVector`
    rhs: np.array()
        Squeezed vector of length mesh.nodeCount() * mesh.dimensions()

    """
    if not isinstance(f, list):
        pg.error("Create Force Vector need list of attribute f with an entry for each dimension.")

    rhs = np.zeros(mesh.nodeCount() * mesh.dim())

    for i in range(mesh.dim()):
        rhs[i*mesh.nodeCount():(i+1)*mesh.nodeCount()] = createLoadVector(mesh, f[i], userData)

    # rhs.reshape(mesh.nodeCount() * mesh.dim()) #contigous not guarantied
    return rhs


def createLoadVector(mesh, f, userData={}):
    """Create right hand side vector based on the given mesh and load values (scalar solution) or force vectors (vector value solution).

    Create right hand side based on the given mesh and load or force
    values.

    TODO
    ----
    - Callable for vector valued problems
    - Callable called dynamic on demand

    Parameters
    ----------
    f: float, array, callable(cell, [userData]), [f_x, f_y, f_z]

        - float will be assumend as constant for all cells
        like rhs = rhs(np.ones(mesh.cellCount() * f),
        - array of length mesh.cellCount() will be processed as load value for each cell: rhs = rhs(f),
        - array of length mesh.nodeCount() will be assumed to be allready processed correct: rhs = f
        - callable is evaluated on once for each cell and need to return a load value for each cell and can have  optional a userData dictionary: `f_cell = f(cell, [userData={}])`
        rhs = rhs(f(c, userData) for c in mesh.cells())
        - list with length of mesh.dimension() of float or array entries will create a squeezed rhs for vector valued problems
        rhs = squeeze([rhs(f[0]), rhs(f[1]), rhs(f[2])])

    Returns
    -------
    rhs: pg.Vector(mesh.nodeCount())

        Right hand side load vector for scalar values or sqeeuzed vector values.
    """
    ### fix for the lazy
    if isinstance(f, int):
        f = float(f)

    ### f is list [fx, fy, [fz]] for vector problems
    if isinstance(f, list):
        if len(f) == mesh.dim():
            return createForceVector(mesh, f, userData=userData)

    ### f is list of array [f_0, f_1, ..., f_n] for scalar problems
    if isinstance(f, list) or hasattr(f, 'ndim'):
        if isinstance(f, list):
            rhs = np.zeros((len(f), mesh.nodeCount()))
            for i, fi in enumerate(f):
                userData['i'] = i
                rhs[i] = createLoadVector(mesh, fi, userData)
            return rhs

        elif f.ndim == 2:
            ### assume rhs [n, nNodes] array is already a valid
            if len(f[0]) == mesh.nodeCount():
                return f

    rhs = pg.Vector(mesh.nodeCount(), 0)

    fArray = None

    if hasattr(f, '__len__'):
        if len(f) == mesh.cellCount():
            # scalar values for each cell
            fArray = f
        elif len(f) == mesh.nodeCount():
            # scalar values for each node
            fArray = f
        elif len(f) == mesh.nodeCount() * mesh.dim():
            # vector values for each node
            # maybe just for speccial cases with allready processed rhs
            return f

    elif hasattr(f, '__call__') and not isinstance(f, pg.Vector):
        fArray = pg.Vector(mesh.cellCount())
        for c in mesh.cells():
            if userData is not None and userData.keys():
                fArray[c.id()] = f(c, userData)
            else:
                fArray[c.id()] = f(c)

    if fArray is None:
        fArray = cellValues(mesh, f, userData=userData)

    if len(fArray) == mesh.cellCount():
        b_l = pg.matrix.ElementMatrix()

        for c in mesh.cells():
            if fArray[c.id()] != 0.0:
                b_l.u(c)
                rhs.add(b_l, fArray[c.id()])

#            print("test reference solution:")
#            rhsRef = pg.Vector(mesh.nodeCount(), 0)
#            for c in mesh.cells():
#                b_l.u(c)
#                for i, idx in enumerate(b_l.idx()):
#                    rhsRef[idx] += b_l.row(0)[i] * fArray[c.id()]
#            np.testing.assert_allclose(rhs, rhsRef)
#            print("Remove revtest in assembleLoadVector after check")

    elif len(fArray) == mesh.nodeCount():
        ### nodal values
        fA = pg.Vector(fArray)
        b_l = pg.matrix.ElementMatrix()
        for c in mesh.cells():
            b_l.u(c)
                # rhs.addVal(b_l.row(0) * fArray[b_l.idx()], b_l.idx())
            rhs.add(b_l, fA)

        # print("test reference solution:")
        # rhsRef = pg.Vector(mesh.nodeCount(), 0)
        # for c in mesh.cells():
        #     b_l.u(c)
        #     for i, idx in enumerate(b_l.idx()):
        #         rhsRef[idx] += b_l.row(0)[i] * fA[idx]

        # np.testing.assert_allclose(rhs, rhsRef)
        # print("Remove revtest in assembleLoadVector after check",
        #       sum(rhs), sum(rhsRef))

            # rhs = pg.Vector(fArray)
    else:
        raise Exception("Load vector have the wrong size: " +
                        str(len(fArray)))

    return rhs


def _assembleUDirichlet(mat, rhs, uDirIndex, uDirichlet):
    """This should be moved directly into the core"""

    if rhs is not None:
        uDir = pg.Vector(mat.rows(), 0.0)
        uDir.setVal(uDirichlet, uDirIndex)
        rhs -= mat * uDir

    for i in uDirIndex:
        mat.cleanRow(i)
        mat.cleanCol(i)
        mat.setVal(i, i, 1.0)

    if rhs is not None:
        rhs[uDirIndex] = uDirichlet
        #rhs.setVal(uDirichlet, uDirIndex)


def assembleDirichletBC(mat, boundaryPairs, rhs=None, time=0.0, userData={},
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
    mat: :gimliapi:`GIMLI::RSparseMatrix`
        System matrix of the system equation.

    boundaryPair: list()
        List of pairs [:gimliapi:`GIMLI::Boundary`, h].
        The value :math:`h` will assigned to the nodes of the boundaries.
        Later assignment overwrites prior.

        :math:`h` need to be a scalar value (float or int) or
        a value generator function that will be executed at runtime.
        See :py:mod:`pygimli.solver.solver.parseArgToBoundaries`
        and :ref:`tut:modelling_bc` for example syntax,

    nodePairs: list()
        List of pairs [ nodeID, uD ].
        The value uD will assigned to the nodes given there ids.
        This node value settings will overwrite any prior settings due to
        boundaryPair.

    rhs: :gimliapi:`GIMLI::RVector`
        Right hand side vector of the system equation will bet set to
        :math:`u_{\text{D}}`

    time: float
        Will be forwarded to value generator.

    userData: class
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
            if isinstance(boundary, pg.core.Node):
                n = boundary
                uDirNodes.append(n)
                uDirVal[n.id()] = uD
            else:
                for i, n in enumerate(boundary.nodes()):
                    uDirNodes.append(n)
                    uDirVal[n.id()] = uD

                    ###  don't the case this could be usefull for scalar fields
                    # if hasattr(uD, '__iter__'):
                    #     uDirVal[n.id()] = uD[i]
                    # else:

    if len(uDirNodes) == 0 and nodePairs is None:
        return

    uniqueNodes = unique(uDirNodes)

    uDirichlet = []
    uDirIndex = []

    for i, n in enumerate(uniqueNodes):
        ud = uDirVal[n.id()]
        if isinstance(ud, list):
            # vector valued problem
            if mat.size() % len(ud) != 0:
                print(mat.size(), len(ud))
                pg.error("Matrix size missmatch for vector valued problem")
            else:
                dof = mat.size() // len(ud)
                for d in range(len(ud)):
                    if ud[d] is not None:
                        uDirIndex.append(n.id() +  d * dof)
                        uDirichlet.append(ud[d])

        else:
            uDirIndex.append(n.id())
            uDirichlet.append(uDirVal[n.id()])

    if nodePairs is not None:
        #print("nodePairs", nodePairs)

        if len(nodePairs) == 2 and isinstance(nodePairs[0], int):
            # assume a single Node [NodeId, val]
            nodePairs = [nodePairs]

        for i, [n, val] in enumerate(nodePairs):
            uDirIndex.append(n)
            if hasattr(val, '__call__'):
                pg.error("callable node pairs need to be implement.")
            uDirichlet.append(val)

    _assembleUDirichlet(mat, rhs, uDirIndex, uDirichlet)


def assembleNeumannBC(rhs, boundaryPairs, a=None, time=0.0, userData={}):
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

    rhs : :py:mod:`Vector`
        Right hand side vector of length node count.

    boundaryPair : list()
        List of pairs [ :gimliapi:`GIMLI::Boundary`, g ].
        The value :math:`g` will assigned to the nodes of the boundaries.
        Later assignment overwrites prior.

        :math:`g` need to be a scalar value (float or int) or
        a value generator function that will be executed at run time.

        See :py:mod:`pygimli.solver.solver.parseArgToBoundaries`
        and :ref:`tut:modelling_bc` for example syntax,

    a : iterable
        Per cell values to scale the Neumann part regarding weak formulation.

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

    Se = pg.matrix.ElementMatrix()

    for pair in boundaryPairs:
        boundary = pair[0]
        val = pair[1]
        g = generateBoundaryValue(boundary, val, time, userData)

        if a is not None:
            try:
                g *= a[boundary.leftCell().id()]
            except BaseException as e:
                print(boundary.leftCell())
                print(boundary.leftCell().id())
                print(a)
                pg.warn('Insufficient cell information.')

        if g != 0.0 and g is not None:
            Se.u(boundary)
            if isinstance(rhs, pg.Vector):
                rhs.add(Se, g)
            else:
                for i, j in enumerate(Se.ids()):
                    rhs[j] += Se.row(0)[i] * g


def assembleRobinBC(mat, boundaryPairs, rhs=None, time=0.0, userData={}):
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
        * fixme, test me

    Parameters
    ----------
    mat: :gimliapi:`GIMLI::SparseMatrix`
        System matrix of the system equation.
    boundaryPair: list()
        List of pairs [:gimliapi:`GIMLI::Boundary`, :math:`\alpha`].
        The value :math:`\alpha` will assigned to the nodes of the boundaries.
        Later assignment overwrites prior.

        :math:`\alpha` needs to be a scalar value (float or int) or
        a value generator (callable) which will be executed at run time.
        See :py:mod:`pygimli.solver.solver.parseArgToBoundaries`
        and :ref:`tut:modelling_bc` for example syntax,
    time: float
        Will be forwarded to value generator.
    userData: dict
        Will be forwarded to value generator.
    """
    if not hasattr(boundaryPairs, '__getitem__'):
        raise BaseException("Boundary pairs need to be a list of "
                            "[boundary, value]")

    Sp = pg.matrix.ElementMatrix()
    Sq = pg.matrix.ElementMatrix()

    #if isinstance(rhs, np.ndarray):
        #rhs = pg.Vector(rhs)

    for pair in boundaryPairs:
        boundary = pair[0]
        val = pair[1]
        # BC = a u + b * du/dn = g
        # solve
        # either du/dn = g/b for a=0 (Neumann)
        # or u = g/a for b = 0 (Dirichlet)
        ### p = gamma / beta
        p = generateBoundaryValue(boundary, val, time, userData)
        #### p = gamma / alpha
        #p = 20.
        #q = -41.0
        q = None

        if p != 0.0 and p is not None:
            #print(p)
            Sp.u2(boundary)
            mat.add(Sp, p)
            #Sp *= p
            #S += Sp

            if q is not None:
                Sq.u(boundary)
                rhs.add(Sq, -p*q)


def assembleBC_(bc, mesh, mat, rhs, a, time=None, userData={}):
    r"""Shortcut to apply all boundary conditions.

    This is a helper function for the solver call.
    Shortcut to apply all boundary conditions will only forward to
    appropriate assemble functions.
    """
    ## we can't iterate because we want the following fixed order
    bct = dict(bc)
    if 'Neumann' in bct:
        assembleNeumannBC(rhs, parseArgToBoundaries(bct.pop('Neumann'), mesh),
                          a=a, time=time, userData=userData)
    if 'Robin' in bct:
        assembleRobinBC(mat, parseArgToBoundaries(bct.pop('Robin'), mesh),
                        rhs=rhs, time=time, userData=userData)
    if 'Dirichlet' in bct:
        assembleDirichletBC(mat,
                            parseArgToBoundaries(bct.pop('Dirichlet'), mesh),
                            rhs=rhs, time=time, userData=userData)
    if 'Node' in bct:
        assembleDirichletBC(mat, [], nodePairs=bct.pop('Node'),
                            rhs=rhs, time=time, userData=userData)

    if len(bct.keys()) > 0:
        pg.warn("Unknown boundary condition[s]" + \
                       str(bct.keys()) + " will be ignored")


def createStiffnessMatrix(mesh, a=None, isVector=False):
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

    a : iterable of type float, int, complex, RMatrix, CMatrix
        Per cell values., e.g., physical parameter. Length of a need to be mesh.cellCount().
        If None given default is 1.

    isVector : bool [False]
        We want to solve for vector valued problems. Resulting SparseMatrix is a SparseMapMatrix and have the dimension
        (nNodes * nDims, nNodes * nDims) with nNodes = mesh.nodeCount() and
        nDims = mesh.dimension().

    Returns
    -------
    A : :gimliapi:`GIMLI::[C]SparseMatrix` | [C]SparseMapMatrix
        Stiffness matrix, with real or complex values.
    """
    if mesh.cellCount() == 0:
        print(mesh)
        raise Exception("Mesh invalid")

    if a is None:
        a = pg.Vector(mesh.cellCount(), 1.0)

    A = None

    if isVector is False:

        if isinstance(a[0], float) or isinstance(a[0], np.float64) and 0:
            A = pg.matrix.SparseMatrix()
            A.fillStiffnessMatrix(mesh, a)
            return A

        dof = 0
        nDof  = mesh.nodeCount()
    else:
        dof = mesh.nodeCount()
        nDof  = mesh.nodeCount() * mesh.dimension()

    #### if vector or scalar(Complex)
    if np.array(a[0]).dtype == np.complex:
        isComplex = True
        A = pg.matrix.CSparseMapMatrix(nDof, nDof)
    else:
        isComplex = False
        A = pg.matrix.SparseMapMatrix(nDof, nDof)

    al = pg.core.ElementMatrix(dof=dof)

    if len(a) != mesh.cellCount():
        pg.error('Number of cell values need to match cell count')

    for c in mesh.cells():
        if isComplex is True:
            #al.gradU2(c, 1.0)
            al.ux2uy2uz2(c)
            A.add(al, scale=a[c.id()])
        else:
            al.gradU2(c, a[c.id()])
            A.add(al)

    if isComplex is True:
        return pg.matrix.CSparseMatrix(A)
    return pg.matrix.SparseMatrix(A)


def createMassMatrix(mesh, b=None):
    r"""Create the mass matrix.

    Calculates the Mass matrix (Finite element identity matrix) the given mesh.

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
        b = pg.Vector(mesh.cellCount(), 1.0)
    elif not hasattr(b, '__iter__'):
        b = pg.Vector(mesh.cellCount(), b)

    B = pg.matrix.SparseMatrix()
    B.fillMassMatrix(mesh, b)
    return B

    # create matrix structure regarding the mesh
    # B.buildSparsityPattern(mesh)
    # define a local element matrix
    # B_l = pg.matrix.ElementMatrix()
    # for c in mesh.cells():
    #    B_l.u2(c)
    # # check if b[i] == B*b
    #    B_l *= b[c.id()]
    #    B += B_l
    # return B


def _feNorm(u, mat):
    """Create a norm within a Finite Element space.

    Create the Finite Element Norm with a preassembled system matrix.
    """
    return np.sqrt(pg.math.dot(u, mat.mult(u)))


def normL2(u, mat=None, mesh=None):
    r"""Create Lebesgue (L2) norm for finite element space.

    Find the L2 Norm for a solution for the finite element space. :math:`u` exact solution
    :math:`{\bf M}` Mass matrix, i.e., Finite element identity matrix.

    .. math::

        L2(f(x)) = || f(x) ||_{L^2} & = (\int |f(x)|^2 \d x)^{1/2} \\
                                    & \approx h (\sum |f(x)|^2 )^{1/2} \\
        L2(u) = || u ||_{L^2} & = (\int |u|^2 \d x)^{1/2} \\
                              & \approx (\sum M (u)) ^{1/2} \\
        e_{L2_rel} = \frac{L2(u)}{L2(u)} & =
                               \frac{(\sum M(u))^{1/2}}{(\sum M u)^{1/2}}

    The error for any approximated solution :math:`u_h` correlates to the L2
    norm of 'L2Norm(u - u_h, M)'. If you like relative values, you can also
    normalize this error with 'L2Norm(u - u_h, M) / L2Norm(u, M)*100'.

    Parameters
    ----------
    u : iterable
        Node based value to compute the L2 norm for.

    mat : Matrix
        Mass element matrix.

    mesh : :gimliapi:`GIMLI::Mesh`
        Mesh with the FE space to generate M if necessary.

    Returns
    -------
    ret : float
        :math:`L2(u)` norm.

    """
    if isinstance(mat, pg.Mesh):
        mesh = mat
        mat = None

    if mat is None and mesh is not None:
        mat = createMassMatrix(mesh)

    if mat is None:
        # M is Identity matrix
        return np.sqrt(pg.math.dot(u, u))

    return _feNorm(u, mat)


def normH1(u, mat=None, mesh=None):
    r"""Create (H1) norm for the finite element space.

    Parameters
    ----------
    u : iterable
        Node based value to compute the H1 norm for.

    mat : Matrix
        Stiffness matrix.

    mesh : :gimliapi:`GIMLI::Mesh`
        Mesh with the FE space to generate S if necessary.

    Returns
    -------
    ret : float
        :math:`H1(u)` norm.

    """
    if isinstance(mat, pg.Mesh):
        mesh = mat
        mat = None

    if mat is None and mesh is not None:
        mat = pg.solver.createStiffnessMatrix(mesh)

    if mat is None:
        raise Exception("Need Stiffness matrix or a mesh here to calculate H1-Norm")

    return _feNorm(u, mat)


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


def solveFiniteElements(mesh, a=1.0, b=None, f=0.0, bc=None,
                        times=None, userData={},
                        verbose=False, **kwargs):
    r"""Solve partial differential equation with Finite Elements.

    This is a syntactic sugar convenience function for using the Finite Element
    functionality of the library core to solve partial differential equation (PDE) that match the following form:

    .. math::

        \frac{\partial u}{\partial t} & = \nabla\cdot(a \nabla u)
        + b u + f(\mathbf{r},t)~~|~~\Omega_{\text{Mesh}}\\
        u & = h~~|~~\Gamma_{\text{Dirichlet}}\\
        \frac{\partial u}{\partial \mathbf{n}} & = g~~|~~\Gamma_{\text{Neumann}}\\
        \alpha u + \beta\frac{\partial u}{\partial \mathbf{n}} & = \gamma~~|~~\Gamma_{\text{Robin}}

    for the scalar :math:`u(\mathbf{r}, t)` or
    vectoral :math:`\mathbf(u)(\mathbf{r}, t)` solution at each node of a given mesh.
    The Domain :math:`\Omega` and the Boundary :math:`\Gamma` are defined
    through the mesh with appropriate boundary marker.

    Note, to ensure vector solution either set vector forces or at least on vector component boundary condition.

    TODO:

        * unsteady ub and dub
        * 'Infinity' Boundary condition (u vanishes at infinity)
        * 'Cauchy' Boundary condition
        (guaranties u and du on same boundary, will never work here because the problem becomes ill posed and would need some inverse strategy to solve.)
        * Example for
            * elastic parameter
            * anisotropic (float/complex)
            * dynamic boundary conditions
            * dynamic load vector
            * nonlinearity

    Parameters
    ----------
    mesh: :gimliapi:`GIMLI::Mesh`
        Mesh represents spatial discretization of the calculation domain
    a: value | array | callable(cell, userData)
        Cell values of type float or complex can be scalar, anisotropic matrix or elastic tensor.
    b: value | array | callable(cell, userData) [None]
        Cell values. None means the term is unused.
    f: value | array(cells) | array(nodes) | callable(args, kwargs)
        force values, for vector fields use (n x dim) values.
    bc: dict()
        Dictionary of boundary conditions.
        Current supported boundary conditions are by dictionary keys:
        'Dirichlet', 'Neumann', 'Robin'.

        The dictionary can contain multiple "key: Arg"
        Arg will be parsed by :py:mod:`pygimli.solver.solver.parseArgPairToBoundaryArray`

        If the dictionary key is 'Node' then fixed values for single node
        indices can by be given. e.g., bc={'Node': [nodeID, value]}.
        Note this is only a shortcut for
        bc={'Dirichlet': [mesh.node(nodeID), value]}.
    times: array [None]
        Solve as time dependent problem for the given times.

    Other Parameters
    ----------------
    **kwargs
        u0: value | array | callable(pos, userData)
            Node values
        theta: float [1]
            - :math:`theta = 0` means explicit Euler, maybe stable for
            :math:`\Delta t \quad\text{near}\quad h`
            - :math:`theta = 0.5`, Crank-Nicolson, maybe instable
            - :math:`theta = 1`, implicit Euler

            If unsure choose :math:`\theta = 0.5 + \epsilon`, which is probably stable.
        dynamic: bool [False]
            Boundary conditions for time depending problems will be considered
            dynamic for each time step.
        stats: bool
            Give some statistics.
        progress: bool
            Give some calculation progress.
        assembleOnly: bool
            Stops after matrix asssemblation.
            Returns the system matrix A and the rhs vector.
        ws: dict
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
        bc = {}

    workSpace = kwargs.pop('ws', dict())
    debug = kwargs.pop('debug', False)
    stats = kwargs.pop('stats', False)

    mesh.createNeighborInfos()
    if verbose:
        print("Mesh: ", str(mesh))

    ## scalar solution default
    vectorValues = False
    dof = mesh.nodeCount()

    ## check if force vector is a vector
    rhs = createLoadVector(mesh, f, userData=userData)

    if len(rhs) > dof:
        if verbose:
            print("Solve vector valued.")
        vectorValues = True
        dof = mesh.nodeCount() * mesh.dimension()

    ## check if some boundary conditions are vector valued



    swatch = pg.core.Stopwatch(True)

    # check for material parameter
    #a = parseArgToArray(a, nDof=mesh.cellCount(), mesh=mesh, userData=userData)
    a = cellValues(mesh, a, userData=userData)
    S = createStiffnessMatrix(mesh, a, isVector=vectorValues)
    M = None

    if b is not None and b != 0:
        b = cellValues(mesh, b, userData=userData)
        M = createMassMatrix(mesh, b)
        #pg.warn("check me")
        A = S - M
    else:
        A = S

    if times is None:
        assembleBC_(bc, mesh, A, rhs, a, time=None, userData=userData)

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
                rhs = pg.math.toComplex(rhs)
            else:
                if u is None:
                    u = pg.Vector(rhs.size(), 0.0)

        assembleTime = swatch.duration(True)

        if stats:
            stats.assembleTime = assembleTime

        if verbose:
            print("Assemblation time: ", assembleTime)

        workSpace['S'] = S
        workSpace['M'] = M
        workSpace['A'] = A
        workSpace['rhs'] = rhs

        if 'assembleOnly' in kwargs:
            return A, rhs

        solver = pg.core.LinSolver(False)
        solver.setMatrix(A, 0)

        if singleForce:
            u = solver.solve(rhs)
        else:
            for i, r in enumerate(rhs):
                u[i] = solver.solve(r)

        solverTime = swatch.duration(True)
        if verbose:
            if stats:
                stats.solverTime = solverTime
            print("Solving time: ", solverTime)

        if len(kwargs.keys()) > 0:
            pg.warn("Unused arguments", *kwargs)
        return u

    else: # times given

        pg.solver.checkCFL(times, mesh, max(a))

        if debug:
            print("start TL", swatch.duration())

        M = createMassMatrix(mesh)
        F = createLoadVector(mesh, f)

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
            assembleBC_(bc, mesh, S, F, a, time=0.0, userData=userData)
            return crankNicolson(times, theta, S, M, F, u0=u0, progress=progress)

        rhs = np.zeros((len(times), dof))
        # no time dependency for rhs so far ... TODO
        rhs[:] = F  # this is slow: optimize

        if debug:
            print("rhs", swatch.duration())
        U = np.zeros((len(times), dof))
        U[0, :] = u0

        # init state
        u = pg.Vector(dof, 0.0)

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
            br = (M + (dt * (theta - 1.)) * S) * U[n - 1] + \
                dt * ((1.0 - theta) * rhs[n - 1] + theta * rhs[n])

            # print ('a',swatch.duration(True))
            # br = M * U[n - 1] - (A * U[n - 1]) * (dt*(1.0 - theta)) + \
            # dt * ((1.0 - theta) * rhs[n - 1] + theta * rhs[n])

            # print ('br',swatch.duration(True))

            # br = M * U[n - 1] - (dt*(1.0 - theta)) * A * U[n - 1] + \
            # dt * ((1.0 - theta) * rhs[n - 1] + theta * rhs[n])
            # print ('c',swatch.duration(True))

            measure += swatch.duration()

            A = M + S * dt * theta

            assembleBC_(bc, mesh, A, br, a, time=times[n], userData=userData)

            if 'assembleOnly' in kwargs:
                return A, br

            # u = S/b
            t_prep = swatch.duration(True)
            solver = pg.core.LinSolver(A, verbose)
            solver.solve(br, u)

            if 'plotTimeStep' in kwargs:
                kwargs['plotTimeStep'](u, times[n])

            U[n, :] = np.asarray(u)

            if progress:
                progress.update(n, ' t_prep: {0}ms t_step {1}s'.format(
                    pg.pf(t_prep*1000),
                    pg.pf(swatch.duration())))
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
            pg.warn("Courant-Friedrichs-Lewy Number:", c,
                    "but should be lower 1 to ensure movement inside a cell "
                    "per timestep. ("
                    "vMax =", vMax,
                    "dt =", dt,
                    "dx =", dx,
                    "dt <", dx/vMax,
                    " | N > ", int((times[-1]-times[0])/(dx/vMax))+1, ")")
    return c

def crankNicolson(times, theta, matS, matI, f, u0=None, progress=None):
    """
        S = constant over time
        f = constant over time
    """
    if len(times) < 2:
        raise BaseException("We need at least 2 times for "
                            "Crank-Nicolsen time discretization." + str(len(times)))
    # sw = pg.core.Stopwatch(True)

    if u0 is None:
        u0 = np.zeros(len(f))

    u = np.zeros((len(times), len(f)))
    u[0, :] = u0
    dt = float(times[1] - times[0])

    rhs = np.zeros((len(times), len(f)))
    rhs[:] = f

    timeAssemble = []
    timeSolve = []

    timeMeasure = False
    if progress:
        timeMeasure = True

    A = matI
    if theta > 0:
        A = matI + matS * (dt * theta)

    solver = pg.core.LinSolver(A, verbose=False)

    St = matI - matS * dt # cache what is possible the theta=0

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
            b = St * u[n-1] + dt * rhs[n-1]
        else:
            b = matI * u[n-1] + matS.mult(dt * (theta - 1.) * u[n-1]) + \
                dt * ((1.0 - theta) * rhs[n-1] + theta * rhs[n])

#        pg.toc()
#        print(np.linalg.norm(b-b1))
        #np.testing.assert_allclose(bRef, b)

        if timeMeasure:
            timeAssemble.append(pg.dur())
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
        self.resU = None
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
        self.resU = deepcopy(u0)

        if isinstance(self.resU, list):
            for r in self.resU:
                r *= 0.0
        else:
            self.resU *= 0.0

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
            # classical 4 step Runge-Kutta rk4
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

                if isinstance(self.resU, list):
                    for i in range(len(self.resU)):
                        self.resU[i] = self.rk4a[jRK] * self.resU[i] + \
                                       self.dt * rhs[i]

                        self.u[i] += self.rk4b[jRK] * self.resU[i]
                else:
                    self.resU = self.rk4a[jRK] * self.resU + self.dt * rhs
                    self.u += self.rk4b[jRK] * self.resU

        self.time += self.dt
        return self.u


if __name__ == "__main__":
    pass
