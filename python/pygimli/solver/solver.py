#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pygimli as pg

from pygimli.utils import unique
import numpy as np

import copy

def parseArgToArray(arg, ndof, mesh=None, userData=None):
    """
    Parse array related arguments to create a valid value array.

    Parameters
    ----------
    arg : float | int | iterable | callable
        The target array value that will be converted to an array.
        
        If arg is a callable with it must fulfill:
        
        >>> arg(MeshEntity, userData=None)
        
        Where MeshEntity is one of
        :gimliapi:`GIMLI::Cell` , 
        :gimliapi:`GIMLI::Node` or 
        :gimliapi:`GIMLI::Boundary`       
        depeding on ndof, where ndof is mesh.cellCount(), 
        mesh.nodeCount() or mesh.boundaryCount(),
        respectively.
        
    ndof : int | [int]
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
  
    if not hasattr(ndof, '__len__'):
        nDofs = [ndof]
    else:
        nDofs = ndof
    
    try:
        return pg.RVector(nDofs[0], float(arg))
    except:
        pass

    if hasattr(arg, '__len__'):
        for n in nDofs:
            if len(arg) == n:
                return arg
        
        raise Exception("Array 'arg' has the wrong size: " + 
                        str(len(arg)) + " != " +  str(dof))
    elif hasattr(arg, '__call__'):
        ret = pg.RVector(nDofs[0], 0.0)
        
        if not mesh:
            raise Exception("Please provide a mesh for the callable"
                            "argument to parse ")
        
        if nDofs[0] == mesh.nodeCount():
            for n in mesh.nodes():
                if userData:
                    ret[n.id()] = arg(n.pos(), userData)
                else:
                    ret[n.id()] = arg(n.pos())
        elif nDofs[0] == mesh.cellCount():
            for c in mesh.cells():
                if userData:
                    ret[c.id()] = arg(c, userData)
                else:
                    ret[c.id()] = arg(c)
        elif nDofs[0] == mesh.boundaryCount():
            for b in mesh.boundaries():
                if userData:
                    ret[b.id()] = arg(b, userData)
                else:
                    ret[b.id()] = arg(b)
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
        
        If arg is a callable with it must fulfill:
        
        >>> arg(:gimliapi:`GIMLI::Boundary`, time=0.0, userData=None)

        and should return an appropriate value.
    
    Returns
    -------S
    val : float or list of ..
    """
    if hasattr(boundary, '__len__'):
        vals = np.zeros(len(boundary))
        for i, b in enumerate(boundary):
            vals[i] = generateBoundaryValue(b, arg[i], time, userData)
        return vals
    val = 0.
    
    if hasattr(arg, '__call__'):
        kwargs = dict()
        args = [boundary]
        
        if time != 0.0:
            args.append(time)
        if userData:
            kwargs['userData'] = userData
        val = arg(*args, **kwargs)
    elif hasattr(arg, '__len__'):
        val = generateBoundaryValue(boundary, arg[boundary.id()], userData)
    else:
        try:
            val = float(arg)
        except ValueError:
            raise(arg, " cannot be converted to a float")
    return val

def parseArgPairToBoundaryArray(pair, mesh):
    """
    Parse boundary related pair argument to 
    [ :gimliapi:`GIMLI::Boundary`, value|callable ] list.
    
    Parameters
    ----------
    
    pair : tuple
    
        - [marker, arg]
        - [[boundary,...], arg]
        
        arg will be parsed by
        :py:mod:`pygimli.solver.solver.generateBoundaryValue` 
        and distributed to each boundary.
        Callable functions will be executed at runtime. 
       
    mesh : :gimliapi:`GIMLI::Mesh`
        Used to find boundaries by marker
         
    Returns
    -------
    
    boundaries : list()
        [ :gimliapi:`GIMLI::Boundary`, value|callable ]
    """
    boundaries = []
    bounds = []
    if type(pair[0]) == int:
        bounds = mesh.findBoundaryByMarker(pair[0])
    elif type(pair[0]) == pg.stdVectorBounds:            
        bounds = pair[0]
        
    for b in bounds:
        val = None
        if hasattr(pair[1], '__call__'):
            # don't execute the callable here in init, 
            # we want to call them at runtime
            val = pair[1]
        else:
            val = generateBoundaryValue(b, pair[1])
        boundaries.append([b, val])
   
    return boundaries
                
def parseArgToBoundaries(args, mesh):
    """
    Parse boundary related arguments to create a valid boundary value list: 
    [ :gimliapi:`GIMLI::Boundary`, value|callable ]
    
    Parameters
    ----------
    
    args : pair | [pair, ...]
        see :py:mod:`pygimli.solver.solver.parseArgPairToBoundaryArray` 
        
    mesh : :gimliapi:`GIMLI::Mesh`
        Used to find boundaries by marker
        
    Returns
    -------
    
    boundaries : list()
        [ :gimliapi:`GIMLI::Boundary`, value|callable ]
    """
    boundaries = []
    if type(args) == list:
        if len(args) == 2:
            try:
                #[[,], [,]]
                if len(args[0]) == 2 and len(args[1]) == 2:
                    boundaries += parseArgPairToBoundaryArray(args[0], mesh)
                    boundaries += parseArgPairToBoundaryArray(args[1], mesh)
                else:
                    boundaries += parseArgPairToBoundaryArray(args, mesh)
            except:
                #[,]
                boundaries += parseArgPairToBoundaryArray(args, mesh)
        else:
            #[[,], [,], ...]
            for arg in args:
                boundaries += parseArgPairToBoundaryArray(arg, mesh)
    
    return boundaries 

def parseMapToCellArray(attributeMap, mesh):
    """
    Parse a value map to cell attributes.
    
    Parameters
    ----------
    mesh : :gimliapi:`GIMLI::Mesh`
        For each cell of mesh a value will be returned.
    
    attributeMap : list | dict
        List of pairs [[int, float]], or dictionary with integer keys
        
    Returns
    -------
    atts : array
        Array of length mesh.cellCount()
    """
    
    atts = pg.RVector(mesh.cellCount(), 0.0)
    
    if type(attributeMap) is dict():
        raise Exception("Please implement me!")
    elif hasattr(attributeMap, '__len__'):
        for pair in attributeMap:
            if hasattr(pair, '__len__'):
                idx = pg.find(mesh.cellMarker() == pair[0]) 
                if len(idx) == 0:
                    print("Warning! parseMapToCellArray: cannot find marker "+\
                        str(pair[0]) + " within mesh.")
                else:
                    atts[idx] = float(pair[1])
            else:
                raise Exception("Please provide a list of [int, value] pairs!" + 
                                str(pair))
    else:
        print("attributeMap: ", attributeMap)
        raise Exception("Cannot interpret attributeMap!")
            
    return atts   
   
def fillEmptyToCellArray(mesh, vals):
    """
    Prolongate empty cell values to complete cell attributes.
    
    Its is possible that you have zero values that need to be filled with
    appropriate attributes. 
    This function tries to fill the empty values successive prolongation of the 
    non zeros. 
    
    Parameters
    ----------
    mesh : :gimliapi:`GIMLI::Mesh`
        For each cell of mesh a value will be returned.
    
    vals : array
        Array of size cellCount().
    
    Returns
    -------
    atts : array
        Array of length mesh.cellCount()
    """
    atts = pg.RVector(mesh.cellCount(), 0.0)
    oldAtts = mesh.cellAttributes()
    mesh.setCellAttributes(vals)
    
    #std::vector< Cell * > 
    #empties = []
    
    mesh.fillEmptyCells(mesh.findCellByAttribute(0.0), background=-1)
    atts = mesh.cellAttributes()
    mesh.setCellAttributes(oldAtts)
    return atts;
       
       

def divergence(mesh, F=lambda r:r, order=1):
    """ 
    MOVE THIS to a better place
    
    Parameters
    ----------
    
    Returns
    -------
    """
    div = 0
    directionCheck = False
    
    if mesh.cellCount() > 0:
        directionCheck = True
        
    for b in mesh.boundaries():
               
        if directionCheck:
            if b.leftCell() is None and b.rightCell() is None:
                print(b.id(), b.leftCell(), b.rightCell())  
                sw = g.Stopwatch(True)
                mesh.createNeighbourInfos()
                print("NeighbourInfos()", sw.duration(True))
                ##return gauss(grid, F)
              
            if not b.leftCell() is None and not b.rightCell() is None: continue
               
        tmpdiv = 0
        shape = b.shape()
        
        if order == 1:
            tmpdiv = shape.norm().dot(F(shape.center())) * shape.domainSize()
        else:
            weights = g.IntegrationRules.instance().weights(shape, order)
            abscissa = g.IntegrationRules.instance().abscissa(shape, order)
                    
            for i, p in enumerate(abscissa): 
                rPos = shape.xyz(p)
                tmpdiv += shape.norm().dot(F(rPos)) * \
                    weights[i] * shape.domainSize()
       
        
        if directionCheck and b.leftCell() is None:
            tmpdiv *= -1
            # raise Exception("invalid mesh: left is None .. every 
            # boundary need leftCell")
                    
        div += tmpdiv
        
    return div

def triDiagToeplitz(dom, a, l, r, start=0, end=-1):
    """ WHATSTHIS? """
    A = pg.RSparseMapMatrix(dom, dom)
    
    if end == -1: end = dom
    
    for i in range(start, end):
        A.addVal(i, i, a)
        if i > start: 
            A.addVal(i, i - 1, l)  
            
        if i < end-1: 
            A.addVal(i, i + 1, r)  
    return A
        
def identity(dom, start=0, end=-1):
    """ WHATSTHIS? """
    A = pg.RSparseMapMatrix(dom, dom)
    
    if end == -1: end = dom
    
    for i in range(start, end):
        A.addVal(i, i, 1)
    return A
    
def showSparseMatrix(A):
    """ helper function """
    S = A
    #S = pg.RSparseMatrix(A)
    rows = S.vecRowIdx()
    cols = S.vecColPtr()
    vals = S.vecVals()

    for i in range(S.rows()):
        for j in range(cols[i], cols[i + 1]):
            print(i, rows[j], vals[j])

def linsolve(A, b, verbose=False):
    r"""
    Direct solution after :math:`\textbf{x}` using cholmod:
    
    .. math::
        \textbf{A}\textbf{x} = \textbf{b}
        
    If :math:`\textbf{A}` is symmetric, sparse and positive definite.
            
    Parameters
    ----------
    A : :gimliapi:`GIMLI::RSparseMatrix` | :gimliapi:`GIMLI::RSparseMapMatrix` |
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
    x = pg.RVector(len(b), .0 )
    
    if type(A) == pg.RSparseMapMatrix:
        S = pg.RSparseMatrix(A)
        solver = pg.LinSolver(S, verbose=verbose)
        solver.solve(b, x)
    elif type(A) == np.ndarray:
        return np.linalg.solve(A, b)
    else:    
        solver = pg.LinSolver(A, verbose=verbose)
        solver.solve(b, x)
        
    return x
    
def assembleForceVector(mesh, f, userData=None):
    """

    Create right hand side vector based on the given mesh and force values.
        
        
    Parameters
    ----------
    f: float, array, callable(cell, [userData])
    
        Force Values
        float -> ones(mesh.nodeCount()) * vals, 
        for each node [0 .. mesh.nodeCount()]
        for each cell [0 .. mesh.cellCount()]
    """

    rhs = pg.RVector(mesh.nodeCount(), 0)
        
    b_l = pg.ElementMatrix()
    
    if hasattr(f, '__call__') and type(f) is not pg.RVector:
        for c in mesh.cells():
            if userData is not None:
                f(c, rhs, userData)
            else:
                f(c, rhs)
    else:
        
        if type(f) == float or type(f) == int:
            fArray = pg.RVector(mesh.cellCount(), f)
        else:
            fArray = f
            
        if len(fArray) == mesh.cellCount():
            for c in mesh.cells():
                b_l.u(c)
                for i, idx in enumerate(b_l.idx()):
                    rhs[idx] += b_l.row(0)[i] * fArray[c.id()]
        elif len(fArray) == mesh.nodeCount():
            rhs = pg.RVector(fArray)
        else:
            raise Exception("Forcevector have the wrong size: " + \
                str(len(fArray)))

    return rhs
    
def assembleNeumannBC(S,
                      boundaryPairs, time=0.0, 
                      userData=None, verbose=False):
    r"""
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
        a value generator function that will be executed at runtime.
        See :py:mod:`pygimli.solver.solver.parseArgToBoundaries` 
                
        See tutorial section for an example,
        e.g., Modelling with Boundary Conditions
        
    time : float
        Will be forwarded to value generator.
        
    userData : class
        Will be forwarded to value generator.
    """

    Se = pg.ElementMatrix()

    if not hasattr(boundaryPairs, '__getitem__'):
        raise("Boundary pairs need to be a list of [boundary, value]")
    
    for pair in boundaryPairs:
        boundary = pair[0]
        val = pair[1]
        g = generateBoundaryValue(boundary, val, time, userData)
                
        if g is not 0.0:
            Se.u2(boundary)
            Se *= g
            S += Se

def assembleUDirichlet_(S, rhs, uDirIndex, uDirchlet):
    """ This should be moved directly into gimli """
    udirTmp = pg.RVector(S.rows(), 0.0)
    udirTmp.setVal(uDirchlet, uDirIndex)
    
    if rhs is not None:
        rhs -= S * udirTmp
        
    for i in uDirIndex:
    
        S.cleanRow(i)
        S.cleanCol(i)
        S.setVal(i, i, 1.0)
    
    if rhs is not None:
        rhs.setVal(uDirchlet, uDirIndex)
    
def assembleDirichletBC(S, boundaryPairs, rhs, time=0.0,
                        userData=None, verbose=False):
    r"""
    Apply Dirichlet boundary condition to the system matrix S and rhs vector.

    .. math::
        u(\arr{r}, t) = u_{\text{D}} \quad\text{with}\quad \arr{r}
        \quad\text{on}\quad \partial\Omega
    
    Parameters
    ----------
    S : :gimliapi:`GIMLI::RSparseMatrix`
        System matrix of the system equation.
      
    boundaryPair : list()
        List of pairs [ :gimliapi:`GIMLI::Boundary`, uD ].
        The value uD will assigned to the nodes of the boundaries.
        Later assignment overwrites prior.
        
        :math:`u_{\text{D}}` need to be a scalar value (float or int) or
        a value generator function that will be executed at runtime.
        See :py:mod:`pygimli.solver.solver.parseArgToBoundaries` 
                
        See tutorial section for an example,
        e.g., Modelling with Boundary Conditions
        
    rhs : :gimliapi:`GIMLI::RVector`
        Right hand side vector of the system equation will bet set to
        :math:`u_{\text{D}}`
        
    time : float
        Will be forwarded to value generator.
        
    userData : class
        Will be forwarded to value generator.
    """

    if not hasattr(boundaryPairs, '__getitem__'):
        raise("Boundary pairs need to be a list of [boundary, value]")
    
    uDirNodes = []
    uDirVal = dict()
    
    for pair in boundaryPairs:
        boundary = pair[0]
        val = pair[1]
        uD = generateBoundaryValue(boundary, val, time, userData)
    
        for n in boundary.nodes():
            uDirNodes.append(n)
            uDirVal[n.id()] = uD
                            
    if len(uDirNodes) == 0:
        return 
        
    uniqueNodes = unique(uDirNodes) 
        
    uDirchlet = pg.RVector(len(uniqueNodes))
    uDirIndex = []
    
    for i, n in enumerate(uniqueNodes):
        uDirIndex.append(n.id())
        uDirchlet[i] = uDirVal[n.id()]
    
    assembleUDirichlet_(S, rhs, uDirIndex, uDirchlet)       

def createStiffnessMatrix(mesh, a=None):
    """
    Assemble the stiffness matrix.
    
    Calculates the scaled stiffness matrix for the given mesh scaled with the per cell 
    values a.
    
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
    
    if type(a[0]) is float or type(a[0]) is np.float64:
        
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
        #A_l *= a[c.id()] 
        #A += A_l
        A.add(A_l, a[c.id()])
    
    return A

def createMassMatrix(mesh, b=None):
    """
    TODO remove b .. not necessary .. b should be scaled in final equation not here 
    Assemble mass element matrix.

    Calculates the mass element matrix for the given mesh scaled with the 
    per cell values b.
    
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
    
    #need callable here
    if b is None:
        b = pg.RVector(mesh.cellCount(), 1.0)
        
    B = pg.RSparseMatrix()
    B.fillMassMatrix(mesh, b)
    return B

    # create matrix structure regarding the mesh
    B.buildSparsityPattern(mesh)
    
    # define a local element matrix 
    B_l = pg.ElementMatrix()
    for c in mesh.cells():
        B_l.u2(c)
        # check if b[i] == B*b
        B_l *= b[c.id()]
        B += B_l
    
    return B

    

def solvePoisson(mesh, a=1.0, b=0.0, f=0.0, times=None, userData=None,
                 verbose=False, stats=None,  **kwargs):
    """
    WRITEME short
    
    WRITEME long
    
    Parameters
    ----------
    
    Returns
    -------
    """
    return solveFiniteElements(mesh, a, b, f, times, userData, verbose, stats,
                               **kwargs)

def solve(mesh, a=1.0, b=0.0, f=0.0, times=None, userData=None,
          verbose=False, stats=None, **kwargs):
    """
    WRITEME short
    
    WRITEME long
    
    Parameters
    ----------
    
    Returns
    -------
    """
    return solveFiniteElements(mesh, a, b, f, times, userData, verbose, stats,
                               **kwargs)

def solveFiniteElements(mesh, a=1.0, b=0.0, f=0.0, times=None, userData=None,
                        verbose=False, stats=None, **kwargs):
    """
    WRITEME short
    
    WRITEME long
    
    Parameters
    ----------
    mesh : :gimliapi:`GIMLI::Mesh`
        Mesh represents spatial discretization of the calculation domain
    
    a   : value | array | callable(cell, userData)
        Cell values
        
    b   : value | array | callable(cell, userData)
        Cell values
        
    u0 : value | array | callable(pos, userData)
        Node values
        
    f : value | array(cells) | array(nodes) | callable(args, kwargs)
        force values 
    
    times : array [None]
        solve as time dependent problem for the given times
        
    theta : float [1]
        - `theta` = 0, explicit Euler, maybe stable for
        - `theta` = 0.5, Crank-Nicolsen, maybe instable 
        - `theta` = 1, implicit Euler

        .. math:: \\Delta t \\quad\\text{near}\\quad h

        Time dependent equation is stable for:
        .. math:: 0.5 <= \\theta <= 1.0
        
        If unsure choose 0.5 + epsilon, which is probably be stable.
     
    progress : bool
        Give some calculation progress.

    Returns
    -------
    
    u : array
        Returns the solution u either 1,n array for stationary problems or 
        for m,n array for m time steps
        
    See Also
    --------
        
    other solver TODO
    """
    debug = kwargs.pop('debug', False)
        
    if verbose:
        print("Mesh: ", str(mesh))

    dof = mesh.nodeCount()
        
    swatch = pg.Stopwatch(True)
    swatch2 = pg.Stopwatch(True)
    
    # check for material parameter
    a = parseArgToArray(a, ndof=mesh.cellCount(), mesh=mesh, userData=userData)
    b = parseArgToArray(b, ndof=mesh.cellCount(), mesh=mesh, userData=userData)
        
    if debug: print("2: ", swatch2.duration(True))
    # assemble the stiffness matrix
    A = createStiffnessMatrix(mesh, a)

    if debug: print("3: ", swatch2.duration(True))
    M = createMassMatrix(mesh, b)
        
    if debug: print("4: ", swatch2.duration(True))
    S = A + M
    
    if debug: print("5: ", swatch2.duration(True))
    if times is None:
        
        rhs = assembleForceVector(mesh, f, userData=userData)
        
        if debug: print("6a: ", swatch2.duration(True))

        if 'duBoundary' in kwargs:
            print(userData)
            assembleNeumannBC(S,
                              parseArgToBoundaries(kwargs['duBoundary'], mesh),
                              time=0.0,
                              userData=userData,
                              verbose=False)

        if debug: print("6b: ", swatch2.duration(True))
        if 'uBoundary' in kwargs:
            assembleDirichletBC(S,
                                parseArgToBoundaries(kwargs['uBoundary'], mesh),
                                rhs, time=0.0,
                                userData=userData,
                                verbose=False)
            
        if debug: print("6c: ", swatch2.duration(True))
        if 'uDirichlet' in kwargs:
            raise("use uBoundary instead")
            assembleUDirichlet_(S, rhs,
                                kwargs['uDirichlet'][0],
                                kwargs['uDirichlet'][1])
            
        u = None    
        
        if type(a[0]) is complex:
            u = pg.CVector(rhs.size(), 0.0)
            rhs = pg.toComplex(rhs)
        else:
            u = pg.RVector(rhs.size(), 0.0)
            
        if debug: print("7: ", swatch2.duration(True))
        
        assembleTime = swatch.duration(True)
        if stats:
            stats.assembleTime = assembleTime
                
        if verbose:
            print(("Asssemblation time: ", assembleTime))

        #showSparseMatrix(S)
        
        solver = pg.LinSolver(True)
        solver.setMatrix(S, 0)
        
        u = solver.solve(rhs)
        
        solverTime = swatch.duration(True)
        if verbose:
            if stats:
                stats.solverTime = solverTime
            print(("Solving time: ", solverTime))
            
        return u
        
    else:
        if debug: print("start TL", swatch.duration())

        M = createMassMatrix(mesh)
        F = assembleForceVector(mesh, f)

        if 'u0' in kwargs:
            u0 = parseArgToArray(kwargs['u0'], dof, mesh, userData)
            
        theta = kwargs.pop('theta', 1)
        
        if not 'duBoundary' in kwargs:
            A = createStiffnessMatrix(mesh, a)
            
            if 'uBoundary' in kwargs:
                assembleDirichletBC(A,
                                    parseArgToBoundaries(kwargs['uBoundary'],
                                                         mesh),
                                    rhs=F)
            
            return crankNicolson(times, theta, A, M, F, u0=u0, verbose=verbose)
        
        rhs = np.zeros((len(times), dof))
        # rhs kann zeitabhängig sein ..wird hier nicht berücksichtigt
        rhs[:] = F # this is slow: optimize
        
        if debug: print("rhs", swatch.duration())
        U = np.zeros((len(times), dof))
        U[0, :] = u0
        
        #init state
        u = pg.RVector(dof, 0.0)
        
        
        if debug: print("u0", swatch.duration())
             
        measure = 0.
        for n in range(1, len(times)):
            swatch.reset()
            
            dt = times[n] - times[n - 1]
            
            # previous timestep
            #print "i: ", i, dt, U[i - 1]
            
            if 'duBoundary' in kwargs:
                # aufschreiben und checken ob neumann auf A oder auf S mit skaliertem val*dt angewendet wird
                A = createStiffnessMatrix(mesh, a)
                assembleNeumannBC(A,
                                  parseArgToBoundaries(kwargs['duBoundary'],
                                                       mesh),
                                  time=times[n],
                                  userData=userData,
                                  verbose=False)

            swatch.reset()
            # (A + a*B)u is fastest, followed by A*u + (B*u)*a and finally A*u + a*B*u and 
            b = (M + (dt * (theta - 1.)) * A ) * U[n - 1] + \
                dt * ((1.0 - theta) * rhs[n - 1] + theta * rhs[n])
            
    
            #print ('a',swatch.duration(True))
            #b = M * U[n - 1] - (A * U[n - 1]) * (dt*(1.0 - theta)) + \
                #dt * ((1.0 - theta) * rhs[n - 1] + theta * rhs[n])
            
            #print ('b',swatch.duration(True))
            
            #b = M * U[n - 1] - (dt*(1.0 - theta)) * A * U[n - 1] + \
                #dt * ((1.0 - theta) * rhs[n - 1] + theta * rhs[n])
            #print ('c',swatch.duration(True))
                        
            measure += swatch.duration()
            
            S = M + A * dt * theta
                
            if 'uBoundary' in kwargs:
                assembleDirichletBC(S,
                                    parseArgToBoundaries(kwargs['uBoundary'],
                                                         mesh),
                                    rhs=b,
                                    time=times[n],
                                    userData=userData,
                                    verbose=verbose)
                
            #u = S/b
            t_prep = swatch.duration(True)
            solver = pg.LinSolver(S, verbose)
            solver.solve(b, u)
            
            if 'plotTimeStep' in kwargs:
                kwargs['plotTimeStep'](u, times[n])
                    
            U[n,:] = np.asarray(u)

            if 'progress' in kwargs:
                if kwargs['progress']:
                    print(("\t" + str(n) +"/" + str(len(times)-1) + 
                          ": " + str(t_prep) +"/" + str(swatch.duration())))
        
        if debug: print("Measure("+str(len(times))+"): ",
                        measure, measure/len(times))
        return U
    
def crankNicolson(times, theta, S, I, f, u0=None, verbose=0):
    """
        S = const over time
        f = const over time
        
    """
    sw = pg.Stopwatch(True)
    swi = pg.Stopwatch(True)
    
    import matplotlib.pyplot as plt
    import numpy as np
    import time
    
    if u0 is None:
        u0 = np.zeros(len(f))
        
    u = np.zeros((len(times), len(f)))
    u[0, :] = u0
    dt = (times[1] - times[0])
            
    rhs = np.zeros((len(times), len(f)))
        
    rhs[:] = f
    
    A = (I + dt * theta * S)
    solver = pg.LinSolver(A, True)
    
    timeIter1 = np.zeros(len(times))
    timeIter2 = np.zeros(len(times))
    for n in range(1, len(times)):
        #if verbose:
            #print(n)
        tic = time.time()
        b = (I + (dt * (theta - 1.)) * S ) * u[n - 1] + \
            dt * ((1.0 - theta) * rhs[n - 1] + theta * rhs[n])
        
        timeIter1[n-1] = time.time()-tic
        
        tic = time.time()
        u[n, : ] = solver.solve(b)
        timeIter2[n-1] = time.time()-tic
    
        #A = (I + dt * theta * S)
        #u[n, : ] = linsolve(A, b)
            
    
    
    #plt.figure()
    #plt.plot(timeIter1)
    #plt.plot(timeIter2)
    #plt.figure()
    if verbose and (n%verbose == 0):
        print("timesteps: ", len(times), sw.duration(), "s ()", np.mean(timeIter1))
        
    return u