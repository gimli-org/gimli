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
                        len(arg) + " != " +  dof)
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

def generateBoundaryValue(arg, boundary, time=0.0, userData=None):
    """
    Generate a value for the given Boundary.
    
    Parameters
    ----------
    arg : convertible | iterable | callable
        
        - convertible into float 
        - iterable of minimum length = boundary.id()
        - callable generator function
        
        If arg is a callable with it must fulfill:
        
        >>> arg(:gimliapi:`GIMLI::Boundary`, time=0.0, userData=None)

        and should return an appropriate value.
        
    b : :gimliapi:`GIMLI::Boundary`
        The related boundary.
    
    Returns
    -------
    val : float
    """
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
        val = generateBoundaryValue(arg[boundary.id()], boundary, userData)
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
            val = generateBoundaryValue(pair[1], b)
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
            print(i,rows[j],vals[j])

def linsolve(A, b, verbose=False):
    r"""
    Direct solution after :math:`\textbf{x}` using cholmod:
    
    .. math::
        \textbf{A}\textbf{x} = \textbf{b}
        
    If :math:`\textbf{A}` is symmetric, sparse and positive definite.
        
    
    Parameters
    ----------
    A : :gimliapi:`GIMLI::RSparseMatrix` | :gimliapi:`GIMLI::RSparseMapMatrix`
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
        g = generateBoundaryValue(val, boundary, time, userData)
                
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
        uD = generateBoundaryValue(val, boundary, time, userData)
    
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
    return solveFEM(mesh, a, b, f, times, userData, verbose, stats,
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
    return solveFEM(mesh, a, b, f, times, userData, verbose, stats,
                     **kwargs)

def solveFEM(mesh, a=1.0, b=0.0, f=0.0, times=None, userData=None,
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
        
    theta : float [None]
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
        if type(a[0]) is float:
            u = pg.RVector(rhs.size(), 0.0)
        else:
            u = pg.CVector(rhs.size(), 0.0)
            rhs = pg.toComplex(rhs)
            
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
        M = createMassMatrix(mesh, pg.RVector(mesh.cellCount(), 1.0))

        rhs = np.zeros((len(times), dof))
        # rhs kann zeitabhängig sein ..wird hier nicht berücksichtigt
        rhs[:] = assembleForceVector(mesh, f) # this is slow: optimize
        
        if debug: print("rhs", swatch.duration())
        U = np.zeros((len(times), dof))
        #init state
        u = pg.RVector(dof, 0.0)
        
        if 'u0' in kwargs:
            U[0, :] = parseArgToArray(kwargs['u0'], dof, mesh, userData)
            
        theta = 1.0
        if 'theta' in kwargs:
             theta = float(kwargs['theta'])
             
        if debug: print("u0", swatch.duration())
             
        measure = 0.
        for n in range(1, len(times)):
            swatch.reset()
            
            dt = times[n] - times[n - 1]
            
            #u[n] = u[n-1] + dt * theta * L(u[n]) + dt * (1-theta) * L(u[n-1])
            
            
            
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
            b = (M - (dt*(1.0 - theta)) * A) * U[n - 1] + \
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
