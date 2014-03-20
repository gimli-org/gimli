# -*- coding: utf-8 -*-

import pygimli as pg

from pygimli.utils import unique
import numpy as np

import copy

def parseArgToArray(arg, ndof, mesh, userData=None):
    """
        What is this
    """
  
    try:
        return pg.RVector(ndof, float(arg))
    except:
        pass
    
    if hasattr(arg, '__len__'):
        if len(arg) == ndof:
            return arg;
        raise Exception("Array 'arg' has the wrong size: " + 
                        len(arg) + " != " +  dof)
    elif hasattr(arg, '__call__'):
        ret = pg.RVector(ndof, 0.0)
        
        if ndof == mesh.nodeCount():
            for n in mesh.nodes():
                if userData:
                    ret[n.id()] = arg(n.pos(), userData)
                else:
                    ret[n.id()] = arg(n.pos())
        elif ndof == mesh.cellCount():
            for c in mesh.cells():
                if userData:
                    ret[c.id()] = arg(c, userData)
                else:
                    ret[c.id()] = arg(c)
        else:
            raise Exception("Cannot parse callable argument " + str(ndof) + 
                            " nodes: " + str(mesh.nodeCount()) + 
                            " cells: " + str(mesh.cellCount()))
            
        return ret
    raise Exception("Cannot parse argument type " + str(type(arg)))
#def parseArgToArray(...)

    
def assembleForceVector(mesh, vals, userData=None):
    
    """Check size of vals for the right hand side .. loading vector"""
    try:
        """ Check size of given f equals the requested rhs vector. 
            So we use it directly.
        """
        if len(vals) == mesh.nodeCount():
            rhs = pg.RVector(vals)
            return rhs
        else:
            print(Exception("f does not fit directly so we generate the RHS"))
    except:
        pass
            
    rhs = pg.RVector(mesh.nodeCount(), 0)
    
    if vals == 0.:
        return rhs
    
    b_l = pg.DElementMatrix()
    
    for c in mesh.cells():
        if hasattr(vals, '__call__'):
            if userData is not None:
                vals(c, rhs, userData)
            else:
                vals(c, rhs)
        else:
            b_l.u(c)
            
            for i, idx in enumerate(b_l.idx()):
                rhs[idx] += b_l.row(0)[i] * vals

    return rhs
# def assembleForceVector()
    
    
def assembleUDirichlet_(S, rhs, uDirIndex, uDirchlet):
    """
        this should be moved directly into gimli
    """
    udirTmp = pg.RVector(S.rows(), 0.0)
    udirTmp.setVal(uDirchlet, uDirIndex)
    
    rhs -= S * udirTmp
        
    for i in uDirIndex:
    
        S.cleanRow(i)
        S.cleanCol(i)
        S.setVal(i, i, 1.0)
    
    rhs.setVal(uDirchlet, uDirIndex)
#def assembleUDirichlet_(...)
       
def assembleNeumannBC(S,
                      boundaryPair, rhs, time=0.0, 
                      userData=None, verbose=False):
    """
    Apply Neumann condition and apply them to system matrix S.
         
    .. math::
        \\frac{\\partial\\textbf{u}}{\\partial\\textbf{n}} = vals
        
    Parameters
    ----------

    S   : pg.DSparseMatrix()
        System matrix of the system equation.
    boundaryPair    : tuple
        Pair of [list_of_boundaris, value], the value will assigned to 
        the nodes of the boundaries. 
        Value can be a scalar (float or int) or a function, which will be called
        with the boundary and a optional time, that return the value for u.
    rhs :unused
        *For compatibility only*
    
    """

    Se = pg.DElementMatrix()

    if type(boundaryPair) == tuple or len(boundaryPair) == 2:
        
        #if verbose:
            #print("Setting " + str(len(boundaryPair[0])) + " bounds to du/dn = " + str(boundaryPair[1]));
    
        for b in boundaryPair[0]:
            val = None
            
            if type(boundaryPair[1]) == float or type(boundaryPair[1]) == int:
                val = boundaryPair[1]
            else:
                kwargs = dict()
                args = [b]
                if time != 0.0:
                    args.append(time)
                if userData:
                    kwargs['userData'] = userData
                    
                val = boundaryPair[1](*args, **kwargs)
            
            if val is not 0.0:
                Se.u2(b)
                Se *= val
                S += Se;
        
#def assembleNeumannBC(...)


def assembleDirichletBC(S, boundaryPair, rhs=None, time=0.0,
                        userData=None, verbose=False):
    """
    Apply Dirichlet boundary condition and apply them to system matrix S.

    .. math::
        u = vals` on boundaries
    
    Parameters
    ----------
    S : pg.DSparseMatrix()
        System matrix of the system equation.
    boundaryPair : tuple
        Pair of [list_of_boundaris, value], the value will assigned to 
        the nodes of the boundaries. 
        Value can be a scalar (float or int) or a function, which will be called
        with the boundary and a optional time, that return the value for u.
    rhs : pg.RVector
        Right hand side vector of the system equation
    time : optional
        time, if some time depended given boundary function

    """
           
    uDirNodes = []
    uDirVal = dict()
    
    if type(boundaryPair) == tuple or len(boundaryPair) == 2:
        
        for b in boundaryPair[0]:
            for n in b.nodes():
                
                uVal = None
                
                if type(boundaryPair[1]) == float or \
                   type(boundaryPair[1]) == int :
                    uVal = boundaryPair[1]
                else:
                    kwargs = dict()
                    args = [n.pos()]
                    if time != 0.0:
                        args.append(time)
                    if userData:
                        kwargs['userData'] = userData
                    
                    uVal = boundaryPair[1](*args, **kwargs)

                if uVal is not None:
                
                    uDirNodes.append(n)
                    #print b.marker(), n.id(), uVal
                    uDirVal[n.id()] = uVal
                else:
                    raise Exception("cannot find dirichlet value for node ", n)
    
    else:
        raise Exception("cannot interpret boundaries sequence:" + str(type(uDirchlet)))
    
    if len(uDirNodes) == 0:
        return 
        
    uniqueNodes = unique(uDirNodes) 
        
    uDirchlet = pg.RVector(len(uniqueNodes))
    uDirIndex = []
    
    for i, n in enumerate(uniqueNodes):
        uDirIndex.append(n.id())
        uDirchlet[i] = uDirVal[n.id()]
    
    #if verbose:
        #print("Setting " + str(len(uDirIndex)) + " nodes to u = " + str(uDirchlet[0]));
    
    assembleUDirichlet_(S, rhs, uDirIndex, uDirchlet)
        
    
def assembleBoundaryConditions(mesh, S, rhs, boundArgs, assembler,
                               time=0.0, 
                               userData=None,
                               verbose=False):
    """
    boundArgs can be:
    [boundaries, value]
    [[boundaries, value],[boundaries, value]]
    
    boundaries can be list of bounds or marker
    value can be float, int or a callable(boundary)
    
    """
    boundaries = list()

    if type(boundArgs[0]) == pg.stdVectorBounds:
        boundaries.append(boundArgs)
    elif type(boundArgs[0]) == int:
        #FIXME
        boundaries.append(copy.deepcopy(boundArgs))
    else:
        #FIXME
        for i in boundArgs:
            boundaries.append(copy.copy(i))
        

    for bound in boundaries:
        
        if type(bound) == list or len(bound) == 2:
            if type(bound[0]) == int:
                bound[0] = mesh.findBoundaryByMarker(bound[0])

            assembler(S, boundaryPair=bound, rhs=rhs, time=time,
                      userData=userData,
                      verbose=verbose)    
            
        else:
            raise Exception("cannot interpret boundaries sequence: " +
                            str(type(bound)))
    
#def assembleBoundaryConditions(...)

def createStiffnessMatrix(mesh, a):
    A = pg.DSparseMatrix()
    A.fillStiffnessMatrix(mesh, a)
    return A
    
    # create matrix structure regarding the mesh
    A.buildSparsityPattern(mesh)

    # define a local element matrix 
    A_l = pg.DElementMatrix()
    for c in mesh.cells():
        A_l.ux2uy2uz2(c)
        A_l *= a[c.id()] 
        A += A_l
    
    return A

def createMassMatrix(mesh, b):
    B = pg.DSparseMatrix()
    B.fillMassMatrix(mesh, b)
    return B

    # create matrix structure regarding the mesh
    B.buildSparsityPattern(mesh)
    
    # define a local element matrix 
    B_l = pg.DElementMatrix()
    for c in mesh.cells():
        B_l.u2(c)
        # check if b[i] == B*b
        B_l *= b[c.id()]
        B += B_l
    
    return B

    

def solvePoisson(mesh, a=1.0, b=0.0, f=0.0, times=None, userData=None,
                 verbose=False, stats=None, *args, **kwargs):
    """
    TODO
    
    The value of :math:`\omega` is larger than 5.

    Variable names are displayed in typewriter font, obtained by using \\mathtt{var}:

    We square the input parameter `a` to obtain
    :math:`\mathtt{alpha}^2`.
    

    Parameters
    ----------
    
    a   : value | array | callable(cell, userData)
        Cell values
    b   : value | array | callable(cell, userData)
        Cell values
    u0      : value | array | callable(pos, userData)
        Node values
    f       : value | array(cells) | array(nodes) | callable(args, kwargs)
        force values 
    theta   : float
        - `theta` = 1, implicit Euler
        - `theta` = 0.5, Crank-Nicolsen, maybe instable 
        - `theta` = 0, explicit Euler, maybe stable for

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
    if verbose:
        print(("Mesh: ", str(mesh)))

    dof = mesh.nodeCount()
        
    swatch = pg.Stopwatch(True)
    swatch2 = pg.Stopwatch(True)
    
    # check for material parameter
    a = parseArgToArray(a, mesh.cellCount(), mesh, userData)
    b = parseArgToArray(b, mesh.cellCount(), mesh, userData)
        
    print("2: ", swatch2.duration(True))
    # assemble the stiffness matrix
    A = createStiffnessMatrix(mesh, a)
    
    print("3: ", swatch2.duration(True))
    M = createMassMatrix(mesh, b)
        
    print("4: ", swatch2.duration(True))
    S = A + M
    
    print("5: ", swatch2.duration(True))
    if times == None:
        
        rhs = assembleForceVector(mesh, f, userData=userData)
        
        print("6a: ", swatch2.duration(True))
        if 'duBoundary' in kwargs:
            assembleBoundaryConditions(mesh, S, rhs, kwargs['duBoundary'], 
                                       assembleNeumannBC,
                                       userData=userData, 
                                       verbose=verbose)

        print("6b: ", swatch2.duration(True))
        if 'uBoundary' in kwargs:
            assembleBoundaryConditions(mesh, S, rhs, kwargs['uBoundary'],
                                       assembleDirichletBC, 
                                       userData=userData,
                                       verbose=verbose)

        print("6c: ", swatch2.duration(True))
        if 'uDirichlet' in kwargs:
            assembleUDirichlet_(S, rhs,
                                kwargs['uDirichlet'][0],
                                kwargs['uDirichlet'][1])
            
            
            
        u = pg.RVector(rhs.size(), 0.0)
        print("7: ", swatch2.duration(True))
        
        assembleTime = swatch.duration(True)
        if stats:
            stats.assembleTime = assembleTime
                
        if verbose:
            print(("Asssemblation time: ", assembleTime))

        
        solver = pg.LinSolver(S, verbose=False)
        solver.solve(rhs, u)

        solverTime = swatch.duration(True)
        if verbose:
            if stats:
                stats.solverTime = solverTime
            print(("Solving time: ", solverTime))
            
        return u
        
    else:
        print("start TL", swatch.duration())
        M = createMassMatrix(mesh, pg.RVector(mesh.cellCount(), 1.0))

        rhs = np.zeros((len(times), dof))
        # rhs kann zeitabhängig sein ..wird hier nicht berücksichtigt
        rhs[:] = assembleForceVector(mesh, f) # this is slow: optimize
        
        print("rhs", swatch.duration())
        U = np.zeros((len(times), dof))
        #init state
        u = pg.RVector(dof, 0.0)
        
        if 'u0' in kwargs:
            U[0, :] = parseArgToArray(kwargs['u0'], dof, mesh, userData)
            
        theta = 1.0
        if 'theta' in kwargs:
             theta = float(kwargs['theta'])
             
        print("u0", swatch.duration())
             
        measure = 0.
        for n in range(1, len(times)):
            swatch.reset()
            
            dt = times[n] - times[n - 1]
            
            # previous timestep
            #print "i: ", i, dt, U[i - 1]
            
            if 'duBoundary' in kwargs:
                # aufschreiben und checken ob neumann auf A oder auf S mit skaliertem val*dt angewendet wird
                A = createStiffnessMatrix(mesh, a)
                assembleBoundaryConditions(mesh, A, None, kwargs['duBoundary'],
                                           assembleNeumannBC,
                                           time=times[n],
                                           verbose=verbose)

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
                assembleBoundaryConditions(mesh, S, b, kwargs['uBoundary'], 
                                           assembleDirichletBC,
                                           time=times[n],
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
        print( "Measure:", measure/len(times))
        return U
# def solvePoisson(..):