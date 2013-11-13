# -*- coding: utf-8 -*-

import pygimli as g

from pygimli.utils import unique
import numpy as np

import copy

def parseArgToArray(arg, ndof, mesh, userData=None):
    """
        What is this
    """
  
    try:
        return g.RVector(ndof, float(arg))
    except:
        pass
    
    if hasattr(arg, '__len__'):
        if len(arg) == ndof:
            return arg;
        raise Exception("Array 'arg' has the wrong size: " + 
                        len(arg) + " != " +  dof)
    elif hasattr(arg, '__call__'):
        ret = g.RVector(ndof, 0.0)
        
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
    rhs = g.RVector(mesh.nodeCount(), 0)
    
    if vals == 0.:
        return rhs
    b_l = g.DElementMatrix()
    
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
    udirTmp = g.RVector(S.rows(), 0.0)
    udirTmp.setVal(uDirchlet, uDirIndex)
    
    rhs -= S * udirTmp
        
    for i in uDirIndex:
    
        S.cleanRow(i)
        S.cleanCol(i)
        S.setVal(i, i, 1.0)
    
    rhs.setVal(uDirchlet, uDirIndex)
#def assembleUDirichlet_(...)
       
       
def assembleNeumannBC(S, boundaryPair, rhs, time=0.0, 
                      userData=None, verbose=False):
    """
    Apply Neumann condition and apply them to system matrix S.
    ..math::`frac{\partial u}{\partial \vec{u}} = vals`
        
    Parameters
    ----------
    S : g.DSparseMatrix()
        System matrix of the system equation.
    boundaryPair: tuple
        Pair of [list_of_boundaris, value], the value will assigned to 
        the nodes of the boundaries. 
        Value can be a scalar (float or int) or a function, which will be called
        with the boundary and a optional time, that return the value for u.
    rhs: unused
        *For compatibility only*
    """

    Se = g.DElementMatrix()

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
                #Se *= val * 0.0134228187919
                S += Se;
                
#def assembleNeumannBC(...)


def assembleDirichletBC(S, boundaryPair, rhs, time=0.0,
                        userData=None, verbose=False):
    """
    Apply Dirichlet boundary condition and apply them to system matrix S.

    ..math::`u = vals` on boundaries
    
    Parameters
    ----------
    S : g.DSparseMatrix()
        System matrix of the system equation.
    boundaryPair: tuple
        Pair of [list_of_boundaris, value], the value will assigned to 
        the nodes of the boundaries. 
        Value can be a scalar (float or int) or a function, which will be called
        with the boundary and a optional time, that return the value for u.
    rhs : g.RVector
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
                
                if type(boundaryPair[1]) == float or type(boundaryPair[1]) == int:
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
        
    uDirchlet = g.RVector(len(uniqueNodes))
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

    if type(boundArgs[0]) == g.stdVectorBounds:
        boundaries.append(boundArgs)
    elif type(boundArgs[0]) == int:
        boundaries.append(copy.deepcopy(boundArgs))
    else:
        boundaries = copy.deepcopy(boundArgs)
        

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
    A = g.DSparseMatrix()
    A.fillStiffnessMatrix(mesh, a)
    return A
    
    # create matrix structure regarding the mesh
    A.buildSparsityPattern(mesh)

    # define a local element matrix 
    A_l = g.DElementMatrix()
    for c in mesh.cells():
        A_l.ux2uy2uz2(c)
        A_l *= a[c.id()] 
        A += A_l
    
    return A

def createMassMatrix(mesh, b):
    B = g.DSparseMatrix()
    B.fillMassMatrix(mesh, b)
    return B

    # create matrix structure regarding the mesh
    B.buildSparsityPattern(mesh)
    
    # define a local element matrix 
    B_l = g.DElementMatrix()
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
    Parameters
    ----------
        a       : value|array|callable(cell, userData)
        b       : value|array|callable(cell, userData)
        u0      : value|array|callable(pos, userData)
        f       : value|array|callable(??????)
        theta   : float
            heat equation is stable for 0.5 <= theta <= 1.0
            theta = 1, implicit Euler
            theta = 0.5, Crank-Nicolsen, maybe instable 
            theta = 0, explicit Euler, maybe stable for dT near h 
            if unsure choose 0.5+epsilon, this may be stable 
        progress : bool
    
    
    """
    if verbose:
        print("Mesh: ", str(mesh))

    dof = mesh.nodeCount()
        
    swatch = g.Stopwatch(True)
    swatch2 = g.Stopwatch(True)
    
    # check for material parameter
    a = parseArgToArray(a, mesh.cellCount(), mesh, userData)
    b = parseArgToArray(b, mesh.cellCount(), mesh, userData)
        
    print "2: ", swatch2.duration(True)
    # assemble the stiffness matrix
    A = createStiffnessMatrix(mesh, a)
    print "3: ", swatch2.duration(True)
    M = createMassMatrix(mesh, b)
        
    print "4: ", swatch2.duration(True)
    S = A + M
    
    print "5: ", swatch2.duration(True)
    if times == None:
        rhs = assembleForceVector(mesh, f, userData=userData)
        
        print "6: ", swatch2.duration(True)
        if 'duBoundary' in kwargs:
            assembleBoundaryConditions(mesh, S, rhs, kwargs['duBoundary'], 
                                       assembleNeumannBC,
                                       userData=userData, 
                                       verbose=verbose)

        if 'uBoundary' in kwargs:
            assembleBoundaryConditions(mesh, S, rhs, kwargs['uBoundary'],
                                       assembleDirichletBC, 
                                       userData=userData,
                                       verbose=verbose)

        if 'uDirichlet' in kwargs:
            assembleUDirichlet_(S, rhs,
                                kwargs['uDirichlet'][0],
                                kwargs['uDirichlet'][1])
            
            
            
        u = g.RVector(rhs.size(), 0.0)
        print "7: ", swatch2.duration(True)
        
        assembleTime = swatch.duration(True)
        if stats:
            stats.assembleTime = assembleTime
                
        if verbose:
            print("Asssemblation time: ", assembleTime)

        solver = g.LinSolver(S, verbose)
        solver.solve(rhs, u)

        solverTime = swatch.duration(True)
        if verbose:
            if stats:
                stats.solverTime = solverTime
            print("Solving time: ", solverTime)
            
        return u
        
    else:
        M = createMassMatrix(mesh, g.RVector(mesh.cellCount(), 1.0))

        rhs = np.zeros((len(times), dof))
        rhs[:] = assembleForceVector(mesh, f) # this is slow: optimize
        
        U = np.zeros((len(times), dof))
        #init state
        u = g.RVector(dof, 0.0)
        
        if 'u0' in kwargs:
            U[0, :] = parseArgToArray(kwargs['u0'], dof, mesh, userData)
            
        theta = 1.0
        if 'theta' in kwargs:
             theta = float(kwargs['theta'])
             
        for i in range(1, len(times)):
            swatch.reset()
            
            dt = times[i] - times[i - 1]
            
            # previous timestep
            #print "i: ", i, dt, U[i - 1]
            
            if 'duBoundary' in kwargs:
                # aufschreiben und checken ob neumann auf A oder auf S mit skaliertem val*dt angewendet wird
                A = createStiffnessMatrix(mesh, a)
                assembleBoundaryConditions(mesh, A, None, kwargs['duBoundary'],
                                           assembleNeumannBC,
                                           time=times[i],
                                           verbose=verbose)

            
            b = (M + A * (dt*(theta - 1.0))) * U[i - 1] + \
                rhs[i - 1] * (dt*(1.0 - theta)) + \
                rhs[i] * dt * theta
            
            S = M + A * dt * theta
                
            if 'uBoundary' in kwargs:
                assembleBoundaryConditions(mesh, S, b, kwargs['uBoundary'], 
                                           assembleDirichletBC,
                                           time=times[i],
                                           verbose=verbose)
                
            #u = S/b
            t_prep = swatch.duration(True)
            solver = g.LinSolver(S, verbose)
            solver.solve(b, u)
            
            if 'plotTimeStep' in kwargs:
                kwargs['plotTimeStep'](u, times[i])
                    
            U[i,:] = np.asarray(u)

            if 'progress' in kwargs:
                if kwargs['progress']:
                    print("\t" + str(i) +"/" + str(len(times)-1) + 
                          ": " + str(t_prep) +"/" + str(swatch.duration()))
                    
        return U
# def solvePoisson(..):