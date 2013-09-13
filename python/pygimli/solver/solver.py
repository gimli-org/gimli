# -*- coding: utf-8 -*-

import pygimli as g

from pygimli.utils import unique
import numpy as np

def assembleForceVector(mesh, vals, scale=1., userData=None):
    rhs = g.RVector(mesh.nodeCount(), 0)
    
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
                rhs[idx] += b_l.row(0)[i] * vals * scale
            
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
        
    
def assembleBoundaryConditions(mesh, S, rhs, boundArgs, assembler, time=0.0, 
                               userData=None,
                               verbose=False):
    """
    boundArgs can be:
    [boundaries, value]
    [[boundaries, value],[boundaries, value]]
    
    boundaries can be list of bounds or marker
    value can be float, int or a callable(boundary)
    
    """
    boundaries = []

    if type(boundArgs[0]) == g.stdVectorBounds or \
       type(boundArgs[0]) == int:
        boundaries.append(boundArgs)
    else:
        boundaries = boundArgs
                
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

def solvePoisson(mesh, a=1.0, b=0.0, f=1.0, times=None, userData=None,
                 verbose=False, *args, **kwargs):
    """
        TODO
    Parameters
    ----------
    
    """
    if verbose:
        print("Mesh: ", str(mesh))

    dof = mesh.nodeCount()
        
    swatch = g.Stopwatch(True)
    
    # check for material parameter
    if type(a) == float:
        a = g.RVector(mesh.cellCount(), a)
    if type(a) == int:
        a = g.RVector(mesh.cellCount(), float(a))
    elif len(a) != mesh.cellCount():
        raise Exception("Material array 'a' has the wrong size: " + 
                        len(a) + " != " +  mesh.cellCount())
    
    # check for material parameter
    if type(b) == float:
        b = g.RVector(mesh.cellCount(), b)
    if type(b) == int:
        b = g.RVector(mesh.cellCount(), float(b))
    elif len(b) != mesh.cellCount():
        raise Exception("Diffusion parameter array 'b' has the wrong size: " + 
                        len(b) + " != " +  mesh.cellCount())
        
    # assemble the stiffness matrix
    A = createStiffnessMatrix(mesh, a)
    B = createMassMatrix(mesh, b)
        
    S = A + B
    
    if times == None:
    
        rhs = assembleForceVector(mesh, f, userData=userData)
        
        if 'duBoundary' in kwargs:
            assembleBoundaryConditions(mesh, S, rhs, kwargs['duBoundary'], 
                                       assembleNeumannBC, userData=userData, 
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
        
        if verbose:
            print("Asssemblation time: ", swatch.duration(True))

        solver = g.LinSolver(S, verbose)
        solver.solve(rhs, u)
        
        if verbose:
            print("Solving time: ", swatch.duration(True))
            
        return u
        
    else:
        ################## better solve this with recursion
        U = np.zeros((len(times), dof))
        #init state
        u = g.RVector(dof, 0.0)
        
        if 'u0' in kwargs:
            U[0,:] = kwargs['u0']
        
        
        for i in range(1, len(times)):
            dt = times[i] - times[i - 1]
            rhs = assembleForceVector(mesh, f, dt)

            # previous timestep
            #print "i: ", i, dt, U[i-1]
            b = rhs + B * U[i - 1]
            
            #S = A * dt + B
            if 'duBoundary' in kwargs:
                # aufschreiben und checken ob neumann auf A oder auf S mit skaliertem val*dt angewendet wird
                A = createStiffnessMatrix(mesh, a)
                assembleBoundaryConditions(mesh, A, None, kwargs['duBoundary'],
                                           assembleNeumannBC,
                                           time=times[i],
                                           verbose=verbose)
            S = A * dt + B
                
            if 'uBoundary' in kwargs:
                assembleBoundaryConditions(mesh, S, b, kwargs['uBoundary'], 
                                           assembleDirichletBC,
                                           time=times[i],
                                           verbose=verbose)
                
            #u = S/b
            solver = g.LinSolver(S, verbose)
            solver.solve(b, u)
            
            if 'plotTimeStep' in kwargs:
                kwargs['plotTimeStep'](u, times[i])
                    
            U[i,:] = np.asarray(u)
    
        return U
# def solvePoisson(..):