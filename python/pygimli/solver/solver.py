# -*- coding: utf-8 -*-

import pygimli as g

from pygimli.utils import unique
import numpy as np

def assembleForceVector(mesh, vals, scale = 1.):
    rhs = g.RVector(mesh.nodeCount(), 0)
    
    b_l = g.DElementMatrix()
    
    for c in mesh.cells():
        b_l.u(c)
            
        for i, idx in enumerate(b_l.idx()):
            rhs[ idx ] += b_l.row(0)[ i ] * vals * scale
            
    return rhs
# def assembleForceVector()
    
def assembleNeumannBC(S, boundaries, vals):

    Se = g.DElementMatrix()
    
    for b in boundaries:
        Se.u2(b)
        Se *= vals
        S += Se;
                
#def assembleNeumannBC(...)


def assembleDirichletBC(S, boundaryPair, rhs, time=0.0):
    """
    Create Dirichlet boundary condition and apply them to system matrix S.

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
        print boundaryPair[0]
        for b in boundaryPair[0]:
            uVal = None
            
            if type(boundaryPair[1]) == float or type(boundaryPair[1]) == int:
                uVal = boundaryPair[1]
            else:
                if time != 0.0:
                    uVal = boundaryPair[1](b, time)
                else:
                    uVal = boundaryPair[1](b)
                                    
            if uVal is not None:
                for n in b.nodes():
                    uDirNodes.append(n)
                        #print b.marker(), n.id(), uVal
                    uDirVal[n.id()] = uVal
    
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
    
    udirTmp = g.RVector(S.rows(), 0.0)
    
    udirTmp.setVal(uDirchlet, uDirIndex)
    
    rhs -= S * udirTmp
        
    for i in uDirIndex:
    
        S.cleanRow(i)
        S.cleanCol(i)
        S.setVal(i, i, 1.0)
    
    rhs.setVal(uDirchlet, uDirIndex)
    
    
def solvePoisson(mesh, a=1.0, timeSteps=None, f=1, verbose=False,
                 *args, **kwargs):
    """
        TODO
    Parameters
    ----------
    
    """
    if verbose:
        print(mesh)

    dof = mesh.nodeCount()
        
    swatch = g.Stopwatch(True)
    
    # define an empty stiffness matrix
    A = g.DSparseMatrix()
    B = g.DSparseMatrix()
    
    # create matrix structure regarding the mesh
    A.buildSparsityPattern(mesh)
    B.buildSparsityPattern(mesh)

    # define a local element matrix 
    A_l = g.DElementMatrix()
    B_l = g.DElementMatrix()

    # check for material parameter
    if type(a) == float:
        a = g.RVector(mesh.cellCount(), a)
    elif len(a) != mesh.cellCount():
        raise Exception("Material array 'a' has the wrong size: " + len(a) + " != " +  mesh.cellCount())
    
    # assemble the stiffness matrix
        
    for c in mesh.cells():
        A_l.ux2uy2uz2(c)
        A_l *= a[ c.id() ] 
        A += A_l
    
        B_l.u2(c)
        B += B_l
    
    
    #print rhs
        #for i in range(S.size()):
            #for j in range(S.vecColPtr()[ i ], S.vecColPtr()[ i + 1 ]):
                    #print i, S.vecVals()[ j ]
                    ##print i, S.rowIdx()[ j ], S.vecVals()[ j ]

    if timeSteps == None:
    
        rhs = assembleForceVector(mesh, f)
        
        if 'neumann' in kwargs:
            gb = kwargs[ 'neumann' ]
            assembleNeumannBC(A, boundaries = gb[ 0 ], vals = gb[ 1 ])
        
        if 'uBoundary' in kwargs:
            
            uBound = kwargs['uBoundary']
            
            bounds = []
            # [bounds, u]
            print uBound, type(uBound)
            if type(uBound[0]) == g.stdVectorBounds:
                
                bounds.append(uBound)
            else:
                bounds = uBound
                
            print bounds
            # [[bounds1, u1], [bounds2, u2]]
            for bound in bounds:

                if bound == tuple or len(bound) == 2:
                    if type(bound[0]) == int:
                        bound[0] = mesh.findBoundaryByMarker(bound[0])

                    assembleDirichletBC(A, boundaryPair=bound, rhs=rhs)    
                else:
                    raise Exception("cannot interpret boundaries sequence: " +
                                    str(type(bound)))
            
            
        u = g.RVector(rhs.size(), 0.0)
        
        if verbose:
            print("asssemblation takes: ", swatch.duration(True))

        solver = g.LinSolver(A, verbose)
        solver.solve(rhs, u)
        
        if verbose:
            print("lin solving takes: ", swatch.duration(True))
            
        return u
        
    else:
        U = g.RMatrix(0, dof)
        #init state
        u = g.RVector(dof, 0.0)
        U.push_back(u)
    
        for i in range(1, len(timeSteps)):
            dt = timeSteps[ i ] - timeSteps[ i - 1 ]
            b = assembleForceVector(mesh, f, dt)

  #% Neumann conditions
  #for j = 1 : size(neumann,1)
    #b(neumann(j,:)) = b(neumann(j,:)) + ...
    #norm(coordinates(neumann(j,1),:)-coordinates(neumann(j,2),:))* ...
    #dt*g(sum(coordinates(neumann(j,:),:))/2,n*dt)/2; 
  #end
  
            # previous timestep
  
            #print "i: ", i, U[i-1]
            b = b + B * U[ i-1 ]
  
            
            S = A * dt + B
            
            if 'uBoundary' in kwargs:
                assembleDirichletBC(S, kwargs[ 'uBoundary' ], b, timeSteps[ i-1 ])    
                             
            solver = g.LinSolver(S, verbose)
            #print "b", b
            solver.solve(b, u)
            
        
            U.push_back(u)
    
        return U
# def solvePoisson(..):