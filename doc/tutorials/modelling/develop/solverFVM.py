#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

"""
import pygimli as pg
import pygimli.solver as solver
from pygimli.viewer import show
from pygimli.meshtools import createMesh

import matplotlib.pyplot as plt
import numpy as np

class WorkSpace:
    pass

def boundaryToCellDistances(mesh):
    return np.array(list(map(lambda b__: boundaryToCellDistancesBound(b__), mesh.boundaries())))

def cellDataToBoundaryData(mesh, v):
    if len(v) != mesh.cellCount():
        raise Exception("len(v) != mesh.cellCount():", len(v), mesh.cellCount())
    return np.array(list(map(lambda b__: cellToFace(b__, v), mesh.boundaries())))

def boundaryNormals(mesh):
    gB = np.zeros((mesh.boundaryCount(), 3))
    for i, b in enumerate(mesh.boundaries()):
        gB[i] = b.norm()
    return gB

def boundaryToCellDistancesBound(b):
    leftCell = b.leftCell() 
    rightCell = b.rightCell()
    df1 = 0.
    df2 = 0.
    if leftCell: 
        df1 = b.center().distance(leftCell.center())
    if rightCell: 
        df2 = b.center().distance(rightCell.center())
    return df1+df2

def cellDataToBoundaryGrad(mesh, v, vGrad):
    """
    """
    if len(v) != mesh.cellCount() or len(vGrad) != mesh.cellCount():
        raise
    gB = mesh.cellDataToBoundaryGradient(v,vGrad)
    return np.vstack([pg.x(gB), pg.y(gB), pg.z(gB)]).T

    gB = np.zeros((mesh.boundaryCount(), 3))
    
    for b in mesh.boundaries():
        leftCell = b.leftCell() 
        rightCell = b.rightCell()
        gr = pg.RVector3(0.0, 0.0, 0.0)
        t = (b.node(1).pos() - b.node(0).pos()).norm()
        
        if leftCell and rightCell:
            df1 = b.center().distance(leftCell.center())
            df2 = b.center().distance(rightCell.center())
        
            gr = b.norm()*(v[rightCell.id()] - v[leftCell.id()]) / (df1+df2)
                        
            grL = t * t.dot(vGrad[leftCell.id()])
            grR = t * t.dot(vGrad[rightCell.id()])
            
            gr += (grL+grR) * 0.5
            
        elif leftCell:
            gr = t * t.dot(vGrad[leftCell.id()])
                        
        gB[b.id(), 0] = gr[0]
        gB[b.id(), 1] = gr[1]
        gB[b.id(), 2] = gr[2]
    return gB

def divergence(mesh, V, span=None):
    div = mesh.divergence(V)
    return div
      
    
    swatch = pg.Stopwatch(True)
    ret = np.zeros((mesh.cellCount(), 3))
    
    for b in mesh.boundaries():
        leftCell = b.leftCell() 
        rightCell = b.rightCell()
        
        #flow = b.norm() * b.size()
        vec = mesh.boundarySizedNormals()[b.id()] * V[b.id()]
        
        #flow = mesh.boundaryFlow()[b.id()]
        if b.leftCell():
            ret[leftCell.id(),:] += vec.array()
            #ret[leftCell.id(),:] += vec
            #ret[leftCell.id(), 0] += vec[0]
            #ret[leftCell.id(), 1] += vec[1]
            #ret[leftCell.id(), 2] += vec[2]
            #ret[leftCell.id()] += vec.array()
            #ret[leftCell.id()] += vec.array()
            
        if b.rightCell():
            ret[rightCell.id(), :] -= vec.array()
            #ret[rightCell.id(), :] -= vec
            #ret[rightCell.id(), 0] -= vec[0]
            #ret[rightCell.id(), 1] -= vec[1]
            #ret[rightCell.id(), 2] -= vec[2]
            
    print('     C', swatch.duration(True))
    
    ret[:,0] /= mesh.cellSizes()
    ret[:,1] /= mesh.cellSizes()
    ret[:,2] /= mesh.cellSizes()
    print('     D', swatch.duration(True))
    
    if type(V[0]) == float:
        return ret
    
    return ret[:,0] + ret[:,1] +ret[:,2]
    

def cellToFace(boundary, vec, harmonic=False):
    """
    DEPRECATED
        Interpolate cell values to face value by weighted arithmetic or harmonic mean.
    """
    leftCell = boundary.leftCell() 
    rightCell = boundary.rightCell()
    df1 = 0.
    df2 = 0.
    u1 = 0.
    u2 = 0.
    if leftCell: 
        df1 = boundary.center().distance(leftCell.center())
        u1 = vec[leftCell.id()]
    if rightCell: 
        df2 = boundary.center().distance(rightCell.center())
        u2 = vec[rightCell.id()]
    uFace  = 0.0
    d12 = (df1 + df2)
    
    if leftCell and rightCell: 
        if harmonic:
            # harmonic mean
            uFace = (u1 * u2) / ((u2-u1)*df2/d12 + u1)
        else:
            # arithmetic mean
            # check left vs. right 
            uFace = (u1 - u2) * df2/d12 + u2
            
    elif leftCell:
        uFace = u1# / df1
    elif rightCell: 
        uFace = u2# / df2
    return uFace

def cellToFaceArithmetic(boundary, AMM):
    leftCell = boundary.leftCell() 
    rightCell = boundary.rightCell()
    df1 = 0.
    df2 = 0.
    harmonic = False
    if leftCell: 
        df1 = boundary.center().distance(leftCell.center())
    if rightCell: 
        df2 = boundary.center().distance(rightCell.center())
    d12 = (df1 + df2)
    
    if leftCell and rightCell: 
        if harmonic:
            pass
            # harmonic mean
            #uFace = (u1 * u2) / ((u2-u1)*df2/d12 + u1)
        else:
            # arithmetic mean
            # check left vs. right 
            AMM.addVal(boundary.id(), leftCell.id(), df2/d12)
            AMM.addVal(boundary.id(), rightCell.id(), -df2/d12 + 1.0)
            #uFace = (u1 - u2) * df2/d12 + u2
            
    elif leftCell:
        AMM.addVal(boundary.id(), leftCell.id(), 1.0)
    elif rightCell: 
        AMM.addVal(boundary.id(), rightCell.id(), 1.0)
        
def cellDataToBoundaryDataMatrix(mesh):
    AMM = pg.RSparseMapMatrix(mesh.boundaryCount(), mesh.cellCount())
        
    for b in mesh.boundaries():
        cellToFaceArithmetic(b, AMM)
        
    return AMM

def cellDataToCellGrad2(mesh, v):
    if len(v) != mesh.cellCount():
        raise
    
    vN = pg.cellDataToPointData(mesh, v)
    gC = np.zeros((mesh.cellCount(), 3))

    for c in mesh.cells():
        gr = c.grad(c.center(), vN)
        gC[c.id(), 0] = gr[0]
        gC[c.id(), 1] = gr[1]
        gC[c.id(), 2] = gr[2]
    return gC

def cellDataToCellGrad(mesh, v, CtB):
    if len(v) != mesh.cellCount():
        print (len(v), mesh.cellCount())
        raise 
    div = mesh.boundaryDataToCellGradient(CtB * v)
    return np.vstack([pg.x(div), pg.y(div), pg.z(div)]).T
        
       
    vF = cellDataToBoundaryData(mesh, v)
    gC = np.zeros((mesh.cellCount(), 3))

    for b in mesh.boundaries():
        
        leftCell = b.leftCell()
        rightCell = b.rightCell()
        vec = b.norm() * vF[b.id()] * b.size()
            
        if leftCell:
            gC[leftCell.id(), 0] += vec[0]
            gC[leftCell.id(), 1] += vec[1]
            gC[leftCell.id(), 2] += vec[2]
        if rightCell:
            gC[rightCell.id(), 0] -= vec[0]
            gC[rightCell.id(), 1] -= vec[1]
            gC[rightCell.id(), 2] -= vec[2]
        
    gC[:,0] /= mesh.cellSizes()
    gC[:,1] /= mesh.cellSizes()
    gC[:,2] /= mesh.cellSizes()  
        
    return gC

def findVelocity(v, cellCount, b, c, nc=None):
    """
        Find velocity for boundary b
    """
    vel = [0.0, 0.0, 0.0]
    if hasattr(v, '__len__'):
        if len(v) == cellCount:
            if nc:
                vel = (v[c.id()] + v[nc.id()])/2.0
            else:
                vel = v[c.id()]
        else:
            vel = c.vec(b.center(), v)
            
    return vel
                    
def findDiffusion(mesh, a, b, c, nc=None):
    D = 0
    
    if nc:
        if len(a) == mesh.boundaryCount():
            D = a[b.id()] / nc.center().distance(c.center()) * b.size()
        else:
            D = a[c.id()] / nc.center().distance(c.center()) * b.size()
        # Diffusion part
        # Interface harmonic median
        #D = 1./(c.center().distance(b.center())/a[c.id()] + 
        #nc.center().distance(b.center())/a[nc.id()]) * b.size():
    else:        
        if len(a) == mesh.boundaryCount():
            D = a[b.id()] / b.center().distance(c.center()) * b.size()                
        else:
            D = a[c.id()] / b.center().distance(c.center()) * b.size()                
    return D
                    
def diffusionConvectionKernel(mesh, a, f, uDirBounds=[], fn=None, v=0, u0=0,
                              scheme='CDS', sparse=False):
    """
        Peclet Number - ratio between convection/diffusion
        
        Advection .. forced convection
    """
    
    AScheme = None
    if scheme == 'CDS':
        # CDS - central differences scheme .. maybe irregular for Peclet-number |F/D| > 2
        # diffusion dominant
        # Error of order 2
        AScheme = lambda peclet_: 1.0 - 0.5 * abs(peclet_)    
    elif scheme == 'UDS':
        # UDS - upwind scheme  
        # Convection dominant
        # Error of order 1
        AScheme = lambda peclet_: 1.0
    elif scheme == 'HS':
        #HS - hybrid scheme. 
        #Diffusion dominant for Peclet-number |(F/D)| < 2
        #Convection dominant else
        AScheme = lambda peclet_: max(0.0, 1.0 - 0.5 * abs(peclet_))
    elif scheme == 'PS':
        #PS - power-law scheme. 
        #Identical to HS for Peclet-number |(F/D)| > 10 and near to ES else
        AScheme = lambda peclet_: max(0.0, (1.0 - 0.1 * abs(peclet_))**5.0)
    elif scheme == 'ES':
        # ES - exponential scheme  
        # Only stationary one-dimensional but exact solution
        AScheme = lambda peclet_: (peclet_) / (np.exp(abs(peclet_))-1.0) if peclet_ != 0.0 else 1
    else:
        raise
        
    useHalfBoundaries = False
    S = None
    dof = mesh.cellCount()
    if useHalfBoundaries:
        dof = mesh.cellCount() + len(uDirBounds)

    if sparse:
        S = pg.RSparseMapMatrix(dof, dof, 0)
    else:
        S = np.zeros((dof, dof))
            
    rhsBoundaryScales = np.zeros(dof)
    
    swatch = pg.Stopwatch(True)
    
    for c in mesh.cells():
        intPre = 0.
        
        for bi in range(c.boundaryCount()):
            b = pg.findBoundary(c.boundaryNodes(bi))
            
            nc = b.leftCell()
            if nc == c: nc = b.rightCell()
              
            if nc or 1:
                n = b.norm(c)
                vel = findVelocity(v, mesh.cellCount(), b, c, nc)
                F = n.dot(vel) * b.size()

                D = findDiffusion(mesh, a, b, c, nc)
                                
                aB = D * AScheme(F / D) + max(-F, 0.0)
                
                if nc:
                    # no boundary
                    if sparse:
                        S.addVal(c.id(), nc.id(), -aB)
                        S.addVal(c.id(), c.id(),  +aB)
                    else:
                        S[c.id(), nc.id()] -= aB
                        S[c.id(), c.id()] += aB
                
                elif not useHalfBoundaries:
                    if b.id() in uDirBounds.keys():
                        
                        if sparse:
                            S.addVal(c.id(), c.id(), aB)
                        else:
                            S[c.id(), c.id()] += aB    
                    
                        rhsBoundaryScales[c.id()] += uDirBounds[b.id()] * aB
        
        if fn != None:
            if sparse:
                S.addVal(c.id(), c.id(), -fn[c.id()] * c.shape().domainSize())
            else:
                S[c.id(), c.id()] -= fn[c.id()] * c.shape().domainSize()
        
    if useHalfBoundaries:        
        for i, [b, val] in enumerate(uDirBounds):
            bIdx = mesh.cellCount() + i
        
            c = b.leftCell()
            if not c:
                c = b.rightCell()
        
            if c:
                n = b.norm(c)
                vel = findVelocity(v, mesh.cellCount(), b, c, nc=None)
                F = n.dot(vel) * b.size()

                D = findDiffusion(mesh, a, b, c)
                aB = D * AScheme(F / D) + max(-F, 0.0)
                        
                        
                if useHalfBoundaries:
                    if sparse:
                        S.setVal(c.id(), c.id(), 1.)
                        S.addVal(c.id(), bIdx, -aB)
                    else:
                        S[bIdx, bIdx] = 1.
                        S[c.id(), bIdx] -= aB
                
                    rhsBoundaryScales[bIdx] = aB
                        
    return S, rhsBoundaryScales

def solveFiniteVolume(mesh, a=1.0, f=0.0, fn=0.0, v=0.0, u0=None,
                      uBoundary=None,
                      times=None,
                      theta=1.0,
                      uL=None, relax=1.0,
                      ws=None, scheme='CDS'):
    """
    """
    #The Workspace if to hold temporary data or preserve matrix rebuild
    swatch = pg.Stopwatch(True)
    sparse=True
    
    workspace = WorkSpace()
    if ws:
        workspace = ws
    
    if theta != 1.0:
        raise(AttributeError("theta != 1 not yet implemented"))
    
    a = solver.parseArgToArray(a, [mesh.cellCount(), mesh.boundaryCount()])
    f = solver.parseArgToArray(f, mesh.cellCount())
    fn = solver.parseArgToArray(fn, mesh.cellCount())
    
    boundsDirichlet = dict()
    
    if not hasattr(workspace, 'S'):
        if uBoundary is not None:
            for bPair in uBoundary:
                marker = bPair[0]
                val = bPair[1]
                bounds = mesh.findBoundaryByMarker(marker)
                if len(bounds) == 0:
                    raise Exception("No boundaries with marker %d found" % marker)
                for b in bounds: 
                    boundsDirichlet[b.id()] = val
 
        workspace.S, workspace.rhsBCScales = diffusionConvectionKernel(mesh=mesh,
                                                                       a=a,
                                                                       f=f,
                                                                       uDirBounds=boundsDirichlet,
                                                                       u0=u0,
                                                                       fn=fn,
                                                                       v=v,
                                                                       scheme=scheme,
                                                                       sparse=sparse)
        dof = len(workspace.rhsBCScales)
        
        workspace.uDir = np.zeros(dof)

        if u0 is not None:
            workspace.uDir = np.array(u0)
        
        if len(boundsDirichlet):
            for key, val in boundsDirichlet.items():
                workspace.uDir[mesh.boundary(key).leftCell().id()] = val
        
        workspace.ap = np.zeros(dof)
    
        print('FVM: WS:', swatch.duration(True))
        # for nonlinears
        
        if uL is not None:
            for i in range(dof):
                v = 0.0
                if sparse:
                    v = workspace.S.getVal(i,i) / relax
                    workspace.S.setVal(i,i,v)
                    #workspace.S[i, i] /= relax
                    #workspace.ap[i] = workspace.S[i, i]
                else:
                    v = workspace.S[i, i] / relax
                    workspace.S[i, i] = v
                    
                workspace.ap[i] = v
    
        if sparse:
            Sm = pg.RSparseMatrix(workspace.S)
            # hold Sm until we have reference counting, loosing Sm here will kill LinSolver later            
            workspace.Sm = Sm
            workspace.solver = pg.LinSolver(Sm)
        
        #print('FVM: Fact:', swatch.duration(True))
        #solver.showSparseMatrix(pg.RSparseMatrix(AMM)) 
    #workspace.rhs = rhs
    
    # ... if not hasattr(workspace, 'S'):

    workspace.rhs = np.zeros(len(workspace.rhsBCScales))
    workspace.rhs[0:mesh.cellCount()] = f * mesh.cellSizes()
        
    if len(workspace.uDir):
        workspace.rhs += workspace.rhsBCScales
        #workspace.rhs += workspace.uDir * workspace.rhsBCScales

    #print(workspace.S)
    #print(workspace.uDir)
    #print(workspace.rhsBCScales)
    #print(workspace.rhs)
    # for nonlinears
    
    if uL is not None:
        workspace.rhs += (1. - relax) * workspace.ap * uL
    # print('FVM: Prep:', swatch.duration(True))
    
    if not hasattr(times, '__len__'):
        
        u = None
        if sparse:
            u = workspace.solver.solve(workspace.rhs)
        else:
            u = np.linalg.solve(workspace.S, workspace.rhs)
        return u[0:mesh.cellCount():1]
    else:
        u = np.zeros((len(times), len(rhs)))
        u[0, 0:len(u0)] = u0
        dt = np.zeros(len(rhs))
        dt[0:mesh.cellCount()] = mesh.cellSizes() / (times[1] - times[0])
        I = np.diag(np.ones(len(S)))
    
        for n in range(1, len(times)):
            b = rhs + u[n - 1] * dt
            A = (I * dt + S)
            ut = np.linalg.solve(A, b)
            u[n, 0:mesh.cellCount()] = ut[0:mesh.cellCount():1]
        return u[:,0:mesh.cellCount()]
    
def createFVPostProzessMesh(mesh, bounds=None):
    """
        Create a mesh suitable for node based post processing of cell centered Finite Volume solutions.
        This is something like cellDataToPointData with extra Dirichlet points but without smoothing due to interpolation.
    """
    def isBoundary(b):
        return b.rightCell() is None and b.leftCell() is not None
  
    if bounds is None:
        bounds = []
        boundsIdx = []
        for b in mesh.boundaries():
            if isBoundary(b):
                bounds.append(b)
                boundsIdx.append(b.id())
                
    poly2 = pg.Mesh(2)
    for p in mesh.cellCenters(): poly2.createNode(p) 
    
    boundSort = []
    boundSort.append(bounds[0])
    boundSort[-1].tag()
    boundSortIdx = [0]
    lastB = boundSort[-1]
    lastN = lastB.node(1)

    while len(boundSort) != len(bounds):
        newB = None
        newN = None
        for b in lastN.boundSet():
            if not b.tagged() and isBoundary(b):
                boundSort.append(b)
                b.tag()
                newB = b
        
                if newB.node(0) == lastN:
                    newN = newB.node(1)
                else:
                    newN = newB.node(0)
    
        lastN = newN
        lastB = newB
        
        if not newB:
            raise

        boundSortIdx.append(boundsIdx.index(newB.id()))
        
    bNodes = list(map(lambda b_ : poly2.createNode(b_.center()), boundSort))
    
    for i in range(len(bNodes)):
        poly2.createEdge(bNodes[i], bNodes[(i + 1)%len(bNodes)])

    mesh2 = createMesh(poly2, switches='-pezY')
    return mesh2, boundSortIdx
#def createFVPostProzessMesh(...)
    

def applyBoundaryValues(uB, mesh, uBBC):
    for marker, val in uBBC:
        for b in mesh.findBoundaryByMarker(marker):
            uB[b.id()] = val

def __d(name, v, showAll=False):
    print(name, np.mean(v), min(v), max(v))
    if showAll:
        print(v)
        
def solveStokes_NEEDNAME(mesh, velBoundary, preBoundary=[], 
                         viscosity=1, pre0=None, vel0=None,
                         tol=1e-4, maxIter=1000,
                         verbose=1):    
    """
    """
    velocityRelaxation = 0.5
    pressureRelaxation = 0.8
    pressureCoeff = None
    preCNorm = []
    divVNorm = []
    
    velBoundaryX = [[marker, vel[0]] for marker, vel in velBoundary]
    velBoundaryY = [[marker, vel[1]] for marker, vel in velBoundary]
    
    class WS:
        pass

    wsux = WS()
    wsuy = WS()
    wsp = WS()

    pressure = None
    if pre0 is None:
        pressure = np.zeros(mesh.cellCount())
    else:
        pressure = np.array(pre0)

    velocity = None
    if vel0 is None:
        velocity = np.zeros((mesh.cellCount(), mesh.dimension()))
    else:
        velocity = np.array(vel0)

    mesh.createNeighbourInfos()

    CtB = mesh.cellToBoundaryInterpolation()
    
    controlVolumes = CtB * mesh.cellSizes()
    
    for i in range(maxIter):
        pressureGrad = cellDataToCellGrad(mesh, pressure, CtB)    
        #__d('vx', pressureGrad[:,0])
        
        velocity[:,0] = solveFiniteVolume(mesh,
                                          a=viscosity,
                                          f=-pressureGrad[:,0], 
                                          uBoundary=velBoundaryX,
                                          uL=velocity[:,0],
                                          relax=velocityRelaxation, 
                                          ws=wsux)
        #for s in wsux.S:
            #print(s)
        #__d('rhs', wsux.rhs, 1)
        #__d('ux', velocity[:,0])
        
        velocity[:,1] = solveFiniteVolume(mesh, 
                                          a=viscosity,
                                          f=-pressureGrad[:,1], 
                                          uBoundary=velBoundaryY,
                                          uL=velocity[:,1],
                                          relax=velocityRelaxation,
                                          ws=wsuy)
        
        ap = np.array(wsux.ap)
        apF = CtB * ap
        uxF = CtB * velocity[:,0]
        uyF = CtB * velocity[:,1]

        applyBoundaryValues(uxF, mesh, velBoundaryX)
        applyBoundaryValues(uyF, mesh, velBoundaryY)
                
        pxF = CtB * pressureGrad[:,0]
        pyF = CtB * pressureGrad[:,1]
     
        pF2 = cellDataToBoundaryGrad(mesh, pressure, pressureGrad)

        velXF = uxF + controlVolumes / apF * (pxF - pF2[:, 0])
        velYF = uyF + controlVolumes / apF * (pyF - pF2[:, 1])
    
        applyBoundaryValues(velXF, mesh, velBoundaryX)
        applyBoundaryValues(velYF, mesh, velBoundaryY)
    
        if pressureCoeff is None:
            pressureCoeff = 1./ apF * mesh.boundarySizes() * boundaryToCellDistances(mesh)
                    
        div = -divergence(mesh, np.vstack([velXF, velYF]).T)
        
        pressureCorrection = solveFiniteVolume(mesh,
                                               a=pressureCoeff,
                                               f=div, 
                                               uBoundary=preBoundary,
                                               ws=wsp)
        
        pressure += pressureCorrection * pressureRelaxation
       
        pressureCorrectionGrad = cellDataToCellGrad(mesh, pressureCorrection, CtB)
    
        velocity[:, 0] -= pressureCorrectionGrad[:, 0] / ap * mesh.cellSizes()
        velocity[:, 1] -= pressureCorrectionGrad[:, 1] / ap * mesh.cellSizes()
     
        preCNorm.append(pg.norm(pressureCorrection))
        divVNorm.append(pg.norm(div))
     
        #__d('div', div)
        #if ( i == 1):
            #sd
        if verbose:
            #print("\rIter: " + str(i) + " div V=" + str(divVNorm[-1]) + " " +  str(preCNorm[-1]), end='')
            print("\r" + str(i) + " div V=" + str(divVNorm[-1]))
     
        if i == maxIter or divVNorm[-1] < tol:
            break
            
    if verbose:
        print(str(i)+ ": " + str(preCNorm[-1]))
    return velocity, pressure, preCNorm, divVNorm
    
if __name__ == '__main__':

    N=21 # 21 reference
    maxIter=11 # 11 reference
    Nx=N
    Ny=N

    x = np.linspace(-1.0, 1.0, Nx+1)
    y = np.linspace(-1.0, 1.0, Ny+1)
    grid = pg.createGrid(x=x, y=y)

    a = pg.RVector(grid.cellCount(), 1.0)

    b7 = grid.findBoundaryByMarker(1)[0]
    for b in grid.findBoundaryByMarker(1):
        if b.center()[1] < b.center()[1]:
            b7 = b
    b7.setMarker(7)

    swatch = pg.Stopwatch(True)
    velBoundary=[[1, [0.0, 0.0]],
                 [2, [0.0, 0.0]],
                 [3, [1.0, 0.0]],
                 [4, [0.0, 0.0]],
                 [7, [0.0, 0.0]]]

    preBoundary=[[7, 0.0]]

    vel, pres, pCNorm, divVNorm = solveStokes_NEEDNAME(grid, velBoundary, preBoundary,
                                             viscosity=a,
                                             maxIter=maxIter,
                                             verbose=1)
   
    print("time", len(pCNorm), swatch.duration(True) )
    referencesolution=1.2889506342694153
    referencesolutionDivV=0.029187181920161752
    print("divNorm: ", divVNorm[-1])
    print("to reference: ", divVNorm[-1]-referencesolutionDivV)
    
    fig = plt.figure()
    ax1 = fig.add_subplot(1, 3, 1)
    ax2 = fig.add_subplot(1, 3, 2)
    ax3 = fig.add_subplot(1, 3, 3)

    show(grid, data=pg.cellDataToPointData(grid, pres),
        logScale=False, showLater=True, colorBar=True, axes=ax1, cbar='b2r')
    show(grid, 
        data=pg.logTransDropTol(pg.cellDataToPointData(grid, vel[:,0]), 1e-2),
        logScale=False, showLater=True, colorBar=True, axes=ax2)
    show(grid,
        data=pg.logTransDropTol(pg.cellDataToPointData(grid, vel[:,1]), 1e-2),
        logScale=False, showLater=True, colorBar=True, axes=ax3)
     
    show(grid, data=vel, axes=ax1)
    show(grid, showLater=True, axes=ax1)
        
    plt.figure()
    plt.semilogy(pCNorm, label='norm')
    plt.legend()

    plt.ioff()
    plt.show()
