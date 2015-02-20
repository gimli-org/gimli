#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

"""

import pygimli as pg
import pygimli.solver as solver
from pygimli.viewer import showMesh
from pygimli.mplviewer import drawMesh, drawModel, drawField
from pygimli.meshtools import createMesh
from solverFVM import solveFiniteVolume, createFVPostProzessMesh, diffusionConvectionKernel, WorkSpace

import matplotlib.pyplot as plt
import numpy as np

def diffusionConvectionKernelGrid(mesh, a, f, v, uDir, scheme='CDS'):
    """
        Peclet Number - ratio between convection/diffusion
        
        Advection .. forced convection
    """
    bounds = mesh.findBoundaryByMarker(1,99)
    S = np.zeros((mesh.cellCount() + len(bounds),
                  mesh.cellCount() + len(bounds)))
    rhs = np.zeros(len(S))
    
    dx = mesh.cell(0).shape().domainSize()
    N = mesh.cellCount()
    
    print(scheme)
    
    AScheme = None
    if scheme == 'CDS':
        # CDS - central differences scheme .. maybe irregular for Peclet-number |F/D| > 2
        # Diffusion dominant
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
        
    for i in range(1, N+1):
    
        De = a[i - 1] / dx
        Dw = a[i - 1] / dx
    
        Fe = v[i - 1]
        Fw = v[i - 1]

        aE = De * AScheme(Fe/De) + max(-Fe, 0.0)
        aW = Dw * AScheme(Fw/Dw) + max( Fw, 0.0)
        aP = aE + aW + (Fe - Fw)
    
        S[i, i + 1] = -aE
        S[i, i - 1] = -aW
        S[i, i] = aP
    
        rhs[i] = f[i - 1] * dx

    #dirichlet    
    aE = a[0] / (dx) - v[0] / 2.
    aW = a[0] / (dx) + v[0] / 2.
    S[0, 0] += 1.
    S[1, 0] += -aE
    S[1, 1] += aE

    S[N+1, N+1] += 1.
    S[N, N+1] += -aW
    S[N, N] += aW

    rhs[0] += uDir[0]
    rhs[N+1] += uDir[1]

    return S, rhs
    
x = np.linspace(.0, 1.0, 11)
dx = x[1]-x[0]

grid = pg.createGrid(x=x)
N = grid.cellCount()

# force vector per cell
f = pg.RVector(grid.cellCount(), 0.0)
# diffusions coefficient
a = pg.RVector(grid.cellCount(), 2.1)
# velocity per cell [x-direction]
v = pg.RVector(grid.cellCount(), 10.1)

print('Peclet-number:', v[0]/(a[0] / dx))
        
ud0=10
udN=0

def uAna(
    u = (1 - np.exp(-v * x / D))/(1. - np.exp(-v * L / D))


for scheme in ['CDS', 'UDS', 'ES', 'HS', 'PS']:
    S, rhs = diffusionConvectionKernelGrid(grid, a, f, v, uDir=[ud0, udN], scheme=scheme)
    u = np.linalg.solve(S, rhs)
    plt.plot(pg.x(grid.cellCenter()), u[1:N+1], 'o-', label='FVM-'+scheme)

plt.legend()
S, rhs = diffusionConvectionKernelGrid(grid, a, f, v, uDir=[ud0, udN], scheme='HS')
print(S)
print(rhs)

#plt.show()
vC = np.vstack((v,np.zeros(len(v)))).T
S, rhs = diffusionConvectionKernel(grid, a=a, f=f, vel=vC, 
                                   uBoundaries=pg.solver.parseArgToBoundaries([[1,ud0],[2,udN]], grid),
                                   scheme='HS')
print(S)
print(rhs)
u = np.linalg.solve(S, rhs)
plt.plot(pg.x(grid.cellCenter()), u[0:N], 'o-', label='FVM-M (mesh)')

ws = WorkSpace()
uFV = solveFiniteVolume(grid, a=a, f=f, vel=vC,
                        uBoundary=[[1,ud0],[2,udN]], ws=ws, scheme='HS')
plt.plot(pg.x(grid.cellCenter()), uFV, 'o-', label='FVM-2 (mesh)')

#solver.showSparseMatrix(ws.Sm)
print(ws.rhs)

#uFV = solveFiniteVolume(grid, a=a, f=f,
                      #uBoundary=[[1,ud0],[2,udN]])
#plt.plot(pg.x(grid.cellCenter()), uFV, 'o-', label='FVM-2 (mesh)')
    
##fem reference
#uFEM = solver.solvePoisson(grid, a=a, f=f,                        
                           #uBoundary=[[1,ud0],[2,udN]])
#plt.plot(pg.x(grid.positions()), uFEM, 'x-', label='FEM')

plt.legend()

plt.show()
#drawMesh(ax, grid)