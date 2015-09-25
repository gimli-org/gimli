#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

"""

import pygimli as pg
import numpy as np

from pygimli.meshtools import polytools as pt
from pygimli.meshtools import createMesh
import pygimli.solver as solver
    
def createModel(maxArea=0.2, nx=80, grid=True):

    mesh=None
    dens=None
    if grid:
        nx = nx
        ny = nx/2
        mesh = pg.createGrid(x=np.linspace(-10, 10, nx),
                            y=np.linspace(-10, 0, ny))
   
        density = mesh.cellAttributes()*0.0 + 1.0
        xstart = (nx-1)*int(ny/1.3) + int(nx/2)
        yoff =nx-1
        for i in range(int(ny/20)):
            density[xstart-nx*0.35+i*yoff: xstart+nx*0.35+i*yoff] *= 2
    else:
        layer1 = pt.createRectangle(start=[-10, 0], end=[10, -10], 
                                 marker=1, boundaryMarker=[1, 3, 2, 4])
        
        block = pt.createRectangle(start=[-10+0.35*10, -1.8918918918918912], 
                                end=[10-0.35*10, -2.2972972972972965], 
                                    marker=2)
        mesh = createMesh([layer1, block], quality=32, area=maxArea,
                          smooth=[1,10])
    
        density = mesh.cellAttributes()*0.0 + 1.0
        density[mesh.cellMarker()==2] = 2
    
    return mesh, density

def createModel2(maxArea=0.2, nx=20, grid=True):
    mesh=None
    if grid:
        nx = nx
        ny = nx
        mesh = pg.createGrid(x=np.linspace(-10, 10, nx),
                             y=np.linspace(0, 20, ny))
   
    else:
        layer1 = pt.createRectangle(start=[-10, 20], end=[10, 10], 
                                    marker=2, boundaryMarker=[1, -1, 2, 4])
        layer2 = pt.createRectangle(start=[-10, 10], end=[10, 0], 
                                    marker=1, boundaryMarker=[1, 3, 2, -1])
        
        mesh = createMesh([layer1, layer2], quality=32, area=maxArea,
                          smooth=[1,10])
        
    density = mesh.cellAttributes()*0.0 + 1.0
    
    for c in mesh.cells():
        if c.center().y() > 10:
            density[c.id()] = 2.0
    
    
    return mesh, density


def calc(out, *args):
    mesh, density = args[0]
    
    print(mesh)
            
    velBoundary=[[1, [0.0,  0.0]],
                [2, [0.0,  0.0]],
                [3, [0.0,  0.0]],
                [4, [0.0,  0.0]]]
    preBoundary=[[1, 0.0],
                [2, 0.0],
                [3, 0.0],
                [4, 0.0],
                ]


    densMatrix = pg.RMatrix()
    vels = []

    swatch = pg.Stopwatch(True)
    class WS():
        pass

    wsfv = WS()

    ax,_ = pg.show(mesh, density)

    nSteps = 300
    dt = 0.01
    dtSteps = 10

    meshC = pg.createGrid(x=np.linspace(-10, 10, 41),
                        y=np.linspace(-10, 0, 21))

    vel=None
    pre=None
    for i in range(nSteps):
        print(i, 'dens', min(density), max(density))
        
        densMatrix.push_back(density)
        vel, pre, pCNorm, divVNorm = solver.solveStokes(mesh, 
                                                        velBoundary=velBoundary,
                                                        preBoundary=preBoundary,
                                                        viscosity=0.1,
                                                        density=density,
                                                        pre0 = pre,
                                                        vel0 = vel,
                                                        f=[0,-9.81],
                                                        maxIter=100,
                                                        tol=1e-4,
                                                        verbose=0,
                                                        ws=wsfv)
        vels.append(vel)
        
        print("stokes:" , swatch.duration(True))
        dens2 = solver.solveFiniteVolume(mesh, a=1./5000, u0=density, vel=vel,
                                        times=np.linspace(0, dt, dtSteps), 
                                        #uBoundary=[4, 0],
                                        scheme='PS', verbose=0)
        print("Convekt:" , swatch.duration(True))
        density=dens2[-1]
        
        
        
        ax.clear()
        pg.show(mesh, density, axes=ax)
        #pg.show(mesh, vel, coarseMesh=meshC, axes=ax, color='white')

    mesh.save(out)
    meshC.save(out+'C')
    densMatrix.save(out+'density')
    np.save(out+'velo.bmat', vels)
    createAnimation(out)
    createAnimation(out, stream=True)

def createAnimation(out, stream=False):
    mesh = pg.load(out + '.bms')
    meshC = pg.load(out + 'C.bms')
    densMatrix = pg.load(out + 'density.bmat')
    vels = np.load(out+'velo.bmat.npy')
    
    if stream:
        pg.mplviewer.saveAnimation(mesh, densMatrix, out + '-stream',
                               plc=None, label='',
                               cMin=2, cMax=3, logScale=True, cmap=None,
                               vData=vels, coarseMesh=meshC, color='white'
                               )
    else:
        pg.mplviewer.saveAnimation(mesh, densMatrix, out,
                               plc=None, label='',
                               cMin=2, cMax=3, logScale=True, cmap=None)


if __name__ == "__main__":
    
    calc('two-grid', createModel2(nx=250, grid=True))
    calc('two-mesh', createModel2(maxArea=0.005, grid=False))
    
    calc('omega-grid', createModel(nx=200, grid=True))
    calc('omega-mesh', createModel(maxArea=0.008, grid=False))
    



    #createAnimation(out)
    #createAnimation(out, stream=True)
    
    pg.wait()