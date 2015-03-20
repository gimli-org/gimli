#!/usr/bin/env python

"""
Test multi
"""

import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation

from pygimli.solver import parseArgToArray
from pygimli.physics.seismics import ricker, solvePressureWave, drawSeismogramm
from pygimli.meshtools import *
import numpy as np

def createCacheName(base, mesh, times):
    return '_seis-cache-' + base + '-' + str(mesh.nodeCount()) + '-' + str(len(times))

def velocityVp(porosity, vMatrix=5000, vFluid=1442, S=1,
               mesh=None):
    r"""
    """
    porosity = parseArgToArray(porosity, mesh.cellCount(), mesh)
    vMatrix = parseArgToArray(vMatrix, mesh.cellCount(), mesh)
    vFluid = parseArgToArray(vFluid, mesh.cellCount(), mesh)
    S = parseArgToArray(S, mesh.cellCount(), mesh)
    
    vAir = 343.0
    
    vel = 1./(np.array((1.-porosity)/vMatrix) + \
              porosity * S / vFluid + \
              porosity * (1.-S)/vAir)
    return vel


def calcSeismics(meshIn, vP):
    
    #meshSeis = meshIn.createH2()
    meshSeis = meshIn
    meshSeis = appendTriangleBoundary(meshSeis, 
                                      xbound=25, ybound=22.0, marker=1,
                                      quality=32.0, area=0.3, smooth=True, 
                                      markerBoundary=1,
                                      isSubSurface=False, verbose=False)
    print(meshSeis)
    meshSeis = meshSeis.createH2()
    meshSeis = meshSeis.createH2()
    #meshSeis = meshSeis.createP2()
    meshSeis.smooth(1, 1, 1, 4)
    vP = pg.interpolate(meshIn, vP, meshSeis.cellCenters())
    
    mesh = meshSeis
    vP = pg.solver.fillEmptyToCellArray(mesh, vP)

    print(mesh)
    #ax, cbar = pg.show(mesh, data=vP)
    #pg.show(mesh, axes=ax)
     
    geophPointsX = np.arange(-19, 19.1, 1)
    geophPoints = np.vstack((geophPointsX, np.zeros(len(geophPointsX)))).T
    sourcePos=geophPoints[4]
     
    c = mesh.findCell(sourcePos)
    h1 = pg.findBoundary(c.boundaryNodes(0)).size()
    h2 = pg.findBoundary(c.boundaryNodes(1)).size()
    h3 = pg.findBoundary(c.boundaryNodes(2)).size()
    print([h1, h2, h3])
    h = pg.median([h1, h2, h3])
    
    #h = pg.median(mesh.boundarySizes())
    dt = 0.5 * h /max(vP)
    cfl = max(vP)*dt/h
    print("Courant-Friedrich-Lewy-Zahl:", cfl)
    
    tmax = 50./min(vP)
    times = np.arange(0.0, tmax, dt)
        
    
    
    solutionName = createCacheName('seis', mesh, times)
    try:   
        U = None
        #u = pg.load(solutionName + '.bmat')
        uI = pg.load(solutionName + 'I.bmat')
    except Exception as e:
        print(e)
        f0 = 1./dt*0.2 
        print("h:", round(h,2),
              "dt:", round(dt,5),
              "1/dt:", round(1/dt,1),
              "f0", round(f0,2),
              "Wavelength: ", round(max(vP)/f0, 2), " m")
        
        uSource = ricker(f0, times, t0=1./f0)
    
        plt.figure()
        plt.plot(times, uSource, '-*')
        plt.show(block=0)
        plt.pause(0.01)
        u = solvePressureWave(mesh, vP, times, sourcePos=sourcePos,
                            uSource=uSource, verbose=10)
        
        u.save(solutionName)
        uI = pg.RMatrix()
        print("interpolate node to cell data ... ")
        pg.interpolate(mesh, u, mesh.cellCenters(), uI)
        print("... done")
        uI.save(solutionName+'I')
    
    #nodes = [mesh.findNearestNode(p) for p in geophPoints]
    
    #fig = plt.figure()
    #axs = fig.add_subplot(1,1,1)
    #drawSeismogramm(axs, mesh, u, nodes, dt, i=None)
    #plt.show()
    
    dpi=92
    scale=1
    fig = plt.figure(facecolor='white',
                     figsize=(scale*800/dpi, scale*490/dpi), dpi=dpi)  
    ax = fig.add_subplot(1,1,1)
    gci = pg.mplviewer.drawModel(ax, mesh, data=uI[0],
                                 cMin=-1, cMax=1, cmap='bwr')
    pg.mplviewer.drawMeshBoundaries(ax, meshIn, hideMesh=1)
    ax.set_xlim((-21, 21))
    ax.set_ylim((-16, 0))
    plt.tight_layout()
    #ax, cbar = pg.show(mesh, data=vP)
    #pg.showNow()
    #ax = fig.add_subplot(1,1,1)
    def animate(i):
        i=i*5
        if i >len(uI)-1:
            return
        print("Frame:", i, "/", len(uI))
        ui = uI[i]
        ui = ui / max(pg.abs(ui))
        ui = pg.logDropTol(ui, 1e-2)
        cMax = max(pg.abs(ui))
        
        pg.mplviewer.setMappableData(gci, 
                                    ui,
                                    cMin=-cMax, cMax=cMax,
                                    logScale=False
                                    )
        
        #plt.pause(0.001)
    
    anim = animation.FuncAnimation(fig, animate,
                                   frames=int(len(uI)/5),
                                   interval=0.001, repeat=0)#, blit=True)
    out = 'seis'
    anim.save(out + ".mp4", writer=None, fps=20, dpi=dpi, codec=None,
              bitrate=24*1024, extra_args=None, metadata=None,
              extra_anim=None, savefig_kwargs=None)
    try:
        print("create frames ... ")
        os.system('mkdir -p anim-' + out)
        os.system('ffmpeg -i ' + out + '.mp4 anim-' + out + '/movie%d.jpg')
    except:
        pass   