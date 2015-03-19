#!/usr/bin/env python

"""
Test multi
"""

from pygimli.solver import parseArgToArray
import numpy as np

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


def calcSeismics(mesh, vP):
    
    meshSeis = appendTriangleBoundary(mesh, 
                                      xbound=50, ybound=50, marker=1,
                                      quality=33.0, area=0.5, smooth=True, 
                                      markerBoundary=1,
                                      isSubSurface=False, verbose=False)
    print(meshSeis)
    meshSeis = meshSeis.createH2()
    vP = pg.interpolate(mesh, vP, meshSeis.cellCenters())
    #mesh = meshSeis.createH2()
    mesh = meshSeis
    vP = pg.solver.fillEmptyToCellArray(mesh, vP)
    print(mesh)
    ax, cbar = pg.show(mesh, data=vP)
    pg.show(mesh, axes=ax)
    pg.showNow()
    
    h = pg.median(mesh.boundarySizes())
    dt = 0.25 * h /max(vP)
    tmax = 50./min(vP)
    times = np.arange(0.0, tmax, dt)
        
    geophPointsX = np.arange(-19, 19.1, 1)
    geophPoints = np.vstack((geophPointsX, np.zeros(len(geophPointsX)))).T
    
    solutionName = createCacheName('seis', mesh, times)
    try:   
        u = pg.load(solutionName + '.bmat')
        uI = pg.load(solutionName + 'I.bmat')
    except Exception as e:
        print(e)
        print("Building .... ")
        f0 = 1./dt*0.2
        
        print("h:", round(h,2), "dt:", round(dt,5), "1/dt:", round(1/dt,1), "f0", round(f0,2))
        
        uSource = ricker(f0, times, t0=1./f0)
    
        plt.plot(times, uSource, '-*')
        plt.show()
        u = solvePressureWave(mesh, vP, times, sourcePos=geophPoints[38],
                            uSource=uSource, verbose=10)
        
        u.save(solutionName)
        uI = pg.RMatrix()
        pg.interpolate(mesh, u, mesh.cellCenters(), uI)
        uI.save(solutionName+'I')
        
        
    nodes = [mesh.findNearestNode(p) for p in geophPoints]
    
    fig = plt.figure()
    axs = fig.add_subplot(1,1,1)
    drawSeismogramm(axs, mesh, u, nodes, dt, i=None)
    plt.show()
        
    fig = plt.figure()
    fig, ax = fig, ax = plt.subplots()
    gci = pg.mplviewer.drawModel(ax, mesh, data=uI[0], cmap='b2r')
        
    #ax, cbar = pg.show(mesh, data=vP)
    #pg.showNow()
    #ax = fig.add_subplot(1,1,1)
    for i in range(1, len(u), 5):
        ui = uI[i]
        ui = ui / max(pg.abs(ui))
        ui = pg.logDropTol(ui, 1e-2)
        cMax = max(pg.abs(ui))
        
        pg.mplviewer.setMappableData(gci, 
                                    ui,
                                    cMin=-cMax, cMax=cMax,
                                    logScale=False
                                    )
        #ax.clear()
        ##pg.show(mesh, axes=ax)
        #drawField(ax, mesh, data=ui,
                  #cMin=-cMax, cMax=cMax,
                  #cmap='RdBu')
        #ax.set_xlim((-25, 25))
        #ax.set_ylim((-25, 0))
        plt.pause(0.001)
    
