#!/usr/bin/env python

"""
Test multi
"""

import sys
import time

import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation

import pygimli as pg
from pygimli.viewer import *
from pygimli.solver import *
from pygimli.meshtools import *

from pygimli.physics.seismics import ricker, wiggle, solvePressureWave
from pygimli.physics.gravimetry import solveGravimetry

import pybert as pb
import pybert.dataview

def createCacheName(base, mesh, times):
    return 'cache-' + base + '-' + str(mesh.nodeCount()) + '-' + str(len(times))

class ERT():
    
    def __init__(self, verbose=False):
        
        self.fop = pb.DCSRMultiElectrodeModelling(verbose=verbose)
        #self.fop.regionManager().region(1).setBackground(True)
        #self.fop.createRefinedForwardMesh(refine=True, pRefine=False)

        self.inv = pg.RInversion(verbose=verbose, dosave=False)

        self.datTrans = pg.RTransLog()
        self.modTrans = pg.RTransLog()

        self.inv.setMaxIter(10)
        self.inv.setTransData(self.datTrans)
        self.inv.setTransModel(self.datTrans)
        #self.inv.setError(data('err'))
        #self.inv.setModel(pg.RVector(fop.regionManager().parameterCount(),
                          #pg.median(data('rhoa'))))
        
        self.inv.setLambda(5)
        self.schemeMg = pb.dataview.dataview.DataSchemeManager()
         
    def show(self, data, values=None, axes=None,
             cMin=None, cMax=None, colorBar=1, **kwargs):
        
        if isinstance(data, pg.DataContainer):
            scheme = kwargs.pop('scheme', 'unknown')
            pseudoscheme = getattr(pb.dataview.Pseudotype, scheme)
            
            if values is None:
                values = data('rhoa')
            
            gci = pb.dataview.dataview.drawDataAsMatrix(axes, data, 
                                                        values,
                                                        pseudotype=pseudoscheme)
            gci.set_clim((cMin, cMax))
            if colorBar:
                pg.mplviewer.colorbar.createColorbar(gci, nLevs=5,
                                                     cMin=cMin,  cMax=cMax,
                                     label='Apparent resistivity in $\Omega m$',
                                     **kwargs)
        else:
            print(type(data))
            raise
        return gci
        
    def createData(self, sensors, scheme):
        scheme = self.schemeMg.scheme(scheme)
        scheme.setInverse(False)
        scheme.addInverse(False)
        return scheme.create(sensorList=sensors)
        
    def simulate(self, mesh, resistivity, data):
        self.fop.setMesh(mesh)
        self.fop.setData(data)
        res = self.fop.response(resistivity)
        data = pb.DataContainerERT(self.fop.data())
        data.set('rhoa', res)
        return data
        

def createTestWorld(maxArea=0.2, verbose=0):
    boundary = []
    boundary.append([-20.0,   0.0]) #0
    boundary.append([-20.0,  -2.0]) #1
    boundary.append([-20.0, -12.0]) #2
    boundary.append([-20.0, -20.0]) #3
    boundary.append([ 20.0, -20.0]) #4
    boundary.append([ 20.0, -12.0]) #5
    boundary.append([ 20.0,  -2.0]) #6
    boundary.append([ 20.0,   0.0]) #7
       
    poly = pg.Mesh(2)
    nodes = [poly.createNode(b) for b in boundary]

    boundaryMarker = [1, 2, 1, 3, 1, 4, 1, 1]
    for i in range(len(boundary)):
        poly.createEdge(nodes[i], nodes[(i+1)%len(boundary)], boundaryMarker[i]) # dirichlet (inflow)
    
    poly.createEdge(nodes[1], nodes[6], 5)
    poly.createEdge(nodes[2], nodes[5], 6)
        
    boundary = []
    boundary.append([-6.0, -4.0]) #0
    boundary.append([-6.0, -8.0]) #1
    boundary.append([ 6.0, -8.0]) #2
    boundary.append([ 6.0, -4.0]) #3
    
    nodes = [poly.createNode(b) for b in boundary]
    for i in range(len(boundary)):
        poly.createEdge(nodes[i], nodes[(i+1)%len(boundary)], 10) # dirichlet (inflow)
  
    poly.addRegionMarker((-19, -1), 1)
    poly.addRegionMarker((-19, -19), 1)
    poly.addRegionMarker((-19, -10), 2)
    poly.addRegionMarker((  0, -7), 3)
    

    mesh = createMesh(poly, quality=33, area=maxArea, smooth=[1,10], verbose=verbose)
    return mesh

def calcVelocity(mesh, viscosity):
    velBoundary=[[1, [0.0, 0.0]],
                 [2, [1.0, 0.0]],
                 [4, [1.0, 0.0]],
                 [3, [0.0, 0.0]]
                ]
    preBoundary=[[1, 0.0],
                 [3, 0.0]]

    solutionName = createCacheName('vel', mesh, times)
    try:
        #vel = pg.load(solutionName + '.bmat')
        vel = np.load(solutionName + '.bmat.npy')
        
    except Exception as e:
        print(e)
        print("Building .... ")
        vel, pres, pCNorm, divVNorm = solveStokes_NEEDNAME(mesh, velBoundary,
                                                           preBoundary,
                                                           viscosity=viscosity,
                                                           maxIter=200,
                                                           verbose=1)
        np.save(solutionName + '.bmat', vel)
    return vel

def calcConcentration(mesh, vel, times):

    f = pg.RVector(mesh.cellCount(), 0.0)
    sourceCell=mesh.findCell([-18., -6])
    f[sourceCell.id()]=0.4

    Peclet = 50.0
   
    solutionName = createCacheName('conc', mesh, times)
    
    try:
        #vel = pg.load(solutionName + '.bmat')
        conc = np.load(solutionName + '.bmat.npy')
        
    except Exception as e:
        print(e)
        print("Building .... ")
        uMesh1 = solveFiniteVolume(mesh, a=1./Peclet, f=f, vel=vel, times=times, 
                          uBoundary=[2, 0],
                          scheme='PS', verbose=10)
        uMesh2 = solveFiniteVolume(mesh, a=1./Peclet, f=0, vel=vel, times=times, 
                          uBoundary=[2, 0], u0=uMesh1[-1],
                          scheme='PS', verbose=10)
        conc = np.vstack((uMesh1, uMesh2[1:]))
        np.save(solutionName + '.bmat', conc)

    return conc

def calcGravKernel(mesh):
    gravPointsX = np.arange(-19, 19.1, 1)
    gravPoints = np.vstack((gravPointsX, np.zeros(len(gravPointsX)))).T
    solutionName = createCacheName('grav', mesh, times)
    try:
        #vel = pg.load(solutionName + '.bmat')
        Gdg = np.load(solutionName + '.bmat.npy')
        
    except Exception as e:
        print(e)
        print("Building .... ")
        Gdg, Gdgz = solveGravimetry(mesh, None, pnts=gravPoints) 
        np.save(solutionName + '.bmat', Gdg)
        
    return gravPointsX, Gdg
        
def transResistivity(mesh, viscosity, concentration, meshI):
    viscERT = pg.interpolate(mesh, viscosity, meshI.cellCenters())
    viscERT = pg.solver.fillEmptyToCellArray(meshI, viscERT)

    Con = pg.RMatrix(len(concentration), len(concentration[0]))
    for i in range(len(concentration)):
        Con[i] = concentration[i]
    ConI = pg.RMatrix(len(concentration), meshI.cellCount())
    pg.interpolate(mesh, Con, meshI.cellCenters(), ConI)
    
    resistivity = pg.RMatrix(len(concentration), meshI.cellCount())
    for i in range(len(concentration)):
        resistivity[i] = 1./(1./viscERT + pg.abs(12.*ConI[i]))
    return resistivity 

def calcApparentResistivities(mesh, resistivities):
    ert = ERT(verbose=False)
    ertPointsX = [pg.RVector3(x,0) for x in np.arange(-19, 19.1, 1)]
    ertScheme = ert.createData(ertPointsX, scheme="Dipole Dipole (CC-PP)")
    
    solutionName = createCacheName('appRes', mesh, times)+ "-" + str(ertScheme.size())
    
    try:
        resA = np.load(solutionName + '.bmat.npy')
        
    except Exception as e:
        print(e)
        print("Building .... ")
        resA = np.zeros((len(resistivities), ertScheme.size()))
        ertScheme.set('k', pb.geometricFactor(ertScheme))
        ertData = ert.simulate(meshERT_FOP, resistivities[0], ertScheme)
        
        errPerc = 1
        errVolt = 1e-5
        voltage = ertData('rhoa') / ertData('k')
        ertData.set('err', pg.abs(errVolt / voltage) + errPerc / 100.0)
        print('err min:', min(ertData('err'))*100, 'max:', max(ertData('err'))*100)
        
        ertData.save('ert-0.dat', 'a b m n rhoa err k')
        
        for i in range(0, len(resistivities)):
            tic =  time.time()
            resA[i] = ert.fop.response(resistivities[i])
                        
            rand = pg.RVector(len(resA[i]))
            pg.randn(rand)
        
            resA[i] *= (1.0 + rand * ertData('err'))
            
            print(i, "/", len(resistivities), " : ", time.time()-tic, "s",
                  "min:", min(resistivities[i]), "max:", max(resistivities[i]),
                  "min:", min(resA[i]), "max:", max(resA[i]) )

        np.save(solutionName + '.bmat', resA)
    return resA, ert, ertScheme


vis=1

swatch = pg.Stopwatch(True)
pg.showLater(1)

mesh = createTestWorld(maxArea=0.2, verbose=0)
meshC = createTestWorld(maxArea=1, verbose=0)

meshERT_FOP = appendTriangleBoundary(createTestWorld(maxArea=0.1, verbose=0), 
                                    xbound=50, ybound=50, marker=1,
                                    quality=34.0, smooth=False, markerBoundary=1,
                                    isSubSurface=False, verbose=False)
print("meshgen:", swatch.duration(True))

times = np.linspace(0, 40, 50)

viscosity = pg.solver.parseMapToCellArray([[1, 3], [2, 1], [3, 12.0]], mesh)
viscosity *= 10.

dpi=92
fig = plt.figure(facecolor='white', figsize=(2*800/dpi, 2*490/dpi), dpi=dpi)  
axVis = fig.add_subplot(3,2,5)
axVel = fig.add_subplot(3,2,6)

axCon = fig.add_subplot(3,3,5)
axGra = fig.add_subplot(3,3,2)
axRes = fig.add_subplot(3,3,6)
axReA = fig.add_subplot(3,3,3)
   
    
if vis:
    show(mesh, data=viscosity, colorBar=1, 
         orientation='vertical', label='Viscosity', axes=axVis)
    show(mesh, axes=axVis)

vel = calcVelocity(mesh, viscosity)

if vis:
    show(mesh, data=pg.cellDataToPointData(mesh,
                                           np.sqrt(vel[:,0]*vel[:,0] + vel[:,1]*vel[:,1])),
         logScale=0, colorBar=1, orientation='vertical', label='|Velocity| in m/s',
         cbar='b2r', axes=axVel)

    pg.mplviewer.drawMeshBoundaries(axVel, mesh, fitView=True)
    show(mesh, data=vel, axes=axVel, coarseMesh=meshC)
    plt.pause(0.1)
print("vel:", swatch.duration(True))


conc = calcConcentration(mesh, vel, times)


if vis:
    gciCon= pg.mplviewer.drawModel(axCon, mesh, data=conc[1],
                                   cMin=0, cMax=0.03)
    cbar = createColorbar(gciCon, orientation='vertical', label='Concentration')
print("con:", swatch.duration(True))


resistivity = transResistivity(mesh, viscosity, conc, meshERT_FOP)

if vis:
    gciRes = pg.mplviewer.drawModel(axRes, meshERT_FOP, 
                                    data=resistivity[0],
                                    cMin=1, cMax=100)
    cbar = createColorbar(gciRes, orientation='vertical', label='Resistivity')
    axRes.set_xlim((-20, 20))
    axRes.set_ylim((-20, 00))
print("res:", swatch.duration(True))
    
dz = np.zeros(len(conc))
densityBrine = 1.2 # g/cm^3
densityBrine *= 1000. # kg/m^3

gravPointsX, Gdg = calcGravKernel(mesh)
print("grav:", swatch.duration(True))

#meshERT = pg.meshtools.createParaMesh2dGrid(ertPointsX,
                                            #paraDZ=0.5)

resA, ert, ertScheme = calcApparentResistivities(meshERT_FOP, resistivity)
gciARes = ert.show(ertScheme, values=resA[0], axes=axReA, 
                   scheme='DipoleDipole',
                   cMin=5, cMax=35, orientation='vertical')
print("appRes:", swatch.duration(True))

def animate(i):
    #for i in range(1, len(conc)):
    tic = time.time()
    if vis:
        pg.mplviewer.setMappableData(gciCon, conc[i], cMin=0, cMax=0.03,
                                     logScale=False)
        gciCon.set_clim(0, 0.03)
        
    axGra.clear()
    
    dDensity = densityBrine * conc[i]
    print(time.time()-tic); tic = time.time()
    
    # calc gravimetric response
    dz = Gdg.dot(dDensity)[:,2]
    axGra.clear()
    axGra.plot(gravPointsX, dz)
    axGra.plot(gravPointsX, gravPointsX*0.0+0.03, 'v', color='black')
    axGra.set_ylabel('Grav in mGal')
    axGra.set_xlim((-20,20))
    axGra.set_ylim((0,0.003))
    axGra.grid((0,0.003))
    
    print(time.time()-tic); tic = time.time()
    axReA.clear()
    ert.show(ertScheme, values=resA[i], axes=axReA, scheme='DipoleDipole',
             cMin=5, cMax=35, colorBar=0)
     
    if vis:
        pg.mplviewer.setMappableData(gciRes, 
                                     resistivity[i],
                                     cMin=1, cMax=100,
                                     logScale=True)
    print(i, time.time()-tic, 
          "sum:", sum(conc[i]),
          "dsum:", (sum(conc[i])-sum(conc[i-1])),
          )
    plt.pause(0.001)

anim = animation.FuncAnimation(fig, animate,
                               frames=int(len(conc)),
                               interval=1)#, blit=True)

solutionName = createCacheName('all', mesh, times)+ "-" + str(ertScheme.size())
anim.save(solutionName + ".mp4", writer=None, fps=20, dpi=dpi, codec=None,
          bitrate=24*1024, extra_args=None, metadata=None,
          extra_anim=None, savefig_kwargs=None)


pg.showNow()
