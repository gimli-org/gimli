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
    
    """
        ERT Manager. should solve most of the common problems
    """
        
    def __init__(self, verbose=False):
        """
        """
        self.fop = pb.DCSRMultiElectrodeModelling(verbose=verbose)
        #self.fop.regionManager().region(1).setBackground(True)
        #self.fop.createRefinedForwardMesh(refine=True, pRefine=False)

        self.inv = pg.RInversion(verbose=verbose, dosave=False)

        self.datTrans = pg.RTransLog()
        self.modTrans = pg.RTransLog()

        self.inv.setTransData(self.datTrans)
        self.inv.setTransModel(self.datTrans)
        #self.inv.setError(data('err'))
        #self.inv.setModel(pg.RVector(fop.regionManager().parameterCount(),
                          #pg.median(data('rhoa'))))
        
        self.schemeMg = pb.dataview.dataview.DataSchemeManager()
        self.paraMesh = None
         
    def show(self, data, values=None, axes=None,
             cMin=None, cMax=None, colorBar=1, **kwargs):
        """
        """
        gci = None
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
        elif isinstance(data, pg.RVector):
            # assuming this is a model from the last inversion run
            pg.show(self.fop.regionManager().paraDomain(), data, 
                    axes=axes, **kwargs)
            
        else:
            print(type(data))
            raise
        return gci
        
    def invert(self, data, values=None, verbose=0, **kwargs):
        """
        Invert the given data.
        
        A parametric mesh for the inversion will be created if non is given
        before.
            
        Parameters
        ----------
        """
        self.fop.setVerbose(verbose)
        self.inv.setVerbose(verbose)
        self.inv.setMaxIter(kwargs.pop('maxiter', 10))
        self.inv.setLambda(kwargs.pop('lambd', 10))
        
        if self.paraMesh is None:
            self.paraMesh = createParaMesh2dGrid(data.sensorPositions(),
                                                 **kwargs)
            self.setParaMesh(self.paraMesh)
            if verbose:
                print(self.paraMesh)
                #pg.show(self.paraMesh)
                
        err = data('err')
        rhoa = data('rhoa')
                    
        startModel = pg.RVector(self.fop.regionManager().parameterCount(),
                                pg.median(rhoa))

        self.fop.setData(data)
        self.inv.setForwardOperator(self.fop)
        
        # check err here 
        self.inv.setData(rhoa)
        self.inv.setError(err)
        self.inv.setModel(startModel)
        
        model = self.inv.run()
        
        if values is not None:
            
            if isinstance(values, pg.RVector):
                values = [values]
            elif isinstance(values, np.ndarray):
                if values.ndim == 1:
                    values = [values]
            
            allModel = pg.RMatrix(len(values)+1, len(model))
            allModel[0] = model
            self.inv.setVerbose(False)
            for i in range(1, len(values)):
                tic=time.time()
                self.inv.setModel(model)
                self.inv.setReferenceModel(model)
                dData = pg.abs(values[i] / rhoa)
                
                relModel = self.inv.invSubStep(pg.log(dData))
                allModel[i] = model * pg.exp(relModel)
                print(i, "/", len(values), " : ", time.time()-tic, "s")
                
            return allModel
        return model

    def setParaMesh(self, mesh):
        """
        Set the parameter mesh for any inversion.
            
        Parameters
        ----------
        """
        self.fop.setMesh(self.paraMesh)
        self.fop.regionManager().region(1).setBackground(True)
        self.fop.createRefinedForwardMesh(refine=True, pRefine=False)
        
    def createData(self, sensors, scheme):
        """
        Create a empty data file. i.e. a measurement scheme.
            
        Parameters
        ----------
        """
        scheme = self.schemeMg.scheme(scheme)
        scheme.setInverse(False)
        scheme.addInverse(False)
        return scheme.create(sensorList=sensors)
        
    def simulate(self, mesh, resistivity, data):
        """
        Simulation a ERT measurement.
        
        Perform the forward task for a given mesh, 
        a resistivity distribution and a measurement scheme.
            
        Parameters
        ----------
        """
        self.fop.setMesh(mesh)
        self.fop.setData(data)
        res = self.fop.response(resistivity)
        data = pb.DataContainerERT(self.fop.data())
        data.set('rhoa', res)
        return data
        

def createTestWorld1(maxArea=0.2, verbose=0):
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

    boundaryMarker = [1, 2, 3, 4, 5, 6, 7, 8]
    for i in range(len(boundary)):
        poly.createEdge(nodes[i], nodes[(i+1)%len(boundary)], boundaryMarker[i]) # dirichlet (inflow)
    
    poly.createEdge(nodes[1], nodes[6], 9)
    poly.createEdge(nodes[2], nodes[5], 10)
        
    boundary = []
    boundary.append([-6.0, -4.0]) #0
    boundary.append([-6.0, -8.0]) #1
    boundary.append([ 6.0, -8.0]) #2
    boundary.append([ 6.0, -4.0]) #3
    
    nodes = [poly.createNode(b) for b in boundary]
    for i in range(len(boundary)):
        poly.createEdge(nodes[i], nodes[(i+1)%len(boundary)], 11) # dirichlet (inflow)
  
    poly.addRegionMarker((-19, -1),  1)
    poly.addRegionMarker((-19, -10), 2)
    poly.addRegionMarker((-19, -19), 3)
    poly.addRegionMarker((  0, -7),  4)
    
    mesh = createMesh(poly, quality=33, area=maxArea,
                      smooth=[1,10], verbose=verbose)
    
    velBoundary=[[2, [1.0, 0.0]],
                 [6, [1.0, 0.0]],
                 [3, [1./3., 0.0]],
                 [8, [0.0, 0.0]],
                 [4, [0.0, 0.0]]
                ]
    preBoundary=[[1, 0.0],
                 #[3, 0.0],
                 [5, 0.0],
                 [7, 0.0]]
    
    viscosity = pg.solver.parseMapToCellArray([[1, 3.0], 
                                               [2, 1.0], 
                                               [3, 6.0],
                                               [4, 12.0]], mesh)
    # 30 works nice with tol=1e-2 and pressureRelaxation = 0.01
    viscosity *= 10.

    return mesh, velBoundary, preBoundary, viscosity

def createTestWorld2(maxArea=0.2, verbose=0):
    boundary = []
    boundary.append([-20.0,   0.0]) #0
    boundary.append([-20.0,  -1.0]) #1
    boundary.append([-20.0,  -4.0]) #2
    boundary.append([-20.0,  -5.0]) #3
    boundary.append([-20.0,  -8.0]) #4
    boundary.append([-20.0, -12.0]) #5
    boundary.append([ 20.0, -12.0]) #6
    boundary.append([ 20.0,  -8.0]) #7
    boundary.append([ 20.0,  -5.0]) #8
    boundary.append([ 20.0,  -4.0]) #9
    boundary.append([ 20.0,  -1.0]) #10
    boundary.append([ 20.0,   0.0]) #11
       
    poly = pg.Mesh(2)
    nodes = [poly.createNode(b) for b in boundary]

    boundaryMarker = [1, 5, 1, 6, 1, 3, 1, 7, 1, 8, 1, 4]
    for i in range(len(boundary)):
        poly.createEdge(nodes[i], nodes[(i+1)%len(boundary)], boundaryMarker[i]) # dirichlet (inflow)
    
    poly.createEdge(nodes[1], nodes[10], 9)
    poly.createEdge(nodes[4], nodes[7], 9)
    
    boundary = []
    boundary.append([-0.5, -4.0]) #12
    boundary.append([-0.5, -5.0]) #13
    boundary.append([ 0.5, -5.0]) #14
    boundary.append([ 0.5, -4.0]) #15
    
    
    nodes2 = [poly.createNode(b) for b in boundary]
    for i in range(len(boundary)):
        poly.createEdge(nodes2[i], nodes2[(i+1)%len(boundary)], 11)
  
    poly.createEdge(nodes[2], nodes2[0], 10)
    poly.createEdge(nodes[3], nodes2[1], 10)
    poly.createEdge(nodes2[2], nodes[8], 10)
    poly.createEdge(nodes2[3], nodes[9], 10)

    poly.createEdge(nodes2[0], nodes2[1], 10)
    poly.createEdge(nodes2[2], nodes2[3], 10)    
        
    poly.addRegionMarker((-19, -0.1), 1)
    poly.addRegionMarker((-19, -1.1), 3)
    poly.addRegionMarker((-19, -4.1), 5)
    poly.addRegionMarker(( 19, -4.1), 5)
    poly.addRegionMarker((-19, -5.1), 4)
    poly.addRegionMarker((-19, -8.1), 2)
    poly.addRegionMarker((0.0, -5.0), 6, 0.02)
    
    poly.addHoleMarker((-19, -4.1))
    poly.addHoleMarker(( 19, -4.1))
    

    mesh = createMesh(poly, quality=33, area=maxArea, smooth=[1,10], verbose=verbose)
    
    velBoundary=[[1, [0.0, 0.0]],
                 [5, [1.0, 0.0]],
                 [8, [0.0, 0.0]],
                 #[6, [0.0, 0.0]],
                 #[7, [-1, 0.0]],
                 [4, [0.0, 0.0]],
                 [3, [0.0, 0.0]],
                 [10,[0.0, 0.0]]
                ]
    preBoundary=[[6, 0.0],
                 [7, 0.0]
                 ]
    
    viscosity = pg.solver.parseMapToCellArray([[1, 1.0],
                                               [2, 1.0],
                                               [3, 1.0],
                                               [4, 1.0],
                                               [5, 1.0],
                                               [6, 1.0]], 
                                              mesh)
    viscosity *= 10.
    
    print("div(V) = ", pg.solver.divergence(mesh, normMap=velBoundary))
    
    return mesh, velBoundary, preBoundary, viscosity

def calcVelocity(mesh, viscosity, velBoundary, preBoundary):

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
                                                           maxIter=1000,
                                                           tol=1e-4,
                                                           verbose=1)
        np.save(solutionName + '.bmat', vel)
    
        #plt.semilogy(divVNorm)
        #pg.show(mesh, pg.cellDataToPointData(mesh, pres), cmap='b2r')
        
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
    viscI = pg.interpolate(mesh, viscosity, meshI.cellCenters())
    viscI = pg.solver.fillEmptyToCellArray(meshI, viscI)

    Con = pg.RMatrix(len(concentration), len(concentration[0]))
    for i in range(len(concentration)):
        Con[i] = concentration[i]
    ConI = pg.RMatrix(len(concentration), meshI.cellCount())
    pg.interpolate(mesh, Con, meshI.cellCenters(), ConI)
    
    resis = pg.RMatrix(len(concentration), meshI.cellCount())
    for i in range(len(concentration)):
        resis[i] = 1./(1./viscI + pg.abs(12.*ConI[i]))
    return resis

def calcApparentResistivities(mesh, resistivities):
    ert = ERT(verbose=False)
    ertPointsX = [pg.RVector3(x,0) for x in np.arange(-19, 19.1, 1)]
    ertScheme = ert.createData(ertPointsX, scheme="Dipole Dipole (CC-PP)")
    
    solutionName = createCacheName('appRes', mesh, times)+ "-" + str(ertScheme.size())
    
    try:
        resA = np.load(solutionName + '.bmat.npy')
        ertScheme = pb.DataContainerERT(solutionName + '.dat')
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
        ertData.save(solutionName + '.dat', 'a b m n rhoa err k')
        
    return resA, ert, ertData

def calcERT(ert, ertScheme, resA):
    solutionName = createCacheName('ERT', mesh, times)+ "-" + str(ertScheme.size())
    try:   
        ertModels = pg.load(solutionName + '.bmat')
        ertMesh = pg.load(solutionName + '.bms')
    except Exception as e:
        print(e)
        print("Building .... ")
        ertModels = ert.invert(ertScheme, values=resA, maxiter=10, 
                               paraDX=0.5, paraDZ= 0.5, nLayers=20, paraDepth=20,
                               lambd=50, verbose=1)
        ertMesh=ert.fop.regionManager().paraDomain()
        ertModels.save(solutionName + '.bmat')
        ertMesh.save(solutionName)
    return ertModels, ertMesh

vis = 1
mp4 = 1

swatch = pg.Stopwatch(True)
pg.showLater(1)

mesh, velBoundary, preBoundary, viscosity = createTestWorld1(maxArea=0.2,
                                                             verbose=0)
print(mesh)
meshC, tmp, tmp, tmp=createTestWorld1(maxArea=1, verbose=0)
meshD, tmp, tmp, tmp=createTestWorld1(maxArea=0.1, verbose=0)

meshERT_FOP = appendTriangleBoundary(meshD, 
                                    xbound=50, ybound=50, marker=1,
                                    quality=34.0, smooth=False, markerBoundary=1,
                                    isSubSurface=False, verbose=False)
print("meshgen:", swatch.duration(True))

times = np.linspace(0, 40, 50)

dpi=92
fig = None
if mp4:
    fig = plt.figure(facecolor='white', figsize=(2*800/dpi, 2*490/dpi), dpi=dpi)  
else:
    fig = plt.figure() 
    
axVis = fig.add_subplot(3,2,5)
axVel = fig.add_subplot(3,2,6)

axCon = fig.add_subplot(3,3,4)
axGra = fig.add_subplot(3,3,1)
axRes = fig.add_subplot(3,3,5)
axReA = fig.add_subplot(3,3,2)

axERT = fig.add_subplot(3,3,6)
axERR = fig.add_subplot(3,3,3)
   
    
if vis:
    show(mesh, data=viscosity, colorBar=1, 
         orientation='vertical', label='Viscosity', axes=axVis)
    show(mesh, axes=axVis)
    
vel = calcVelocity(mesh, viscosity, velBoundary, preBoundary)

if vis:
    axVel, cbar = show(mesh, data=pg.cellDataToPointData(mesh,
                                           np.sqrt(vel[:,0]*vel[:,0] + vel[:,1]*vel[:,1])),
         logScale=0, colorBar=1,
         orientation='vertical', label='|Velocity| in m/s',
         axes=axVel
         )

    pg.mplviewer.drawMeshBoundaries(axVel, mesh, fitView=True)
    pg.viewer.showBoundaryNorm(mesh, velBoundary, color='red', axes=axVel)
    show(mesh, data=vel, axes=axVel, coarseMesh=meshC)
    plt.pause(0.1)
print("vel:", swatch.duration(True))

conc = calcConcentration(mesh, vel, times)


if vis:
    gciCon= pg.mplviewer.drawModel(axCon, mesh, data=conc[1],
                                   cMin=0, cMax=0.03, 
                                   #cmap='b2r'
                                   )
    cbar = createColorbar(gciCon, orientation='vertical', label='Concentration')
print("con:", swatch.duration(True))


resistivities = transResistivity(mesh, viscosity, conc, meshERT_FOP)

if vis:
    gciRes = pg.mplviewer.drawModel(axRes, meshERT_FOP, 
                                    data=resistivities[0],
                                    cMin=1, cMax=100)
    cbar = createColorbar(gciRes, orientation='vertical', label='Resistivity')
    axRes.set_xlim((-20, 20))
    axRes.set_ylim((-20, 00))
print("res:", swatch.duration(True))
    
gravPointsX, Gdg = calcGravKernel(mesh)
print("grav:", swatch.duration(True))

resA, ert, ertScheme = calcApparentResistivities(meshERT_FOP, resistivities)
gciARes = ert.show(ertScheme, values=resA[0], axes=axReA, 
                   scheme='DipoleDipole',
                   cMin=5, cMax=35, orientation='vertical')
print("appRes:", swatch.duration(True))

ertModels, meshERT = calcERT(ert, ertScheme, resA)

if vis:
    gciERT = pg.mplviewer.drawModel(axERT, meshERT, 
                                    data=ertModels[0],
                                    cMin=1, cMax=100)
    cbar = createColorbar(gciERT, orientation='vertical', label='Resistivity')
    gciERR = pg.mplviewer.drawModel(axERR, meshERT, 
                                    data=ertModels[0]/ertModels[0],
                                    cMin=1/6, cMax=6, cmap='b2r')
    cbar = createColorbar(gciERR, orientation='vertical', label='Ratio')

print("ert:", swatch.duration(True))



dz = np.zeros(len(conc))
densityBrine = 1.2 # g/cm^3
densityBrine *= 1000. # kg/m^3

def animate(i):
    #for i in range(1, len(conc)):
    tic = time.time()
    if vis:
        pg.mplviewer.setMappableData(gciCon, conc[i], logScale=False)
        gciCon.set_clim(0, 0.03)
        
    axGra.clear()
    
    dDensity = densityBrine * conc[i]
    #print(time.time()-tic); tic = time.time()
    
    # calc gravimetric response
    dz = Gdg.dot(dDensity)[:,2]
    
    axGra.clear()
    axGra.plot(gravPointsX, dz)
    axGra.plot(gravPointsX, gravPointsX*0.0, 'v', color='black')
    axGra.set_ylabel('Grav in mGal')
    axGra.set_xlim((-20, 20))
    axGra.set_ylim((0, 0.003))
    axGra.grid()
    
    #print(time.time()-tic); tic = time.time()
    axReA.clear()
    ert.show(ertScheme, values=resA[i], axes=axReA, scheme='DipoleDipole',
             cMin=5, cMax=35, colorBar=0)
     
    if vis:
        pg.mplviewer.setMappableData(gciRes, 
                                     resistivities[i],
                                     cMin=1, cMax=100,
                                     logScale=True)
        pg.mplviewer.setMappableData(gciERT, 
                                     ertModels[i+1],
                                     cMin=1, cMax=100,
                                     logScale=True)
        pg.mplviewer.setMappableData(gciERR, 
                                     ertModels[i+1]/ertModels[0],
                                     cMin=1/6, cMax=6,
                                     logScale=True)
        
    print(i, time.time()-tic, 
          "sum:", sum(conc[i]),
          "dsum:", (sum(conc[i])-sum(conc[i-1])),
          )
    if mp4 or vis:
        plt.pause(0.001)

anim = animation.FuncAnimation(fig, animate,
                               frames=int(len(conc)),
                               interval=1)#, blit=True)
solutionName = createCacheName('all', mesh, times)+ "-" + str(ertScheme.size())

if mp4:
    anim.save(solutionName + ".mp4", writer=None, fps=20, dpi=dpi, codec=None,
          bitrate=24*1024, extra_args=None, metadata=None,
          extra_anim=None, savefig_kwargs=None)

pg.showNow()
