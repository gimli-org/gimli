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

from pygimli.physics.seismics import ricker, solvePressureWave, drawSeismogramm
from pygimli.physics.gravimetry import solveGravimetry

import pybert as pb
import pybert.dataview

def createCacheName(base, mesh, times):
    return 'cache-' + base + '-' + str(mesh.nodeCount()) + '-' + str(len(times))


class MethodManager(object):
    """
        General manager to maintenance a measurement method. 
        
        The method manager holds one instance of a forward operator and a 
        appropriate inversion method to handle simulation and reconstruction of
        common geophysical problems.
        
    """
    def __init__(self, verbose=False):
        self.fop = self.createFOP(verbose)
        self.tD = None
        self.tM = None
        self.inv = self.createInv(verbose)
    
    def __str__(self):
        return self.__repr__()
    
    def __repr__(self):
        """ String representation of the class """
        out = "Abstract MethodManager."
        pass
    
    def createFOP(self, refine=True): 
        """ Create forward operator working on refined mesh """
        raise Exception("Overload me!")

    def createInv(self, verbose=True, dosave=False):
        """ Create inversion instance, data- and model transformations. """
        raise Exception("Overload me!")
    
    # Visualization stuff
    def show(self, data, values=None, axes=None,
             cMin=None, cMax=None, colorBar=1, **kwargs):
        """ Forward the visualization """
        pass
        
    def showData(self, ax=None, response=None):
        """ show data in form of travel time curves """
        pass

    def showMesh(self, ax=None):
        """ show mesh in given axes or in a new figure """
        pass

    def showResult(self, ax=None, cMin=None, cMax=None, logScale=False,
                   **kwargs):
        """ show resulting velocity vector """
        pass

    # Mesh related methods
    def createMesh(self, ):
        """ Create a mesh aka the parametrization """
        pass
    
    def setParaMesh(self, mesh):
        """ Set mesh for the inverse problem """
        pass
    
    # Data related methods
    def createData(self, sensors, scheme):   
        """ Create an empty data set """
        pass
    
    def setData(self, data):
        """ Set data """
        pass

    def importData(self, filename):
        """ Import data """
        pass
        
    def checkData(self):
        """ Check data validity """
        pass

    def estimateError(self, absoluteError=0.001, relativeError=0.001):
        """ estimate error composed of an absolute and a relative part """
        pass
    
    # Work related methods
    def invert(self, data, values=None, verbose=0, **kwargs):
        """ Invert the data and fill the parametrization. """ 
        pass
             
    def simulate(self, mesh, values, data=None):
        """ Run a simulation aka the forward task. """ 
        pass
    
    def saveResult(self, folder=None, size=(16, 10),
                   **kwargs):
        """
        Saves the results in the specified folder.
        """
        pass

class ERT(MethodManager):
    
    """
        ERT Manager. 
        
        Should solve most of the common problems
    """
        
    def __init__(self, verbose=False):
        """
        """
        super(type(self), self).__init__(verbose)
                
        self.schemeMg = pb.dataview.dataview.DataSchemeManager()
        self.paraMesh = None
         
    def createFOP(self, verbose):
        """ Create resistivity modelling forward operator. """
        return pb.DCSRMultiElectrodeModelling(verbose=verbose)
         
    def createInv(self, verbose):
        """ Create resistivity inversion instance. """
        self.tD = pg.RTransLog()
        self.tM = pg.RTransLog()
        inv = pg.RInversion(verbose=verbose, dosave=False)
        inv.setTransData(self.tD)
        inv.setTransModel(self.tM)
        return inv
        
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
    
class Gravimetry(MethodManager):
    """
        General Gravimetry Method Manager
    """
    
    def __init__(self, verbose=False):
        """Default constructor."""
        super(type(Gravimetry), self).__init__(verbose)
        

def createTestWorld1(maxArea=0.2, verbose=0):
    boundary = []
    boundary.append([-20.0,   0.0]) #0
    boundary.append([-20.0,  -2.0]) #1
    boundary.append([-20.0,  -8.0]) #2
    boundary.append([-20.0,  -15.0]) #3
    boundary.append([ 20.0,  -15.0]) #4
    boundary.append([ 20.0,  -8.0]) #5
    boundary.append([ 20.0,  -2.0]) #6
    boundary.append([ 20.0,   0.0]) #7
       
    boundary = boundary + [pg.RVector3(x,0) for x in np.arange(19, -19.1, -1)]   
    boundaryMarker = [1, 2, 3, 4, 5, 6, 7, 8] + [8 for x in range(40)]
    
    poly = pg.Mesh(2)
    nodes = [poly.createNode(b) for b in boundary]
    
    for i in range(len(boundary)):
        poly.createEdge(nodes[i], nodes[(i+1)%len(boundary)], boundaryMarker[i]) # dirichlet (inflow)
    
    poly.createEdge(nodes[1], nodes[6], 9)
    poly.createEdge(nodes[2], nodes[5], 10)
        
    boundary = []
    boundary.append([-6.0, -3.5]) #0
    boundary.append([-6.0, -6.0]) #1
    boundary.append([ 6.0, -6.0]) #2
    boundary.append([ 6.0, -3.5]) #3
    
    nodes = [poly.createNode(b) for b in boundary]
    for i in range(len(boundary)):
        poly.createEdge(nodes[i], nodes[(i+1)%len(boundary)], 11) # dirichlet (inflow)
  
    poly.addRegionMarker((-19, -1),  1)
    poly.addRegionMarker((-19, -3), 2)
    poly.addRegionMarker((-19, -14), 3)
    poly.addRegionMarker((  0, -5),  4, 0.1)
    
    mesh = createMesh(poly, quality=33, area=maxArea,
                      smooth=[1,10], verbose=verbose)
    
    velBoundary=[[2, [1.0, 0.0]],
                 [6, [1.0, 0.0]],
                 [3, [1./3., 0.0]],
                 [8, [0.0, 0.0]],
                 [5, [1./3, 0.0]],
                 [4, [1./3, 0.0]]
                ]
    preBoundary=[[1, 0.0],
                 #[4, 0.0], // bottom
                 #[5, 0.0],
                 [7, 0.0]]
    
    porosity = pg.solver.parseMapToCellArray([[1, 0.3], 
                                              [2, 0.4], 
                                              [3, 0.2],
                                              [4, 0.115]], mesh)
    # 30 works nice with tol=1e-2 and pressureRelaxation = 0.01
    
    return mesh, velBoundary, preBoundary, porosity

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
    poly.addRegionMarker((0.0, -5.0), 6, 0.01)
    
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

    print("div(V) = ", pg.solver.divergence(mesh, normMap=velBoundary))
    
    return mesh, velBoundary, preBoundary, viscosity

        
def density(porosity, rhoMatrix=2510, rhoFluid=1000, S=1,
            mesh=None):
    r"""
    """
    porosity = parseArgToArray(porosity, mesh.cellCount(), mesh)
    rhoMatrix = parseArgToArray(rhoMatrix, mesh.cellCount(), mesh)
    rhoFluid = parseArgToArray(rhoFluid, mesh.cellCount(), mesh)
    S = parseArgToArray(S, mesh.cellCount(), mesh)
    
    dens = np.array(rhoMatrix * (1.-porosity)) + rhoFluid * porosity * S
    return dens
 
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

def permeabiltyEngelhardtPitter(porosity, q=3.5, s=5e-3,
                                mesh=None, meshI=None):
    r"""
    For sand and sandstones
    
    .. math:: 
        k & = 2\cdot 10^7 \frac{\phi^2}{(1-\phi)^2}* \frac{1}{S^2} \\
        S & = q\cdot s \\
        s & = \sum_{i=1}(\frac{P_i}{r_i})
        
    * :math:`\phi` - porosity 0.0 --1.0
    * :math:`S` - in cm^-1  specific surface in cm^2/cm^3
    * :math:`q` - (3 for spheres, > 3 shape differ from sphere)
        3.5 sand
    * :math:`s` - in cm^-1 (s = 1/r for particles with homogeneous radii r)
    * :math:`P_i` - Particle ration with radii :math:`r_i` on 1cm^3 Sample
    
    Returns
    -------
    k :
        in Darcy
    """
    porosity = parseArgToArray(porosity, mesh.cellCount(), mesh)
    q = parseArgToArray(q, mesh.cellCount(), mesh)
    s = parseArgToArray(s, mesh.cellCount(), mesh)
    
    S = q * s
    k = 2e-7 * (porosity**2 / (1.0-porosity)**2) * 1.0/ S**2
    
    if meshI:
        k = pg.interpolate(mesh, k, meshI.cellCenters()) 
        k = pg.solver.fillEmptyToCellArray(meshI, k)
    return k
        
        
def resistivityArchie(rBrine, porosity, a=1.0, m=2.0, S=1.0, n=2.0,
                      mesh=None, meshI=None):
    """
    .. math:: 
        \rho = a\rho_{\text{Brine}}\phi^{-m}\S_w^{-n}
        
    * :math:`\rho` - the electrical conductivity of the fluid saturated rock
    * :math:`\rho_{\text{Brine}}` - electrical conductivity of the brine
    * :math:`\phi` - porosity 0.0 --1.0
    * :math:`a` - tortuosity factor. (common 1)
    * :math:`m` - cementation exponent of the rock
            (usually in the range 1.3 -- 2.5 for sandstones)
    * :math:`n` - is the saturation exponent (usually close to 2)
     
    """
    rB = None
    
    if rBrine.ndim == 1:
        rB = pg.RMatrix(1, len(rBrine))
        rB[0] = parseArgToArray(rBrine, mesh.cellCount(), mesh)
    elif rBrine.ndim == 2:
        rB = pg.RMatrix(len(rBrine), len(rBrine[0]))
        for i in range(len(rBrine)):
            rB[i] = rBrine[i]
     
    porosity = parseArgToArray(porosity, mesh.cellCount(), mesh)
    a = parseArgToArray(a, mesh.cellCount(), mesh)
    m = parseArgToArray(m, mesh.cellCount(), mesh)
    S = parseArgToArray(S, mesh.cellCount(), mesh)
    n = parseArgToArray(n, mesh.cellCount(), mesh)
    
    r = pg.RMatrix(len(rBrine), len(rBrine[0]))
    for i in range(len(r)):
        r[i] = rB[i] * a * porosity**(-m) * S**(-n)
            
    rI = pg.RMatrix(len(r), meshI.cellCount())
    if meshI:
        pg.interpolate(mesh, r, meshI.cellCenters(), rI) 
        
    for i in range(len(rI)):
        rI[i] = pg.solver.fillEmptyToCellArray(meshI, rI[i])
        
    return rI
    
def calcVelocity(mesh, permeabilty, velBoundary, preBoundary):

    viscosity = 1./permeabilty
    solutionName = createCacheName('vel', mesh, times)
    try:
        #vel = pg.load(solutionName + '.bmat')
        vel = np.load(solutionName + '.bmat.npy')
        
    except Exception as e:
        print(e)
        print("Building .... ")
        class WS:
            pass
        ws = WS
        vel, pres, pCNorm, divVNorm = solveStokes_NEEDNAME(mesh, velBoundary,
                                                           preBoundary,
                                                           viscosity=viscosity,
                                                           maxIter=1000,
                                                           tol=1e-4, pRelax=0.1,
                                                           verbose=1, ws=ws)
        np.save(solutionName + '.bmat', vel)
    
        #plt.semilogy(divVNorm)
        #pg.show(mesh, pg.cellDataToPointData(mesh, pres), cmap='b2r')
        #data=pg.cellDataToPointData(mesh,
                                    #np.sqrt(vel[:,0]*vel[:,0] + vel[:,1]*vel[:,1]))
        
        #ax, cbar = pg.show(mesh, ws.div, colorBar=1, label='divergence')
        #ax, cbar = pg.show(mesh, data, colorBar=1, label='Velocity')
        #pg.show(mesh, vel, axes=ax)
        #pg.showNow()
        
    return vel

def calcConcentration(mesh, vel, times):
    r"""
    .. math::
        
    """
    solutionName = createCacheName('conc', mesh, times)
    
    try:
        #vel = pg.load(solutionName + '.bmat')
        conc = np.load(solutionName + '.bmat.npy')
        
    except Exception as e:
        print(e)
        print("Building .... ")
        
        f = pg.RVector(mesh.cellCount(), 0.0)
        sourceCell=mesh.findCell([-18., -6])
        f[sourceCell.id()]=2.0
        print(sourceCell.size())

        Peclet = 50.0


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
        #Gdg, Gdgz = solveGravimetry(mesh, None, pnts=gravPoints, complete=True) 
        Gdg = solveGravimetry(mesh, None, pnts=gravPoints) 
        np.save(solutionName + '.bmat', Gdg)
        
    return gravPointsX, Gdg

def calcApparentResistivities(mesh, resistivities):
    ert = ERT(verbose=False)
    ertPointsX = [pg.RVector3(x,0) for x in np.arange(-19, 19.1, 1)]
    ertScheme = ert.createData(ertPointsX, scheme="Dipole Dipole (CC-PP)")
    
    solutionName = createCacheName('appRes', mesh, times)+ "-" + str(ertScheme.size())
    
    try:
        rhoa = np.load(solutionName + '.bmat.npy')
        ertData = pb.DataContainerERT(solutionName + '.dat')
    except Exception as e:
        print(e)
        print("Building .... ")
        rhoa = np.zeros((len(resistivities), ertScheme.size()))
        ertScheme.set('k', pb.geometricFactor(ertScheme))
        ertData = ert.simulate(meshERT_FOP, resistivities[0], ertScheme)
        
        errPerc = 1
        errVolt = 1e-5
        voltage = ertData('rhoa') / ertData('k')
        ertData.set('err', pg.abs(errVolt / voltage) + errPerc / 100.0)
        print('err min:', min(ertData('err'))*100, 'max:', max(ertData('err'))*100)
        ertData.save(solutionName + '.dat', 'a b m n rhoa err k')
        #sys.exit()
        for i in range(0, len(resistivities)):
            tic =  time.time()
            rhoa[i] = ert.fop.response(resistivities[i])
                        
            rand = pg.RVector(len(rhoa[i]))
            pg.randn(rand)
        
            rhoa[i] *= (1.0 + rand * ertData('err'))
            
            print(i, "/", len(resistivities), " : ", time.time()-tic, "s",
                  "min:", min(resistivities[i]), "max:", max(resistivities[i]),
                  "min:", min(rhoa[i]), "max:", max(rhoa[i]) )

        np.save(solutionName + '.bmat', rhoa)
        
        
    return rhoa, ert, ertData

def calcERT(ert, ertScheme, rhoa):
    solutionName = createCacheName('ERT', mesh, times)+ "-" + str(ertScheme.size())
    try:   
        ertModels = pg.load(solutionName + '.bmat')
        ertMesh = pg.load(solutionName + '.bms')
    except Exception as e:
        print(e)
        print("Building .... ")
        ertModels = ert.invert(ertScheme, values=rhoa, maxiter=10, 
                               paraDX=0.5, paraDZ= 0.5, nLayers=20, paraDepth=20,
                               lambd=50, verbose=1)
        ertMesh=ert.fop.regionManager().paraDomain()
        ertModels.save(solutionName + '.bmat')
        ertMesh.save(solutionName)
    return ertModels, ertMesh

def calcSeismics(mesh, vP):
    
    meshSeis = appendTriangleBoundary(mesh, 
                                      xbound=50, ybound=50, marker=1,
                                      quality=33.0, area=0.5, smooth=True, 
                                      markerBoundary=1,
                                      isSubSurface=False, verbose=False)
    print(meshSeis)
    meshSeis = meshSeis.createP2()
    vP = pg.interpolate(mesh, vP, meshSeis.cellCenters())
    #mesh = meshSeis.createH2()
    mesh = meshSeis
    vP = pg.solver.fillEmptyToCellArray(mesh, vP)
    print(mesh)
    ax, cbar = pg.show(mesh, data=vP)
    pg.show(mesh, axes=ax)
    pg.showNow()
    
    h = pg.median(mesh.boundarySizes())
    dt = 0.5 * h /max(vP)
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
        f0 = 1./dt*0.1
        
        print("h:", round(h,2), "dt:", round(dt,5), "1/dt:", round(1/dt,1), "f0", round(f0,2))
        
        uSource = ricker(f0, times, t0=1./f0)
    
        plt.plot(times, uSource)
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
    
    
    
vis = 1
mp4 = 1

swatch = pg.Stopwatch(True)
pg.showLater(1)

modelFkt = createTestWorld1

mesh, velBoundary, preBoundary, porosity = modelFkt(maxArea=0.2,
                                                    verbose=0)
print("meshgen:", swatch.duration(True))

times = np.linspace(0, 50, 51)
densityBrine = 1.2 * 1000 # kg/m^3
 
permeabilty = permeabiltyEngelhardtPitter(porosity, mesh=mesh)    

vP = velocityVp(porosity, mesh=mesh)

calcSeismics(mesh, vP)
print("Seismics:", swatch.duration(True))

vel = calcVelocity(mesh, permeabilty, velBoundary, preBoundary)
print("vel:", swatch.duration(True))

conc = calcConcentration(mesh, vel, times)
print("con:", swatch.duration(True))

dens0 = density(porosity, rhoMatrix=2510, rhoFluid=1000, mesh=mesh)
densC = density(porosity, rhoMatrix=2510, 
                rhoFluid=1000 *(1.-conc) + densityBrine * conc, mesh=mesh)
dDens = densC - dens0


gravPointsX, Gdg = calcGravKernel(mesh)
print("grav:", swatch.duration(True))



meshD, tmp, tmp, tmp=modelFkt(maxArea=0.1, verbose=0)
meshERT_FOP = appendTriangleBoundary(meshD, 
                                    xbound=50, ybound=50, marker=1,
                                    quality=34.0, smooth=False, markerBoundary=1,
                                    isSubSurface=False, verbose=False)
resistivities = resistivityArchie(rBrine=1./(1./20. + abs(1.*conc)),
                                  porosity=porosity, S=1, 
                                  mesh=mesh, meshI=meshERT_FOP)
print("res:", swatch.duration(True))

rhoa, ert, ertScheme = calcApparentResistivities(meshERT_FOP, resistivities)
print("rhoa:", swatch.duration(True))

ertModels, meshERT = calcERT(ert, ertScheme, rhoa)
print("ert:", swatch.duration(True))

vP = velocityVp(porosity, mesh=mesh)
print("vp:", swatch.duration(True))




dpi=92
fig = None

orientation = 'horizontal'

if vis:

    if mp4:
        fig = plt.figure(facecolor='white', figsize=(2*800/dpi, 2*490/dpi), dpi=dpi)  
    else:
        fig = plt.figure() 
    
    axPor = fig.add_subplot(4,4,1)
    axPer = fig.add_subplot(4,4,2)
    axDen = fig.add_subplot(4,4,3)
    axVp  = fig.add_subplot(4,4,4)
    
    axVel = fig.add_subplot(4,4,5)
    axCon = fig.add_subplot(4,4,6)
    axDDe = fig.add_subplot(4,4,7)
    axRes = fig.add_subplot(4,4,8)
    
    axGra = fig.add_subplot(4,4,11)
    axReA = fig.add_subplot(4,4,12)
    
    axERT = fig.add_subplot(4,4,15)
    axERR = fig.add_subplot(4,4,16)
        
    # ** Porosity **
    show(mesh, data=porosity, colorBar=1, 
         orientation=orientation, label='Porosity', axes=axPor)
    show(mesh, axes=axPor)
    
    # ** Permeabilty **
    show(mesh, data=permeabilty, colorBar=1, 
         orientation=orientation, label='Permeabilty', axes=axPer)
    show(mesh, axes=axPer)
    
    # ** Density **
    show(mesh, data=dens0, colorBar=1, 
         orientation=orientation, label='Density in kg/m$^3$', axes=axDen)
    show(mesh, axes=axPer)
    
    # ** Velocity abs **
    axVel, cbar = show(mesh, data=pg.cellDataToPointData(mesh,
                                np.sqrt(vel[:,0]*vel[:,0] + vel[:,1]*vel[:,1])),
                       logScale=0, colorBar=1,
                       orientation=orientation, label='|Velocity| in m/s',
                       axes=axVel
                       )

    # ** Velocity vector **
    pg.mplviewer.drawMeshBoundaries(axVel, mesh, fitView=True)
    pg.viewer.showBoundaryNorm(mesh, velBoundary, color='red', axes=axVel)
    
    meshC, tmp, tmp, tmp=modelFkt(maxArea=1, verbose=0)
    show(mesh, data=vel, axes=axVel, coarseMesh=meshC)
    
    # ** vP **
    show(mesh, data=vP, colorBar=1, 
         orientation=orientation, label='Vp m/s', axes=axVp)
    show(mesh, axes=axVp)
    
    # Prepare time lapse figures
    # ** Concentration **
    gciCon= pg.mplviewer.drawModel(axCon, mesh, data=conc[1],
                                   cMin=0, cMax=0.1, 
                                   logScale=False
                                   #cmap='b2r'
                                   )
    cbar = createColorbar(gciCon, orientation=orientation, label='Concentration')
    
    gciDDe = pg.mplviewer.drawModel(axDDe, mesh, data=dDens[1],
                                    cMin=0, cMax=20, 
                                    #cmap='b2r'
                                   )
    cbar = createColorbar(gciDDe, orientation=orientation, label='Delta density in kg/m$^3$')
    
    # ** Resistivity (model) **
    gciRes = pg.mplviewer.drawModel(axRes, meshERT_FOP, 
                                    data=resistivities[0],
                                    cMin=20, cMax=700,
                                    )
    cbar = createColorbar(gciRes, orientation=orientation, label='Resistivity')
    axRes.set_xlim((-20, 20))
    axRes.set_ylim((-14, 00))
    
    # ** Apparent resistivities (data) **
    gciARes = ert.show(ertScheme, values=rhoa[0], axes=axReA, 
                       scheme='DipoleDipole',
                       #cMin=100, cMax=300, 
                       orientation=orientation)

    # ** ERT (model) **
    gciERT = pg.mplviewer.drawModel(axERT, meshERT, 
                                    data=ertModels[0],
                                    cMin=20, cMax=700)
    cbar = createColorbar(gciERT, orientation=orientation, label='Resistivity')
    # ** ERT ratio (model) **
    gciERR = pg.mplviewer.drawModel(axERR, meshERT, 
                                    data=ertModels[0]/ertModels[0],
                                    cMin=1/4, cMax=4, cmap='b2r')
    cbar = createColorbar(gciERR, orientation=orientation, label='Ratio')


def animate(i):
    tic = time.time()
        
    axGra.clear()
    dz = Gdg.dot(dDens[i])
    
    axGra.clear()
    axGra.plot(gravPointsX, dz)
    axGra.plot(gravPointsX, gravPointsX*0.0, 'v', color='black')
    axGra.set_ylabel('Grav in mGal')
    axGra.set_xlim((-20, 20))
    axGra.set_ylim((0, 0.003))
    axGra.grid()
    
    axReA.clear()
    ert.show(ertScheme, values=rhoa[i], axes=axReA, scheme='DipoleDipole',
             #cMin=100, cMax=300, 
             colorBar=0)
     
    if vis:
        
        pg.mplviewer.setMappableData(gciCon, conc[i], 
                                     #cMin=0, cMax=0.03,
                                     logScale=False)
        gciCon.set_clim(0, 0.1)
        pg.mplviewer.setMappableData(gciDDe, dDens[i], 
                                     cMin=0, cMax=20,
                                     logScale=False)
        pg.mplviewer.setMappableData(gciRes, 
                                     resistivities[i],
                                     cMin=20, cMax=700,
                                     logScale=True)
        pg.mplviewer.setMappableData(gciERT, 
                                     ertModels[i+1],
                                     cMin=20, cMax=700,
                                     logScale=True)
        pg.mplviewer.setMappableData(gciERR, 
                                     ertModels[i+1]/ertModels[0],
                                     cMin=1/4, cMax=4,
                                     logScale=True)
        
    print(i, round(time.time()-tic, 2),
          "t=", round(i *(times[1]-times[0]),1),
          "dt:", round(times[1]-times[0],1),
          "sum mass:", round(sum(dDens[i]*mesh.cellSizes()),1), "kg/m "
          "sum:", sum(conc[i]),
          "dsum:", (sum(conc[i])-sum(conc[i-1])),
          )
    if mp4 or vis:
        pass
        plt.pause(0.001)

#animate(50)
plt.show()
for i in range(1, len(times)*2-1):
    animate(i)

#anim = animation.FuncAnimation(fig, animate,
                               #frames=int(len(conc)),
                               #interval=1)#, blit=True)

#solutionName = createCacheName('all', mesh, times)+ "-" + str(ertScheme.size())

#if mp4:
    #anim.save(solutionName + ".mp4", writer=None, fps=20, dpi=dpi, codec=None,
          #bitrate=24*1024, extra_args=None, metadata=None,
          #extra_anim=None, savefig_kwargs=None)

pg.showNow()
