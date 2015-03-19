#!/usr/bin/env python

"""
    Fluid flow stuff .. calc 
    
    - static velocity field
    - calc nonsteady advection/diffusion
"""

import matplotlib.pyplot as plt
import numpy as np

import pygimli as pg
from pygimli import physics
from pygimli.solver import parseArgToArray
from pygimli.solver import *

def createCacheName(base, mesh, timeSteps=[]):
    return 'cache-' + base + "-" + str(mesh.nodeCount()) + '-' + str(len(timeSteps))

def permeabiltyEngelhardtPitter(poro, q=3.5, s=5e-3,
                                mesh=None, meshI=None):
    r"""
    For sand and sandstones
    
    .. math:: 
        k & = 2\cdot 10^7 \frac{\phi^2}{(1-\phi)^2}* \frac{1}{S^2} \\
        S & = q\cdot s \\
        s & = \sum_{i=1}(\frac{P_i}{r_i})
        
    * :math:`\phi` - poro 0.0 --1.0
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
    poro = parseArgToArray(poro, mesh.cellCount(), mesh)
    q = parseArgToArray(q, mesh.cellCount(), mesh)
    s = parseArgToArray(s, mesh.cellCount(), mesh)
    
    S = q * s
    k = 2e7 * (poro**2 / (1.0-poro)**2) * 1.0/S**2 * physics.constants.Darcy
        
    if meshI:
        k = pg.interpolate(mesh, k, meshI.cellCenters()) 
        k = pg.solver.fillEmptyToCellArray(meshI, k)
    return k
        
def hydraulicConductivity(perm, visc=1.0, dens=1.0,
                          mesh=None, meshI=None):
    perm = parseArgToArray(perm, mesh.cellCount(), mesh)
    visc = parseArgToArray(visc, mesh.cellCount(), mesh)
    dens = parseArgToArray(dens, mesh.cellCount(), mesh)
    
    k = perm * dens/visc * pg.physics.constants.g
    
    if meshI:
        k = pg.interpolate(mesh, k, meshI.cellCenters()) 
        k = pg.solver.fillEmptyToCellArray(meshI, k)
    return k
    

def calcStokesVelocity(mesh, visc, velBoundary, preBoundary):
   
    solutionName = createCacheName('vel', mesh)
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
                                                           viscosity=visc,
                                                           maxIter=1000,
                                                           tol=1e-4, pRelax=0.1,
                                                           verbose=1, ws=ws)
        np.save(solutionName + '.bmat', vel)
       
    return vel

def calcConcentration(mesh, vel, times, injectPos, peclet=50, scale=1):
    r"""
    .. math::
        
    """
    solutionName = createCacheName('conc', mesh, times)
    
    try:
        conc = np.load(solutionName + '.bmat.npy')
        
    except Exception as e:
        print(e)
        print("Building .... ")
        
        f = pg.RVector(mesh.cellCount(), 0.0)
        sourceCell=mesh.findCell(injectPos)
        
        f[sourceCell.id()] = scale
        
        print(sourceCell.size())

        uMesh1 = solveFiniteVolume(mesh, a=1./peclet, f=f, vel=vel, times=times, 
                          uBoundary=[2, 0],
                          scheme='PS', verbose=10)
        uMesh2 = solveFiniteVolume(mesh, a=1./peclet, f=0, vel=vel, times=times, 
                          uBoundary=[2, 0], u0=uMesh1[-1],
                          scheme='PS', verbose=10)
        conc = np.vstack((uMesh1, uMesh2[1:]))
        np.save(solutionName + '.bmat', conc)

    return conc


if __name__ == "__main__":
    pass