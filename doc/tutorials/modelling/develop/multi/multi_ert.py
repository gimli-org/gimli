#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pygimli as pg
import pybert as pb
import numpy as np

def createCacheName(base, mesh=None):
    nc = ''
    if mesh:
        nc = str(mesh.nodeCount())
    return 'cache-' + base + "-" + nc

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
     
    porosity = pg.solver.parseArgToArray(porosity, mesh.cellCount(), mesh)
    a = pg.solver.parseArgToArray(a, mesh.cellCount(), mesh)
    m = pg.solver.parseArgToArray(m, mesh.cellCount(), mesh)
    S = pg.solver.parseArgToArray(S, mesh.cellCount(), mesh)
    n = pg.solver.parseArgToArray(n, mesh.cellCount(), mesh)
    
    r = pg.RMatrix(len(rBrine), len(rBrine[0]))
    for i in range(len(r)):
        r[i] = rB[i] * a * porosity**(-m) * S**(-n)
            
    rI = pg.RMatrix(len(r), meshI.cellCount())
    if meshI:
        pg.interpolate(mesh, r, meshI.cellCenters(), rI) 
        
    for i in range(len(rI)):
        rI[i] = pg.solver.fillEmptyToCellArray(meshI, rI[i])
        
    return rI


def simulateERTData(saturation, meshSat, cache=False, verbose=0):
    swatch = pg.Stopwatch(True)
    if verbose:
        print("res 1:", swatch.duration(True))
    
    ertScheme = pb.DataContainerERT('20dd.shm')
    
    meshERT = pg.meshtools.createParaMesh(ertScheme, quality=34,
                                          paraMaxCellSize=0.1, 
                                          boundaryMaxCellSize=10)
    
    conductivityBack = 1./100.
    conductivityBrine = 100.
    conductivity = saturation * conductivityBrine + conductivityBack
    
    resis = resistivityArchie(rBrine=1./conductivity,
                              porosity=0.3, S=1.0, 
                              mesh=meshSat, meshI=meshERT)
    if verbose:
        print("res 2:", swatch.duration(True))
    #pg.show(meshERT, resis[-1], colorBar=1)

    ert = pb.manager.Resistivity(verbose=False)
    
    solutionName = createCacheName('appRes', meshERT) + \
            "-" + str(ertScheme.size()) + \
            "-" + str(len(saturation)) + \
            "-" + str(len(saturation[0]))   
    
    try:
        if cache:
            rhoa = np.load(solutionName + '.bmat.npy')
            ertData = pb.DataContainerERT(solutionName + '.dat')
        else:
            raise Exception("no cacheing")
    except Exception as e:
        if verbose:
            print(e)
            print("Building .... ")
        rhoa = np.zeros((len(resis), ertScheme.size()))
        err = np.zeros((len(resis), ertScheme.size()))
        ertScheme.set('k', pb.geometricFactor(ertScheme))
        
        ertData = ert.simulate(meshERT, resis[0], ertScheme)
        #ertData.estimateError(errPerc=1, errVolt=1e-5, verbose=1) would be nice??
        errPerc = 1
        errVolt = 1e-5
        voltage = ertData('rhoa') / ertData('k')
        ertData.set('err', pg.abs(errVolt / voltage) + errPerc / 100.0)
        if verbose:
            print('err min:', min(ertData('err'))*100, 'max:', max(ertData('err'))*100)
        ertData.save(solutionName + '.dat', 'a b m n rhoa err k')
        
        #sys.exit()
        for i in range(0, len(resis)):
            pg.tic()
            rhoa[i] = ert.fop.response(resis[i])
                        
            rand = pg.RVector(len(rhoa[i]))
            pg.randn(rand)
            err[i] = ertData('err')
            rhoa[i] *= (1.0 + rand * ertData('err'))
            
            if verbose:
                print(i, "/", len(resis), " : ", pg.dur(), "s",
                  "min r:", min(resis[i]), "max r:", max(resis[i]),
                  "min r_a:", min(rhoa[i]), "max r_a:", max(rhoa[i]) )

        np.save(solutionName + '.bmat', rhoa)
        
        
    return meshERT, ertData, resis, rhoa, err


if __name__ == '__main__':
    pass
pg.wait()

