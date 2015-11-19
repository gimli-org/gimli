#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pygimli as pg
import numpy as np

def createCacheName(base, mesh, timeSteps=[], dt=None, peclet=0):
    return 'cache-' + base + "-" + str(mesh.nodeCount()) \
            + '-' + str(len(timeSteps)) \
            + '-' + str(dt) \
            + '-' + str(peclet)


def calcSaturation(mesh, vel, times, peclet=5, scale=1, cache=False, verbose=False):
    r"""
    .. math::
        
    """
    injectPos=[1.01, -1.31]
    
    f = pg.RVector(mesh.cellCount(), 0.0)
    sourceCell=mesh.findCell(injectPos)
    
    f[sourceCell.id()] = scale*sourceCell.size()
    
    
    uMesh1 = pg.solver.solveFiniteVolume(mesh, a=1./peclet, f=f, vel=vel,
                                         times=times, 
                                         uB=[1, 0],
                                         scheme='PS',
                                         verbose=0)#verbose*len(times)/10)
    #uMesh2 = solveFiniteVolume(mesh, a=1./peclet, f=0, vel=vel, times=times, 
                      #uBoundary=[2, 0], u0=uMesh1[-1],
                      #scheme='PS', verbose=10)
    saturation = uMesh1
    #saturation = np.vstack((uMesh1, uMesh2[1:]))

    return saturation

if __name__ == '__main__':
    pass
pg.wait()

