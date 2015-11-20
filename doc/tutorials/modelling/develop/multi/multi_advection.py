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
    injectPos=[0.21, -1.61]
    
    f = pg.RVector(mesh.cellCount(), 0.0)
    sourceCell=mesh.findCell(injectPos)
    
    f[sourceCell.id()] = scale*sourceCell.size()
    t = times[:len(times)/2]
    uMesh1 = pg.solver.solveFiniteVolume(mesh, a=1./peclet, f=f, vel=vel,
                                         times=t, 
                                         uB=[1, 0],
                                         scheme='PS',
                                         verbose=0)
    t = times[:len(times)/2+1]
    
    uMesh2 = pg.solver.solveFiniteVolume(mesh, a=1./peclet, f=0, vel=vel,
                                         times=t, 
                                         uB=[1, 0],
                                         u0=uMesh1[-1],
                                         scheme='PS',
                                         verbose=0)

    #saturation = uMesh1
    saturation = np.vstack((uMesh1, uMesh2[1:]))

    return saturation

if __name__ == '__main__':
    pass
pg.wait()

