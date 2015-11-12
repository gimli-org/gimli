#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pygimli as pg
import numpy as np

def createCacheName(base, mesh, timeSteps=[], dt=None, peclet=0):
    return 'cache-' + base + "-" + str(mesh.nodeCount()) \
            + '-' + str(len(timeSteps)) \
            + '-' + str(dt) \
            + '-' + str(peclet)


def calcSaturation(mesh, vel, times, injectPos, peclet=5, scale=1, cache=False, verbose=False):
    r"""
    .. math::
        
    """
    solutionName = createCacheName('saturation', mesh, times, times[1]-times[0],
                                   peclet)
    
    try:
        if cache:
            saturation = np.load(solutionName + '.bmat.npy')
        else:
            raise Exception("no cacheing")
        
    except Exception as e:
        if verbose:
            print(e)
            print("Building .... ")
        
        f = pg.RVector(mesh.cellCount(), 0.0)
        sourceCell=mesh.findCell(injectPos)
        
        f[sourceCell.id()] = scale*sourceCell.size()
        
        uMesh1 = pg.solver.solveFiniteVolume(mesh, a=1./peclet, f=f, vel=vel,
                                             times=times, 
                                             uB=[1, 0],
                                             scheme='PS',
                                             verbose=verbose*len(times)/10)
        #uMesh2 = solveFiniteVolume(mesh, a=1./peclet, f=0, vel=vel, times=times, 
                          #uBoundary=[2, 0], u0=uMesh1[-1],
                          #scheme='PS', verbose=10)
        saturation = uMesh1
        #saturation = np.vstack((uMesh1, uMesh2[1:]))
        np.save(solutionName + '.bmat', saturation)

    return saturation


if __name__ == '__main__':
    pass
pg.wait()

