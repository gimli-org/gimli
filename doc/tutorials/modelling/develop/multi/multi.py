#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pygimli as pg
import numpy as np

from multi_darcy_flow import darcyFlow
from multi_advection import calcSaturation
from multi_ert import simulateERTData


class MulitFOP(pg.ModellingBase):
    def __init__(self, verbose):
        pg.ModellingBase.__init__(self, verbose)

        self.regionManager().setParameterCount(3)
        self.setStartModel([0.001, 0.001, 0.001])
                
        self.dataERT = None
        self.timesAdvection = np.linspace(1, 500, 10)
        
    def response(self, par):
        print('create response for: ' + str(par))
        model = par[0:3]
        print(".. darcy step")
        mesh, vel, p, k = darcyFlow(model=model, verbose=0)
  
        print(".. advection step")
        saturation = calcSaturation(mesh, vel.T, self.timesAdvection, 
                                    injectPos=[2., -2.], peclet=5000, verbose=0)
        
        print(".. dc step")
        meshERT, dataERT, r, rhoa, err = simulateERTData(saturation[::4], mesh, verbose=0)
        self.dataERT = dataERT
        self.ertError = err

        return rhoa.flatten()
    
    
    
def test(model):
    
    mesh, vel, p, k = darcyFlow(model=model)

    times = np.linspace(1, 500, 10)
   
    saturation = calcSaturation(mesh, vel.T, times, 
                                injectPos=[2., -2.], peclet=5000)
    
    meshERT, dataERT, resistivities, rhoa, err = simulateERTData(saturation[::4], mesh)
        
        
    ax, _ = pg.show(mesh)
    pg.show(mesh, p, label='Pressure in ??', axes=ax)
    ax, _ = pg.show(mesh)
    pg.show(mesh, k, label='Permeabilty in ??', axes=ax)
    
    ax, _ = pg.show(mesh)
    pg.show(mesh, np.sqrt(vel[0]**2+vel[1]**2), label='velocity in m/s', axes=ax)
    pg.show(mesh, vel, axes=ax)

    pg.show(mesh, saturation[-1]+1e-3, label='brine saturation')
    pg.show(meshERT, resistivities[-1], label='resistivity in Ohm m')


if __name__ == '__main__':
    
    #test([0.0001, 0.1, 0.001])
    
    fop = MulitFOP(verbose=1)
    ert = fop.response([0.0001, 0.1, 0.001])

    err = fop.ertError.flatten()
        
    inv = pg.RInversion(ert, fop, verbose=1, dosave=1)
        
    tD = pg.RTransLog()
    tM = pg.RTransLog()
    inv.setTransData(tD)
    inv.setTransModel(tM)
    
    inv.setRelativeError(err)
    inv.setLambda(0)
    coeff = inv.run()
    print(coeff)

# actual inversion run yielding coefficient model



pg.wait()

