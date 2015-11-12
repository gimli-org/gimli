#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pygimli as pg
import numpy as np

from multi_darcy_flow import darcyFlow
from multi_advection import calcSaturation
from multi_ert import simulateERTData


class WorkSpace():
    pass

class MulitFOP(pg.ModellingBase):
    def __init__(self, tMax=5000, tSteps=10, peclet=500, verbose=False):
        pg.ModellingBase.__init__(self, verbose=verbose)

        self.tSteps = tSteps
        self.regionManager().setParameterCount(3)
        self.setStartModel([1e-7, 0.01, 0.0001])
            
        self.ws = WorkSpace()
        self.timesAdvection = np.linspace(1, tMax, tSteps)
        self.peclet = peclet
            
    def response(self, par):
        ws = self.ws
        print('create response for: ' + str(par))
        model = par[0:3]
        print(".. darcy step")
        ws.mesh, ws.vel, ws.p, ws.k = darcyFlow(model=model, verbose=0)
  
        print(".. advection step")
        ws.saturation = calcSaturation(ws.mesh, ws.vel.T,
                                       self.timesAdvection, 
                                       injectPos=[2., -2.],
                                       peclet=self.peclet, verbose=0)
        
        print(".. dc step")
        sampleTime = [0, self.tSteps/2, self.tSteps-1]
        ws.meshERT, ws.dataERT, ws.res, ws.rhoa, ws.err = simulateERTData(ws.saturation[sampleTime], ws.mesh, verbose=1)

        return ws.rhoa.flatten()
    
    
    
def test(model, tMax=5000, tSteps=40, peclet=50000):
    fop = MulitFOP(verbose=1, tMax=tMax, tSteps=tSteps, peclet=peclet)
    
    ert = fop.response(model)

    ws = fop.ws
    mesh = ws.mesh
    ax, _ = pg.show(mesh)
    pg.show(mesh, fop.ws.p, label='Pressure in ??', axes=ax)
    ax, _ = pg.show(mesh)
    pg.show(mesh, fop.ws.k, label='Permeabilty in ??', axes=ax)
    
    ax, _ = pg.show(mesh)
    pg.show(mesh, np.sqrt(fop.ws.vel[0]**2+fop.ws.vel[1]**2), label='velocity in m/s', axes=ax)
    pg.show(mesh, fop.ws.vel, axes=ax)

    pg.show(mesh, fop.ws.saturation[-1]+1e-3, label='brine saturation')
    pg.show(ws.meshERT, fop.ws.res[-1], label='resistivity in Ohm m')

    pg.wait()
if __name__ == '__main__':
    
    #test(model=[1e-7, 0.01, 0.0001], tMax=5000, tSteps=40, peclet=500000)
    
    fop = MulitFOP(verbose=1, tMax=5000, tSteps=40, peclet=500000)
    ert = fop.response([1e-7, 0.01, 0.001])

    err = fop.ws.ertError.flatten()
        
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

