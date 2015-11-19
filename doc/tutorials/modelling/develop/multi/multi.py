#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pygimli as pg
import pybert as pb
import numpy as np

from multi_darcy_flow import darcyFlow
from multi_advection import calcSaturation
from multi_ert import simulateERTData, showERTData

from multiprocessing import Process, JoinableQueue
import copy

from math import ceil
import time
class WorkSpace():
    pass

def responseProc(model, peclet, timesAdvection, tSteps, i, j, outQueue=None, ws=None):
    
    if outQueue is not None: print("proc:", i)
            
    if ws is None:
        ws = WorkSpace()

    if outQueue is None: print(".. darcy step")

    ws.mesh, ws.vel, ws.p, ws.k = darcyFlow(model, verbose=0)

    #pg.show(model, model.cellAttributes(), label='cell cellAttributes'); pg.wait()
    #pg.show(ws.mesh, ws.vel.T, label='velocity'); pg.wait()
    #pg.show(ws.mesh, ws.k, label='Permeabilty'); pg.wait()
    
    if outQueue is None: print(".. saturation step")        
        
    ws.sat = calcSaturation(ws.mesh, ws.vel.T, 
                            timesAdvection, 
                            peclet=peclet, 
                            verbose=0)
    iProc = i
    
    if outQueue is None: print(".. ert step")
        
        
    sampleTime = [0, tSteps/4, tSteps/2, 3*tSteps/4, tSteps-1]
    #sampleTime = [0, tSteps-1]
    ws.meshERT, ws.dataERT, ws.res, ws.rhoa, ws.err = simulateERTData(ws.sat[sampleTime],
                                                                      ws.mesh, i=iProc, j=j, verbose=0)

    if outQueue:
        #print(i, ws.rhoa.shape, np.min(ws.rhoa.flatten()), np.max(ws.rhoa.flatten()))
        outQueue.put([i, ws.rhoa])
        
    else:
        return ws.rhoa
    
class MulitFOP(pg.ModellingBase):
    def __init__(self, verbose=False, **kwargs):

        pg.ModellingBase.__init__(self, verbose=verbose)
        self.init(kwargs.pop('mesh', None), 
                  kwargs.pop('tMax', 50000),
                  kwargs.pop('tSteps', 10),
                  kwargs.pop('peclet', 50000))
        self.iter = 0
        
    def init(self, mesh, tMax, tSteps, peclet):
    
        self.parMesh = pg.Mesh(mesh)
        self.setMesh(mesh)
        self.createRefinedForwardMesh(refine=False, pRefine=False)
        
        self.tMax = tMax
        self.tSteps = tSteps
        self.timesAdvection = np.linspace(1, tMax, tSteps)
        self.peclet = peclet
        self._J = pg.RMatrix()
        self.setJacobian(self._J)   
        self.ws = WorkSpace()
        
    def createJacobian(self, model):
        print("Create Jacobian for", str(model))
        nModel = len(model)
        resp = self.response(model)
        nData = len(resp)
        
        self._J.resize(nData, nModel)
        self._J *= 0.
        
        dModel = pg.RVector(len(model))
        rhoaJ = np.zeros((len(model), len(resp)))

        fak = 1.05
        self.iter += 1

        pg.tic()
        output = JoinableQueue()
        
        nProcs = float(pg.numberOfCPU())
                        
        for pCount in range(int(np.ceil(nModel/nProcs))):
            procs = []
            print(pCount*nProcs, "/" ,nModel)
            for i in range(int(pCount*nProcs), int((pCount+1)*nProcs)):

                if i < nModel:
                    modelChange = pg.RVector(model)
                    modelChange[i] *= fak
                    dModel[i] = modelChange[i]-model[i]
                
                    self.mapModel(modelChange)
                    modelMesh = pg.Mesh(self.mesh())
            
                    #rhoaJ[i] = responseProc(modelMesh, self.peclet, self.timesAdvection, self.tSteps, i, self.iter).flatten()
            
                    procs.append(Process(target=responseProc, args=(modelMesh, 
                                                                    self.peclet, self.timesAdvection, self.tSteps,
                                                                    i, self.iter, output)))

            for i, p in enumerate(procs):
                p.start()
                
            for i, p in enumerate(procs):
                p.join()
                iModel, rhoa = output.get()
                #print(iModel, rhoa.shape, np.min(rhoa.flatten()), np.max(rhoa.flatten()))
                rhoaJ[iModel] = rhoa.flatten()
        pg.toc()

        self._J *= 0.0
        for i in range(nModel):
            dData = rhoaJ[i] - resp
            self._J.setCol(i, dData/dModel[i])
            
        #print('#'*40 + 'Jac:')
        #for i in range(nModel):
            #print(i, sum(pg.abs(self._J.col(i))))
        
        #exit()
    def response(self, par):
        model = par
        if self.mesh():
            if len(par) == self.mesh().cellCount():
                self.mesh().setCellAttributes(par)
            else:
                self.mapModel(par)
        model = self.mesh()
        
        rhoa = responseProc(model, self.peclet, self.timesAdvection, self.tSteps,
                     0, 0, None, self.ws)
                
        return rhoa.flatten()   
        
    
def simulateSynth(tMax=5000, tSteps=40, peclet=50000, show=False, load=False):
    
    paraMesh = pg.createGrid(x=[0, 10], y=[-5, -3.5, -0.5, 0])
    paraMesh.cell(0).setMarker(0) # bottom
    paraMesh.cell(1).setMarker(1) # center
    paraMesh.cell(2).setMarker(2) # top
    
    paraMesh = paraMesh.createH2()
    paraMesh = paraMesh.createH2()
    paraMesh = paraMesh.createH2()
    paraMesh = paraMesh.createH2()
    paraMesh = paraMesh.createH2()
    
    for c in paraMesh.cells():
        if c.center()[0] > 3 and c.center()[0] < 7 and \
            c.center()[1] > -2.5 and c.center()[1] < -1.0:
            c.setMarker(3)
    
    model=[0.001, 0.01, 1e-7, 0.001]
    
    #pg.show(paraMesh, np.array(model)[paraMesh.cellMarker()], label='Permeabilty Model')
        
    fop = MulitFOP(mesh=paraMesh, tMax=tMax, tSteps=tSteps, peclet=peclet, verbose=1)

    if load:
        rhoa = np.load('synthRhoa.npy') 
        err = np.load('synthErr.npy') 
        scheme = pb.DataContainerERT('synth.shm')
    else:
        
        rhoa = fop.response(np.array(model)[paraMesh.cellMarker()])

        err = fop.ws.err
        rand = pg.RVector(len(rhoa)); pg.randn(rand)
        rhoa *= (1.0 + rand * err.flatten())
        rhoa = rhoa.reshape(err.shape)
        scheme = fop.ws.dataERT
    
        np.save('synthRhoa', rhoa) 
        np.save('synthErr', err) 
        scheme.save('synth.shm', 'a b m n')
        model=[0.001, 0.01, 1e-7]
        
    if show:
        ws = fop.ws
        mesh = ws.mesh
        ax, _ = pg.show(mesh)
        pg.show(mesh, fop.ws.p, label='Pressure in ??', axes=ax)
        ax, _ = pg.show(mesh)
        pg.show(mesh, fop.ws.k, label='Permeabilty in ??', axes=ax)
        
        ax, _ = pg.show(mesh)
        pg.show(mesh, np.sqrt(fop.ws.vel[0]**2+fop.ws.vel[1]**2), label='velocity in m/s', axes=ax)
        pg.show(mesh, fop.ws.vel, axes=ax)

        pg.show(mesh, fop.ws.sat[-1]+1e-3, label='brine saturation')
        pg.show(ws.meshERT, fop.ws.res[-1], label='resistivity in Ohm m')

        showERTData(scheme, rhoa)
        
        pg.wait()
    
    return rhoa, err, fop
    
if __name__ == '__main__':
    
    rhoa, err, fop = simulateSynth(tMax=5000, tSteps=10, peclet=500000, show=0, load=1)
    
    #### create paramesh and startmodel
    paraMesh = pg.createGrid(x=[0, 10], y=[-5, -3.5, -0.5, 0])
    paraMesh.cell(0).setMarker(0) # bottom
    paraMesh.cell(1).setMarker(1) # center
    paraMesh.cell(2).setMarker(2) # top
    paraMesh = paraMesh.createH2()
    paraMesh = paraMesh.createH2()
    paraMesh = paraMesh.createH2()
    
    fop.setMesh(paraMesh)
    fop.regionManager().region(0).setFixValue(0.001)
    fop.regionManager().region(2).setFixValue(1e-7)
    
    fop.createRefinedForwardMesh(refine=False, pRefine=False)
    
    paraDomain = fop.regionManager().paraDomain()
    
    startModel = pg.RVector(fop.regionManager().parameterCount(), 0.007)
    fop.setStartModel(startModel) 
    
    ##### create paramesh and startmodel

    #rhoa2 = fop.response(startModel)
    #print('norm: ', np.linalg.norm(np.log10(rhoa)-np.log10(rhoa2), ord=2))
    #print('rms: ' , pg.rms(np.log10(rhoa)-np.log10(rhoa2)))
    #print('chi: ' , pg.rms( (np.log10(rhoa)-np.log10(rhoa2)) / err ))
    #print('chi²(lin): ' , pg.utils.chi2(rhoa, rhoa2, err))
    #print('chi²(log): ' , pg.utils.chi2(rhoa, rhoa2, err, pg.RTransLog()))
      
    #pg.show(paraDomain, startModel[paraDomain.cellMarker()], label='Permeabilty START', cMin=1e-7, cMax=0.02)    
    #pg.wait()
    
    inv = pg.RInversion(rhoa.flatten(), fop, verbose=1, dosave=1)
     
    tD = pg.RTransLog()
    tM = pg.RTransLogLU(1e-8,1)
    inv.setTransData(tD)
    inv.setTransModel(tM)
    
    inv.setRelativeError(err.flatten())
    #inv.setLambda(0)
    inv.setMaxIter(20)
    inv.setLineSearch(True)
    inv.setLambda(100)
    #inv.setMarquardtScheme(0.8)
    
    coeff = inv.run()
    pg.show(paraDomain, coeff[paraDomain.cellMarker()], label='Permeabilty INV')    
    coeff.save("permModel.vector")

# actual inversion run yielding coefficient model



pg.wait()

