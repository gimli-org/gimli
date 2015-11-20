#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pygimli as pg
import pybert as pb
import numpy as np

from multi_darcy_flow import darcyFlow
from multi_advection import calcSaturation
from multi_ert import simulateERTData, showERTData

from multiprocessing import Process, Array
import copy

from math import ceil
import time
class WorkSpace():
    pass

def responseProc(model, peclet, timesAdvection, tSteps,  i, j, halfInjection=0, resShm=None, ws=None):
    
    if resShm is not None: print("proc:", i)
            
    if ws is None:
        ws = WorkSpace()

    if resShm is None: print(".. darcy step")

    ws.mesh, ws.vel, ws.p, ws.k = darcyFlow(model, p0=0.25, verbose=0)

    #pg.show(model, model.cellAttributes(), label='cell cellAttributes'); pg.wait()
    #pg.show(ws.mesh, ws.vel.T, label='velocity'); pg.wait()
    #pg.show(ws.mesh, ws.k, label='Permeabilty'); pg.wait()
    
    if resShm is None: print(".. saturation step")        
        
    ws.sat = calcSaturation(ws.mesh, ws.vel.T, 
                            timesAdvection, 
                            peclet=peclet, 
                            halfInjection=halfInjection,
                            verbose=0)
    iProc = i
    
    if resShm is None: print(".. ert step")
        
        
    sampleTime = pg.IndexArray(np.floor(np.linspace(0, len(timesAdvection)-1, tSteps)))
    #print(sampleTime)
              
    #sampleTime = [0, tSteps-1]
    ws.meshERT, ws.dataERT, ws.res, ws.rhoa, ws.err = simulateERTData(ws.sat[sampleTime],
                                                                      ws.mesh, i=iProc, j=j, verbose=0)

    if resShm:
        rhoa = ws.rhoa.flatten()
        
        for j in range(len(rhoa)):
            resShm[j] = rhoa[j]
    else:
        print((np.min(ws.rhoa), np.max(ws.rhoa)))
        return ws.rhoa
    
class MulitFOP(pg.ModellingBase):
    def __init__(self, verbose=False, **kwargs):

        pg.ModellingBase.__init__(self, verbose=verbose)
        self.init(kwargs.pop('mesh', None), 
                  kwargs.pop('tMax', 50000),
                  kwargs.pop('satSteps', 50),
                  kwargs.pop('ertSteps', 5),
                  kwargs.pop('peclet', 50000))
        self.iter = 0
        
    def init(self, mesh, tMax, satSteps, ertSteps, peclet):
    
        self.parMesh = pg.Mesh(mesh)
        self.setMesh(mesh)
        self.createRefinedForwardMesh(refine=False, pRefine=False)
        
        self.tMax = tMax
        self.satSteps = satSteps
        self.ertSteps = ertSteps
        self.timesAdvection = np.linspace(1, tMax, satSteps)
        self.peclet = peclet
        self.halfInjection = 0
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
        nProcs = float(pg.numberOfCPU())
            
        rhoaShm = []
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
            
                    #rhoaJ[i] = responseProc(modelMesh, self.peclet, self.timesAdvection, self.ertSteps, i, self.iter).flatten()

                    rhoaShm.append(Array('d', len(resp)))
                    procs.append(Process(target=responseProc, args=(modelMesh, 
                                                                    self.peclet, self.timesAdvection, self.ertSteps,
                                                                    i, self.iter, self.halfInjection, rhoaShm[i])))

            for i, p in enumerate(procs):
                p.start()
                
            for i, p in enumerate(procs):
                p.join()
                
        pg.toc()

        self._J *= 0.0
        for i in range(nModel):
            dData = rhoaShm[i] - resp
            self._J.setCol(i, dData/dModel[i])
            
        self._J.save('jacobi' + str(self.iter))
        print('#'*40 + 'Jac:')
        self.cov = np.zeros(nModel)
        for i in range(nModel):
            self.cov[i] = sum(pg.abs(self._J.col(i))) 
            
            
        print(self.cov)
        #exit()
    def response(self, par):
        model = par
        if self.mesh():
            if len(par) == self.mesh().cellCount():
                self.mesh().setCellAttributes(par)
            else:
                self.mapModel(par)
        model = self.mesh()
        
        rhoa = responseProc(model, self.peclet, self.timesAdvection, self.ertSteps,
                     0, 0, self.halfInjection, None, self.ws)
                
        return rhoa.flatten()   
        
    
def simulateSynth(model, tMax=5000, satSteps=50, ertSteps=10, peclet=500000, halfInjection=0, show=False, load=False):
    
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
            
    #pg.show(paraMesh, np.array(model)[paraMesh.cellMarker()], label='Permeabilty Model')
    #pg.show(paraMesh, paraMesh.cellMarker(), label='marker')
    #pg.wait()
    
    fop = MulitFOP(mesh=paraMesh, tMax=tMax, 
                   satSteps=satSteps,
                   ertSteps=ertSteps, 
                   peclet=peclet, verbose=1)
    
    fop.halfInjection = halfInjection
    
    if load:
        rhoa = np.load('synthRhoa.npy') 
        err = np.load('synthErr.npy') 
        scheme = pb.DataContainerERT('synth.shm')
    else:
        
        rhoa = fop.response(pg.RVector(model)[paraMesh.cellMarker()])

        err = fop.ws.err
        rand = pg.RVector(len(rhoa)); pg.randn(rand)
        rhoa *= (1.0 + rand * err.flatten())
        rhoa = rhoa.reshape(err.shape)
        scheme = fop.ws.dataERT
    
        np.save('synthRhoa', rhoa) 
        np.save('synthErr', err) 
        scheme.save('synth.shm', 'a b m n')
        
    if show:
        ws = fop.ws
        mesh = ws.mesh
        ax = [pg.plt.subplot(2,3,i+1) for i in range(6)]
        
        pg.show(mesh, axes=ax[0])
        pg.show(mesh, fop.ws.p, label='Pressure in ??', axes=ax[0])
        
        pg.show(mesh, fop.ws.k, label='Permeabilty in ??', axes=ax[1])
        
        pg.show(mesh, np.sqrt(fop.ws.vel[0]**2+fop.ws.vel[1]**2), label='Velocity in m/s', axes=ax[2])
        pg.show(mesh, fop.ws.vel, axes=ax[2])

        pg.show(mesh, fop.ws.sat[-1]+1e-3, label='Brine concentration for t=' + str(int(tMax/3600)) + "h",  axes=ax[3])
        pg.show(ws.meshERT, fop.ws.res[-1], label='Resistivity in Ohm m', axes=ax[4])
        ax[4].set_xlim([-1, 11])
        ax[4].set_ylim([-6, 0])
        showERTData(scheme, rhoa, axes=ax[5])
        
        pg.wait()
    
    return rhoa, err, fop
    
if __name__ == '__main__':
    paraRefine = 2
    show = 1
    load = 0
    halfInjection=1
        
    rhoa, err, fop = simulateSynth(model=[1e-4, 5e-3, 1e-8, 2.5e-3],
                                   tMax=72000, satSteps=50, ertSteps=10,
                                   peclet=5e6, halfInjection=halfInjection,
                                   show=show, load=load)
    paraMesh = pg.createGrid(x=[0, 5, 10], y=[-5, -3.5, -0.5, 0])
    paraMesh.cell(0).setMarker(0) # bottom
    paraMesh.cell(1).setMarker(0) # bottom
    paraMesh.cell(2).setMarker(1) # center
    paraMesh.cell(3).setMarker(1) # center
    paraMesh.cell(4).setMarker(2) # top
    paraMesh.cell(5).setMarker(2) # top
    
    
    for i in range(paraRefine):
        paraMesh = paraMesh.createH2()
        
    fop.setMesh(paraMesh)
    fop.regionManager().region(0).setFixValue(1e-4)
    fop.regionManager().region(2).setFixValue(1e-8)
    
    fop.createRefinedForwardMesh(refine=False, pRefine=False)
    
    startModel = pg.RVector(fop.regionManager().parameterCount(), 4e-3)

    fop.setStartModel(startModel) 
    
    inv = pg.RInversion(rhoa.flatten(), fop, verbose=1, dosave=1)
     
    tD = pg.RTransLog()
    tM = pg.RTransLogLU(1e-9,1)
    inv.setTransData(tD)
    inv.setTransModel(tM)
    
    inv.setRelativeError(err.flatten())
    #inv.setLambda(0)
    inv.setMaxIter(20)
    inv.setLineSearch(True)
    inv.setLambda(100)
    #inv.setMarquardtScheme(0.8)
    
    coeff = inv.run()
    
    pd = fop.regionManager().paraDomain()
    pd.save("permModel")
    coeff.save("permModel.vector")
    cov = pg.RVector(np.log10(fop.cov))
    cov.save("permModelCov.vector")

    #pg.show(paraDomain, coeff[paraDomain.cellMarker()], label='Permeabilty INV')    
# actual inversion run yielding coefficient model



pg.wait()

