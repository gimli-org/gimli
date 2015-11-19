#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pygimli as pg
import pybert as pb
import numpy as np

from multi_darcy_flow import darcyFlow
from multi_advection import calcSaturation
from multi_ert import simulateERTData, showERTData

from multiprocessing import Process, Queue

class WorkSpace():
    pass

class MulitFOP(pg.ModellingBase):
    def __init__(self, paraMesh=None, tMax=5000, tSteps=10, peclet=500, verbose=False):
        pg.ModellingBase.__init__(self, verbose=verbose)

        self.tSteps = tSteps
        if paraMesh:
            self.setMesh(paraMesh)
            self.createRefinedForwardMesh(refine=False, pRefine=False)
        else:
            self.regionManager().setParameterCount(3)
        
        self.ws = WorkSpace()
        self.timesAdvection = np.linspace(1, tMax, tSteps)
        self.peclet = peclet
        self._J = pg.RMatrix()
        self.setJacobian(self._J)   
        self.iter = 0
        
    def createJacobian(self, model):
        self.iter += 1  
        print("Create Jacobian for", str(model))
        nModel = len(model)
        resp = self.response(model)
        nData = len(resp)
        
        self._J.resize(nData, nModel)
        self._J *= 0.
        
        respChange = pg.RMatrix(nModel, nData)
        dModel = pg.RVector(len(model))

        fak = 1.05

        pg.tic()
        output = Queue()
        procs = []
        for i in range(nModel):
            #self.createResponseProc(model, fak, i, output) # fails reprodu
            procs.append(Process(target=self.createResponseProc, args=(model, fak, i, output)))

        for p in procs:
            p.start()
                        
        for p in procs:
            p.join()
        pg.toc()

        for i in range(nModel):
            Imodel, respChange, dModel = output.get()
            dData = respChange - resp
            self._J.setCol(Imodel, dData/dModel)
        
        print('#'*40 + 'Jac:')
        for i in range(nModel):
            print(i, sum(pg.abs(self._J.col(i))))
            
    
    def createResponseProc(self, par, fak, i, output):
        print("proc:", i)
        modelChange = pg.RVector(par)
        modelChange[i] *= fak
        dModel = modelChange[i]-par[i]
                
        model = par
        if self.mesh():
            meshIn = pg.Mesh(self.mesh())
            meshIn.setCellAttributes(modelChange[meshIn.cellMarker()])
            model = meshIn
        
        mesh, vel, p, k = darcyFlow(model, verbose=0)
        
        sat = calcSaturation(mesh, vel.T, 
                             self.timesAdvection, 
                             peclet=self.peclet, 
                             verbose=0)
        
        sampleTime = [0, self.tSteps/4, self.tSteps/2, 3*self.tSteps/4, self.tSteps-1]
        meshERT, dataERT, res, rhoa, err = simulateERTData(sat[sampleTime], mesh, verbose=0)

        output.put([i, rhoa.flatten(), dModel])
                                
                
    def response(self, par):
        ws = self.ws
        print('create response for: ' + str(par))
        model = par

        if self.mesh():
            self.mapModel(model)
            model = self.mesh()
            
            #if len(par) == model.cellCount():
                #model.setCellAttributes(par)
        
        
        print(self.mesh())
        print(self.mesh().cellMarker())
        print(self.mesh().cellAttributes())
        print(min(model.cellAttributes()), max(model.cellAttributes()))
        pg.show(model, model.cellAttributes(), label='Permeabilty (model)')
        #pg.wait()
        print(".. darcy step")
        ws.mesh, ws.vel, ws.p, ws.k = darcyFlow(model, verbose=0)
        #pg.show(ws.mesh, ws.k, label='Permeabilty iter:' + str(self.iter))
        #pg.wait()
  
        print(".. advection step")
        ws.saturation = calcSaturation(ws.mesh, ws.vel.T,
                                       self.timesAdvection, 
                                       peclet=self.peclet, verbose=0)
        
        print(".. dc step")
        sampleTime = [0, self.tSteps/4, self.tSteps/2, 3*self.tSteps/4, self.tSteps-1]
        ws.meshERT, ws.dataERT, ws.res, ws.rhoa, ws.err = simulateERTData(ws.saturation[sampleTime], ws.mesh, verbose=0)

        print(".. resp min/max", min(ws.rhoa.flat), max(ws.rhoa.flat))
        return ws.rhoa.flatten()
    
    
def simulateSynth(tMax=5000, tSteps=40, peclet=50000, show=True, load=False):
    
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
        
    fop = MulitFOP(paraMesh, tMax=tMax, tSteps=tSteps, peclet=peclet, verbose=1)

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

        pg.show(mesh, fop.ws.saturation[-1]+1e-3, label='brine saturation')
        pg.show(ws.meshERT, fop.ws.res[-1], label='resistivity in Ohm m')

        showERTData(scheme, rhoa)
        
        pg.wait()
    
    return rhoa, err, fop
    
if __name__ == '__main__':
    
    rhoa, err, fop = simulateSynth(tMax=5000, tSteps=40, peclet=500000, show=0, load=1)
    
    #### create paramesh and startmodel
    paraMesh = pg.createGrid(x=[0, 10], y=[-5, -3.5, -0.5, 0])
    paraMesh.cell(0).setMarker(0) # bottom
    paraMesh.cell(1).setMarker(1) # center
    paraMesh.cell(2).setMarker(2) # top
    paraMesh = paraMesh.createH2()
    #paraMesh = paraMesh.createH2()
    
    fop.setMesh(paraMesh)
    fop.regionManager().region(0).setFixValue(0.001)
    fop.regionManager().region(2).setFixValue(1e-7)
    
    fop.createRefinedForwardMesh(refine=False, pRefine=False)
    
    paraDomain = fop.regionManager().paraDomain()
    
    startModel = pg.RVector(fop.regionManager().parameterCount(), 0.007) # else
    #startModel[0] = 0.001 # bottom
    #startModel[-1] = 1e-7 # top
    fop.setStartModel(startModel) 
    
    ##### create paramesh and startmodel
    #0.000609525030457, 0.00695938025714, 2.15533816606e-07
    #[0.000609525030503, 0.0069593802455, 2.15534036123e-07]



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
    inv.setLineSearch(True)
    inv.setLambda(100)
    #inv.setMarquardtScheme(0.8)
    
    coeff = inv.run()
    pg.show(paraDomain, coeff[paraDomain.cellMarker()], label='Permeabilty INV', cMin=1e-7, cMax=0.02)    
    print(coeff)

# actual inversion run yielding coefficient model



pg.wait()

