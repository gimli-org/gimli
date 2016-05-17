#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""    
    Special meta forward operator for modelling with petrophysical relations
    
"""

import pygimli as pg
import numpy as np

from pygimli.physics import MethodManager


class PetroModelling(pg.ModellingBase):
    """
    Combine petrophysical relation m(p) with modelling class f(p)
    
    Combine petrophysical relation m(p) with modelling class f(p) to invert 
    for m (or any inversion transformation) instead of p.
    """
    def __init__(self, fop, trans, mesh=None, verbose=False):
        """ save forward class and transformation, create Jacobian matrix """
        super().__init__(verbose=verbose)
        self.fop = fop
        self.trans = trans  # class defining m(p)
        #self.setData(self.fop.data())
        
        if mesh is not None:
            self.setMesh(mesh)

    def setData(self, data):
        pg.ModellingBase.setData(self, data)
        self.fop.setData(data)
            
    def setMesh(self, mesh):
        
        if mesh is None and self.fop.mesh() is None:
            raise StandardException("Please provide a mesh for this forward operator")
        
        if mesh is not None:
            self.fop.setMesh(mesh)
            self.fop.createRefinedForwardMesh(refine=False)

        #self.setMesh(f.mesh(), ignoreRegionManager=True) # not really nessary
        self.setRegionManager(self.fop.regionManagerRef())

        self.nModel = self.regionManager().parameterCount()
        
        self.J = pg.RMultRMatrix(self.fop.jacobian())
        self.setJacobian(self.J)
        
    def response(self, model):
        """ use inverse transformation to get p(m) and compute response """
        tModel = self.trans.trans(model)
        ret = self.fop.response(tModel)
        return ret

    def createJacobian(self, model):
        """ fill the individual jacobian matrices"""
        par = self.trans.trans(model)
        self.fop.createJacobian(par)
        self.J.r = self.trans.deriv(model)  # set inner derivative


class PetroJointModelling(pg.ModellingBase):
    """ Cumulative (joint) forward operator for petrophysical inversions """
    def __init__(self, f=None, p=None, mesh=None, verbose=True):
        super().__init__(verbose=verbose)
        
        self.fops = None
        self.J = None
        self.Ji = None
        self.mesh = None

        if f is not None and p is not None:
            self.setFopsAndTrans(f, p)

    def setFopsAndTrans(self, fops, trans):
        self.fops = [PetroModelling(fi, pi, self.mesh) for fi, pi in zip(fops, trans)]

    def setMesh(self, mesh):
        
        self.mesh=mesh
        for i, fi in enumerate(self.fops):
            fi.setMesh(mesh)
            
        self.setRegionManager(self.fops[0].regionManagerRef())
        nModel = self.regionManager().parameterCount()
        
        self.initJacobian()
        
    def setData(self, data):
        for i, fi in enumerate(self.fops):
            fi.setData(data[i])
           
        self.initJacobian()    
             
    def initJacobian(self):
        self.J = pg.RBlockMatrix()
        self.Ji = [self.J.addMatrix(fi.jacobian()) for fi in self.fops]
        nData = 0
        for i, fi in enumerate(self.fops):
            self.J.addMatrixEntry(self.Ji[i], nData, 0)
            nData += fi.data().size()  # update total vector length
           
        self.setJacobian(self.J)
    
    def response(self, model):
        """ cumulative response """
        resp = []
        for f in self.fops:
            resp.extend(f.response(model))
        return resp

    def createJacobian(self, model):
        """ just force creating individual Jacobian matrices """
        for f in self.fops:
            f.createJacobian(model)
            

class InvertJointPetro(MethodManager):

    def __init__(self, managers, trans, verbose=False, debug=False, **kwargs):
        MethodManager.__init__(self, verbose=verbose, debug=debug, **kwargs)
        
        self.managers = managers
        self.trans = trans
        self.fops = []
        self.dataVals = pg.RVector(0)
        self.dataErrs = pg.RVector(0)
        self.model = pg.RVector(0)
        
        self.tD = pg.RTransCumulative()
        self.tM = managers[0].tM
        
        for i, mgr in enumerate(self.managers):
            fop = mgr.createFOP(verbose)
            fop.setVerbose(verbose=verbose)
            self.fops.append(fop)
            
        self.fop.setFopsAndTrans(self.fops, self.trans)
                        
    @staticmethod
    def createFOP(verbose=False):
        """Create forward operator
        """
        fop = PetroJointModelling()
        
        return fop
     
    def createInv(self, fop, verbose=True, doSave=False):
        
        inv = pg.RInversion(verbose, doSave)
        inv.setForwardOperator(fop)

        return inv
    
    def setData(self, data):
        """
        """
        if type(data) is list:
            if len(data) == len(self.managers):
                self.tD.clear()
                self.dataVals.clear()
                self.dataErrs.clear()
                
                self.fop.setData(data)
                
                for i, mgr in enumerate(self.managers):
                    t = mgr.tD
                    self.tD.add(t, data[i].size())
    
                    self.dataVals = pg.cat(self.dataVals, data[i](mgr.dataToken()))
            
                    if mgr.errIsAbsolute:
                        self.dataErrs = pg.cat(self.dataErrs, data[i]('err')/data[i](mgr.dataToken()))
                    else:
                        self.dataErrs = pg.cat(self.dataErrs, data[i]('err'))
                
                self.data = data
                
                self.inv.setTransData(self.tD)
                self.inv.setTransModel(self.tM)
            else:
                raise BaseException("To few datacontainer given")
            
    
    def setMesh(self, mesh):
        self.fop.setMesh(mesh)
            
    def invert(self, data, mesh, lam=20, limits=None):
        
        self.setData(data)
        self.setMesh(mesh)

        if limits is not None:
            if hasattr(self.tM, 'setLowerBound'):
                if self.verbose:
                    print('Lower limit set to', limits[0])
                self.tM.setLowerBound(limits[0])
            if hasattr(self.tM, 'setUpperBound'):
                if self.verbose:
                    print('Upper limit set to', limits[1])
                self.tM.setUpperBound(limits[1])
        
        nModel = self.fop.regionManager().parameterCount()
        
        startModel = pg.RVector(nModel, 0.0)
        
        for i in range(len(self.managers)):
            startModel += pg.RVector(nModel,
                                     pg.median(self.trans[i].inv(self.managers[i].createApparentData(self.data[i]))))
        startModel /= len(self.managers)
        self.fop.setStartModel(startModel)
    
        
        self.inv.setData(self.dataVals)
        self.inv.setRelativeError(self.dataErrs)
        self.inv.setLambda(lam)

        self.model = self.inv.start()
        self.model = self.model(self.fop.regionManager().paraDomain().cellMarkers())
        return self.model
        
    def showModel(self, **showkwargs):
        if len(showkwargs):
            pg.show(self.fop.regionManager().paraDomain(), 
                    self.model, **showkwargs)

class InvertPetro(InvertJointPetro):

    def __init__(self, manager, trans,  **kwargs):
        InvertJointPetro.__init__(self, [manager], [trans], **kwargs)
        
    def invert(self, data, **kwargs):
        return InvertJointPetro.invert(self, [data], **kwargs)
          

if __name__ == "__main__":
    pass

