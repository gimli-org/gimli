#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""    
    Special meta forward operator for modelling with petrophysical relations
    
"""

import pygimli as pg
import numpy as np

class PetroModelling(pg.ModellingBase):
    """
    Combine petrophysical relation m(p) with modelling class f(p)
    
    Combine petrophysical relation m(p) with modelling class f(p) to invert 
    for m (or any inversion transformation) instead of p.
    """
    def __init__(self, f, trans, mesh=None, verbose=False):
        """ save forward class and transformation, create Jacobian matrix """
        super().__init__(verbose=verbose)
        self.f = f
        self.trans = trans  # class defining m(p)
        self.setData(self.f.data())
        
        if mesh is None and f.mesh() is None:
            raise StandardException("Please provide a mesh for this forward operator")
            
        if mesh is not None:
            self.f.setMesh(mesh)
            self.f.createRefinedForwardMesh(refine=False)
       
        #self.setMesh(f.mesh(), ignoreRegionManager=True) # not really nessary
        self.setRegionManager(self.f.regionManagerRef())

        self.nModel = self.regionManager().parameterCount()
        
        self.J = pg.RMultRMatrix(self.f.jacobian())
        self.setJacobian(self.J)

    def response(self, model):
        """ use inverse transformation to get p(m) and compute response """
        tModel = self.trans.trans(model)
        ret = self.f.response(tModel)
        return ret

    def createJacobian(self, model):
        """ fill the individual jacobian matrices"""
        par = self.trans.trans(model)
        self.f.createJacobian(par)
        self.J.r = self.trans.deriv(model)  # set inner derivative

####################
# we need this?? Ca
####################
#class PetroInvModelling(pg.ModellingBase):
    #"""
    #class for combining petrophysical relation m(p) with modelling class f(p)
    #to invert for m (or any inversion transformation) instead of p
    #"""
    #def __init__(self, f, trans, mesh=None, verbose=False):
        #""" save forward class and transformation, create Jacobian matrix """
        #super().__init__(verbose=verbose)
        #self.f = f
        #self.trans = trans  # class defining m(p)
        #if mesh is None:
            #self.setMesh(f.mesh())  # this part is not perfectly clean
        #else:
            #self.setMesh(mesh)
        #self.createRefinedForwardMesh(False)
        #self.regionManager().setConstraintType(1)
        #self.nModel = f.regionManager().parameterCount()
        #tmpModel = pg.RVector(self.nModel, 1)  # for now a ones model
        #self.J = pg.RMultRMatrix(self.f.jacobian(), tmpModel)
        #self.setJacobian(self.J)

    #def response(self, model):
        #""" use inverse transformation to get p(m) and compute response """
        #return self.f.response(self.trans.invTrans(model))

    #def createJacobian(self, model):
        #""" fill the individual jacobian matrices"""
        #par = self.trans.invTrans(model)
        #self.f.createJacobian(par)
        #self.J.r = 1./self.trans.deriv(par)  # set inner derivative


class PetroJointModelling(pg.ModellingBase):
    """ Cumulative (joint) forward operator for petrophysical inversions """
    def __init__(self, f, p, mesh, verbose=True):
        super().__init__(verbose=verbose)

        self.fP = [PetroModelling(fi, pi, mesh) for fi, pi in zip(f, p)]

        self.setRegionManager(self.fP[0].regionManagerRef())
        
        nModel = self.regionManager().parameterCount()
        
        self.J = pg.RBlockMatrix()
        self.Ji = [self.J.addMatrix(fi.jacobian()) for fi in self.fP]
        nData = 0
        for i, fi in enumerate(self.fP):
            self.J.addMatrixEntry(self.Ji[i], nData, 0)
            nData += fi.data().size()  # update total vector length
           
        self.setJacobian(self.J)

    def response(self, model):
        """ cumulative response """
        resp = []
        for f in self.fP:
            resp.extend(f.response(model))
        return resp

    def createJacobian(self, model):
        """ just force creating individual Jacobian matrices """
        for f in self.fP:
            f.createJacobian(model)
            


def invertPetro(manager, data, trans=pg.RTrans(), lam=20, mesh=None,
                limits=None, verbose=False, **showkwargs):
    """
        Simplistic petrophysical inversion framework. 
        
        For testing purposes currently .. maybe better in pg.inversion
    """
    fop = manager.createFOP(verbose)
    fop.setVerbose(verbose=verbose)
    fop.setData(data)
    
    dataVal = data(manager.dataToken())
    dataErr = data('err')

    if mesh is None:
        # CR: don't happy with this
        manager.setMesh(data)
        
        mesh = manager.createMesh()
        
    f = PetroModelling(fop, trans, mesh)
    tD = manager.tD
    tM = manager.tM
    
    if limits is not None:
        if hasattr(tM, 'setLowerBound'):
            if verbose:
                print('Lower limit set to', limits[0])
            tM.setLowerBound(limits[0])
        if hasattr(tM, 'setUpperBound'):
            if verbose:
                print('Upper limit set to', limits[1])
            tM.setUpperBound(limits[1])
                
    inv = pg.RInversion(dataVal, f, tD, tM, verbose=verbose)

    if manager.errIsAbsolute:
        inv.setAbsoluteError(dataErr)
    else:
        inv.setRelativeError(dataErr)

    startModel = pg.RVector(fop.regionManager().parameterCount(),
                            pg.median(trans.inv(manager.createApparentData(data))))
    inv.setModel(startModel)
    inv.setLambda(lam)

    model = inv.run()

    if len(showkwargs):
        pg.show(fop.regionManager().paraDomain(), model, **showkwargs)

    # we need to return f and/or inv .. because fop is temporary here without 
    # an increased python reference counter.
    # fop is set to the inv object but the python-internal reference counter 
    # is not increased so fop and all content will be else removed here
    return model, f, inv

if __name__ == "__main__":
    pass

