#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""    
    Special meta forward operator for modelling with petrophysical relations
    
"""

import pygimli as pg
import numpy as np

from pygimli.physics import MethodManager

#class PetroModelling(pg.ModellingBase):
    #"""
    #Combine petrophysical relation m(p) with modelling class f(p)
    
    #Combine petrophysical relation m(p) with modelling class f(p) to invert 
    #for m (or any inversion transformation) instead of p.
    #"""
    #def __init__(self, f, trans, mesh=None, verbose=False):
        #""" save forward class and transformation, create Jacobian matrix """
        #super().__init__(verbose=verbose)
        #self.f = f
        #self.trans = trans  # class defining m(p)
        #self.setData(self.f.data())
        
        #if mesh is not None:
            #self.setMesh(mesh)

    #def setData(self, data):
        #self.f.setData()
            
    #def setMesh(self, mesh):
        
        #if mesh is None and f.mesh() is None:
            #raise StandardException("Please provide a mesh for this forward operator")
        
        #if mesh is not None:
            #self.f.setMesh(mesh)
            #self.f.createRefinedForwardMesh(refine=False)

        ##self.setMesh(f.mesh(), ignoreRegionManager=True) # not really nessary
        #self.setRegionManager(self.f.regionManagerRef())

        #self.nModel = self.regionManager().parameterCount()
        
        #self.J = pg.RMultRMatrix(self.f.jacobian())
        #self.setJacobian(self.J)
        
    #def response(self, model):
        #""" use inverse transformation to get p(m) and compute response """
        #tModel = self.trans.trans(model)
        #ret = self.f.response(tModel)
        #return ret

    #def createJacobian(self, model):
        #""" fill the individual jacobian matrices"""
        #par = self.trans.trans(model)
        #self.f.createJacobian(par)
        #self.J.r = self.trans.deriv(model)  # set inner derivative

#####################
## we need this?? Ca
#####################
##class PetroInvModelling(pg.ModellingBase):
    ##"""
    ##class for combining petrophysical relation m(p) with modelling class f(p)
    ##to invert for m (or any inversion transformation) instead of p
    ##"""
    ##def __init__(self, f, trans, mesh=None, verbose=False):
        ##""" save forward class and transformation, create Jacobian matrix """
        ##super().__init__(verbose=verbose)
        ##self.f = f
        ##self.trans = trans  # class defining m(p)
        ##if mesh is None:
            ##self.setMesh(f.mesh())  # this part is not perfectly clean
        ##else:
            ##self.setMesh(mesh)
        ##self.createRefinedForwardMesh(False)
        ##self.regionManager().setConstraintType(1)
        ##self.nModel = f.regionManager().parameterCount()
        ##tmpModel = pg.RVector(self.nModel, 1)  # for now a ones model
        ##self.J = pg.RMultRMatrix(self.f.jacobian(), tmpModel)
        ##self.setJacobian(self.J)

    ##def response(self, model):
        ##""" use inverse transformation to get p(m) and compute response """
        ##return self.f.response(self.trans.invTrans(model))

    ##def createJacobian(self, model):
        ##""" fill the individual jacobian matrices"""
        ##par = self.trans.invTrans(model)
        ##self.f.createJacobian(par)
        ##self.J.r = 1./self.trans.deriv(par)  # set inner derivative


#class PetroJointModelling(pg.ModellingBase):
    #""" Cumulative (joint) forward operator for petrophysical inversions """
    #def __init__(self, f=None, p=None, mesh=None, verbose=True):
        #super().__init__(verbose=verbose)
        
        #self.fP = None
        #self.J = None
        #self.Ji = None
        #self.mesh = None
        #if f is not None and p is not None:
            #self.setFopsAndTrans(f, p)

    #def setMesh(self, mesh):
        
        #self.mesh=mesh
        #for f in self.fp:
            #f.setMesh(mesh)
        
    #def setData(self, data):
        #for i, f in enumerate(self.fp):
            #f.setData(data[i])
                
    #def setFopsAndTrans(self, fops, trans):
        #self.fP = [PetroModelling(fi, pi, self.mesh) for fi, pi in zip(fops, trans)]
        #self.setRegionManager(self.fP[0].regionManagerRef())
        
        #nModel = self.regionManager().parameterCount()
        
        #self.J = pg.RBlockMatrix()
        #self.Ji = [self.J.addMatrix(fi.jacobian()) for fi in self.fP]
        #nData = 0
        #for i, fi in enumerate(self.fP):
            #self.J.addMatrixEntry(self.Ji[i], nData, 0)
            #nData += fi.data().size()  # update total vector length
           
        #self.setJacobian(self.J)


    #def response(self, model):
        #""" cumulative response """
        #resp = []
        #for f in self.fP:
            #resp.extend(f.response(model))
        #return resp

    #def createJacobian(self, model):
        #""" just force creating individual Jacobian matrices """
        #for f in self.fP:
            #f.createJacobian(model)
            


#def invertPetro0(manager, data, trans=pg.RTrans(), lam=20, mesh=None,
                #limits=None, verbose=False, **showkwargs):
    #"""
        #Simplistic petrophysical inversion framework. 
        
        #For testing purposes currently .. maybe better in pg.inversion
    #"""
    #fop = manager.createFOP(verbose)
    #fop.setVerbose(verbose=verbose)
    #fop.setData(data)
    
    #dataVal = data(manager.dataToken())
    #dataErr = data('err')

    #if mesh is None:
        ## CR: don't happy with this
        #manager.setMesh(data)
        
        #mesh = manager.createMesh()
        
    #f = PetroModelling(fop, trans, mesh)
    #tD = manager.tD
    #tM = manager.tM
    
    #if limits is not None:
        #if hasattr(tM, 'setLowerBound'):
            #if verbose:
                #print('Lower limit set to', limits[0])
            #tM.setLowerBound(limits[0])
        #if hasattr(tM, 'setUpperBound'):
            #if verbose:
                #print('Upper limit set to', limits[1])
            #tM.setUpperBound(limits[1])
                
    #inv = pg.RInversion(dataVal, f, tD, tM, verbose=verbose)

    #if manager.errIsAbsolute:
        #inv.setAbsoluteError(dataErr)
    #else:
        #inv.setRelativeError(dataErr)

    #startModel = pg.RVector(fop.regionManager().parameterCount(),
                            #pg.median(trans.inv(manager.createApparentData(data))))
    #inv.setModel(startModel)
    #inv.setLambda(lam)

    #model = inv.run()

    #if len(showkwargs):
        #pg.show(fop.regionManager().paraDomain(), model, **showkwargs)

    ## we need to return f and/or inv .. because fop is temporary here without 
    ## an increased python reference counter.
    ## fop is set to the inv object but the python-internal reference counter 
    ## is not increased so fop and all content will be else removed here
    #return model, f, inv


#class InvertJointPetro(MethodManager):

    #def __init__(self, managers, trans, verbose=False, debug=False, **kwargs):
        #MethodManager.__init__(self, verbose=verbose, debug=debug, **kwargs)
        
        #self.managers = managers
        #self.trans = trans
        #self.fops = []
        #self.dataVal = pg.RVector(0)
        #self.dataErr = pg.RVector(0)
        
        #self.tD = pg.RTransCumulative()
        #self.tM = managers[0].tM
        #for i, mgr in enumerate(self.managers):
            #fop = mgr.createFOP(verbose)
            #fop.setVerbose(verbose=verbose)
            #self.fops.append(fop)
            
        #self.fop.setFopsAndTrans(self.fops, self.trans)
        
                        
    #@staticmethod
    #def createFOP(verbose=False):
        #"""Create forward operator
        #"""
        #fop = PetroJointModelling()
        
        #return fop
     
    #def createInv(self, fop, verbose=True, doSave=False):
        
        #inv = pg.RInversion(verbose, doSave)
        #inv.setForwardOperator(fop)

        #return inv
    
    #def setData(self, data):
        #"""
        #"""
        #if type(data) is list:
            #if len(data) == len(self.managers):
                #self.tD.clear()
                #self.dataVal.clear()
                #self.dataErr.clear()
                
                #for i, mgr in enumerate(self.managers):
                    #t = mgr.tD
                    #self.tD.add(t, data[i].size())
    
                    #self.dataVal = pg.cat(self.dataVal, data[i](mgr.dataToken()))
            
                    #if mgr.errIsAbsolute:
                        #self.dataErr = pg.cat(self.dataErr, data[i]('err')/data[i](mgr.dataToken()))
                    #else:
                        #self.dataErr = pg.cat(self.dataErr, data[i]('err'))
                
                #self.data=data
                
                #self.inv.setTransData(self.tD)
                #self.inv.setTransModel(self.tM)
            #else:
                #raise BaseException("To few datacontainer given")
            
    
    #def setMesh(self, mesh):
        #self.fop.setMesh(mesh)
            
    #def invert(self, data, mesh, lam=20, limits=None):
        
        #self.setData(data)
        #self.setMesh(mesh)

        #if limits is not None:
            #if hasattr(self.tM, 'setLowerBound'):
                #if verbose:
                    #print('Lower limit set to', limits[0])
                #self.tM.setLowerBound(limits[0])
            #if hasattr(self.tM, 'setUpperBound'):
                #if verbose:
                    #print('Upper limit set to', limits[1])
                #self.tM.setUpperBound(limits[1])
                
        #nModel = self.fop.regionManager().parameterCount()
        
        #self.inv.setData(self.dataVals)
        #self.inv.setRelativeError(self.dataErr)
        
        #startModel = pg.RVector(nModel, 0.0)
        
        #for i in range(len(self.managers)):
            #startModel += pg.RVector(nModel,
                                    #pg.median(self.trans[i].inv(self.managers[i].createApparentData(self.data[i]))))
        #startModel /= len(self.managers)
        
    
        #self.inv.setModel(startModel)
        #self.inv.setLambda(lam)

        #model = self.inv.run()
        #self.inv.echoStatus()
        #return model
        
    #def showModel(self, **showkwargs):
        #if len(showkwargs):
            #pg.show(fop.regionManager().paraDomain(), model, **showkwargs)

#class InvertPetro(InvertJointPetro):

    #def __init__(self, manager, trans,  **kwargs):
        #InvertJointPetro.__init__(self, [manager], [trans], **kwargs)
        
    #def invert(self, data, **kwargs):
        #InvertJointPetro.invert([data], **kwargs)
          

if __name__ == "__main__":
    pass

