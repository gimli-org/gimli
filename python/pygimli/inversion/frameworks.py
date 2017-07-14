# -*- coding: utf-8 -*-
"""pyGIMLi - Inversion Frameworks.

These are basic inversion frameworks that usually need a forward operator to run.
"""
import numpy as np

import pygimli as pg

class Modelling(pg.ModellingBase):
    """Abstract Forward Operator.

    Abstract Forward Operator that use one or more different ModellingBase classes.
    Can be seen as some kind of proxy Forward Operator.

    """
    def __init__(self, **kwargs):

        fop = kwargs.pop('fop', None)
        super(Modelling, self).__init__(self, **kwargs)

        if fop is not None:
            self.setForwardOperator(fop)

    def setForwardOperator(self, fop):
        self.fop = fop

    def setMesh(self, mesh):
        if self.fop != None:
            self.fop.setMesh(mesh)
            self.setRegionManager(self.fop.regionManagerRef())

    def setData(self, data):
        self.fop.setData(data)


class MeshInversion():
    def __init__(self, **kwargs):
        self.fop = None
        self.dataVals = None
        self.dataErrs = None
        self.tM = None
        self.tD = None
        self.inv = pg.Inversion(**kwargs)

    def setForwardOperator(self, fop):
        self.fop = fop
        self.inv.setForwardOperator(fop)

    def setMesh(self, mesh):
        self.fop.setMesh(mesh)

    def setData(self, data):
        self.fop.setData(data)
        self.dataVals = None
        self.dataErrs = None
        raise Exception('implementme')

    def invert(self, data=None, mesh=None, lam=20, **kwargs):
        if data is not None:
            self.setData(data)

        if mesh is not None:
            self.setMesh(mesh)

        startModel = kwargs.pop('startModel', None)
        nModel = self.fop.regionManager().parameterCount()
        startModel = pg.Vector(nModel, startModel)
        self.fop.setStartModel(startModel)

        zWeight = kwargs.pop('zWeight', 1.0)
        self.fop.regionManager().setZWeight(zWeight)
        self.inv.setData(self.dataVals)
        self.inv.setRelativeError(self.dataErrs)
        self.inv.setLambda(lam)
        print('#'*100)
        print(self.fop.jac, self.fop.jac.rows(), self.fop.jac.rows())
        self.mod = self.inv.run()
        self.mod = self.mod(self.fop.regionManager().paraDomain().cellMarkers())
        return self.mod
