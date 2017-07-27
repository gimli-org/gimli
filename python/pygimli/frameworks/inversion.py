# -*- coding: utf-8 -*-
"""pyGIMLi - Inversion Frameworks.

These are basic inversion frameworks that usually need a forward operator to run.
"""
import numpy as np

import pygimli as pg

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

        self.mod = self.inv.run()
        self.mod = self.mod(self.fop.regionManager().paraDomain().cellMarkers())
        return self.mod
