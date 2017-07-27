# -*- coding: utf-8 -*-
"""pyGIMLi - Inversion Frameworks.

These are basic modelling proxies.
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

