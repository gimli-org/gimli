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
        super().__init__(**kwargs)
        self.fop = None

        if fop is not None:
            self.setForwardOperator(fop)

    def setForwardOperator(self, fop):
        self.fop = fop

    def setMesh(self, mesh):
        print("setMesh", mesh)
        if self.fop is not None:
            self.fop.setMesh(mesh)
            self.setRegionManager(self.fop.regionManagerRef())
        else:
            super().setMesh(mesh)

    def createStartModel(dataValues, **kwargs):
        """ Create Starting model.

        Create Starting model based on current data values and additional args.
        """
        raise Exception("Implement me in derived classes")

    # Mandatory api
    def setData(self, data):
        self.fop.setData(data)

#class Block1DModelling(Modelling):
    #def __init__(self, **kwargs):
        #super().__init__(**kwargs)

    #def createStartModel(dataValues, nLayers, **kwargs):
        #raise Exception("Implement me in derived classes", self)

#class MeshModelling(Modelling):

    #def __init__(self, **kwargs):
        #super(MeshModelling, self).__init__(self, **kwargs)

