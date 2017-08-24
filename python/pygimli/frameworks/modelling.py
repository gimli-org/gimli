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
        super(Modelling, self).__init__(**kwargs)
        self.__regionProperties = {}
        self.__transModel = pg.RTransLog()

        self.fop = None

        if fop is not None:
            self.setForwardOperator(fop)

    @property
    def transModel(self):
        self._applyRegionProperties()
        if self.regionManager().haveLocalTrans():
            return self.regionManager().transModel()
        return self.__transModel

    def setForwardOperator(self, fop):
        self.fop = fop

    def setMesh(self, mesh):
        if self.fop is not None:
            self.fop.setMesh(mesh)
            self.setRegionManager(self.fop.regionManagerRef())
        else:
            super(Modelling, self).setMesh(mesh)

    def setRegionProperties(self, region, startModel=None, limits=None, trans=None):
        """
        """
        if region not in self.__regionProperties:
            self.__regionProperties[region] = {'startModel': 0,
                                               'limits': [0, 0],
                                               'trans': 'Log',
                                              }

        if startModel is not None:
            self.__regionProperties[region]['startModel'] = startModel
        if limits is not None:
            self.__regionProperties[region]['limits'] = limits
        if trans is not None:
            self.__regionProperties[region]['trans'] = trans

    def _applyRegionProperties(self):
        """
        """
        for rID, vals in self.__regionProperties.items():

            self.regionManager().region(rID).setStartModel(vals['startModel'])

            self.regionManager().region(rID).setModelTransStr_(vals['trans'])

            if vals['limits'][0] > 0:
                self.regionManager().region(rID).setLowerBound(vals['limits'][0])
            if vals['limits'][1] > 0:
                self.regionManager().region(rID).setUpperBound(vals['limits'][1])


    def createStartModel(dataValues, **kwargs):
        """ Create Starting model.

        Create Starting model based on current data values and additional args.
        """
        raise Exception("Implement me in derived classes")

    # Mandatory api
    def setData(self, data):
        #raise Exception("Needed? Implement me in derived classes")
        self.fop.setData(data)

    def estimateError(self, data, **kwargs):
        """Create
        """
        raise Exception("Needed?? Implement me in derived classes")
        #data = data * (pg.randn(len(data)) * errPerc / 100. + 1.)
        #return data

    def setDataBasis(self, **kwargs):
        """Set Data basis, e.g., DataContainer, times, coordinates."""
        data = kwargs.pop('data', None)
        if isinstance(data, pg.DataContainer) and self.fop is not None:
            self.fop.setData(data)


class Block1DModelling(Modelling):
    """
    """
    def __init__(self, nBlocks=1, **kwargs):
        super(Block1DModelling, self).__init__(**kwargs)
        self.__withMultiThread = True
        self.__nBlocks = nBlocks

    def setLayers(self, nLayers):

        if nLayers < 2:
            raise Exception("Number of layers need to be at least 2")

        mesh = pg.createMesh1DBlock(nLayers, self.__nBlocks)
        self.setMesh(mesh)
        #self.setStartModel(pg.RVector(0))

        for i in range(self.__nBlocks + 1):
            self.setRegionProperties(i, trans='log')

        if self.__withMultiThread:
            self.setMultiThreadJacobian(2*nLayers - 1)


        #raise Exception("Implement me in derived classes", self)

#class MeshModelling(Modelling):

    #def __init__(self, **kwargs):
        #super(MeshModelling, self).__init__(self, **kwargs)

