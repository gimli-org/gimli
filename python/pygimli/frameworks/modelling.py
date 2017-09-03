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
        self._regionProperties = {}
        self._transModel = pg.RTransLog()
        self.fop = None

        self.data = None # dataContainer

        if fop is not None:
            self.setForwardOperator(fop)

    @property
    def transModel(self):
        self._applyRegionProperties()
        if self.regionManager().haveLocalTrans():
            return self.regionManager().transModel()
        return self._transModel

    def regionManager(self):
        """
        """
        # init RM if necessary
        super(Modelling, self).regionManager()
        # set all local properties
        self._applyRegionProperties()
        return super(Modelling, self).regionManager()

    def setForwardOperator(self, fop):
        self.fop = fop

    def setMesh(self, mesh, ignoreRegionManager=False):

        if (not ignoreRegionManager):
            self.clearRegionProperties()

        if self.fop is not None:
            print("Modelling:setMesh", self.fop)
            self.fop.setMesh(mesh, ignoreRegionManager)

            if (not ignoreRegionManager):
                self.setRegionManager(self.fop.regionManagerRef())
        else:
            super(Modelling, self).setMesh(mesh, ignoreRegionManager)

    def clearRegionProperties(self):
        self._regionProperties = {}

    def setRegionProperties(self, region, startModel=None, limits=None, trans=None):
        """
        """
        if region not in self._regionProperties:
            self._regionProperties[region] = {'startModel': 0,
                                               'limits': [0, 0],
                                               'trans': 'Log',
                                              }

        if startModel is not None:
            self._regionProperties[region]['startModel'] = startModel
        if limits is not None:
            self._regionProperties[region]['limits'] = limits
        if trans is not None:
            self._regionProperties[region]['trans'] = trans

    def _applyRegionProperties(self):
        """
        """
        RM = super(Modelling, self).regionManager()
        for rID, vals in self._regionProperties.items():

            RM.region(rID).setStartModel(vals['startModel'])

            RM.region(rID).setModelTransStr_(vals['trans'])

            if vals['limits'][0] > 0:
                RM.region(rID).setLowerBound(vals['limits'][0])
            if vals['limits'][1] > 0:
                RM.region(rID).setUpperBound(vals['limits'][1])

    def createStartModel(self, dataValues, **kwargs):
        """ Create Starting model.

        Create Starting model based on current data values and additional args.
        """
        raise Exception("Implement me in derived classes")

    # Mandatory api
    def setData(self, data):
        #raise Exception("Needed? Implement me in derived classes")
        if self.fop is not None:
            self.fop.setData(data)

    def setDataContainer(self, data):
        if self.fop is not None:
            self.fop.setData(data)
        else:
            super(Modelling, self).setData(data)
            self.data = data

    def setDataBasis(self, **kwargs):
        """Set Data basis, e.g., DataContainer, times, coordinates."""
        data = kwargs.pop('dataContainer', None)
        if isinstance(data, pg.DataContainer):
            self.setDataContainer(data)

    def estimateError(self, data, **kwargs):
        """Create
        """
        raise Exception("Needed?? Implement me in derived classes")
        #data = data * (pg.randn(len(data)) * errPerc / 100. + 1.)
        #return data

    def drawModel(self, ax, model):
        """
        """
        print(ax, model)
        raise Exception("No yet implemented")

    def drawData(self, ax, data, err=None, label=None):
        """
        """
        print(ax, data, err, label)
        raise Exception("No yet implemented")


class Block1DModelling(Modelling):
    """
    """
    def __init__(self, nBlocks=1, **kwargs):
        super(Block1DModelling, self).__init__(**kwargs)
        self._withMultiThread = True
        self._nBlocks = nBlocks

    def setLayers(self, nLayers):
        if nLayers < 2:
            raise Exception("Number of layers need to be at least 2")

        mesh = pg.createMesh1DBlock(nLayers, self._nBlocks)
        self.setMesh(mesh)
        #self.setStartModel(pg.RVector(0))

        for i in range(self._nBlocks + 1):
            self.setRegionProperties(i, trans='log')

        if self._withMultiThread:
            self.setMultiThreadJacobian(2*nLayers - 1)

        self._applyRegionProperties()

    def drawModel(self, ax, model):
        pg.mplviewer.drawModel1D(ax=ax,
                                 model=model,
                                 plot='loglog',
                                 xlabel='Model parameter')
        return ax

    def drawData(self, ax, data, err=None, label=None):
        nData = len(data)
        yVals = range(nData)
        ax.loglog(data, yVals, 'rx-')
        if err is not None:
            ax.errorbar(data, yVals,
                        xerr=err*data,
                        linewidth=1, color='red', linestyle='-')

        ax.set_ylim(max(yVals), min(yVals))
        ax.set_xlabel('Data')
        ax.set_ylabel('Data Number')
        return ax


class MeshModelling(Modelling):
    """
    """
    def __init__(self, **kwargs):
        super(MeshModelling, self).__init__(**kwargs)

    @property
    def paraDomain(self):
        return self.regionManager().paraDomain()

    def setMesh(self, mesh, ignoreRegionManager=False):

        super(MeshModelling, self).setMesh(mesh, ignoreRegionManager)

    def drawModel(self, ax, model):
        pg.mplviewer.drawModel(ax=ax,
                               mesh=self.paraDomain,
                               data=model,
                               label='Model parameter')
        return ax


class PetroModelling(Modelling):
    """
    Combine petrophysical relation m(p) with modeling class f(p).

    Combine petrophysical relation m(p) with modeling class f(p) to invert
    for m (or any inversion transformation) instead of p.
    """
    def __init__(self, fop, trans, **kwargs):
        """Save forward class and transformation, create Jacobian matrix."""
        mesh = kwargs.pop('mesh', None)

        super(PetroModelling, self).__init__(**kwargs)
        self.fop = fop      # class defining f(p)
        self.trans = trans  # class defining m(p)
        print("Petro_init:", self.fop)
        print(self.fop.regionManagerRef())
        #self.setRegionManager(self.fop.regionManagerRef())

        if mesh is not None:
            self.setMesh(mesh)

        self._jac = pg.MultRightMatrix(self.fop.jacobian())
        self.setJacobian(self._jac)

        #TODO global TransModel will break RegionConcept
        self._transModel = pg.RTransLogLU()

    def drawModel(self, ax, model):
        self.fop.drawModel(ax, model)

    def response(self, model):
        """Use inverse transformation to get p(m) and compute response."""
        tModel = self.trans(model)
        ret = self.fop.response(tModel)
        return ret

    def createJacobian(self, model):
        """Fill the individual jacobian matrices."""
        self.fop.createJacobian(self.trans(model))
        self._jac.r = self.trans.deriv(model)  # set inner derivative

    def setDataBasis(self, **kwargs):
        self.fop.setDataBasis(**kwargs)

    #def setMesh(self, mesh):
        #"""TODO."""
        #if mesh is None and self.fop.mesh() is None:
            #raise BaseException("Please provide a mesh for "
                                #"this forward operator")

        #if mesh is not None:
            #self.fop.setMesh(mesh)
            #self.fop.createRefinedForwardMesh(refine=False)

        ## self.setMesh(f.mesh(), ignoreRegionManager=True) # not really nessary
        #self.setRegionManager(self.fop.regionManagerRef())


