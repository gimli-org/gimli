#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Class for managing first arrival travel time inversions"""
import numpy as np

import pygimli as pg
from pygimli.frameworks import MeshModelling
from pygimli.manager import MeshMethodManager

from . raplot import drawTravelTimeData


class TravelTimeDijkstraModelling(MeshModelling):
    def __init__(self, **kwargs):
        self.dijkstra = pg.TravelTimeDijkstraModelling()
        
        super(TravelTimeDijkstraModelling, self).__init__(**kwargs)

        # expose 
        self.jacobian = self.dijkstra.jacobian
        self.createJacobian = self.dijkstra.createJacobian

        self.setJacobian(self.dijkstra.jacobian())
 
    def regionManagerRef(self):
        # necessary because core dijkstra use its own RM
        return self.dijkstra.regionManagerRef()
     
    def setMesh(self, mesh, ignoreRegionManager=False):
        """
        """
        # necessary because core dijkstra use its own mesh management
        self.dijkstra.setMesh(mesh)
        super(TravelTimeDijkstraModelling, self).setMesh(mesh, ignoreRegionManager)

    def setData(self, data):
        self.dijkstra.setData(data)
        
    def createStartModel(self, t):
        pg.p(t)
        return self.dijkstra.createDefaultStartModel()
        
    def response(self, par):
        return self.dijkstra.response(par)

    def drawData(self, ax, data, err=None, label=None):
        """
        """
        print("showData", ax, data, err, label)
        return self.showVA(data, ax=ax)

    def showVA(self, data, t=None, name='va', pseudosection=False,
               squeeze=True, full=True, ax=None):
        """Show apparent velocity as image plot.

        TODO showXXX commands need to return ax and cbar .. if there is one

        """
        if data is None:
            data = self.dataContainer

        if ax is None:
            fig, ax = plt.subplots()

        if t is None:
            t = data('t')

        px = pg.x(data.sensorPositions())
        gx = np.array([px[int(g)] for g in data("g")])
        sx = np.array([px[int(s)] for s in data("s")])
        offset = self.getOffset(data=data, full=full)
        va = offset / t

        if pseudosection:
            midpoint = (gx + sx) / 2
            pg.mplviewer.dataview.plotVecMatrix(midpoint, offset, va, squeeze=True, ax=ax,
                          label='Apparent slowness [s/m]')
        else:
            pg.mplviewer.dataview.plotVecMatrix(gx, sx, va, squeeze=squeeze, ax=ax,
                          label='Apparent velocity [m/s]')
        return va

    def getOffset(self, data=None, full=False):
        """Return vector of offsets (in m) between shot and receiver."""

        if data is None:
            data = self.dataContainer

        if full:
            pos = data.sensorPositions()
            s, g = data('s'), data('g')
            nd = data.size()
            off = [pos[int(s[i])].distance(pos[int(g[i])]) for i in range(nd)]
            return np.absolute(off)
        else:
            px = pg.x(data.sensorPositions())
            gx = np.array([px[int(g)] for g in data("g")])
            sx = np.array([px[int(s)] for s in data("s")])
            return np.absolute(gx - sx)

    def getMidpoint(self, data=None):
        """Return vector of offsets (in m) between shot and receiver."""
        if data is None:
            data = self.dataContainer

        px = pg.x(data.sensorPositions())
        gx = np.array([px[int(g)] for g in data("g")])
        sx = np.array([px[int(s)] for s in data("s")])
        return (gx + sx) / 2


class TravelTimeManager(MeshMethodManager):
    """Manager for refraction seismics (traveltime tomography)

    TODO Document main members and use default MethodeManager interface
    e.g., self.inv, self.fop, self.paraDomain, self.mesh, self.data
    """

    def __init__(self, **kwargs):
        """Init function with optional data load"""
        self._useFMM = False

        super(TravelTimeManager, self).__init__(**kwargs)

        self._dataToken = 't'



        # check if needed

        self.figs = {}
        self.axs = {}

        self.errIsAbsolute = True

        # should be forwarded so it can be accessed from outside
        self.mesh = None
        self.poly = None
        self.error = None
        self.velocity = None
        self.response = None
        self.__dict__.update(**kwargs)
        # self.start = []
        self.pd = None



        data = kwargs.pop('data', None)
        if isinstance(data, str):
            self.loadData(data)
        elif isinstance(data, pg.DataContainer):
            self.setDataContainer(data)
            self.basename = kwargs.pop('name', 'new')

#        if self.dataContainer is not None:
#            self.createMesh()
#        self.fop = self.createFOP(verbose=self.verbose)
#        self.inv = self.createInv(self.fop,
#                                  verbose=self.verbose, doSave=self.doSave)

#    def __str__(self):  # no need to overwrite with identical content
#        """string representation of the class"""
#        return self.__repr__()
#
    def __repr__(self):  # to be moved to Mesh/Data Method manager
        """String representation of the class for the print function"""
        out = type(self).__name__ + " object"
        if hasattr(self, 'dataContainer'):
            out += "\n" + self.dataContainer.__str__()
        if hasattr(self, 'mesh'):
            out += "\n" + self.mesh.__str__()
        return out

    @property
    def useFMM(self):
        return self.__useFMM

    @useFMM.setter
    def useFMM(self, u):
        self.__useFMM = u
        self.initForwardOperator()

    def relErrVals(self, data):
        """Return pure data values from a given DataContainer."""
        return pg.abs(data('err')) / pg.abs(self.dataVals(data))

    def createForwardOperator(self, **kwargs):
        """Create default forward operator for Traveltime modelling.

        Your want your Manager use a special forward operator use
        Set self.useFMM = True to forces Fast Marching Method,
        otherwise Dijkstra is used.
        """
        useFMM = kwargs.pop('useFMM', False)

        if self._useFMM or useFMM:
            from .FMModelling import TravelTimeFMM
            fop = TravelTimeFMM(**kwargs)
        else:
            fop = TravelTimeDijkstraModelling(**kwargs)

        return fop

    def simulate(self, mesh, slowness, scheme, **kwargs):
        """
        Simulate an Traveltime measurement.

        Perform the forward task for a given mesh,
        a slowness distribution (per cell) and return data
        (Traveltime) for a measurement scheme.
        This is a static method since it does not interfere with the Managers
        inversion approaches.

        Parameters
        ----------
        mesh : :gimliapi:`GIMLI::Mesh`
            Mesh to calculate for.

        slowness : array(mesh.cellCount()) | array(N, mesh.cellCount())
            slowness distribution for the given mesh cells can be:

            * a single array of len mesh.cellCount()
            * a matrix of N slowness distributions of len mesh.cellCount()
            * a res map as [[marker0, res0], [marker1, res1], ...]

        scheme : :gimliapi:`GIMLI::DataContainer`
            data measurement scheme

        **kwargs :
            * noisify : add normal distributed noise based on scheme('err')
                IMPLEMENTME

        Returns
        -------
        t : array(N, data.size()) | DataContainer
            The resulting simulated travel time values.
            Either one column array or matrix in case of slowness matrix.
            A DataContainer is return if noisify set to True.

        """

        fop = self.createForwardOperator()

        fop.setData(scheme)
        fop.setMesh(mesh, ignoreRegionManager=True)

        if len(slowness) == mesh.cellCount():
            if max(slowness) > 1.:
                print('Warning: slowness values larger than 1 (' +
                      str(max(slowness)) + ').. assuming that are velocity '
                      'values .. building reciprocity')
                t = fop.response(1./slowness)
            else:
                t = fop.response(slowness)
        else:
            print(mesh)
            print("slowness: ", slowness)
            raise BaseException("Simulate called with wrong slowness array.")

        ret = pg.DataContainer(scheme)
        ret.set('t', t)

        noiseLevel = kwargs.pop('noiseLevel', 0)

        if noiseLevel > 0:
            if not ret.allNonZero('err'):
                ret.set('t', t)
                ret.set('err', pg.physics.Refraction.estimateError(
                        ret, absoluteError=kwargs.pop('noiseAbs', 1e-4),
                        relativeError=noiseLevel))

            if self.verbose:
                print("Data error estimates (min:max) ",
                      min(ret('err')), ":", max(ret('err')))

            t += pg.randn(ret.size()) * ret('err')
            ret.set('t', t)

        if kwargs.pop('returnArray', False):
            return t

        return ret


    def invert(self, data=None, err=None, mesh=None, **kwargs):
        """Invert measured data.
        """

        #ensure data and error sizes here
        dataVal = None
        if isinstance(data, pg.DataContainer):
            self.fop.setDataBasis(dataContainer=data)

            dataVals = self.dataVals(data)
            errVals = self.relErrVals(data)

        else:
            dataVal = data

        self.transData = pg.RTransLog()
        self.inv.transData = self.transData

        slowness = super(TravelTimeManager, self).invert(dataVals=dataVals,
                                                         errVals=errVals,
                                                         mesh=mesh,
                                                         **kwargs)
        self.model = 1./slowness
        return self.model

    def loadData(self, filename, **kwargs):
        """"""
        dat = pg.DataContainer(filename, 's g')
        self.fop.setDataContainer(dat)
        return dat

if __name__ == '__main__':
    pg.wait()
