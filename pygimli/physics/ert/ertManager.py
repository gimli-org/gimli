#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""2.5D non-optimized totalfield forward operator for ERT.

Please use the BERT package for more advanced forward operator
https://gitlab.com/resistivity-net/bert
"""

import os.path
import numpy as np
import matplotlib.pyplot as plt

import pygimli as pg
from pygimli.frameworks import MeshMethodManager
from .ertModelling import ERTModelling, ERTModellingReference
from .ert import createInversionMesh, createGeometricFactors, estimateError
from pygimli.utils import getSavePath


class ERTManager(MeshMethodManager):
    """ERT Manager.

    Method Manager for Electrical Resistivity Tomography (ERT)

    Todo
    ----
        * 3d
        * 3dtopo
        * complex on/off
        * closed geometry
        * transdim
        * singularity removal
        * ERT specific inversion options:
            * ...
    """

    def __init__(self, data=None, **kwargs):
        """Create ERT Manager instance.

        Parameters
        ----------
        data: :gimliapi:`GIMLI::DataContainerERT` | str
            You can initialize the Manager with data or give them a dataset
            when calling the inversion.

        Other Parameters
        ----------------
        * useBert: bool [True]
            Use Bert forward operator instead of the reference implementation.
        * sr: bool [True]
            Calculate with singularity removal technique.
            Recommended but needs the primary potential.
            For flat earth cases the primary potential will be calculated
            analytical. For domains with topography the primary potential
            will be calculated numerical using a p2 refined mesh or
            you provide primary potentials with setPrimPot.
        """
        self.useBert = kwargs.pop('useBert', True)
        self.sr = kwargs.pop('sr', True)

        super().__init__(data=data, **kwargs)
        self.inv.dataTrans = pg.trans.TransLogLU()

    def setSingularityRemoval(self, sr=True):
        """Turn singularity removal on or off."""
        self.reinitForwardOperator(sr=True)

    def createForwardOperator(self, **kwargs):
        """Create and choose forward operator."""
        verbose = kwargs.pop('verbose', False)
        self.useBert = kwargs.pop('useBert', self.useBert)
        self.sr = kwargs.pop('sr', self.sr)
        if self.useBert:
            pg.verbose('Create ERTModelling FOP')
            fop = ERTModelling(sr=self.sr, verbose=verbose)
        else:
            pg.verbose('Create ERTModellingReference FOP')
            fop = ERTModellingReference(**kwargs)

        return fop

    def load(self, fileName):
        """Load ERT data.

        Forwarded to :py:mod:`pygimli.physics.ert.load`

        Parameters
        ----------
        fileName: str
            Filename for the data.

        Returns
        -------
        data: :gimliapi:`GIMLI::DataContainerERT`
        """
        self.data = pg.physics.ert.load(fileName)
        return self.data

    def createMesh(self, data=None, **kwargs):
        """Create default inversion mesh.

        Forwarded to :py:mod:`pygimli.physics.ert.createInversionMesh`
        """
        d = data or self.data

        if d is None:
            pg.critical('Please provide a data file for mesh generation')

        return createInversionMesh(d, **kwargs)

    def setPrimPot(self, pot):
        """Set primary potential from external is not supported anymore."""
        pg.critical("Not implemented.")

    def simulate(self, mesh, scheme, res, **kwargs):
        """Simulate an ERT measurement.

        Perform the forward task for a given mesh, resistivity distribution &
        measuring scheme and return data (apparent resistivity) or potentials.

        For complex resistivity, the apparent resistivities is complex as well.

        The forward operator itself only calculates potential values for the
        electrodes in the given data scheme.
        To calculate apparent resistivities, geometric factors (k) are needed.
        If there are no values k in the DataContainerERT scheme, the function
        tries to calculate them, either analytically or numerically by using a
        p2-refined version of the given mesh.

        TODO
        ----
        * 2D + Complex + SR

        Args
        ----
        mesh : :gimliapi:`GIMLI::Mesh`
            2D or 3D Mesh to calculate for.

        res : float, array(mesh.cellCount()) | array(N, mesh.cellCount()) |
              list
            Resistivity distribution for the given mesh cells can be:
            . float for homogeneous resistivity (e.g. 1.0)
            . single array of length mesh.cellCount()
            . matrix of N resistivity distributions of length mesh.cellCount()
            . resistivity map as [[regionMarker0, res0],
                                  [regionMarker0, res1], ...]

        scheme : :gimliapi:`GIMLI::DataContainerERT`
            Data measurement scheme.

        Keyword Args
        ------------
        verbose: bool[False]
            Be verbose. Will override class settings.
        calcOnly: bool [False]
            Use fop.calculate instead of fop.response. Useful if you want
            to force the calculation of impedances for homogeneous models.
            No noise handling. Solution is put as token 'u' in the returned
            DataContainerERT.
        noiseLevel: float [0.0]
            add normally distributed noise based on
            scheme['err'] or on noiseLevel if error>0 is not contained
        noiseAbs: float [0.0]
            Absolute voltage error in V
        returnArray: bool [False]
            Returns an array of apparent resistivities instead of
            a DataContainerERT
        returnFields: bool [False]
            Returns a matrix of all potential values (per mesh nodes)
            for each injection electrodes.

        Returns
        -------
        DataContainerERT | array(data.size()) | array(N, data.size()) |
        array(N, mesh.nodeCount()):
            Data container with resulting apparent resistivity data and
            errors (if noiseLevel or noiseAbs is set).
            Optional returns a Matrix of rhoa values
            (for returnArray==True forces noiseLevel=0).
            In case of a complex valued resistivity model, phase values are
            returned in the DataContainerERT (see example below), or as an
            additionally returned array.

        Examples
        --------
        # >>> from pygimli.physics import ert
        # >>> import pygimli as pg
        # >>> import pygimli.meshtools as mt
        # >>> world = mt.createWorld(start=[-50, 0], end=[50, -50],
        # ...                        layers=[-1, -5], worldMarker=True)
        # >>> scheme = ert.createData(
        # ...                     elecs=pg.utils.grange(start=-10, end=10, n=21),
        # ...                     schemeName='dd')
        # >>> for pos in scheme.sensorPositions():
        # ...     _= world.createNode(pos)
        # ...     _= world.createNode(pos + [0.0, -0.1])
        # >>> mesh = mt.createMesh(world, quality=34)
        # >>> rhomap = [
        # ...    [1, 100. + 0j],
        # ...    [2, 50. + 0j],
        # ...    [3, 10.+ 0j],
        # ... ]
        # >>> data = ert.simulate(mesh, res=rhomap, scheme=scheme, verbose=True)
        # >>> rhoa = data.get('rhoa').array()
        # >>> phia = data.get('phia').array()
        """
        verbose = kwargs.pop('verbose', self.verbose)
        calcOnly = kwargs.pop('calcOnly', False)
        returnFields = kwargs.pop("returnFields", False)
        returnArray = kwargs.pop('returnArray', False)
        noiseLevel = kwargs.pop('noiseLevel', 0.0)
        noiseAbs = kwargs.pop('noiseAbs', 1e-4)
        seed = kwargs.pop('seed', None)
        sr = kwargs.pop('sr', self.sr)

        # segfaults with self.fop (test & fix)
        fop = self.createForwardOperator(useBert=self.useBert,
                                         sr=sr, verbose=verbose)
        fop.data = scheme
        fop.setMesh(mesh, ignoreRegionManager=True)

        rhoa = None
        phia = None

        isArrayData = False
        # parse the given res into mesh-cell-sized array
        if isinstance(res, int) or isinstance(res, float):
            res = np.ones(mesh.cellCount()) * float(res)
        elif isinstance(res, complex):
            res = np.ones(mesh.cellCount()) * res
        elif hasattr(res[0], '__iter__'):  # ndim == 2
            if len(res[0]) == 2:  # res seems to be a res map
                # check if there are markers in the mesh that are not defined
                # the rhomap. better signal here before it results in errors
                meshMarkers = list(set(mesh.cellMarkers()))
                mapMarkers = [m[0] for m in res]
                if any([mark not in mapMarkers for mark in meshMarkers]):
                    left = [m for m in meshMarkers if m not in mapMarkers]
                    pg.critical("Mesh contains markers without assigned "
                                "resistivities {}. Please fix given "
                                "rhomap.".format(left))
                res = pg.solver.parseArgToArray(res, mesh.cellCount(), mesh)
            else:  # probably nData x nCells array
                # better check for array data here
                isArrayData = True

        if isinstance(res[0], np.complex) or isinstance(res, pg.CVector):
            pg.info("Complex resistivity values found.")
            fop.setComplex(True)
        else:
            fop.setComplex(False)

        if not scheme.allNonZero('k') and not calcOnly:
            if verbose:
                pg.info('Calculate geometric factors.')
            scheme.set('k', fop.calcGeometricFactor(scheme))

        ret = pg.DataContainerERT(scheme)
        # just to be sure that we don't work with artifacts
        ret['u'] *= 0.0
        ret['i'] *= 0.0
        ret['r'] *= 0.0

        if isArrayData:
            rhoa = np.zeros((len(res), scheme.size()))
            for i, r in enumerate(res):
                rhoa[i] = fop.response(r)
                if verbose:
                    print(i, "/", len(res), " : ", pg.dur(), "s",
                          "min r:", min(r), "max r:", max(r),
                          "min r_a:", min(rhoa[i]), "max r_a:", max(rhoa[i]))
        else:  # res is single resistivity array
            if len(res) == mesh.cellCount():

                if calcOnly:
                    fop.mapERTModel(res, 0)

                    dMap = pg.core.DataMap()
                    fop.calculate(dMap)
                    if fop.complex():
                        pg.critical('Implement me')
                    else:
                        ret["u"] = dMap.data(scheme)
                        ret["i"] = np.ones(ret.size())

                    if returnFields:
                        return pg.Matrix(fop.solution())
                    return ret
                else:
                    if fop.complex():
                        res = pg.utils.squeezeComplex(res)

                    resp = fop.response(res)

                    if fop.complex():
                        rhoa, phia = pg.utils.toPolar(resp)
                    else:
                        rhoa = resp
            else:
                print(mesh)
                print("res: ", res)
                raise BaseException(
                    "Simulate called with wrong resistivity array.")

        if not isArrayData:
            ret['rhoa'] = rhoa

            if phia is not None:
                ret.set('phia', phia)
        else:
            ret.set('rhoa', rhoa[0])
            if phia is not None:
                ret.set('phia', phia[0])

        if returnFields:
            return pg.Matrix(fop.solution())

        if noiseLevel > 0:  # if errors in data noiseLevel=1 just triggers
            if not ret.allNonZero('err'):
                # 1A  and #100µV
                ret.set('err', self.estimateError(ret,
                                                  relativeError=noiseLevel,
                                                  absoluteUError=noiseAbs,
                                                  absoluteCurrent=1))
                print("Data error estimate (min:max) ",
                      min(ret('err')), ":", max(ret('err')))

            rhoa *= 1. + pg.randn(ret.size(), seed=seed) * ret('err')
            ret.set('rhoa', rhoa)

            ipError = None
            if phia is not None:
                if scheme.allNonZero('iperr'):
                    ipError = scheme('iperr')
                else:
                    # np.abs(self.data("phia") +TOLERANCE) * 1e-4absoluteError
                    if noiseLevel > 0.5:
                        noiseLevel /= 100.

                    if 'phiErr' in kwargs:
                        ipError = np.ones(ret.size()) * kwargs.pop('phiErr') \
                            / 1000
                    else:
                        ipError = abs(ret["phia"]) * noiseLevel

                    if verbose:
                        print("Data IP abs error estimate (min:max) ",
                              min(ipError), ":", max(ipError))

                phia += pg.randn(ret.size(), seed=seed) * ipError
                ret['iperr'] = ipError
                ret['phia'] = phia

        # check what needs to be setup and returned

        if returnArray:
            if phia is not None:
                return rhoa, phia
            else:
                return rhoa

        return ret

    def checkData(self, data):
        """Return data from container.

        THINKABOUT: Data will be changed, or should the manager keep a copy?
        """
        if isinstance(data, pg.DataContainer):

            if not data.allNonZero('k'):
                pg.warn("Data file contains no geometric factors (token='k').")
                data['k'] = createGeometricFactors(data, verbose=True)

            if self.fop.complex():
                if not data.haveData('rhoa'):
                    pg.critical('Datacontainer have no "rhoa" values.')
                if not data.haveData('ip'):
                    pg.critical('Datacontainer have no "ip" values.')

                # pg.warn('check sign of phases')
                rhoa = data['rhoa']
                phia = -data['ip']/1000  # 'ip' is defined for neg mrad.
                # we should think about some 'phia' in rad

                return pg.utils.squeezeComplex(pg.utils.toComplex(rhoa, phia))

            else:
                if not data.haveData('rhoa'):

                    if data.allNonZero('r'):
                        pg.info("Creating apparent resistivies from "
                                "impedences rhoa = r * k")
                        data['rhoa'] = data['r'] * data['k']
                    elif data.allNonZero('u') and data.allNonZero('i'):
                        pg.info("Creating apparent resistivies from "
                                "voltage and currrent rhoa = u/i * k")
                        data['rhoa'] = data['u']/data['i'] * data['k']
                    else:
                        pg.critical("Datacontainer have neither: "
                                    "apparent resistivies 'rhoa', "
                                    "or impedances 'r', "
                                    "or voltage 'u' along with current 'i'.")

                return data['rhoa']

        return data

    def checkErrors(self, err, dataVals):
        """Return relative error.

        Default we assume 'err' are relative vales.
        """
        if isinstance(err, pg.DataContainer):
            rae = None

            if not err.allNonZero('err'):
                pg.warn("Datacontainer have no 'err' values. "
                        "Fallback of 1mV + 3% using "
                        "ERTManager.estimateError(...) ")
                rae = self.estimateError(err, absoluteError=0.001,
                                         relativeError=0.03)
            else:
                rae = err['err']

            if self.fop.complex():

                ipe = None

                if err.haveData('iperr'):
                    amp, phi = pg.utils.toPolar(dataVals)
                    # assuming ipErr are absolute dPhi in mrad
                    ipe = err['iperr'] / abs((phi*1000))
                else:
                    pg.warn("Datacontainer have no 'iperr' values. "
                            "Fallback set to 0.01")
                    ipe = np.ones(err.size()) * 0.01

                # pg._y("err", min(rae), max(rae), rae)
                # pg._y("iperr", min(ipe), max(ipe), ipe)
                return pg.cat(rae, ipe)

        return rae  # not set if err is no DataContainer (else missing)

    def estimateError(self, data=None, **kwargs):
        """Estimate error composed of an absolute and a relative part.

        Parameters
        ----------
        absoluteError : float [0.001]
            Absolute data error in Ohm m. Need 'rhoa' values in data.

        relativeError : float [0.03]
            relative error level in %/100

        absoluteUError : float [0.001]
            Absolute potential error in V. Need 'u' values in data. Or
            calculate them from 'rhoa', 'k' and absoluteCurrent if no 'i'
            is given

        absoluteCurrent : float [0.1]
            Current level in A for reconstruction for absolute potential V

        Returns
        -------
        error : Array
        """
        if data is None:  #
            error = estimateError(self.data, **kwargs)
            self.data["err"] = error
        else:  # the old way: better use ert.estimateError directly
            error = estimateError(data, **kwargs)

        return error

    def coverage(self):
        """Coverage vector considering the logarithmic transformation."""
        covTrans = pg.core.coverageDCtrans(self.fop.jacobian(),
                                           1.0 / self.inv.response,
                                           1.0 / self.inv.model)

        paramSizes = np.zeros(len(self.inv.model))
        for c in self.fop.paraDomain.cells():
            paramSizes[c.marker()] += c.size()

        return np.log10(covTrans / paramSizes)

    def standardizedCoverage(self, threshhold=0.01):
        """Return standardized coverage vector (0|1) using thresholding."""
        return 1.0*(abs(self.coverage()) > threshhold)

    def saveResult(self, folder=None, size=(16, 10), **kwargs):
        """Save all results in the specified folder.

        Saved items are:
            Inverted profile
            Resistivity vector
            Coverage vector
            Standardized coverage vector
            Mesh (bms and vtk with results)
        """
        subfolder = self.__class__.__name__
        path = getSavePath(folder, subfolder)

        pg.info('Saving resistivity data to: {}'.format(path))

        np.savetxt(path + '/resistivity.vector',
                   self.model)
        np.savetxt(path + '/resistivity-cov.vector',
                   self.coverage())
        np.savetxt(path + '/resistivity-scov.vector',
                   self.standardizedCoverage())

        m = pg.Mesh(self.paraDomain)
        m['Resistivity'] = self.paraModel(self.model)
        m['Resistivity (log10)'] = np.log10(m['Resistivity'])
        m['Coverage'] = self.coverage()
        m['S_Coverage'] = self.standardizedCoverage()
        m.exportVTK(os.path.join(path, 'resistivity'))
        m.saveBinaryV2(os.path.join(path, 'resistivity-pd'))
        self.fop.mesh().save(os.path.join(path, 'resistivity-mesh'))

        if self.paraDomain.dim() == 2:
            fig, ax = plt.subplots(figsize=size)
            self.showResult(ax=ax, coverage=self.coverage(), **kwargs)
            fig.savefig(path + '/resistivity.pdf', bbox_inches="tight")
            return path, fig, ax
        return path


if __name__ == "__main__":
    pass
