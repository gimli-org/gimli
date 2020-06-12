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
from pygimli.frameworks import MeshModelling, MeshMethodManager
from pygimli.utils import getSavePath
from .visualization import showERTData

from pygimli import pf


def simulate(mesh, scheme, res, sr=True, useBert=True,
             verbose=False, **kwargs):
    """ERT forward calculation.

    Convenience function to use the ERT modelling operator
    if you like static functions.

    See :py:mod:`pygimli.ert.ERTManager.simulate` for description
    of the arguments.

    Parameters
    ----------
    mesh: :gimliapi:`GIMLI::Mesh` | str
        Modelling domain. Mesh can be a file name here.
    scheme: :gimliapi:`GIMLI::DataContainerERT` | str
        Data configuration. Scheme can be a file name here.
    res: see :py:mod:`pygimli.ert.ERTManager.simulate`
        Resistivity distribution.
    sr: bool [True]
        Use singularity removal technique.
    useBert: bool [True]
        Use Bert forward operator instead of the reference implementation.
    **kwargs:
        Forwarded to :py:mod:`pygimli.ert.ERTManager.simulate`
    """
    ert = ERTManager(useBert=useBert, sr=sr, verbose=verbose)

    if isinstance(mesh, str):
        mesh = pg.load(mesh)

    if isinstance(scheme, str):
        scheme = pg.physics.ert.load(scheme)

    return ert.simulate(mesh=mesh, res=res, scheme=scheme,
                        verbose=verbose, **kwargs)


@pg.cache
def createGeometricFactors(scheme, numerical=None, mesh=None,
                           h2=True, p2=True, verbose=False):
    """Create geometric factors for a given data scheme.

    Create geometric factors for a data scheme with and without topography.
    Calculation will be done analytical (only for half space geometry)
    or numerical.

    This function caches the result depending on scheme, mesh and pg.version()

    Parameters
    ----------
    scheme: :gimliapi:`GIMLI::DataContainerERT`
        Datacontainer of the scheme.
    numerical: bool | None [False]
        If numerical is None, False is assumed, we try to guess topography
        and warn if we think we found them.
        If set to True or False, numerical calculation will used respectively.
    mesh: :gimliapi:`GIMLI::Mesh` | str
        Mesh for numerical calculation. If not given, analytical geometric
        factors for halfspace earth are guessed or a default mesh will be
        created. The mesh will be h and p refined. If given topo is set to
        True. If the numerical effort is to high or the accuracy to low
        you should consider to calculate the factors manual.
    h2: bool [True]
        Default refinement to achieve high accuracy calculation
    p2: bool [True]
        Default refinement to achieve high accuracy calculation
    verbose: bool
        Give some output.
    """
    if numerical is None:
        numerical = False
        if min(pg.z(scheme)) != max(pg.z(scheme)):
            verbose = True
            pg.warn('Sensor z-coordinates not equal. Is there topography?')

    if numerical is False and mesh is None:
        if verbose:
            pg.info('Calculate analytical flat earth geometric factors.')

        return pg.core.geometricFactors(scheme, forceFlatEarth=True)

    if mesh is None:
        mesh = createInversionMesh(scheme)

    if verbose:
        pg.info('mesh', mesh)

    if h2 is True:
        m = mesh.createH2()
        if verbose:
            pg.info('h2 refine', m)

    if p2 is True:
        m = m.createP2()
        if verbose:
            pg.info('p2 refine', m)

    if verbose:
        pg.info('Calculate numerical geometric factors.')

    d = simulate(m, res=1.0, scheme=scheme, sr=False, useBert=True,
                 calcOnly=True, verbose=True)
    return 1./d['u']


def createInversionMesh(data, **kwargs):
    """Create default mesh for ERT inversion.

    Parameters
    ----------
    data: :gimliapi:`GIMLI::DataContainerERT`
        Data Container needs at least sensors to define the geometry of the
        mesh.

    Other Parameters
    ----------------
        Forwarded to :py:mod:`pygimli.meshtools.createParaMesh`

    Returns
    -------
    mesh: :gimliapi:`GIMLI::Mesh`
        Inversion mesh with default marker (1 for background,
        2 parametric domain)
    """
    mesh = pg.meshtools.createParaMesh(data.sensors(), **kwargs)
    return mesh


class ERTModellingBase(MeshModelling):
    """Modelling base class for ERT modelling."""

    def __init__(self, **kwargs):
        super(ERTModellingBase, self).__init__(**kwargs)

    def drawData(self, ax, data=None, **kwargs):
        """Draw data in given axe."""
        kwargs['label'] = kwargs.pop('label', pg.unit('res'))
        kwargs['cMap'] = kwargs.pop('cMap', pg.utils.cMap('res'))

        if hasattr(data, '__iter__'):
            vals = data
            data = self.data
        elif data is None:
            data = self.data

        vals = kwargs.pop('vals', data['rhoa'])

        return showERTData(data, vals=vals, ax=ax, **kwargs)

    def drawModel(self, ax, model, **kwargs):
        """Draw the para domain with option model values."""
        kwargs.setdefault('label', pg.unit('res'))
        kwargs.setdefault('cMap', pg.utils.cMap('res'))

        return super(ERTModellingBase, self).drawModel(ax=ax, model=model,
                                                       logScale=True,
                                                       **kwargs)


class ERTModelling(ERTModellingBase):
    """Forward operator for Electrical Resistivty Tomography.

    Note
    ----
    Convention for complex resistiviy inversion:
    We want to use logarithm transformation for the imaginary part of model
    so we need the startmodel to have positive imaginary parts.
    The sign is flipped back to physical correct assumption before we call
    the response function.
    The Jacobian is calculated with negative imaginary parts and will
    be a conjugated complex block matrix for further calulations.
    """

    def __init__(self, sr=True, verbose=False):
        super(ERTModelling, self).__init__()

        # don't use DC*fop or its regionmanager directly
        #
        self._core = None
        if sr:
            self._core = pg.core.DCSRMultiElectrodeModelling(verbose=verbose)
        else:
            self._core = pg.core.DCMultiElectrodeModelling(verbose=verbose)

        self._core.initJacobian()
        self.setJacobian(self._core.jacobian())

        # called from the ERTManager .. needed?
        self.solution = self._core.solution
        self.setComplex = self._core.setComplex
        self.complex = self._core.complex
        self.calculate = self._core.calculate
        self.calcGeometricFactor = self._core.calcGeometricFactor
        self.mapERTModel = self._core.mapERTModel

        self._conjImag = False  # the imaginary parts are flipped for log trans

    def setDefaultBackground(self):
        """Set the default background behaviour."""
        if self.complex():
            self.regionManager().addRegion(3, self._baseMesh, 2)

        regionIds = self.regionManager().regionIdxs()
        pg.info("Found {} regions.".format(len(regionIds)))
        if len(regionIds) > 1:
            bk = pg.sort(regionIds)[0]
            pg.info("Region with smallest marker ({0}) "
                    "set to background".format(bk))
            self.setRegionProperties(bk, background=True)

    def createStartModel(self, dataVals):
        """Create Starting model for ERT inversion."""
        if self.complex():
            dataC = pg.utils.toComplex(dataVals)
            nModel = self.regionManager().parameterCount() // 2
            smRe = np.ones(nModel) * np.median(np.median(dataC.real))
            smIm = np.ones(nModel) * np.median(np.median(dataC.imag))

            if min(smIm) < 0:
                # we want positive phase model
                sm = smRe - 1j * smIm
                pg.info("Model imaginary part being flipped to positive.")
                self._conjImag = True
            else:
                sm = smRe + 1j * smIm

            return pg.utils.squeezeComplex(sm)  # complex impedance
        else:
            return super(ERTModelling, self).createStartModel(dataVals)

    def flipImagPart(self, v):
        """Flip imaginary port (convention)."""
        z = pg.utils.toComplex(v)
        pg.warn('pre min/max={0} / {1} im: {2} / {3}'.format(
            pf(min(z.real)), pf(max(z.real)),
            pf(min(z.imag)), pf(max(z.imag))))

        v = pg.utils.squeezeComplex(pg.utils.toComplex(v),
                                    conj=self._conjImag)

        z = pg.utils.toComplex(v)
        pg.warn('pos min/max={0} / {1} im: {2} / {3}'.format(
            pf(min(z.real)), pf(max(z.real)),
            pf(min(z.imag)), pf(max(z.imag))))
        return v

    def response(self, mod):
        """Forward response (apparent resistivity)."""
        # ensure the mesh is initialized
        self.mesh()
        if self.complex() and self._conjImag:
            pg.warn('flip imaginary part for response calc')
            mod = self.flipImagPart(mod)

        resp = self._core.response(mod)

        if self.complex() and self._conjImag:
            pg.warn('backflip imaginary part after response calc')
            resp = self.flipImagPart(resp)

        return resp

    def createJacobian(self, mod):
        """Compute Jacobian matrix and store but not return."""
        # ensure the mesh is initialized
        self.mesh()
        if self.complex():
            if self._conjImag:
                pg.warn("Flipping imaginary part for jacobian calc")
                mod = self.flipImagPart(mod)

            self._core.createJacobian(mod)
            self._J = pg.utils.squeezeComplex(self._core.jacobian(),
                                              conj=self._conjImag
                                              )
            self.setJacobian(self._J)
            # pg._r("create Jacobian", self, self._J)
            return self._J

        return self._core.createJacobian(mod)

    def setDataPost(self, data):
        """"""
        self._core.setData(data)

    def setMeshPost(self, mesh):
        """"""
        self._core.setMesh(mesh, ignoreRegionManager=True)


class ERTModellingReference(ERTModellingBase):
    """Reference implementation for 2.5D Electrical Resistivity Tomography."""

    def __init__(self, **kwargs):
        super(ERTModelling, self).__init__()

        self.subPotentials = None
        self.lastResponse = None

        # only for mixed boundary hack since this need to know resistivies.
        self.resistivity = None

        # abscissa k and weight for 2.5 inverse cos-transform
        self.k = None
        self.w = None

    def response(self, model):
        """Solve forward task and return apparent resistivity for self.mesh."""
        # NOTE TODO can't be MT until mixed boundary condition depends on
        # self.resistivity
        pg.tic()
        if not self.data.allNonZero('k'):
            pg.error('Need valid geometric factors: "k".')
            pg.warn('Fallback "k" values to -sign("rhoa")')
            self.data.set('k', -pg.math.sign(self.data('rhoa')))

        mesh = self.mesh()

        nDof = mesh.nodeCount()
        elecs = self.data.sensorPositions()

        nEle = len(elecs)
        nData = self.data.size()

        self.resistivity = res = self.createMappedModel(model, -1.0)

        if self.verbose:
            print("Calculate response for model:", min(res), max(res))

        rMin = elecs[0].dist(elecs[1]) / 2.0
        rMax = elecs[0].dist(elecs[-1]) * 2.0

        k, w = self.getIntegrationWeights(rMin, rMax)

        self.k = k
        self.w = w

        # pg.show(mesh, res, label='res')
        # pg.wait()

        rhs = self.createRHS(mesh, elecs)

        # store all potential fields
        u = np.zeros((nEle, nDof))
        self.subPotentials = [pg.Matrix(nEle, nDof) for i in range(len(k))]

        for i, ki in enumerate(k):
            uE = pg.solve(mesh, a=1./res, b=-(ki * ki)/res, f=rhs,
                          bc={'Robin': ['*', self.mixedBC]},
                          userData={'sourcePos': elecs, 'k': ki},
                          verbose=False, stats=0, debug=False)
            self.subPotentials[i] = uE
            u += w[i] * uE

        # collect potential matrix,
        # i.e., potential for all electrodes and all injections
        pM = np.zeros((nEle, nEle))

        for i in range(nEle):
            pM[i] = pg.interpolate(mesh, u[i, :], destPos=elecs)

        # collect resistivity values for all 4 pole measurements
        r = np.zeros(nData)

        for i in range(nData):
            iA = int(self.data('a')[i])
            iB = int(self.data('b')[i])
            iM = int(self.data('m')[i])
            iN = int(self.data('n')[i])

            uAB = pM[iA] - pM[iB]
            r[i] = uAB[iM] - uAB[iN]

        self.lastResponse = r * self.data('k')

        if self.verbose:
            print("Resp min/max: {0} {1} {2}s".format(min(self.lastResponse),
                                                      max(self.lastResponse),
                                                      pg.dur()))

        return self.lastResponse

    def createJacobian(self, model):
        """TODO WRITEME."""
        if self.subPotentials is None:
            self.response(model)

        J = self.jacobian()
        J.resize(self.data.size(), self.regionManager().parameterCount())

        cells = self.mesh().findCellByMarker(0, -1)
        Si = pg.matrix.ElementMatrix()
        St = pg.matrix.ElementMatrix()

        u = self.subPotentials

        pg.tic()
        if self.verbose:
            print("Calculate sensitivity matrix for model: ",
                  min(model), max(model))

        Jt = pg.Matrix(self.data.size(),
                       self.regionManager().parameterCount())

        for kIdx, w in enumerate(self.w):
            k = self.k[kIdx]
            w = self.w[kIdx]

            Jt *= 0.
            A = pg.matrix.ElementMatrixMap()

            for i, c in enumerate(cells):
                modelIdx = c.marker()

                # 2.5D
                Si.u2(c)
                Si *= k * k
                Si += St.ux2uy2uz2(c)

                # 3D
                # Si.ux2uy2uz2(c); w = w* 2

                A.add(modelIdx, Si)

            for dataIdx in range(self.data.size()):

                a = int(self.data('a')[dataIdx])
                b = int(self.data('b')[dataIdx])
                m = int(self.data('m')[dataIdx])
                n = int(self.data('n')[dataIdx])
                Jt[dataIdx] = A.mult(u[kIdx][a] - u[kIdx][b],
                                     u[kIdx][m] - u[kIdx][n])

            J += w * Jt

        m2 = model*model
        k = self.data('k')

        for i in range(J.rows()):
            J[i] /= (m2 / k[i])

        if self.verbose:
            sumsens = np.zeros(J.rows())
            for i in range(J.rows()):
                sumsens[i] = pg.sum(J[i])
            print("sens sum: median = ", pg.math.median(sumsens),
                  " min = ", pg.min(sumsens),
                  " max = ", pg.max(sumsens))

    def calcGeometricFactor(self, data):
        """Calculate geometry factors for a given dataset."""
        if pg.y(data.sensorPositions()) == pg.z(data.sensorPositions()):
            k = np.zeros(data.size())
            for i in range(data.size()):
                a = data.sensorPosition(data('a')[i])
                b = data.sensorPosition(data('b')[i])
                m = data.sensorPosition(data('m')[i])
                n = data.sensorPosition(data('n')[i])
                k[i] = 1./(2.*np.pi) * (1./a.dist(m) - 1./a.dist(n) -
                                        1./b.dist(m) + 1./b.dist(n))
            return k
        else:
            raise BaseException("Please use BERT for non-standard "
                                "data sets" + str(data))

    def uAnalytical(self, p, sourcePos, k):
        """
        Calculate analytical potential for homogeneous halfspace.

            For sigma = 1 [S m]
        """
        r1A = (p - sourcePos).abs()
        # Mirror on surface at depth=0
        r2A = (p - pg.RVector3(1.0, -1.0, 1.0) * sourcePos).abs()

        if r1A > 1e-12 and r2A > 1e-12:
            return (pg.math.besselK0(r1A * k) + pg.math.besselK0(r2A * k)) / \
                    (2.0 * np.pi)
        else:
            return 0.

    def getIntegrationWeights(self, rMin, rMax):
        """TODO WRITEME."""
        nGauLegendre = max(int((6.0 * np.log10(rMax / rMin))), 4)
        nGauLaguerre = 4

        k = pg.Vector()
        w = pg.Vector()

        k0 = 1.0 / (2.0 * rMin)
        pg.GaussLegendre(0.0, 1.0, nGauLegendre, k, w)
        kLeg = k0 * k * k
        wLeg = 2.0 * k0 * k * w / np.pi

        pg.GaussLaguerre(nGauLaguerre, k, w)
        kLag = k0 * (k + 1.0)
        wLag = k0 * np.exp(k) * w / np.pi

        return pg.cat(kLeg, kLag), pg.cat(wLeg, wLag)

    def mixedBC(self, boundary, userData):
        """Apply mixed boundary conditions."""
        if boundary.marker() != pg.core.MARKER_BOUND_MIXED:
            return 0

        sourcePos = pg.center(userData['sourcePos'])

        k = userData['k']
        r1 = boundary.center() - sourcePos

        # Mirror on surface at depth=0
        r2 = boundary.center() - pg.RVector3(1.0, -1.0, 1.0) * sourcePos
        r1A = r1.abs()
        r2A = r2.abs()

        rho = 1.

        if self.resistivity is not None:
            rho = self.resistivity[boundary.leftCell().id()]

        n = boundary.norm()

        if r1A > 1e-12 and r2A > 1e-12:
            # see mod-dc-2d example for robin like BC and the negative sign
            if (pg.math.besselK0(r1A * k) + pg.math.besselK0(r2A * k)) > 1e-12:

                return k / rho * (r1.dot(n) / r1A * pg.math.besselK1(r1A * k) +
                                  r2.dot(n) / r2A * pg.math.besselK1(r2A * k))\
                                / (pg.math.besselK0(r1A * k) +
                                   pg.math.besselK0(r2A * k))
            else:
                return 0.
        else:
            return 0.

    def pointSource(self, cell, f, userData):
        r"""
        Define function for the current source term.

        :math:`\delta(x-pos), \int f(x) \delta(x-pos)=f(pos)=N(pos)`
            Right hand side entries will be shape functions(pos)
        """
        i = userData['i']
        sourcePos = userData['sourcePos'][i]

        if cell.shape().isInside(sourcePos):
            f.setVal(cell.N(cell.shape().rst(sourcePos)), cell.ids())

    def createRHS(self, mesh, elecs):
        """Create right-hand-side vector."""
        rhs = np.zeros((len(elecs), mesh.nodeCount()))
        for i, e in enumerate(elecs):
            c = mesh.findCell(e)
            rhs[i][c.ids()] = c.N(c.shape().rst(e))
        return rhs


class ERTManager(MeshMethodManager):
    """ERT Manager.

    Method Manager for Electrical Resistivity Tomography (ERT)

    TODO:
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

        Perform the forward task for a given mesh, a resistivity distribution
        a measuring scheme and return data (apparent resistivity) or potentials.

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
        # TODO: Remove pybert dependencies
        # >>> import pybert as pb
        # >>> import pygimli as pg
        # >>> import pygimli.meshtools as mt
        # >>> world = mt.createWorld(start=[-50, 0], end=[50, -50],
        # ...                        layers=[-1, -5], worldMarker=True)
        # >>> scheme = pb.createData(
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
        # >>> ert = pb.ERTManager()
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

        #segfaults with self.fop (test & fix)
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
                        ipError = np.ones(ret.size()) * kwargs.pop('phiErr') / 1000
                    else:
                        ipError = abs(ret["phia"]) * noiseLevel

                    if verbose:
                        print("Data IP abs error estimate (min:max) ",
                              min(ipError), ":", max(ipError))

                phia += np.randn(ret.size(), seed=seed) * ipError
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

        return rae

    def estimateError(self, data, absoluteError=0.001, relativeError=0.03,
                      absoluteUError=None, absoluteCurrent=0.1):
        """Estimate error composed of an absolute and a relative part.

        This is a static method and will not alter any member of the Manager

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
        if relativeError >= 0.5:
            print("relativeError set to a value > 0.5 .. assuming this "
                  "is a percentage Error level dividing them by 100")
            relativeError /= 100.0

        if absoluteUError is None:
            if not data.allNonZero('rhoa'):
                pg.critical("We need apparent resistivity values "
                            "(rhoa) in the data to estimate a "
                            "data error.")
            error = relativeError + absoluteError / data('rhoa')
        else:
            u = None
            i = absoluteCurrent
            if data.haveData("i"):
                i = data('i')

            if data.haveData("u"):
                u = data('u')
            else:
                if data.haveData("r"):
                    u = data('r') * i
                elif data.haveData("rhoa"):
                    if data.haveData("k"):
                        u = data('rhoa') / data('k') * i
                    else:
                        pg.critical("We need (rhoa) and (k) in the"
                                    "data to estimate data error.")

                else:
                    pg.critical("We need apparent resistivity values "
                                "(rhoa) or impedances (r) "
                                "in the data to estimate data error.")

            error = pg.abs(absoluteUError / u) + relativeError

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
            fig, ax = plt.subplots()
            fig.set_size_inches(size)

            self.showResult(ax=ax, coverage=self.coverage(), **kwargs)
            fig.savefig(path + '/resistivity.pdf')

            return path, fig, ax
        return path


def createERTData(elecs, schemeName='none', **kwargs):
    """Create data scheme for compatibility (advanced version in BERT).

    Parameters
    ----------
    sounding : bool [False]
        Create a 1D VES Schlumberger configuration.
        elecs need to be an array with elecs[0] = mn/2 and elecs[1:] = ab/2.
    """
    if kwargs.pop('sounding', False):
        data = pg.DataContainerERT()
        data.setSensors(pg.cat(-elecs[::-1], elecs))

        nElecs = len(elecs)
        for i in range(nElecs-1):
            data.createFourPointData(i, i, 2*nElecs-i-1, nElecs-1, nElecs)

        return data

    if schemeName != "dd":
        import pybert as pb  # that's bad!!! TODO: remove pybert deps
        return pb.createData(elecs, schemeName, **kwargs)

    isClosed = kwargs.pop('closed', False)

    data = pg.DataContainerERT()
    data.setSensors(elecs)

    nElecs = len(elecs)
    a = []
    b = []
    m = []
    n = []
    eb = 0
    for i in range(nElecs):
        for j in range(eb + 2, nElecs):
            ea = i
            eb = ea + 1
            em = j
            en = em + 1

            if isClosed:
                en = en % nElecs

            if en < nElecs and en != ea:
                a.append(ea)
                b.append(eb)
                m.append(em)
                n.append(en)

    data.resize(len(a))
    data.add('a', a)
    data.add('b', b)
    data.add('m', m)
    data.add('n', n)
    data.set('valid', np.ones(len(a)))

    return data


if __name__ == "__main__":
    pass
