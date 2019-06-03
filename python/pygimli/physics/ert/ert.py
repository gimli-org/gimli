#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""2.5D non-optimized totalfield forward operator for ERT.

Please use the BERT package for more advanced forward operator
https://gitlab.com/resistivity-net/bert
"""

import numpy as np

import pygimli as pg
from pygimli.frameworks import Modelling, MeshModelling
from pygimli.frameworks import MeshMethodManager
from .visualization import showERTData

def simulate(mesh, res, scheme, sr=True, useBert=True,
             verbose=False, **kwargs):
    """Convenience function to use the ERT modelling operator.

    Convenience function to use the ERT modelling operator
    if you like static functions.

    See :py:mod:`pygimli.ert.ERTManager.simulate` for description
    of the arguments.

    Parameters
    ----------
    res : see :py:mod:`pygimli.ert.ERTManager.simulate`
        Resistivity distribution.
    mesh : :gimliapi:`GIMLI::Mesh` | str
        Modelling domain. Mesh can be a file name here.
    scheme : :gimliapi:`GIMLI::DataContainerERT` | str
        Data configuration. Scheme can be a file name here.
    sr : bool [True]
        Use singularity removal technique.
    useBert : bool [True]
        Use Bert forward operator instead of the reference implementation.
    **kwargs:
        Forwarded to :py:mod:`pygimli.ert.ERTManager.simulate`
    """
    ert = ERTManager(useBert=useBert, sr=sr, verbose=verbose)

    if isinstance(mesh, str):
        mesh = pg.load(mesh)

    if isinstance(scheme, str):
        scheme = pb.load(scheme)

    return ert.simulate(mesh, res, scheme, verbose=verbose, **kwargs)


class ERTModellingBase(MeshModelling):
    def __init__(self, **kwargs):
        super(ERTModellingBase, self).__init__(**kwargs)

    def drawData(self, ax, data=None, **kwargs):
        """Draw data in given axe."""

        if hasattr(data, '__iter__'):
            vals = data
            data = self.data
        elif data is None:
            data = self.data

        if vals is None:
            vals = data['rhoa']

        return showERTData(data, vals=vals, ax=ax, **kwargs)

    def drawModel(self, ax, model, **kwargs):
        """Draw the para domain with option model values"""
        ax, cBar = pg.show(mesh=self.paraDomain,
                           data=model,
                           label=kwargs.pop('label', pg.utils.unit('res')),
                           ax=ax,
                           cMap=pg.utils.cMap('res'),
                           **kwargs)
        return ax, cBar


class BertModelling(ERTModellingBase):
    def __init__(self, sr=True, verbose=False):
        """Constructor, optional with data container and mesh."""
        super(BertModelling, self).__init__()

        # don't use DC*fop or its regionmanager directly
        #
        self.bertFop = None
        if sr:
            self.bertFop = pg.DCSRMultiElectrodeModelling(verbose=verbose)
        else:
            self.bertFop = pg.DCMultiElectrodeModelling(verbose=verbose)

        self.bertFop.initJacobian()
        self.setJacobian(self.bertFop.jacobian())

        self.response = self.bertFop.response
        self.jacobian = self.bertFop.jacobian
        self.createJacobian = self.bertFop.createJacobian


        ## called from the ERTManager .. needed?
        self.solution = self.bertFop.solution
        self.setComplex = self.bertFop.setComplex
        self.complex = self.bertFop.complex
        self.calculate = self.bertFop.calculate
        self.calcGeometricFactor = self.bertFop.calcGeometricFactor
        self.mapERTModel = self.bertFop.mapERTModel

    def setDataPost(self, data):
        """"""
        self.bertFop.setData(data)

    def setMeshPost(self, mesh):
        """"""
        self.bertFop.setMesh(mesh, ignoreRegionManager=True)


class ERTModelling(ERTModellingBase):
    """Reference implementation for 2.5D Electrical Resistivity Tomography."""

    def __init__(self, **kwargs):
        """"Constructor, optional with data container and mesh."""
        super(ERTModelling, self).__init__()

        self.subPotentials = None
        self.lastResponse = None

        # only for mixed boundary hack since this need to know resistivies.
        self.resistivity = None

        # abscissa k and weight for 2.5 inverse cos-transform
        self.k = None
        self.w = None

    def response(self, model):
        """Solve forward task.

        Create apparent resistivity values for a given resistivity distribution
        for self.mesh.
        """
        ### NOTE TODO can't be MT until mixed boundary condition depends on
        ### self.resistivity
        pg.tic()
        if not self.data.allNonZero('k'):
            pg.error('Need valid geometric factors: "k".')
            pg.warn('Fallback "k" values to -sign("rhoa")')
            self.data.set('k', -pg.sign(self.data('rhoa')))

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
        self.subPotentials = [pg.RMatrix(nEle, nDof) for i in range(len(k))]

        for i, ki in enumerate(k):
            ws = dict()
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
        Si = pg.ElementMatrix()
        St = pg.ElementMatrix()

        u = self.subPotentials

        pg.tic()
        if self.verbose:
            print("Calculate sensitivity matrix for model: ",
                  min(model), max(model))

        Jt = pg.RMatrix(self.data.size(),
                        self.regionManager().parameterCount())

        for kIdx, w in enumerate(self.w):
            k = self.k[kIdx]
            w = self.w[kIdx]

            Jt *= 0.
            A = pg.ElementMatrixMap()

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
            print("sens sum: median = ", pg.median(sumsens),
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
            return (pg.besselK0(r1A * k) + pg.besselK0(r2A * k)) / \
                    (2.0 * np.pi)
        else:
            return 0.

    def getIntegrationWeights(self, rMin, rMax):
        """TODO WRITEME."""
        nGauLegendre = max(int((6.0 * np.log10(rMax / rMin))), 4)
        nGauLaguerre = 4

        k = pg.RVector()
        w = pg.RVector()

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
        if boundary.marker() != pg.MARKER_BOUND_MIXED:
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
            if (pg.besselK0(r1A * k) + pg.besselK0(r2A * k)) > 1e-12:

                return 1./rho * k * (r1.dot(n) / r1A * pg.besselK1(r1A * k) +
                                     r2.dot(n) / r2A * pg.besselK1(r2A * k)) /\
                                (pg.besselK0(r1A * k) + pg.besselK0(r2A * k))
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
        """TODO WRITEME."""
        rhs = np.zeros((len(elecs), mesh.nodeCount()))
        for i, e in enumerate(elecs):
            c = mesh.findCell(e)
            rhs[i][c.ids()] = c.N(c.shape().rst(e))
        return rhs


class ERTManager(MeshMethodManager):
    """ERT Manager.

    Method Manager for Electrical Resistivity Tomography (ERT)

    TODO:
        * 2d
        * 2dtopo
        * 3d
        * 3dtopo
        * complex on/off
        * closed geometry
        * transdim
        * singularity removal
        * ERT specific inversion options:
            * ...
    """
    def __init__(self, **kwargs):
        """Constructor.

        Parameters
        ----------
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
        kwargs['useBert'] = kwargs.pop('useBert', True)
        kwargs['sr'] = kwargs.pop('sr', True)

        super(ERTManager, self).__init__(**kwargs)
        self.inv.dataTrans = pg.RTransLogLU()

    def setSingularityRemoval(self, sr=True):
        """Turn singularity removal on or off."""
        self.reinitForwardOperator(sr=True)

    def createForwardOperator(self, **kwargs):
        """Create and choose forward operator. """
        useBert = kwargs.pop('useBert', False)
        verbose = kwargs.pop('verbose', False)
        if useBert:
            fop = BertModelling(sr=kwargs.pop('sr', True), verbose=verbose)
        else:
            fop = ERTModelling(**kwargs)

        return fop

    def setPrimPot(self, pot):
        """
        """
        Implementme

    def simulate(self, mesh, res, scheme, **kwargs):
        """Simulate an ERT measurement.

        Perform the forward task for a given mesh, a resistivity distribution
        (per cell), a measurement
        scheme and will return data (apparent resistivity) or potential fields.

        This function can also operate on complex resistivity models, thereby
        computing complex apparent resistivities.

        The forward operator itself only calculate potential values
        for the given scheme file.
        To calculate apparent resistivities, geometric factors (k) are needed.
        If there are no values k in the DataContainerERT scheme, then we will
        try to calculate them, either analytic or by using a p2-refined
        version of the given mesh.

        TODO
        ----
        * 2D + Complex + SR

        Parameters
        ----------
        mesh : :gimliapi:`GIMLI::Mesh`
            2D or 3D Mesh to calculate for.

        res : float, array(mesh.cellCount()) | array(N, mesh.cellCount()) | list
            Resistivity distribution for the given mesh cells can be:
            . float for homogeneous resistivity
            . single array of length mesh.cellCount()
            . matrix of N resistivity distributions of length mesh.cellCount()
            . resistivity map as [[regionMarker0, res0],
                                  [regionMarker0, res1], ...]

        scheme : :bertapi:`Bert::DataContainerERT`
            Data measurement scheme.

        Other Parameters
        ----------------
        verbose: bool[False]
            Be verbose. Will override class settings.
        calcOnly: bool [False]
            Use fop.calculate instead of fop.response. Useful if you want
            to force the calculation of impedances for homogeneous models.
            No noise handling. Solution is put as token 'u' in the returned
            DataContainerERT.
        noiseLevel: float [0.0]
            add normally distributed noise based on
            scheme('err') or on noiseLevel if scheme did not contain 'err'
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
        rhoa : DataContainerERT | array(N, data.size()) | array(N, data.size()),
            array(N, data.size())
                Data container with resulting apparent resistivity data and
                errors (if noiseLevel or noiseAbs is set).
                Optional returns a Matrix of rhoa values
                (for returnArray==True forces noiseLevel=0).
                In case of a complex valued resistivity model, phase values will be
                returned in the DataContainerERT (see example below), or as an
                additional returned array.

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

        fop = self.fop
        fop.data = scheme
        fop.setMesh(mesh, ignoreRegionManager=True)
        fop.verbose = verbose

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

                    dMap = pg.DataMap()
                    fop.calculate(dMap)
                    if fop.complex():
                        pg.critical('Implement me')
                    else:
                        ret.set("u", dMap.data(scheme))
                        ret.set("i", np.ones(ret.size()))

                    if returnFields:
                        return pg.Matrix(fop.solution())
                    return ret
                else:
                    if fop.complex():
                        res = pg.physics.SIP.squeezeComplex(res)
                    
                    resp = fop.response(res)

                    if fop.complex():
                        rhoa, phia = pg.physics.SIP.toPolar(resp)
                        rhoa *= scheme['k']
                        phia *= 1.0 # we don't want to change the sign
                    else:
                        rhoa = resp
            else:
                print(mesh)
                print("res: ", res)
                raise BaseException("Simulate called with wrong resistivity array.")


        if not isArrayData:
            ret.set('rhoa', rhoa)

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

            rhoa *= 1. + pg.randn(ret.size()) * ret('err')
            ret.set('rhoa', rhoa)

            ipError = None
            if phia is not None:
                if scheme.allNonZero('iperr'):
                    ipError = scheme('iperr')
                else:
                    # np.abs(self.data("phia") +TOLERANCE) * 1e-4absoluteError
                    if noiseLevel > 0.5:
                        noiseLevel /= 100.

                    ipError = ret("phia") * noiseLevel

                    if verbose:
                        print("Data IP abs error estimate (min:max) ",
                               min(ipError), ":", max(ipError))

                phia *= (1. + pg.randn(ret.size()) * noiseLevel)
                ret.set('iperr', ipError)
                ret.set('phia', phia)

        # check what needs to be setup and returned

        if returnArray:
            if phia is not None:
                return rhoa, phia
            else:
                return rhoa

        return ret

    def dataCheck(self, data):
        """Return data from container"""
        if isinstance(data, pg.DataContainer):
            if not data.haveData('rhoa'):
                pg.critical('Datacontainer have no "rhoa" values.')
            return data['rhoa']
        return data

    def estimateError(self, data, absoluteError=0.001, relativeError=0.03,
                      absoluteUError=None, absoluteCurrent=0.1):
        """ Estimate error composed of an absolute and a relative part.
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
        """Return coverage vector considering the logarithmic transformation.
        """
        covTrans = pg.coverageDCtrans(self.fop.jacobian(),
                                      1.0 / self.inv.response,
                                      1.0 / self.inv.model)
        return np.log10(covTrans / self.inv.paraDomain.cellSizes())

    def standardizedCoverage(self, threshhold=0.01):
        """Return standardized coverage vector (0|1) using thresholding.
        """
        return 1.0*(abs(self.coverage()) > threshhold)


def createERTData(elecs, schemeName='none', **kwargs):
    """ Simple data creator for compatibility (advanced version in BERT).

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

    if schemeName is not "dd":
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
