#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""2.5D non-optimized totalfield forward operator for ERT.

Please use the BERT package for more advanced forward operator
https://gitlab.com/resistivity-net/bert
"""

import numpy as np

import pygimli as pg
from pygimli.frameworks import MeshModelling
from pygimli.manager import MeshMethodManager


class ERTModelling(MeshModelling):
    """Reference implementation for 2.5D Electrical Resistivity Tomography."""

    def __init__(self, **kwargs):
        """"Constructor, optional with data container and mesh."""
        super().__init__()

        self.subPotentials = None
        self.lastResponse = None

        # only for mixed boundary hack since this need to know resistivies.
        self.resistivity = None

        # abscissa k and weight for 2.5 inverse cos-transform
        self.k = None
        self.w = None

    def createStartModel(self, rhoa):
        sm = pg.RVector(self.regionManager().parameterCount(), pg.median(rhoa))
        return sm

    def response(self, model):
        """Solve forward task.

        Create apparent resistivity values for a given resistivity distribution
        for self.mesh.
        """
        ### NOTE TODO can't be MT until mixed boundary condition depends on
        ### self.resistivity
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

        rhs = self.createRHS(mesh, elecs)

        # store all potential fields
        u = np.zeros((nEle, nDof))
        self.subPotentials = [pg.RMatrix(nEle, nDof) for i in range(len(k))]
        for i, ki in enumerate(k):
            ws = dict()
            # pg.p(ki, min(res), max(res))
            uE = pg.solve(mesh, a=1./res, b=(ki * ki)/res, f=rhs,
                          bc={'Robin': ['*', self.mixedBC]},
                          userData={'sourcePos': elecs, 'k': ki},
                          verbose=False, stats=0, debug=False)
            # pg.p(min(uE.flat), max(uE.flat))
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
            print("Resp: ", min(self.lastResponse), max(self.lastResponse))

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

    def drawData(self, ax, data, err=None, label=None, vals=None, **kwargs):
        """Draw data in given axe."""

        if hasattr(data, '__iter__'):
            vals = data
            data = self.data
        elif data is None:
            data = self.data

        if vals is None:
            vals = data('rhoa')

        #ISSUE
        # why is plotERT data not used instead?
        # CR: plotERT is BERT function and should not used in lib pygimli directly

        mid, sep = midconfERT(data)
        dx = np.median(np.diff(np.unique(mid))) * 2
        ax, _, ymap = pg.mplviewer.patchValMap(
            vals, mid, sep, dx=dx, ax=ax, logScale=True,
            label=r'Apparent resistivity in $\Omega$m')

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


def simulate(mesh, res, scheme, verbose=False, **kwargs):
    """Convenience function of you like static functions.
    
    Parameters
    ----------
    
    """
    sr = kwargs.pop('sr', True)
    ert = ERTManager(verbose=verbose)
    ert.setSingularityRemoval(sr)

    if isinstance(mesh, str):
        mesh = pg.load(mesh)

    if isinstance(scheme, str):
        scheme = pb.load(scheme)

    return ert.simulate(mesh, res, scheme, verbose, **kwargs)


class ERTManager(MeshMethodManager):
    """ERT Manager.
    """
    def __init__(self, **kwargs):
        """Constructor.
        
        Parameters
        ----------
        * useBert: bool [False]
            Use Bert forward operator instead of the reference implementation.
        * sr: bool [True]
            Use singularity removal technique to solve the forward problem
            Currently only possible by the BERT solver.
        """
        kwargs['useBert'] = kwargs.pop('useBert', False)
        kwargs['sr'] = kwargs.pop('sr', True)

        super(ERTManager, self).__init__(**kwargs)
        self._dataToken = 'rhoa'
        self.dataTrans = pg.RTransLog()
        self.inv.dataTrans = self.dataTrans

    def setSingularityRemoval(self, sr=True):
        """Turn on singularity removal."""
        self.reinitForwardOperator(sr=True)

    def createForwardOperator(self, **kwargs):
        """Create and choose forward operator. """
        useBert = kwargs.pop('useBert', False)

        if useBert:
            sr = kwargs.pop('sr', True)
            if sr:
                fop = pb.DCSRMultiElectrodeModelling(verbose=verbose)
            else:
                fop = pb.DCMultiElectrodeModelling(verbose=verbose)
        else:
            fop = ERTModelling(**kwargs)

        return fop

    def simulate(self, mesh, model, scheme, **kwargs):
        """Forward calculation for given mesh, data and resistivity."""
        pg.critical('implementme')

        # use !!self.fop!!
        fop = self.createForwardOperator()

        fop.setDataBasis(dataContainer=scheme)
        fop.setMesh(mesh, ignoreRegionManager=True)

        if not scheme.allNonZero('k'):
            scheme.set('k', pg.RVector(scheme.size(), -1))

        res = model
        rhoa = None
        isArrayData = None

        if hasattr(res[0], '__iter__'):
            isArrayData = True
            rhoa = np.zeros((len(res), scheme.size()))
            for i, r in enumerate(res):
                rhoa[i] = fop.response(r)
        else:
            rhoa = fop.response(res)

        noiseLevel = kwargs.pop('noiseLevel', 0.0)
        if noiseLevel > 0:
            err = kwargs.pop('noiseLevel', 0.03) + kwargs.pop('noiseAbs',
                                                              1e-4) / rhoa
            scheme.set('err', err)
            rhoa *= 1. + pg.randn(scheme.size()) * err

            if not isArrayData:
                scheme.set('rhoa', rhoa)

        if kwargs.pop('returnArray', False):
            return rhoa
        return scheme

    def invert(self, data=None, err=None, **kwargs):
        """Invert measured data.
        """
        #ensure data and error sizes here
        dataVals = None
        if isinstance(data, pg.DataContainer):
            self.fop.setDataSpace(dataContainer=data)
            dataVals = self.dataValues(data)
            errVals = self.errorValues(data, relative=True)
        else:
            dataVal = data
            errVals = err

        return super(ERTManager, self).invert(dataVals=dataVals,
                                              errVals=errVals,
                                              **kwargs)


def createERTData(elecs, schemeName='none', **kwargs):
    """ Simple data creator for compatibility (advanced version in BERT)."""
    if schemeName is not "dd":
        import pybert as pb  # that's bad!!! TODO: remove pybert deps
        return pb.createData(elecs, schemeName, **kwargs)

    if isinstance(elecs, pg.RVector):
        sPos = []
        for e in elecs:
            sPos.append(pg.RVector3(e, 0., 0.))
        elecs = sPos

    isClosed = kwargs.pop('closed', False)

    data = pg.DataContainerERT()
    data.registerSensorIndex('a')
    data.registerSensorIndex('b')
    data.registerSensorIndex('m')
    data.registerSensorIndex('n')
    data.setSensorPositions(elecs)
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























from pygimli.manager import MeshMethodManager0


class ERTModelling0(pg.ModellingBase):
    """Minimal Forward Operator for 2.5D Electrical resistivity Tomography."""

    def __init__(self, mesh=None, data=None, verbose=False):
        """Constructor, optionally with data container and mesh."""
        super(ERTModelling, self).__init__(verbose=verbose)

        self.setVerbose(verbose=verbose)
        self.electrodes = None
        self.subPotentials = None
        self.lastResponse = None
        self.resistivity = None
        self.k = None
        self.w = None
        self.setMesh(mesh)
        self.setData(data)

    def setData(self, data):
        """Set a DataContainer."""
        # pg.ModellingBase.setData(data)
        if data is not None:
            self.electrodes = data.sensorPositions()

        self.data = data

    def setMesh(self, mesh, ignoreRegionManager=True):
        """Set Mesh."""
        if mesh is not None:  # ignore default different from ModBase (False)
            pg.ModellingBase.setMesh(self, mesh, ignoreRegionManager)
            # Landscape complains that ModellingBaseMT does not have setMesh

    def calcGeometricFactors(self, data):
        """Calculate geometry factors for a given dataset."""
        if pg.y(data) == pg.z(data):
            k = np.zeros(data.size())

            for i in range(data.size()):
                a = data.sensorPosition(data('a')[i])
                b = data.sensorPosition(data('b')[i])
                m = data.sensorPosition(data('m')[i])
                n = data.sensorPosition(data('n')[i])
                k[i] = 2.*np.pi * 1./(1./a.dist(m) - 1./a.dist(n) - 
                                      1./b.dist(m) + 1./b.dist(n))
            return k
        else:
            raise BaseException("Please use BERT for non-standard "
                                "data sets" + str(data))

    def uAnalytical(self, p, sourcePos, k):
        """Calculates the analytical solution for the 2.5D geoelectrical problem.
    
        Solves the 2.5D geoelectrical problem for one wave number k.
        It calculates the normalized (for injection current 1 A and sigma=1 S/m) 
        potential at position p for a current injection at position sourcePos.
        Injection at the subsurface is recognized via mirror sources along the
        surface at depth=0.
        
        Parameters
        ----------
        p : pg.Pos
            Position for the sought potential
        sourcePos : pg.Pos
            Current injection position.
        k : float
            Wave number 

        Returns
        -------
        u : float
            Solution u(p)
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
        """Retrieve Gauss-Legende/Laguerre integration weights."""
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
        r"""Define function for the current source term.

        :math:`\delta(x-pos), \int f(x) \delta(x-pos)=f(pos)=N(pos)`
            Right-hand-side entries will be shape functions(pos)
        """
        i = userData['i']
        sourcePos = userData['sourcePos'][i]

        if cell.shape().isInside(sourcePos):
            f.setVal(cell.N(cell.shape().rst(sourcePos)), cell.ids())

    def createRHS(self, mesh, elecs):
        """Create right-hand side."""
        rhs = np.zeros((len(elecs), mesh.nodeCount()))
        for i, e in enumerate(elecs):
            c = mesh.findCell(e)
            rhs[i][c.ids()] = c.N(c.shape().rst(e))
        return rhs

    def response(self, model):
        """Solve forward task.

        Create apparent resistivity values for a given resistivity distribution
        for self.mesh.
        """
        mesh = self.mesh()

        nDof = mesh.nodeCount()
        nEle = len(self.electrodes)
        nData = self.data.size()
 
        self.resistivity = res = self.createMappedModel(model, -1.0)

        if self.verbose():
            print("Calculate response for model:", min(res), max(res))

        rMin = self.electrodes[0].dist(self.electrodes[1]) / 2.0
        rMax = self.electrodes[0].dist(self.electrodes[-1]) * 2.0

        k, w = self.getIntegrationWeights(rMin, rMax)

        self.k = k
        self.w = w

        rhs = self.createRHS(mesh, self.electrodes)

        # store all potential fields
        u = np.zeros((nEle, nDof))
        self.subPotentials = [pg.RMatrix(nEle, nDof) for i in range(len(k))]

        for i, ki in enumerate(k):
            ws = {'u': self.subPotentials[i]}
            uE = pg.solve(mesh, a=1./res, b=-(ki * ki)/res, f=rhs,
                          bc={'Robin': self.mixedBC},
                          userData={'sourcePos': self.electrodes, 'k': ki},
                          verbose=self.verbose(), stats=0, debug=False,
                          ws=ws
                          )
            pg.show(mesh, np.log10(abs(uE[0])) )
            u += w[i] * uE

        # collect potential matrix,
        # i.e., potential for all electrodes and all injections
        pM = np.zeros((nEle, nEle))

        for i in range(nEle):
            pM[i] = pg.interpolate(mesh, u[i, :], destPos=self.electrodes)
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

        if self.verbose():
            print("Resp: ", min(self.lastResponse), max(self.lastResponse))

        return self.lastResponse

    def createJacobian(self, model):
        """Create Jacobian matrix."""
        if self.subPotentials is None:
            self.response(model)

        J = self.jacobian()
        J.resize(self.data.size(), self.regionManager().parameterCount())

        cells = self.mesh().findCellByMarker(0, -1)
        Si = pg.ElementMatrix()
        St = pg.ElementMatrix()

        u = self.subPotentials

        pg.tic()
        if self.verbose():
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

        if self.verbose():
            sumsens = np.zeros(J.rows())
            for i in range(J.rows()):
                sumsens[i] = pg.sum(J[i])
            print("sens sum: median = ", pg.median(sumsens),
                  " min = ", pg.min(sumsens),
                  " max = ", pg.max(sumsens))




class ERTManager0(MeshMethodManager0):
    """Minimalistic ERT Manager to keep compatibility.

    More advanced version comes with BERT.
    """
    def __init__(self, **kwargs):
        """Constructor."""
        super(ERTManager0, self).__init__(**kwargs)
        self.fop = ERTManager0.createFOP()
        self.setDataToken('rhoa')

    def showData(self, data=None, vals=None, ax=None):
        """Show mesh in given axes or in a new figure."""
        if data is None:
            data = self.data

        if vals is None:
            vals = data('rhoa')

        # why is plotERT data not used instead?
        # #c42: because pybert is not yet part of pygimli
        mid, sep = midconfERT(data)
        dx = np.median(np.diff(np.unique(mid)))*2
        ax, _, _ = pg.mplviewer.patchValMap(
            vals, mid, sep, dx=dx, ax=ax, logScale=True,
            label=r'Apparent resistivity ($\Omega$m)')

        return ax

    @staticmethod
    def createFOP(verbose=False):
        """Create forward operator working on refined mesh."""
        fop = ERTModelling0(verbose=verbose)
        return fop

    def createInv(self, fop, verbose=True, dosave=False):
        """Create inversion instance."""
        self.tD = pg.RTransLog()
        self.tM = pg.RTransLogLU()

        inv = pg.RInversion(verbose, dosave)
        inv.setTransData(self.tD)
        inv.setTransModel(self.tM)
        inv.setForwardOperator(fop)
        return inv

    @staticmethod
    def simulate(mesh, res, scheme, verbose=False, **kwargs):
        """Forward calculation vor given mesh, data and resistivity."""
        fop = ERTModelling0(verbose=verbose)
        # fop = ERTManager.createFOP(verbose=verbose)

        fop.setData(scheme)
        fop.setMesh(mesh, ignoreRegionManager=True)

        if not scheme.allNonZero('k'):

            if min(pg.y(scheme)) != max(pg.y(scheme)) or min(pg.z(scheme)) != max(pg.z(scheme)):
                pg.info("Non flat earth topography found. "
                    "We will set geometric factors to -1 to emulate "
                    "electrical impedance tomography (EIT). If you want to "
                    "use ERT will full topography support. "
                    "Please consider the use of pyBERT.")

                scheme.set('k', pg.RVector(scheme.size(), -1))
            else:
                scheme.set('k', fop.calcGeometricFactors(scheme))

        rhoa = None
        isArrayData = None

        if hasattr(res[0], '__iter__'):
            isArrayData = True
            rhoa = np.zeros((len(res), scheme.size()))
            for i, r in enumerate(res):
                rhoa[i] = fop.response(r)
        else:
            rhoa = fop.response(res)

        pg.renameKwarg('noisify', 'noiseLevel', kwargs)

        noiseLevel = kwargs.pop('noiseLevel', 0.0)

        if noiseLevel > 0:
            noiseAbs = kwargs.pop('noiseAbs', 1e-4)
            err = noiseLevel + noiseAbs / rhoa
            scheme.set('err', err)
            if verbose:
                pg.info("Set noise (" + str(noiseLevel*100) + "% + " + str(noiseAbs) + " V) min:",
                      min(err), "max:", max(err))
            rhoa *= 1. + pg.randn(scheme.size()) * err

        if isArrayData is None:
            scheme.set('rhoa', rhoa)

        if kwargs.pop('returnArray', False):
            return rhoa
        return scheme

    def createApparentData(self, data):
        """Create apparent data (what the hack is this?)"""
        return data('rhoa')

    def dataVals(self, data):
        """Return pure data values from a given DataContainer."""
        return data('rhoa')

    def relErrorVals(self, data):
        """Return pure data values from a given DataContainer."""
        return data('err')


def midconfERT(data):
    """Return the midpoint and configuration key for ERT data.

    Return the midpoint and configuration key for ERT data.

    Parameters
    ----------
    data : pybert.DataContainerERT
        data container with sensorPositions and a/b/m/n fields

    Returns
    -------
    mid : np.array of float
        representative midpoint (middle of MN, AM depending on array)
    conf : np.array of float
        configuration/array key consisting of
        1) array type (Wenner-alpha/beta, Schlumberger, PP, PD, DD, MG)
        2) potential dipole length
        3) separation factor
    """
#    xe = np.hstack((pg.x(data.sensorPositions()), np.nan))  # not used anymore
    x0 = data.sensorPosition(0).x()
    dx = [data.sensorPosition(i).distance(data.sensorPosition(i+1)) for i in
          range(data.sensorCount()-1)]
    xe = np.hstack((0., np.cumsum(np.round(dx, 1)), np.nan))
    de = np.median(np.diff(xe[:-1])).round(1)
    ne = np.round(xe/de)
    a, b, m, n = data('a'), data('b'), data('m'), data('n')
    # check if xe[a]/a is better suited (has similar size)
    a = np.array([ne[int(i)] for i in data('a')])
    b = np.array([ne[int(i)] for i in data('b')])
    m = np.array([ne[int(i)] for i in data('m')])
    n = np.array([ne[int(i)] for i in data('n')])
    ab, am, an = np.abs(a-b), np.abs(a-m), np.abs(a-n)
    bm, bn, mn = np.abs(b-m), np.abs(b-n), np.abs(m-n)
    # 2-point (default) 00000
    sep = np.abs(a-m)
    mid = (a+m) / 2 * de + x0
    # 3-point (PD, DP) (now only b==-1 or n==-<1, check also for a and m)
    imn = np.isfinite(n)*np.isnan(b)
    mid[imn] = (m[imn]+n[imn]) / 2 * de + x0
    sep[imn] = np.minimum(am[imn], an[imn]) + 10000 + 100 * (mn[imn]-1) + \
        (np.sign(a[imn]-m[imn])/2+0.5) * 10000
    iab = np.isfinite(b)*np.isnan(n)
    mid[iab] = (a[iab]+b[iab]) / 2 * de + x0  # better 20000 or -10000?
    sep[iab] = np.minimum(am[iab], bm[iab]) + 10000 + 100 * (ab[iab]-1) + \
        (np.sign(a[iab]-n[iab])/2+0.5) * 10000
    #  + 10000*(a-m)
    # 4-point alpha: 30000 (WE) or 4000 (SL)
    iabmn = np.isfinite(a) & np.isfinite(b) & np.isfinite(m) & np.isfinite(n)
    ialfa = np.copy(iabmn)
    ialfa[iabmn] = (ab[iabmn] > mn[iabmn])
    mid[ialfa] = (m[ialfa] + n[ialfa]) / 2 * de + x0
    spac = np.minimum(bn[ialfa], bm[ialfa])
    abmn3 = np.round((3*mn[ialfa]-ab[ialfa])*10000)/10000
    sep[ialfa] = spac + (mn[ialfa]-1)*100*(abmn3 != 0) + \
        30000 + (abmn3 < 0)*10000
    # gradient

    # %% 4-point beta
    ibeta = np.copy(iabmn)
    ibeta[iabmn] = (bm[iabmn] >= mn[iabmn]) & (~ialfa[iabmn])
    mid[ibeta] = (a[ibeta] + b[ibeta] + m[ibeta] + n[ibeta]) / 4 * de + x0
    sep[ibeta] = 50000 + (ab[ibeta]-1) * 100 + np.minimum(
        np.minimum(am[ibeta], an[ibeta]), np.minimum(bm[ibeta], bn[ibeta]))

    # %% 4-point gamma
    return mid, sep


if __name__ == "__main__":
    pass
