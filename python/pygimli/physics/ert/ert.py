#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
2.5D non optimized totalfield FOP for ERT.

Please use the BERT package for more advanced FOP
https://gitlab.com/resistivity-net/bert
"""

import numpy as np

import pygimli as pg

from pygimli.manager import MeshMethodManager


class ERTModelling(pg.ModellingBase):
    """Minimal Forward Operator for 2.5D Electrical resistivity Tomography."""

    def __init__(self, mesh=None, data=None, verbose=False):
        """"Constructor."""
        super().__init__()

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
        if mesh is not None:
            pg.ModellingBase.setMesh(self, mesh)

    def calcGeometricFactor(self, data):
        """Calculate geometry factors for a given dataset."""
        if pg.y(data.sensorPositions()) == pg.z(data.sensorPositions()):
            k = np.zeros(data.size())
            for i in range(data.size()):
                a = data.sensorPosition(data('a')[i])
                b = data.sensorPosition(data('b')[i])
                m = data.sensorPosition(data('m')[i])
                n = data.sensorPosition(data('n')[i])
                k[i] = 1./(2.*np.pi) * \
                       1./(a.dist(m) - a.dist(n) - b.dist(m) + b.dist(b))
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
        """TODO WRITEME."""
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
            uE = pg.solve(mesh, a=1./res, b=(ki * ki)/res, f=rhs,
                          duB=self.mixedBC,
                          userData={'sourcePos': self.electrodes, 'k': ki},
                          verbose=False, stat=0, debug=False,
                          ret=self.subPotentials[i])
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


class ERTManager(MeshMethodManager):
    """Minimalistic ERT Manager to keep compatibility. More advanced version
    comes with BERT.
    """
    def __init__(self, **kwargs):
        """Constructor."""
        super(MeshMethodManager, self).__init__(**kwargs)
        self.setDataToken('rhoa')

    def showData(self, ax=None):
        """Show mesh in given axes or in a new figure."""
        raise Exception('use the BERT Manager')

    @staticmethod
    def createFOP(verbose=False):
        """Create forward operator working on refined mesh."""
        fop = ERTModelling(verbose=verbose)
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
        """"""
        fop = ERTModelling(verbose=verbose)
        #fop = ERTManager.createFOP(verbose=verbose)

        fop.setData(scheme)
        fop.setMesh(mesh, ignoreRegionManager=True)

        if not scheme.allNonZero('k'):
            scheme.set('k', pg.RVector(scheme.size(), -1))

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

            err = kwargs.pop('noiseLevel', 0.03) + kwargs.pop('noiseAbs', 1e-4) / rhoa
            scheme.set('err', err)
            rhoa *= 1. + pg.randn(scheme.size()) * err

            if not isArrayData:
                scheme.set('rhoa', rhoa)

        if kwargs.pop('returnArray', False):
            return rhoa
        return scheme

    def createApparentData(self, data):
        """ what the hack is this?"""
        return data('rhoa')


def createERTData(elecs, schemeName='none', **kwargs):
    """ Simple data creator to keep compatibility. More advanced version
    comes with BERT.
    """
    if schemeName is not "dd":
        import pybert as pb
        return bp.createData(elecs, schemeName, **kwargs)

    if isinstance(elecs, pg.RVector):
        sPos = []
        for e in elecs:
            sPos.append(pg.RVector3(e, 0., 0.))
        elecs = sPos

    isClosed = kwargs.pop('closed', False)

    data = pg.DataContainer()
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
                en = en%nElecs

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

    data.save('tmp.dat')
    return data

if __name__ == "__main__":
    pass
