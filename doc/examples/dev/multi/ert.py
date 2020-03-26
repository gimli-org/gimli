#!/usr/bin/env python

"""Test multi. Is this used anymore?"""

import time

import numpy as np

import pygimli as pg
from pygimli.solver import parseArgToArray
from pygimli.meshtools import createParaMesh2dGrid, appendTriangleBoundary

import pybert as pb


def createCacheName(base, mesh=None):
    nc = ''
    if mesh:
        nc = str(mesh.nodeCount())
    return 'cache-' + base + "-" + nc


class ERT():
    """
        ERT Manager.

        Should solve most of the common problems
    """

    def __init__(self, verbose=False):
        """
        """
        self.fop = self.createFOP(verbose)
        self.tD = None
        self.tM = None
        self.inv = self.createInv(verbose)
        self.schemeMg = pb.dataview.dataview.DataSchemeManager()
        self.paraMesh = None

    def createFOP(self, verbose):
        """ Create resistivity modelling forward operator. """
        return pb.DCSRMultiElectrodeModelling(verbose=verbose)

    def createInv(self, verbose):
        """ Create resistivity inversion instance. """
        self.tD = pg.trans.TransLog()
        self.tM = pg.trans.TransLog()
        inv = pg.Inversion(verbose=verbose, dosave=False)
        inv.setTransData(self.tD)
        inv.setTransModel(self.tM)
        return inv

    def show(self, data, values=None, axes=None,
             cMin=None, cMax=None, colorBar=1, **kwargs):
        """
        """
        gci = None
        if isinstance(data, pg.DataContainer):
            scheme = kwargs.pop('scheme', 'unknown')
            pseudoscheme = getattr(pb.dataview.Pseudotype, scheme)

            if values is None:
                values = data('rhoa')

            gci = pb.dataview.dataview.drawDataAsMatrix(
                    axes, data, values, pseudotype=pseudoscheme)
            gci.set_clim((cMin, cMax))
            if colorBar:
                pg.viewer.mpl.colorbar.createColorbar(
                        gci, nLevs=5, cMin=cMin,  cMax=cMax,
                        label='Apparent resistivity in $\Omega m$', **kwargs)
        elif isinstance(data, pg.Vector):
            # assuming this is a model from the last inversion run
            pg.show(self.fop.regionManager().paraDomain(), data,
                    axes=axes, **kwargs)

        else:
            print(type(data))
            raise
        return gci

    def invert(self, data, values=None, verbose=0, **kwargs):
        """
        Invert the given data.

        A parametric mesh for the inversion will be created if non is given
        before.

        Parameters
        ----------
        """
        self.fop.setVerbose(verbose)
        self.inv.setVerbose(verbose)
        self.inv.setMaxIter(kwargs.pop('maxiter', 10))
        self.inv.setLambda(kwargs.pop('lambd', 10))

        if self.paraMesh is None:
            self.paraMesh = createParaMesh2dGrid(data.sensorPositions(),
                                                 **kwargs)
            self.setParaMesh(self.paraMesh)
            if verbose:
                print(self.paraMesh)
#                pg.show(self.paraMesh)

        err = data('err')
        rhoa = data('rhoa')

        startModel = pg.Vector(self.fop.regionManager().parameterCount(),
                                pg.math.median(rhoa))

        self.fop.setData(data)
        self.inv.setForwardOperator(self.fop)

        # check err here
        self.inv.setData(rhoa)
        self.inv.setError(err)
        self.inv.setModel(startModel)

        model = self.inv.run()

        if values is not None:

            if isinstance(values, pg.Vector):
                values = [values]
            elif isinstance(values, np.ndarray):
                if values.ndim == 1:
                    values = [values]

            allModel = pg.Matrix(len(values), len(model))

            self.inv.setVerbose(False)
            for i in range(len(values)):
                print(i)
                tic = time.time()
                self.inv.setModel(model)
                self.inv.setReferenceModel(model)
                dData = pg.abs(values[i] / rhoa)

                relModel = self.inv.invSubStep(pg.log(dData))
                allModel[i] = model * pg.exp(relModel)
                print(i, "/", len(values), " : ", time.time()-tic,
                      "s min/max: ", min(allModel[i]), max(allModel[i]))

            return allModel
        return model

    def setParaMesh(self, mesh):
        """
        Set the parameter mesh for any inversion.

        Parameters
        ----------
        """
        self.fop.setMesh(mesh)
        self.fop.regionManager().region(1).setBackground(True)
        self.fop.createRefinedForwardMesh(refine=True, pRefine=False)

    def createData(self, sensors, scheme):
        """
        Create a empty data file. i.e. a measurement scheme.

        Parameters
        ----------
        """
        scheme = self.schemeMg.scheme(scheme)
        scheme.setInverse(False)
        scheme.addInverse(False)
        return scheme.create(sensorList=sensors)

    def simulate(self, mesh, resistivity, data):
        """
        Simulation a ERT measurement.

        Perform the forward task for a given mesh,
        a resistivity distribution and a measurement scheme.

        Parameters
        ----------
        """
        self.fop.setMesh(mesh)
        self.fop.setData(data)
        res = self.fop.response(resistivity)
        data = pb.DataContainerERT(self.fop.data())
        data.set('rhoa', res)
        return data


def resistivityArchie(rBrine, porosity, a=1.0, m=2.0, S=1.0, n=2.0,
                      mesh=None, meshI=None):
    """
    .. math::
        \rho = a\rho_{\text{Brine}}\phi^{-m}\S_w^{-n}

    * :math:`\rho` - the electrical conductivity of the fluid saturated rock
    * :math:`\rho_{\text{Brine}}` - electrical conductivity of the brine
    * :math:`\phi` - porosity 0.0 --1.0
    * :math:`a` - tortuosity factor. (common 1)
    * :math:`m` - cementation exponent of the rock
            (usually in the range 1.3 -- 2.5 for sandstones)
    * :math:`n` - is the saturation exponent (usually close to 2)

    """
    rB = None

    if rBrine.ndim == 1:
        rB = pg.Matrix(1, len(rBrine))
        rB[0] = parseArgToArray(rBrine, mesh.cellCount(), mesh)
    elif rBrine.ndim == 2:
        rB = pg.Matrix(len(rBrine), len(rBrine[0]))
        for i in range(len(rBrine)):
            rB[i] = rBrine[i]

    porosity = parseArgToArray(porosity, mesh.cellCount(), mesh)
    a = parseArgToArray(a, mesh.cellCount(), mesh)
    m = parseArgToArray(m, mesh.cellCount(), mesh)
    S = parseArgToArray(S, mesh.cellCount(), mesh)
    n = parseArgToArray(n, mesh.cellCount(), mesh)

    r = pg.Matrix(len(rBrine), len(rBrine[0]))
    for i in range(len(r)):
        r[i] = rB[i] * a * porosity**(-m) * S**(-n)

    rI = pg.Matrix(len(r), meshI.cellCount())
    if meshI:
        pg.interpolate(mesh, r, meshI.cellCenters(), rI)

    for i in range(len(rI)):
        rI[i] = pg.solver.fillEmptyToCellArray(meshI, rI[i])

    return rI


def calcApparentResistivities(mesh, meshERT, poro, rhoBrine):
    ert = ERT(verbose=False)

    meshFOP = appendTriangleBoundary(meshERT,
                                     xbound=50, ybound=50, marker=1,
                                     quality=34.0, smooth=False,
                                     markerBoundary=1,
                                     isSubSurface=False, verbose=False)

    swatch = pg.core.Stopwatch(True)

    print("res 1:", swatch.duration(True))

    resis = resistivityArchie(rBrine=rhoBrine, porosity=poro, S=1.0,
                              mesh=mesh, meshI=meshFOP)

    print("res 2:", swatch.duration(True))

    ertPointsX = [pg.RVector3(x, 0) for x in np.arange(-19, 19.1, 1)]
    ertScheme = ert.createData(ertPointsX, scheme="Dipole Dipole (CC-PP)")

    solutionName = createCacheName('appRes', mesh) + "-" + \
        str(ertScheme.size()) + "-" + str(len(rhoBrine))

    try:
        rhoa = np.load(solutionName + '.bmat.npy')
        ertData = pb.DataContainerERT(solutionName + '.dat')
    except Exception as e:
        print(e)
        print("Building .... ")
        rhoa = np.zeros((len(resis), ertScheme.size()))
        ertScheme.set('k', pb.geometricFactor(ertScheme))
        ertData = ert.simulate(meshFOP, resis[0], ertScheme)

        errPerc = 1
        errVolt = 1e-5
        voltage = ertData('rhoa') / ertData('k')
        ertData.set('err', pg.abs(errVolt / voltage) + errPerc / 100.0)
        print('err min:', min(ertData('err'))*100, 'max:',
              max(ertData('err'))*100)
        ertData.save(solutionName + '.dat', 'a b m n rhoa err k')
        for i in range(0, len(resis)):
            tic = time.time()
            rhoa[i] = ert.fop.response(resis[i])

            rand = pg.Vector(len(rhoa[i]))
            pg.math.randn(rand)

            rhoa[i] *= (1.0 + rand * ertData('err'))

            print(i, "/", len(resis), " : ", time.time()-tic, "s",
                  "min:", min(resis[i]), "max:", max(resis[i]),
                  "min:", min(rhoa[i]), "max:", max(rhoa[i]))

        np.save(solutionName + '.bmat', rhoa)

    return meshFOP, resis, ertData, rhoa


def calcERT(ertScheme, rhoa):
    ert = ERT(verbose=False)

    solutionName = createCacheName('ERT') + "-" + str(ertScheme.size()) + \
        "-" + str(len(rhoa))
    try:
        ertModels = pg.load(solutionName + '.bmat')
        ertMesh = pg.load(solutionName + '.bms')
    except Exception as e:
        print(e)
        print("Building .... ")
        ertModels = ert.invert(ertScheme, values=rhoa, maxiter=10, lambd=50,
                               paraDX=0.5, paraDZ=0.5,
                               nLayers=20, paraDepth=15,
                               verbose=1)
        ertMesh = ert.fop.regionManager().paraDomain()
        ertModels.save(solutionName + '.bmat')
        ertMesh.save(solutionName)

    ertRatioModels = pg.Matrix(ertModels)
    for i in range(len(ertModels)):
        ertRatioModels[i] /= ertModels[0]

    return ertMesh, ertModels, ertRatioModels


if __name__ == "__main__":
    pass
