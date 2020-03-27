#!/usr/bin/env python

import time
import numpy as np
import matplotlib.pyplot as plt

import pygimli as pg
from pygimli.solver import parseArgToArray
from pygimli.physics.gravimetry import solveGravimetry


def density(poro, densMatrix=2510, densFluid=1000, satur=1,
            mesh=None):
    r"""
        densMatrix, densFluid in kg/m^3
    """
    poro = parseArgToArray(poro, mesh.cellCount(), mesh)
    densMatrix = parseArgToArray(densMatrix, mesh.cellCount(), mesh)
    densFluid = parseArgToArray(densFluid, mesh.cellCount(), mesh)
    satur = parseArgToArray(satur, mesh.cellCount(), mesh)

    dens = np.array(densMatrix * (1.-poro)) + densFluid * poro * satur
    return dens


class Gravimetry():
    """
        General Gravimetry Method Manager
    """

    def __init__(self, verbose=False):
        """Default constructor."""
        self.fop = self.createFOP(verbose)
        self.tD = None
        self.tM = None
        self.inv = self.createInv(verbose)

    def createFOP(self, verbose):
        return pg.physics.gravimetry.GravimetryModelling(verbose=verbose)

    def createInv(self, verbose):
        self.tD = pg.trans.Trans()
        self.tM = pg.trans.Trans()
        inv = pg.Inversion(verbose=verbose, dosave=False)
        inv.setTransData(self.tD)
        inv.setTransModel(self.tM)
        return inv

    def setParaMesh(self, mesh):
        """
        Set the parameter mesh for any inversion.

        Parameters
        ----------
        """
        self.fop.setMesh(mesh)
        self.fop.createRefinedForwardMesh(refine=False, pRefine=False)

    def invert(self, sensorPositions, gz, errAbs,
               verbose=0, **kwargs):
        """
        """
        self.fop.setVerbose(verbose)
        self.inv.setMaxIter(kwargs.pop('maxiter', 10))
        self.inv.setLambda(kwargs.pop('lambd', 10))

        self.fop.setSensorPositions(sensorPositions)
        mesh = kwargs.pop('mesh', None)
        if mesh is None:
            raise('implement me')

        self.setParaMesh(mesh)

        startModel = pg.Vector(self.fop.regionManager().parameterCount(), 0.0)

        self.inv.setForwardOperator(self.fop)

        self.fop.regionManager().setConstraintType(10)
        # check err here
        self.inv.setData(gz)
        self.inv.setAbsoluteError(errAbs)
        self.inv.setModel(startModel)

        model = self.inv.run()
        return model
        # tl can start here
        values = model
        if values is not None:

            if isinstance(values, pg.Vector):
                values = [values]
            elif isinstance(values, np.ndarray):
                if values.ndim == 1:
                    values = [values]

            allModel = pg.Matrix(len(values)+1, len(model))
            allModel[0] = model
            self.inv.setVerbose(False)
            for i in range(1, len(values)):
                tic = time.time()
                self.inv.setModel(model)
                self.inv.setReferenceModel(model)
                dData = pg.abs(values[i] - data)

                # relModel = self.inv.invSubStep(pg.log(dData))
                # allModel[i] = model * pg.exp(relModel)

                relModel = self.inv.invSubStep(dData)
                allModel[i] = model + relModel

                print(i, "/", len(values), " : ", time.time()-tic, "s")

            return allModel

        return model

    def simulate(self, mesh, dDensity):
        self.fop.setMesh(mesh)
        # TODO!


def calcInvBlock(mesh, dens, out='gravInv'):

    # extract block delta density
    densBlock = pg.Vector(dens)
    densMarker2 = dens[pg.find(mesh.cellMarker() == 2)[0]]
#    densBlock[(mesh.cellMarker() == 1)|(mesh.cellMarker() == 3)] = densMarker2
    densBlock[pg.find((mesh.cellMarker() == 1) | (mesh.cellMarker() == 3))] = \
        densMarker2
    densBlock -= densMarker2

    # define meausrement positions
    gravPointsX = np.linspace(-20, 20, 41)
    sensorPositions = np.vstack((gravPointsX, np.zeros(len(gravPointsX)))).T

    # solve analytical
    gz = solveGravimetry(mesh, densBlock, pnts=sensorPositions, complete=False)

    # noisyfy
    errAbs = 0.00001
    dzerr = np.random.randn(len(sensorPositions)) * errAbs
    gz = gz + dzerr

    # createParamesh
    paraMesh = pg.createGrid(x=np.linspace(-20, 20, 41),
                             y=np.linspace(-20, 0, 21))

    # init Gravimetry manager (should do meshing, simulation and noisying)
    Grav = Gravimetry(verbose=True)

    model = Grav.invert(sensorPositions, gz, errAbs, verbose=1, mesh=paraMesh)

    fig, ax = plt.subplots()
    ax.plot(pg.x(sensorPositions), gz, label='gz')
    ax.plot(pg.x(sensorPositions), Grav.inv.response(), label='response')
    ax.legend()
    ax.grid()
    ax.set_xlabel('$x$ [m]')
    ax.set_ylabel('$\partial u / \partial z$ [mGal]')
    plt.show(block=False)
    ax.figure.savefig(out, bbox_inches='tight')

    return Grav, densBlock

    #savefig(mesh, plc, densBlock, '$\Delta$ Density [kg$/$m$^3$]')



def simulateGravimetry(mesh, dDens):
    gravPointsX = np.arange(-19, 19.1, 1)
    gravPoints = np.vstack((gravPointsX, np.zeros(len(gravPointsX)))).T
    solutionName = createCacheName('grav', mesh, times)
    try:
        #vel = pg.load(solutionName + '.bmat')
        Gdg = np.load(solutionName + '.bmat.npy')

    except Exception as e:
        print(e)
        print("Building .... ")
        #Gdg, Gdgz = solveGravimetry(mesh, None, pnts=gravPoints, complete=True)
        Gdg = solveGravimetry(mesh, None, pnts=gravPoints)
        np.save(solutionName + '.bmat', Gdg)

    #dz = Gdg.dot(dDens.transpose([1,0])).T

    dz = np.zeros((len(dDens), len(gravPoints)))
    for i in range(len(dDens)):
        dzerr = np.random.randn(len(gravPoints)) * 0.01
        dz[i] = Gdg.dot(dDens[i]) + dzerr

    print(Gdg.shape, dDens.shape, dz.shape)
    return gravPoints, dz

def invertGravimetry(gravPoints, dz):

    dzerr = np.random.randn(len(gravPoints)) * 0.0001
    dz = dz + dzerr

    mesh = pg.createGrid(x=np.linspace(-20, 20, 41),
                         y=np.linspace(-20, 0, 21))

    grav = Gravimetry(verbose=True)

    model = grav.invert(gravPoints, dz, verbose=1, mesh=mesh)

    plt.plot(pg.x(gravPoints), dz)
    plt.plot(pg.x(gravPoints), grav.inv.response())

    paraDomain=grav.fop.regionManager().paraDomain()
    pg.show(paraDomain, model, colorBar=1, hold=1)

    pg.showNow()
    plt.show()
    pass

def animateGravimetry(mesh, dDens, gravPoints, dz):
    dpi=92
    orientation = 'horizontal'
    fig = plt.figure(facecolor='white', figsize=(2*800/dpi, 2*490/dpi), dpi=dpi)

    axDDe = fig.add_subplot(3,1,1)
    axGra = fig.add_subplot(3,1,2)

    #axGra = fig.fig.add_subplot(1,3,2)
        # ** Density **

    gciDDe = pg.viewer.mpl.drawModel(axDDe, mesh, data=dDens[1],
                                    cMin=0, cMax=20,
                                    )
    cbar = createColorbar(gciDDe, orientation=orientation,
                          label='Delta density in kg/m$^3$')

    def ani(i):
        axGra.clear()
        axGra.plot(pg.x(gravPoints), dz[i])
        axGra.plot(pg.x(gravPoints), pg.y(gravPoints), 'v', color='black')
        axGra.set_ylabel('Grav in mGal')
        axGra.set_xlim((-20, 20))
        axGra.set_ylim((0, 0.001))
        axGra.grid()

        pg.viewer.mpl.setMappableData(gciDDe, abs(dDens[i]),
                                     cMin=0, cMax=20,
                                     logScale=False)
    for i in range(len(dz)):
        ani(i)
        plt.pause(0.001)


if __name__ == "__main__":
    pass