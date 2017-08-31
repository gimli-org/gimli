# -*- coding: utf-8 -*-
"""pyGIMLi - Inversion Frameworks.

These are basic inversion frameworks that usually need a forward operator to run.
"""
import numpy as np

import pygimli as pg

class Inversion(object):
    """Basic inversion framework.
    Attributes
    ----------
    verbose : bool
        Give verbose output
    debug : bool
        Give debug output
    """
    def __init__(self, **kwargs):
        self.__verbose = kwargs.pop('verbose', False)
        self.__debug = kwargs.pop('debug', False)

        self.dataVals = None
        self.errorVals = None

        self.transData = pg.RTransLin()

        self.inv = pg.Inversion(self.__verbose, self.__debug)

        self.maxIter = kwargs.pop('maxIter', 20)

        fop = kwargs.pop('fop', None)
        if fop is not None:
            self.setForwardOperator(fop)

        self.inv.setDeltaPhiAbortPercent(0.5)

    @property
    def verbose(self):
        return self.__verbose
    @verbose.setter
    def verbose(self, v):
        self.__verbose = v
        if self.inv is not None:
            self.inv.setVerbose(self.__verbose)

    @property
    def debug(self):
        return self.__debug
    @debug.setter
    def debug(self, v):
        self.__debug = v
        if self.inv is not None:
            self.inv.setDosave(self.__debug)

    @property
    def maxIter(self):
        return self.inv.maxIter()
    @maxIter.setter
    def maxIter(self, v):
        if self.inv is not None:
            self.inv.setMaxIter(v)

    @property
    def response(self):
        if len(self.inv.response()) > 0:
            return self.inv.response()
        else:
            raise Exception("There was no inversion run so there is no response yet")
    @property
    def model(self):
        if len(self.inv.model()) > 0:
            return self.inv.model()
        else:
            raise Exception("There was no inversion run so there is last model")

    # backward compatibility
    @property
    def dataErrs(self):
        return self.errorVals
    @dataErrs.setter
    def dataErrs(self, v):
        self.errorVals = v

    def setDeltaChiStop(self, it):
        self.inv.setDeltaPhiAbortPercent(it)

    def setForwardOperator(self, fop):
        self.fop = fop
        self.inv.setForwardOperator(fop)

    def setData(self, data):
        if isinstance(data, pg.DataContainer):
            raise Exception("should not been here .. its Managers job")
            self.fop.setData(data)
        else:
            self.dataVals = data

    def setError(self, err):
        self.errorVals = err

    def run(self, dataVals, errorVals, **kwargs):
        """Run inversion.
        """
        self.verbose = kwargs.pop('verbose', self.verbose)
        self.maxIter = kwargs.pop('maxIter', self.maxIter)

        showProgress = kwargs.pop('showProgress', False)

        lam = kwargs.pop('lam', 20)

        self.setData(dataVals)
        self.setError(errorVals)

        if self.dataVals is None:
            raise Exception("Inversion framework need data values to run")

        if self.errorVals is None:
            raise Exception("Inversion framework need data error values to run")

        self.inv.setTransModel(self.fop.transModel)
        self.inv.setTransData(self.transData)

        self.inv.setData(self.dataVals)
        self.inv.setRelativeError(self.errorVals)
        self.inv.setLambda(lam)

        maxIter = self.maxIter
        self.maxIter = 1

        if self.verbose:
            print("inv.start()")

        self.inv.start()
        self.maxIter = maxIter

        if showProgress:
            self.showProgress()

        lastChi2 = self.inv.chi2()
        chi2History = [lastChi2]

        for i in range(1, self.maxIter):
            if self.verbose:
                print("inv.iter", i, "...", end='')

            self.inv.oneStep()
            resp = self.inv.response()
            chi2 = self.inv.chi2()

            if showProgress:
                self.showProgress()

            self.inv.setLambda(self.inv.getLambda() * self.inv.lambdaFactor())

            if self.inv.robustData():
                self.inv.robustWeighting()

            if self.inv.blockyModel():
                self.inv.constrainBlocky()

            chi2History.append(chi2)

            if self.verbose:
                print("chi² = ", round(chi2,2), "lam:", self.inv.getLambda())

            if chi2 < 1:
                if self.verbose:
                    print("Abbort criteria reached: chi² < 1")
                break

            if chi2 < lastChi2 and lastChi2/chi2 < (1.0 + self.inv.deltaPhiAbortPercent()/100):
                if self.verbose:
                    print("Abbort criteria reached: dChi²=",
                          round((1-lastChi2/chi2) * 100, 2),
                         "(", self.inv.deltaPhiAbortPercent(), '%)')
                break

            lastChi2 = chi2

        if len(kwargs.keys()) > 0:
            print("Warning! unhandled keyword arguments", kwargs)

        return self.inv.model()

    def showProgress(self):
        """Called if showProgress=True is set for the inversion run."""
        raise Exception("Implement me in derived classes", self)


class MarquardtInversion(Inversion):
    """Marquardt scheme (local damping with decreasing regularization strength
    """
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.inv.setLocalRegularization(True)
        self.inv.stopAtChi1(False)
        self.inv.setLambdaFactor(0.9)

    def run(self, data, error, **kwargs):

        self.fop.regionManager().setConstraintType(0)

        return super(MarquardtInversion, self).run(data, error, **kwargs)


class Block1DInversion(MarquardtInversion):
    def __init__(self, **kwargs):
        self.axs = None
        super().__init__(**kwargs)

    def run(self, dataVals, errVals, nLayer=4, **kwargs):

        if len(self.fop.startModel()) == 0:
            self.fop.createStartModel(dataVals, nLayer)

        return super(Block1DInversion, self).run(dataVals, errVals, **kwargs)

    def showProgress(self):

        if self.axs is None:
            fig, axs = pg.plt.subplots(1, 2)
            self.axs = axs

        ax = self.axs

        ax[0].clear()
        ax[1].clear()
        if hasattr(self.fop, 'drawModel'):
            self.fop.drawModel(ax[0], self.inv.model())
        else:
            pg.mplviewer.drawModel1D(ax=ax[0],
                                     model=self.inv.model(),
                                     plot='loglog',
                                     xlabel='Model parameter')

        if hasattr(self.fop, 'drawData'):
            self.fop.drawData(ax[1], self.dataVals, self.errorVals, label='Data')
            self.fop.drawData(ax[1], self.inv.response(), label='Response')
        else:
            nData = len(self.dataVals)
            yVals = range(nData)
            ax[1].loglog(self.dataVals, yVals, 'rx-')
            ax[1].errorbar(self.dataVals, yVals,
                           xerr=self.errorVals*self.dataVals,
                           linewidth=1, color='red', linestyle='-')
            ax[1].loglog(self.inv.response(), yVals, 'bo-')
            ax[1].set_ylim(max(yVals), min(yVals))
            ax[1].set_xlabel('Data')
            ax[1].set_ylabel('Data Number')

        ax[1].text(0.01, 0.96,
                   "iter: %d, rrms: %.2g, $\chi^2$: %.2g" %
                    (self.inv.iter(), self.inv.relrms(), self.inv.chi2()),
                    transform=ax[1].transAxes)

        pg.plt.pause(0.05)
















class MeshInversion(Inversion):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def setMesh(self, mesh):
        self.fop.setMesh(mesh)

    def invert(self, data=None, mesh=None, lam=20, **kwargs):
        if data is not None:
            self.setData(data)

        if mesh is not None:
            self.setMesh(mesh)

        startModel = kwargs.pop('startModel', None)
        nModel = self.fop.regionManager().parameterCount()
        startModel = pg.Vector(nModel, startModel)
        self.fop.setStartModel(startModel)

        zWeight = kwargs.pop('zWeight', 1.0)
        self.fop.regionManager().setZWeight(zWeight)
        self.inv.setData(self.dataVals)
        self.inv.setRelativeError(self.errorVals)
        self.inv.setLambda(lam)

        self.mod = self.inv.run()
        self.mod = self.mod(self.fop.regionManager().paraDomain().cellMarkers())
        return self.mod
