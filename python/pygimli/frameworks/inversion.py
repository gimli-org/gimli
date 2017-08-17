# -*- coding: utf-8 -*-
"""pyGIMLi - Inversion Frameworks.

These are basic inversion frameworks that usually need a forward operator to run.
"""
import numpy as np

import pygimli as pg

class Inversion(object):
    def __init__(self, **kwargs):
        self.fop = None
        self.dataVals = None
        self.dataErrs = None
        self.transData = pg.RTransLin()
        self.transModel = pg.RTransLog()
        fop = kwargs.pop('fop', None)

        self.inv = pg.Inversion(**kwargs)
        self.inv.setDeltaPhiAbortPercent(0.5)

        if fop is not None:
            self.setForwardOperator(fop)

    def setMaxIter(self, it):
        self.inv.setMaxIter(it)

    def setDeltaChiStop(self, it):
        self.inv.setDeltaPhiAbortPercent(it)

    def setForwardOperator(self, fop):
        self.fop = fop
        self.inv.setForwardOperator(fop)

    def setData(self, data):
        if isinstance(data, pg.DataContainer):
            self.fop.setData(data)
        else:
            self.dataVals = data

    def setError(self, err):
        self.dataErrs = err

    def invert(self, data=None, lam=20, **kwargs):
        pass

    def run(self, data, error, lam=20, **kwargs):

        verbose = kwargs.pop('verbose', False)

        self.inv.setVerbose(verbose)

        showProgress = kwargs.pop('showProgress', False)

        self.setData(data)
        self.setError(error)

        self.inv.setTransData(self.transData)
        self.inv.setTransModel(self.transModel)

        self.inv.setData(self.dataVals)
        self.inv.setRelativeError(self.dataErrs)
        self.inv.setLambda(lam)

        #inv.setOptimizeLambda(1)

        maxIter = self.inv.maxIter()
        self.inv.setMaxIter(1)
        if verbose:
            print("inv.start()")

        self.inv.start()
        self.inv.setMaxIter(maxIter)

        if showProgress:
            self.showProgress()

        lastChi2 = self.inv.chi2()
        chi2History = [lastChi2]
        for i in range(1, self.inv.maxIter()):
            if verbose:
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

            if verbose:
                print("chi² = ", round(chi2,2), "lam:", self.inv.getLambda())

            if chi2 < 1:
                if verbose:
                    print("Abbort criteria reached: chi² < 1")
                break

            if chi2 < lastChi2 and lastChi2/chi2 < (1.0 + self.inv.deltaPhiAbortPercent()/100):
                if verbose:
                    print("Abbort criteria reached: dChi²=",
                          round((1-lastChi2/chi2) * 100, 2),
                         "(", self.inv.deltaPhiAbortPercent(), '%)')
                break

            lastChi2 = chi2

        return self.inv.model()

    def showProgress(self):
        """Called if showProgress=True is set for the inversion run."""
        raise Exception("Implement me in derived classes", self)

class MarquardtInversion(Inversion):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.setLambdaDecrease(0.9)

    def setLambdaDecrease(self, a):
        """Set lambda decrease factor for Marquardt Inversion Scheme"""
        self.inv.setMarquardtScheme(a)


class Block1DInversion(MarquardtInversion):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.axs = None

    def run(self, data, error=None, nLayer=4, **kwargs):

        if len(self.fop.startModel()) == 0:
            self.fop.createStartModel(data, nLayer)

        return super().run(data, error, **kwargs)

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
            self.fop.drawData(ax[1], self.dataVals, self.dataErrs, label='Data')
            self.fop.drawData(ax[1], self.inv.response(), label='Response')
        else:
            nData = len(self.dataVals)
            yVals = range(nData)
            ax[1].loglog(self.dataVals, yVals, 'rx-')
            ax[1].errorbar(self.dataVals, yVals,
                           xerr=self.dataErrs*self.dataVals,
                           linewidth=1, color='red', linestyle='-')
            ax[1].loglog(self.inv.response(), yVals, 'bo-')
            ax[1].set_ylim(max(yVals), min(yVals))
            ax[1].set_xlabel('Data')
            ax[1].set_ylabel('Data Number')

        ax[1].text(0.01, 0.96,
                   "iter: %d, rrms: %.2g, $\chi^2$: %.2g" %
                    (self.inv.iter(), self.inv.relrms(), self.inv.chi2()),
                    transform=ax[1].transAxes)

        pg.plt.pause(0.1)


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
        self.inv.setRelativeError(self.dataErrs)
        self.inv.setLambda(lam)

        self.mod = self.inv.run()
        self.mod = self.mod(self.fop.regionManager().paraDomain().cellMarkers())
        return self.mod
