#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""TODO DOCUMENTME."""
import time
import random
import sys

import numpy as np
import matplotlib.pyplot as plt

import pygimli as pg
from pygimli.physics.sNMR import MRS, MRS1dBlockQTModelling
from pygimli.mplviewer import drawModel1D
from pygimli.utils import iterateBounds


def correctBranches(ab2, mn2, rhoa):
    """Shifts the branches of a dc sounding to generate a matching curve."""
    um = np.unique(mn2)
    for i in range(len(um) - 1):
        r0, r1 = [], []
        ac = np.intersect1d(ab2[mn2 == um[i]], ab2[mn2 == um[i + 1]])
        for a in ac:
            r0.append(rhoa[(ab2 == a) * (mn2 == um[i])][0])
            r1.append(rhoa[(ab2 == a) * (mn2 == um[i + 1])][0])

        if len(r0) > 0:
            fak = np.mean(np.array(r0) / np.array(r1))
            print("branch correction with factor", fak)
            if np.isfinite(fak) and fak > 0.:
                rhoa[mn2 == um[i + 1]] *= fak


class MRSVESBlockModelling(pg.ModellingBase):
    """Joint MRS(QT) and VES forward modelling"""
    def __init__(self, nlay, K, z, t, ab2, mn2, verbose=False):
        self.nlay_ = nlay
        mesh = pg.createMesh1DBlock(nlay, 3)  # 3 paramaters: thk,wc,t2,res
        pg.ModellingBase.__init__(self, mesh)
        self.fMRS = MRS1dBlockQTModelling(nlay, K, z, t)
        self.fVES = pg.DC1dModelling(nlay, ab2, mn2)

    def response(self, par):
        nl = self.nlay_
        mVES = pg.cat(par(0, nl - 1), par(nl * 3 - 1, nl * 4 - 1))  # thk+res
        mMRS = par(0, 3 * nl - 1)  # thk+wc+t2
        return pg.cat(self.fMRS.response(mMRS), self.fVES.response(mVES))


class MRSVES(MRS):
    """Class for joint MRS/VES inversion using LS or GA"""
    def __init__(self, mrsfile, vesfile=None, correctbranch=False):
        """init class and load MRS (*.mrsi) (and VES data *.ves)"""
        MRS.__init__(self, mrsfile)
        self.lowerBound.append(1.)
        self.upperBound.append(3000.)
        self.startval.append(300.)
        if vesfile is not None:
            self.loadVES(vesfile, correctbranch)

    def loadVES(self, vesfile, correctbranches=True):
        """
            Loads vertical electrical sounding (VES)

            Loads VES from file
            with columns AB/2, MN/2 and apparent resistivity
            loadVES( vesfile, correctbranches=True) corrects branches
        """
        self.ab2, self.mn2, self.rhoa = np.loadtxt(vesfile,
                                                   usecols=(0, 1, 2),
                                                   unpack=True)
        self.startval[3] = np.median(self.rhoa)

        if correctbranches:
            correctBranches(self.ab2, self.mn2, self.rhoa)

    def createFOP(self, nlay, verbose=False):
        """creates forward operator for given number of layers"""
        self.nlay = nlay
        self.f = MRSVESBlockModelling(
            nlay,
            self.K,
            self.z,
            self.t,
            self.ab2,
            self.mn2,
            verbose)
        self.trans = []  # save in class to keep them alive
        for i in range(4):  # thk, cw, T2*, res
            self.f.region(i).setStartValue(self.startval[i])
            self.trans.append(
                pg.RTransLogLU(
                    self.lowerBound[i],
                    self.upperBound[i]))
            self.f.region(i).setTransModel(self.trans[-1])

    def createInv(self, nlay, lam=100., errVES=3, verbose=True):
        """Create Marquardt type inversion instance with data transformatio"""
        self.createFOP(nlay)
        self.tMod = pg.RTransLog()
        self.tMRS = pg.RTrans()
        self.tVES = pg.RTransLog()
        self.transData = pg.RTransCumulative()
        self.transData.push_back(self.tMRS, len(self.data))
        self.transData.push_back(self.tVES, len(self.rhoa))
        data = pg.cat(self.data, self.rhoa)
        self.INV = pg.RInversion(data, self.f, self.transData, verbose)
        self.INV.setLambda(lam)
        self.INV.setMarquardtScheme(0.8)
        self.INV.stopAtChi1(False)  # now in MarquardtScheme
        self.INV.setDeltaPhiAbortPercent(0.5)
#        self.INV.setMaxIter(1)
        error = pg.cat(self.error, self.rhoa * errVES / 100.)
        self.INV.setAbsoluteError(error)

    def runInv(self, uncertainty=False):
        """run actual inversion (assumes one was created before)"""
        self.model = np.array(mrsves.INV.run())
        if uncertainty:
            self.modelL, self.modelU = iterateBounds(self.INV,
                                                     dchi2=self.INV.chi2() / 2, change=1.2)

    def plotResult(self, filename=None, nrows=1, figsize=(10, 6)):
        """plot Result as three (time 1 or more) plots"""
        nl = self.nlay
        thk = self.model[:nl - 1]
        wc = self.model[nl - 1:2 * nl - 1]
        t2 = self.model[2 * nl - 1:3 * nl - 1]
        res = self.model[3 * nl - 1:4 * nl - 1]
        fig, ax = plt.subplots(
            nrows=nrows, ncols=3, sharey=(
                nrows == 1), figsize=figsize, squeeze=False)
        drawModel1D(ax[0, 0], thk, wc * 100., xlabel=r'$\theta$ [%]')
        drawModel1D(ax[0,
                       1],
                    thk,
                    t2 * 1e3,
                    plotfunction='semilogx',
                    xlabel=r'$T_2^*$ [ms]')
        drawModel1D(ax[0, 2], thk, res, plotfunction='semilogx')

        if self.modelL is not None and self.modelU is not None:
            thkL, thkU = self.modelL[:nl - 1], self.modelU[:nl - 1]
            wcL, wcU = self.modelL[nl -
                                   1:2 *
                                   nl -
                                   1], self.modelU[nl -
                                                   1:2 *
                                                   nl -
                                                   1]
            t2L, t2U = self.modelL[2 *
                                   nl -
                                   1:3 *
                                   nl -
                                   1], self.modelU[2 *
                                                   nl -
                                                   1:3 *
                                                   nl -
                                                   1]
            resL, resU = self.modelL[3 *
                                     nl -
                                     1:4 *
                                     nl -
                                     1], self.modelU[3 *
                                                     nl -
                                                     1:4 *
                                                     nl -
                                                     1]
            zc = np.cumsum(thk)
            zm = np.hstack((zc - thk / 2, np.sum(thk) + 3.))
            ax[0, 0].errorbar((wc[:-1] + wc[1:]) / 2 * 100., zc,
                              fmt='.', yerr=np.vstack((thk - thkL, thkU - thk)))
            ax[0, 0].errorbar(
                wc * 100., zm, fmt='.', xerr=np.vstack((wc - wcL, wcU - wc)) * 100.)
            ax[0, 1].set_xlim(
                self.lowerBound[1] * 100., self.upperBound[1] * 100.)
            ax[0, 1].errorbar(np.sqrt(
                t2[:-1] * t2[1:]) * 1e3, zc, fmt='.', yerr=np.vstack((thk - thkL, thkU - thk)))
            ax[0, 1].errorbar(
                t2 * 1e3, zm, fmt='.', xerr=np.vstack((t2 - t2L, t2U - t2)) * 1e3)
            ax[0, 1].set_xlim(
                self.lowerBound[2] * 1e3, self.upperBound[2] * 1e3)
            ax[0, 2].errorbar(np.sqrt(res[:-1] * res[1:]), zc,
                              fmt='.', yerr=np.vstack((thk - thkL, thkU - thk)))
            ax[0, 2].errorbar(
                res, zm, fmt='.', xerr=np.vstack((res - resL, resU - res)))
            ax[0, 2].set_xlim(self.lowerBound[3], self.upperBound[3])
            ax[0, 2].set_xlabel(r'$\rho$ [$\Omega$m]')
        if filename is not None:
            fig.savefig(filename, bbox_inches='tight')

        return fig, ax

    def plotResultAndFit(self, filename=None, figsize=(10, 10)):
        """Whats this?"""
        fig, ax = self.plotResult(nrows=2, figsize=figsize)
        clim = self.showCube(ax[1, 0], self.data * 1e9, islog=False)
        resp = np.array(self.INV.response())
        self.showCube(
            ax[1, 1], resp[:len(self.data)] * 1e9, islog=False, clim=clim)
        ax[1, 2].loglog(self.rhoa, self.ab2, 'bx')
        ax[1, 2].loglog(resp[len(self.data):], self.ab2, 'r-')
        ax[1, 2].set_ylim((max(self.ab2), min(self.ab2)))
        ax[1, 2].grid(True)

        if filename is not None:
            fig.savefig(filename, bbox_inches='tight')
        return fig, ax

    def exportResult(self, basename):
        """export result in column file (z,thk,wc,t2,res)"""
        nl = self.nlay
        thk = self.model[:nl - 1]
        wc = self.model[nl - 1:2 * nl - 1]
        t2 = self.model[2 * nl - 1:3 * nl - 1]
        res = self.model[3 * nl - 1:4 * nl - 1]
        z = np.hstack((0., np.cumsum(thk)))
        thk0 = np.hstack((thk, 0.))
        ALL = np.column_stack((z, thk0, wc, t2, res))
        print(ALL)
        np.savetxt(basename + '-result.txt', ALL)
#        if hasattr(self,'modelL') and hasattr(self,'modelU'):

        if self.modelL is not None and self.modelU is not None:
            thkL, thkU = self.modelL[:nl - 1], self.modelU[:nl - 1]
            wcL, wcU = self.modelL[nl -
                                   1:2 *
                                   nl -
                                   1], self.modelU[nl -
                                                   1:2 *
                                                   nl -
                                                   1]
            t2L, t2U = self.modelL[2 *
                                   nl -
                                   1:3 *
                                   nl -
                                   1], self.modelU[2 *
                                                   nl -
                                                   1:3 *
                                                   nl -
                                                   1]
            resL, resU = self.modelL[3 *
                                     nl -
                                     1:4 *
                                     nl -
                                     1], self.modelU[3 *
                                                     nl -
                                                     1:4 *
                                                     nl -
                                                     1]
            thkL0, thkU0 = np.hstack((thkL, 0.)), np.hstack((thkU, 0.))
            ALL = np.column_stack(
                (z,
                 thk0,
                 wc,
                 t2,
                 res,
                 wcL,
                 t2L,
                 resL,
                 wcU,
                 t2U,
                 resU,
                 thkL0,
                 thkU0))
            np.savetxt(basename + '-resultLU.txt', ALL)

    def runEMO(self, nlay=5, pop_size=100, max_generations=100):
        """
            Run evolutionary multi-objective optimization (EMO)

            Run EMO using inspyred for now fixed to NSGA-II algorithm
            after Deb (2002) TODO (cite correctly)
            (non-dominated sorting genetic algorithm)
        """
        import inspyred

        def genMods(individual):
            """generate MRS and VES models from unit vector"""
            model = individual * (self.lUB - self.lLB) + self.lLB
#            model = pg.asvector(individual) * (self.lUB - self.lLB) + self.lLB
            if self.logpar:
                model = pg.exp(model)

            modMRS = model(0, nlay * 3 - 1)
            modVES = pg.cat(
                model(
                    0,
                    nlay -
                    1),
                model(
                    nlay *
                    3 -
                    1,
                    nlay *
                    4 -
                    1))
            return modMRS, modVES

        def mygenerate(random, args):
            """generate a random vector of model size"""
            return [random.random() for i in range(nlay * 4 - 1)]

        @inspyred.ec.evaluators.evaluator
        def datafit(individual, args):
            """return data fits for MRS and VES as Pareto object"""
            modMRS, modVES = genMods(individual)
            MRSmisfit = (self.data - self.f(modMRS)) / self.error
            VESmisfit = (
                np.log(self.rhoa) - np.log(self.fVES(modVES))) / np.log(1.02)
            return inspyred.ec.emo.Pareto(
                [np.mean(MRSmisfit**2), np.mean(VESmisfit**2)])

        self.createFOP(nlay)
        self.fVES = pg.DC1dModelling(nlay, self.ab2, self.mn2)

        lowerBound = pg.cat(
            pg.cat(
                pg.RVector(
                    nlay - 1,
                    self.lowerBound[0]),
                pg.RVector(
                    nlay,
                    self.lowerBound[1])),
            pg.cat(pg.RVector(nlay, self.lowerBound[2]), pg.RVector(nlay, self.lowerBound[3])))
        upperBound = pg.cat(
            pg.cat(
                pg.RVector(
                    nlay - 1,
                    self.upperBound[0]),
                pg.RVector(
                    nlay,
                    self.upperBound[1])),
            pg.cat(pg.RVector(nlay, self.upperBound[2]), pg.RVector(nlay, self.upperBound[3])))
        if self.logpar:
            self.lLB, self.lUB = pg.log(lowerBound), pg.log(
                upperBound)  # ready mapping functions
        else:
            self.lLB, self.lUB = lowerBound, upperBound

        rand = random.Random()
        rand.seed(int(time.time()))
        ea = inspyred.ec.emo.NSGA2(rand)
        ea.variator = [
            inspyred.ec.variators.blend_crossover,
            inspyred.ec.variators.gaussian_mutation]
        ea.terminator = inspyred.ec.terminators.generation_termination
        ea.observer = [
            inspyred.ec.observers.stats_observer,
            inspyred.ec.observers.file_observer]
        self.pop = ea.evolve(evaluator=datafit,
                             generator=mygenerate,
                             maximize=False,
                             bounder=inspyred.ec.Bounder(0., 1.),
                             pop_size=pop_size,
                             max_generations=max_generations)

if __name__ == "__main__":
    from optparse import OptionParser

    parser = OptionParser(
        "usage: %prog [options] basename",
        version="%prog: " +
        pg.__version__)
    parser.add_option(
        "-v",
        "--verbose",
        dest="verbose",
        action="store_true",
        help="be verbose",
        default=False)
    parser.add_option("-n", "--nLayers", dest="nlay",
                            help="number of layers", type="int", default="5")
    parser.add_option(
        "-G",
        "--nsga2",
        dest="ga",
        action="store_true",
        help="use NSGA-II algorithm",
        default=False)
    parser.add_option(
        "-U",
        "--uncertainty",
        dest="uncertainty",
        action="store_true",
        help="compute uncertainty",
        default=False)

    (options, args) = parser.parse_args()

    if len(args) == 0:
        #datafile = 'example.xyz'
        parser.print_help()
        print("Please add a mesh or model name.")
        sys.exit(2)
    else:
        basename = args[0]

    mrsves = MRSVES(basename + '.mrsi', basename + '.ves')
    if options.ga:  # multi-objective optimization using NSGA-II
        mrsves.lowerBound = (1.0, 0.05, 0.02, 1.)  # d, theta, T2*, res
        mrsves.upperBound = (30., 0.45, 0.60, 100.)  # d, theta, T2*, res
        mrsves.startval = (10., 0.30, 0.20, 10.)  # d, theta, T2*, res
        mrsves.logpar = True  # log-normal distribution of parameters
        mrsves.runEMO(5, popsize=300, maxgen=100)
    else:
        mrsves.createInv(options.nlay)
        mrsves.runInv(options.uncertainty)
        mrsves.exportResult(basename)
        fig, ax = mrsves.plotResultAndFit(filename=basename + '-result.pdf')
        plt.show()
