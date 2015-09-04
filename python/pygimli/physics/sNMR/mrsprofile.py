import sys
import glob
import os
import numpy as np
import matplotlib.pyplot as plt
import pygimli as pg
from pygimli.utils import showStitchedModels
from . mrs import MRS
from . modelling import MRS1dBlockQTModelling
from matplotlib.cm import Spectral


class JointMRSModelling(MRS1dBlockQTModelling):
    """MRS Laterally constrained modelling based on BlockMatrices"""
    def __init__(self, mrs, nlay=2, verbose=False):
        """Parameters: FDEM data class and number of layers"""
        super(JointMRSModelling, self).__init__(nlay, mrs[0].K, mrs[0].z,
                                                mrs[0].t, verbose)
        self.mrs = mrs
        for mrsi in mrs:
            mrsi.createFOP(nlay)

    def response(self, model):
        """response function (all responses together)"""
        response = pg.RVector()
        for mrs in self.mrs:
            response = pg.cat(response, mrs.f(model))


class MRSLCI(pg.ModellingBase):
    """MRS Laterally constrained modelling based on BlockMatrices"""
    def __init__(self, profile, nlay=2, verbose=False):
        """Parameters: FDEM data class and number of layers"""
        super(MRSLCI, self).__init__(verbose)
        self.nlay = nlay
        self.nx = len(profile)
        self.np = 3 * nlay - 1
#        self.mesh2d = pg.createMesh2D(npar, self.nx)
#        self.mesh2d = pg.createMesh2D(self.nx, self.np)
        self.mesh2d = pg.createMesh2D(range(self.np+1), range(self.nx+1))
        self.mesh2d.rotate(pg.RVector3(0, 0, -np.pi/2))
        self.setMesh(self.mesh2d)
        self.J = pg.RBlockMatrix()
        self.FOP1d = []
        ipos = 0
        for i, mrs in enumerate(profile):
            mrs.createFOP(nlay, verbose=False)
            self.FOP1d.append(mrs.f)
            n = self.J.addMatrix(self.FOP1d[-1].jacobian())
            self.J.addMatrixEntry(n, ipos, self.np * i)
            ipos += len(mrs.data)

        self.J.recalcMatrixSize()
        print(self.J.rows(), self.J.cols())
        self.setJacobian(self.J)

    def response(self, model):
        """cut-together forward responses of all soundings"""
        modA = np.asarray(model).reshape(self.nx, self.np)
#        modA = np.reshape(model, (self.nx, self.np))
        resp = pg.RVector(0)
        for i, modi in enumerate(modA):
            resp = pg.cat(resp, self.FOP1d[i].response(modi))

        return resp

    def createJacobian(self, model):
        """ fill individual blocks of Block-Jacobian matrix """
        modA = np.asarray(model).reshape(self.nx, self.np)
#        modA = np.reshape(model, (self.nx, self.np))
        for i in range(self.nx):
            self.FOP1d[i].createJacobian(modA[i])


class MRSprofile():
    """manager class for several MRS data along a profile for joint inversion

    Attributes
    ----------
    mrs : list of MRS objects (single soundings)

    x : list of positions for the soundings


    Methods
    -------
    load - load mrs files from a directory

    set X - set x vector

    showData - show MRS data

    independentBlock1dInversion - perform independent 1D block inversion

    block1dInversion - 1D block inversion of all data sets together

    blockLCInversion - 1D block laterally constrained inversion of all data

    printFits - print total misfit (chi^2, rms) and individual values

    showModel - show LCI model
    """
    def __init__(self, filename, x=None, dx=1, x0=0, **kwargs):
        self.mrs = []
        if '*' in filename:
            filename = glob.glob(filename)
        elif os.path.isdir(filename):
            filename = glob.glob(filename+'/*.mrsi')
        self.load(filename, **kwargs)
        if x is not None:
            self.setX(x)
        else:
            self.x = np.arange(len(self.mrs)) * dx + x0

    def __repr__(self):
        return "MRS profile with "+str(len(self.mrs))+" soundings"

    def setX(self, x=None, x0=0, dx=1):
        """define positions for soundings and sort accordingly"""
        if x is None:
            x = np.arange(len(self.mrs)) * dx + x0
            print(x)
        ind = np.argsort(x)
        self.mrs = self.mrs(ind)
        self.x = np.sort(x)

    def load(self, filenames, **kwargs):
        """ load mrs files in a list of (single) MRS handlers
            filename can be a list of mrsi files or a directory to search
            Additional parameters: usereal, mint, maxt (see MRS.load)
        """
        self.mrs = [MRS(filename, **kwargs) for filename in filenames]

    def showData(self, figsize=(15, 10), clim=None):
        """show all data cubes in subplots"""
        from math import sqrt, ceil
        nsond = len(self.mrs)
        nc = ceil(sqrt(nsond*3))
        nr = ceil(nsond/nc)
        fig, ax = plt.subplots(nrows=int(nr), ncols=int(nc), figsize=figsize)
        for i, mrs in enumerate(self.mrs):
            mrs.showCube(ax=ax.flat[i], vec=mrs.data*1e9, islog=False,
                         clim=clim)
        return fig, ax

    def independentBlock1dInversion(self, nlay=2, lam=100, startModel=None):
        """ independent inversion of all soundings
            Parameters: nlay, lam, startModel (see MRS.run parameters)
        """
        self.WMOD, self.TMOD = [], []
        self.RMSvec, self.Chi2vec, self.nData = [], [], []
        for i, mrs in enumerate(self.mrs):
            if hasattr(startModel[0], '__iter__'):
                startvec = startModel[i]
            else:
                startvec = startModel
            mrs.run(nlay=nlay, startvec=startvec, lam=lam, verbose=False)
            mrs.INV.echoStatus()
            thk, wc, t2 = mrs.result()
            self.WMOD.append(np.hstack((thk, wc)))
            self.TMOD.append(np.hstack((thk, t2)))
            self.RMSvec.append(mrs.INV.absrms()*1e9)  # in nV
            self.Chi2vec.append(mrs.INV.chi2())
            self.nData.append(len(mrs.data))

        self.totalChi2 = sum(np.array(self.Chi2vec)*np.array(self.nData)) / \
            sum(self.nData)
        self.totalRMS = np.sqrt(sum(np.array(self.RMSvec)**2 *
                                    np.array(self.nData)) / sum(self.nData))

    def block1dInversion(self, nlay=2, startModel=None, verbose=True):
        """invert all data together by one 1D model (variant 1 - all equal)"""
        self.mrsall = MRS()
        self.mrsall.z = self.mrs[0].z
        self.mrsall.t = self.mrs[0].t
        self.mrsall.data, self.mrsall.error, self.mrsall.q = [], [], []
        self.mrsall.K = np.zeros((0, self.mrs[0].K.shape[1]))
        for mrs in self.mrs:
            self.mrsall.data = np.hstack((self.mrsall.data, mrs.data))
            self.mrsall.error = np.hstack((self.mrsall.error, mrs.error))
            self.mrsall.q = np.hstack((self.mrsall.q, mrs.q))
            self.mrsall.K = np.vstack((self.mrsall.K, mrs.K))

        self.mrsall.run(nlay, stvec=startModel, verbose=verbose)
        if verbose:
            self.mrsall.showResult()
        return self.mrsall.model

    def block1dInversionNew(self, nlay=2, lam=100., verbose=True):
        """invert all data together by a 1D model (more general solution)"""
        data, error = pg.RVector(), pg.RVector()
        for mrs in self.mrs:
            data = pg.cat(data, mrs.data)
            error = pg.cat(error, mrs.error)

        f = JointMRSModelling(self.mrs, nlay)
        mrsobj = self.mrs[0]
        for i in range(3):
            f.region(i).setParameters(mrsobj.startval[i], mrsobj.lowerBound[i],
                                      mrsobj.upperBound[i])

        INV = pg.RInversion(data, f, verbose)
        INV.setLambda(lam)
        INV.setMarquardtScheme(0.8)
#        INV.stopAtChi1(False)  # should be already in MarquardtScheme
        INV.setDeltaPhiAbortPercent(0.5)
        INV.setAbsoluteError(error)
        model = INV.run()
        return model

    def blockLCInversion(self, nlay=2, startModel=None, **kwargs):
        """laterally constrained (piece-wise 1D) block inversion"""
        data, error, self.nData = pg.RVector(), pg.RVector(), []
        for mrs in self.mrs:
            data = pg.cat(data, mrs.data)
            error = pg.cat(error, mrs.error)
            self.nData.append(len(mrs.data))

        fop = MRSLCI(self.mrs, nlay=nlay)
        fop.region(0).setZWeight(kwargs.pop('zWeight', 0))
        fop.region(0).setConstraintType(kwargs.pop('cType', 1))
        transData, transMod = pg.RTrans(), pg.RTransLog()  # LU(1., 500.)
        if startModel is None:
            startModel = self.block1dInversion(nlay, verbose=False)
        model = kwargs.pop('startvec', np.tile(startModel, len(self.mrs)))
        INV = pg.RInversion(data, fop, transData, transMod, True, False)
        INV.setModel(model)
        INV.setReferenceModel(model)
        INV.setAbsoluteError(error)
        INV.setLambda(kwargs.pop('lam', 100))
        INV.setMaxIter(kwargs.pop('maxIter', 20))
        INV.stopAtChi1(False)
        INV.setLambdaFactor(0.9)
        INV.setDeltaPhiAbortPercent(0.1)
        model = INV.run()
        self.WMOD, self.TMOD = [], []
        for par in np.reshape(model, (len(self.mrs), 3*nlay-1)):
                thk = par[0:nlay-1]
                self.WMOD.append(np.hstack((thk, par[nlay-1:2*nlay-1])))
                self.TMOD.append(np.hstack((thk, par[2*nlay-1:3*nlay-1])))

        ind = np.hstack((0, np.cumsum(self.nData)))
        resp = INV.response()
        misfit = data - resp
        emisfit = misfit / error
        misfit *= 1e9
        self.totalChi2 = INV.chi2()
        self.totalRMS = INV.absrms()*1e9
        self.RMSvec, self.Chi2vec = [], []
        for i in range(len(self.mrs)):
            self.RMSvec.append(np.sqrt(np.mean(misfit[ind[i]:ind[i+1]]**2)))
            self.Chi2vec.append(np.mean(emisfit[ind[i]:ind[i+1]]**2))

    def printFits(self):
        """show single fits and total fit"""
        np.set_printoptions(precision=2)
        print("Single RMS [nV]:")
        print(np.array(self.RMSvec))
        print("Single Chi^2:")
        print(np.array(self.Chi2vec))
        print('Total RMS/Chi^2 value:')
        print(np.round(self.totalRMS, 2), np.round(self.totalChi2, 2))

    def showModel(self, showFit=0, cmap=Spectral, figsize=(13, 12),
                  wlim=[None, None], tlim=[None, None]):
        """ show 2d model as stitched 1d models along with fit"""
        fig, ax = plt.subplots(nrows=2+showFit, figsize=figsize, sharex=True)
        showStitchedModels(self.WMOD, x=self.x, ax=ax[-2], islog=False,
                           cmap=cmap, cmin=wlim[0], cmax=wlim[1],
                           title=r'$\theta$ (-)')
        showStitchedModels(self.TMOD, x=self.x, ax=ax[-1], cmap=cmap,
                           cmin=tlim[0], cmax=tlim[1], title=r'$T_2^*$ (s)')
        xl = ax[-1].get_xlim()
        ax[0].set_xlabel('x (m)')
        ax[0].xaxis.set_label_position('top')
        ax[0].xaxis.tick_top()
        ax[0].set_ylabel(r'$\chi^2$ (-)')
        for axi in ax[-2:]:
            axi.set_ylabel('z (m)')
        if showFit > 0:
            ax0b = ax[0].twinx()
            ax0b.set_ylabel('rms (nV)')
            ax[0].plot(self.x, self.Chi2vec, 'rx', label=r'$\chi^2$')
            ax0b.plot(self.x, self.RMSvec, 'bx', label='rms [nV]')
            ax[0].legend(numpoints=1, loc=2)
            ax0b.legend(numpoints=1, loc=1)
            ax[0].set_xlim(xl)

        return fig, ax

if __name__ == "__main__":
    name = sys.argv[-1]
    nlay = 3
    profile = MRSprofile(name)  # directory or list of names
    print(profile)
    fig, ax = profile.showData()  # subplots with data cubes
    model = profile.block1dInversion(nlay=nlay)  # all in one model
    profile.independentBlock1dInversion(nlay=2)  # each separately
    profile.printFits()
    fig, ax = profile.showModel(showFit=1)
    fig.savefig(name+'-Ind-N'+str(nlay)+'.pdf', bbox_inches='tight')
    profile.blockLCInversion(nlay=nlay, lam=100, startModel=model)
    profile.printFits()
    fig, ax = profile.showModel(showFit=1)
    fig.savefig(name+'-LCI-N'+str(nlay)+'.pdf', bbox_inches='tight')
