import sys
import glob
import os
import numpy as np
import matplotlib.pyplot as plt
import pygimli as pg
from pygimli.utils import showStitchedModels
from . mrs import MRS
from . modelling import MRS1dBlockQTModelling


class NDMatrix(pg.RBlockMatrix):  # to be moved to a better place
    """
    Diagonal block (block-Jacobi) matrix derived from pg.BlockMatrix
    (to be moved to a better place at a later stage)
    """
    def __init__(self, num, nrows, ncols):
        super(type(self), self).__init__()  # call inherited init function
        self.Ji = []  # list of individual block matrices
        for i in range(num):
            self.Ji.append(pg.RMatrix())
            self.Ji[-1].resize(nrows, ncols)
            n = self.addMatrix(self.Ji[-1])
            self.addMatrixEntry(n, nrows * i, ncols * i)

        self.recalcMatrixSize()
        print(self.rows(), self.cols())


class JointMRSModelling(MRS1dBlockQTModelling):
    """ MRS Laterally constrained modelling based on BlockMatrices """
    def __init__(self, mrs, nlay=2, verbose=False):
        """ Parameters: FDEM data class and number of layers """
        super(JointMRSModelling, self).__init__(nlay, mrs[0].K, mrs[0].z,
                                                mrs[0].t, verbose)
        self.mrs = mrs
        for mrsi in mrs:
            mrsi.createFOP(nlay)

    def response(self, model):
        """ response function (all responses together) """
        response = pg.RVector()
        for mrs in self.mrs:
            response = pg.cat(response, mrs.f(model))


class MRSLCI(pg.ModellingBase):
    """ MRS Laterally constrained modelling based on BlockMatrices """
    def __init__(self, profile, nlay=2, verbose=False):
        """ Parameters: FDEM data class and number of layers """
        super(MRSLCI, self).__init__(verbose)
        self.nlay = nlay
        self.nx = len(profile)
        npar = 3 * nlay - 1
        self.mesh2d = pg.Mesh()
        self.mesh2d.create2DGrid(range(npar+1), range(self.nx+1))
        self.setMesh(self.mesh2d)

        # self.J = NDMatrix(self.nx, self.nf*2, npar)
        self.J = pg.RBlockMatrix()
        self.FOP1d = []
        ipos = 0
        for i, mrs in enumerate(profile):
            mrs.createFOP(nlay, verbose=False)
            self.FOP1d.append(mrs.f)
            n = self.J.addMatrix(self.FOP1d[-1].jacobian())
            self.J.addMatrixEntry(n, ipos, npar * i)
            ipos += len(mrs.data)

        self.J.recalcMatrixSize()
        print(self.J.rows(), self.J.cols())
        self.setJacobian(self.J)

    def response(self, model):
        """ cut-together forward responses of all soundings """
        modA = np.reshape(model, (self.nx, self.nlay * 3 - 1))
        resp = pg.RVector(0)
        for i, modi in enumerate(modA):
            resp = pg.cat(resp, self.FOP1d[i].response(modi))

        return resp

    def createJacobian(self, model):
        modA = np.reshape(model, (self.nx, self.nlay * 3 - 1))
        for i in range(self.nx):
            self.FOP1d[i].createJacobian(modA[i])


class MRSprofile():
    """ manager class for several MRS data along a profile """
    def __init__(self, filename, x=None):
        self.mrs = []
        if '*' in filename:
            filename = glob.glob(filename)
        elif os.path.isdir(filename):
            filename = glob.glob(filename+'/*.mrsi')
        self.load(filename)
        self.x = range(len(self.mrs))
        if x is not None:
            self.x = x

    def __repr__(self):
        return "MRS profile with "+str(len(self.mrs))+" soundings"

    def load(self, filenames, usereal=True):
        """ load mrs files in a list of MRS handlers """
        self.mrs = [MRS(filename, usereal=usereal) for filename in filenames]

    def independentBlock1dInversion(self, nlay=2, startmodel=None):
        """ independent inversion of all soundings """
        self.WMOD, self.TMOD = [], []
        self.RMSvec, self.Chi2vec, self.nData = [], [], []
        for mrs in self.mrs:
            mrs.run(nlay=nlay, startvec=startmodel, verbose=False)
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

    def showFits(self):
        """ show single fits and total fit """
        np.set_printoptions(precision=2)
        print("Single RMS [nV]:")
        print(np.array(self.RMSvec))
        print("Single Chi^2:")
        print(np.array(self.Chi2vec))
        print('Total RMS/Chi^2 value:')
        print(np.round(self.totalRMS, 2), np.round(self.totalChi2, 2))

    def block1dInversion(self, nlay=2, verbose=True):
        """ invert all data together by a 1D model """
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

        self.mrsall.run(nlay, verbose=verbose)
        if verbose:
            self.mrsall.showResult()
        return self.mrsall.model

    def blockLCInversion(self, nlay=2, startModel=None, lam=100.):
        """ laterally constrained (piece-wise 1D) block inversion """
        data, error = pg.RVector(), pg.RVector()
        for mrs in self.mrs:
            data = pg.cat(data, mrs.data)
            error = pg.cat(error, mrs.error)

        fop = MRSLCI(self.mrs, nlay=nlay)
        fop.region(0).setZWeight(0.0)
        fop.region(0).setConstraintType(1)
        transData, transMod = pg.RTrans(), pg.RTransLog()  # LU(1., 500.)
        if startModel is None:
            startModel = self.block1dInversion(nlay, verbose=False)
        model = np.tile(startModel, len(self.mrs))
        INV = pg.RInversion(data, fop, transData, transMod, True, False)
        INV.setModel(model)
        INV.setReferenceModel(model)
        INV.setAbsoluteError(error)
        INV.setLambda(lam)
        INV.stopAtChi1(False)
        model = INV.run()
        # %%
        self.WMOD, self.TMOD = [], []
        for par in np.reshape(model, (len(self.mrs), nlay*3-1)):
                thk = par[0:nlay-1]
                self.WMOD.append(np.hstack((thk, par[nlay-1:2*nlay-1])))
                self.TMOD.append(np.hstack((thk, par[2*nlay-1:3*nlay-1])))

    def showModel(self):
        # %%
        fig, ax = plt.subplots(nrows=2, figsize=(15, 10))
        showStitchedModels(self.WMOD, ax=ax[0], islog=False,
                           title='water content')
        showStitchedModels(self.TMOD, ax=ax[1], title='decay time')
        return fig, ax

if __name__ == "__main__":
    name = sys.argv[-1]
    nlay = 3
    profile = MRSprofile(name)
    print(profile)
    model = profile.block1dInversion(nlay=nlay)
    profile.independentBlock1dInversion(nlay=2)
    profile.showFits()
    fig, ax = profile.showModel()
    fig.savefig(name+'-Ind-N'+str(nlay)+'.pdf', bbox_inches='tight')
    profile.blockLCInversion(nlay=nlay, lam=100, startModel=model)
    fig, ax = profile.showModel()
    fig.savefig(name+'-LCI-N'+str(nlay)+'.pdf', bbox_inches='tight')
