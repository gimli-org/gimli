#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    Frequency Domain Electromagnetics (FDEM) functions and class
"""

import sys

import pygimli as pg
from pygimli.viewer import show1dmodel, drawModel1D
from pygimli.utils import opt_import

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap


class NDMatrix(pg.RBlockMatrix):

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


def cmapDAERO():
    """standardized colormap from A-AERO projects (purple=0.3 to red=500)"""
    CMY = np.array([
        [127, 255, 31], [111, 255, 47], [95, 255, 63], [79, 255, 79],
        [63, 255, 95], [47, 255, 111], [31, 255, 127], [16, 255, 159],
        [0, 223, 159], [0, 191, 159], [0, 159, 207], [0, 127, 175],
        [0, 95, 175], [0, 63, 175], [0, 47, 175], [0, 31, 191], [0, 0, 255],
        [0, 0, 159], [15, 0, 127], [47, 0, 143], [79, 0, 143], [111, 0, 143],
        [143, 0, 127], [159, 31, 63], [175, 47, 31], [207, 63, 0],
        [223, 111, 0], [231, 135, 0], [239, 159, 0], [255, 191, 47],
        [239, 199, 63], [223, 207, 79], [207, 239, 111]], dtype=float)
    RGB = 1.0 - CMY/255
    return LinearSegmentedColormap.from_list('D-AERO', RGB)


def xfplot(ax, DATA, x, freq, everyx=5, orientation='horizontal', aspect=40):
    """ Plots a matrix according to x and frequencies """
    nt = list(range(0, len(x), everyx))
    im = ax.imshow(DATA.T, interpolation='nearest')
    ax.set_ylim(plt.ylim()[::-1])
    ax.set_xticks(nt)
    ax.set_xticklabels(["%g" % xi for xi in x[nt]])
    ax.set_yticks(list(range(0, len(freq) + 1, 2)))
    ax.set_yticklabels(["%g" % freq[i] for i in range(0, len(freq), 2)])
    ax.set_xlabel('x [m]')
    ax.set_ylabel('f [Hz]')
    plt.colorbar(im, ax=ax, orientation=orientation, aspect=aspect)
    return im


class FDEM2dFOPold(pg.ModellingBase):
    """ old variant of 2D FOP (to be deleted) """
    def __init__(self, data, nlay=2, verbose=False):
        """ constructor with data and (optionally) number of layers """
        pg.ModellingBase.__init__(self, verbose)
        self.nlay = nlay
        self.FOP1d = data.FOP(nlay)
        self.nx = len(data.x)
        self.nf = len(data.freq())
        self.mesh_ = pg.createMesh1D(self.nx, 2 * nlay - 1)
        self.setMesh(self.mesh_)

    def response(self, model):
        """ yields forward model response """
        modA = np.asarray(model).reshape((self.nlay * 2 - 1, self.nx)).T
        resp = pg.RVector(0)
        for modi in modA:
            resp = pg.cat(resp, self.FOP1d.response(modi))

        return resp


class FDEM2dFOP(pg.ModellingBase):
    """ FDEM 2d-LCI modelling class based on BlockMatrices """
    def __init__(self, data, nlay=2, verbose=False):
        """ Parameters: FDEM data class and number of layers """
        super(FDEM2dFOP, self).__init__(verbose)
        self.nlay = nlay
        self.FOP = data.FOP(nlay)
        self.nx = len(data.x)
        self.nf = len(data.freq())
        npar = 2 * nlay - 1
        self.mesh1d = pg.createMesh1D(self.nx, npar)
        self.mesh_ = pg.createMesh1D(self.nx, 2 * nlay - 1)
        self.setMesh(self.mesh_)

        # self.J = NDMatrix(self.nx, self.nf*2, npar)
        self.J = pg.RBlockMatrix()
        self.FOP1d = []
        for i in range(self.nx):
            self.FOP1d.append(pg.FDEM1dModelling(
                nlay, data.freq(), data.coilSpacing, -data.height))
            n = self.J.addMatrix(self.FOP1d[-1].jacobian())
            self.J.addMatrixEntry(n, self.nf * 2 * i, npar * i)

        self.J.recalcMatrixSize()
        print(self.J.rows(), self.J.cols())

    def response(self, model):
        """ cut-together forward responses of all soundings """
        modA = np.asarray(model).reshape((self.nlay * 2 - 1, self.nx)).T
        resp = pg.RVector(0)
        for modi in modA:
            resp = pg.cat(resp, self.FOP.response(modi))

        return resp

    def createJacobian(self, model):
        modA = np.asarray(model).reshape((self.nlay * 2 - 1, self.nx)).T
        for i in range(self.nx):
            self.FOP1d[i].createJacobian(modA[i])


class HEM1dWithElevation(pg.ModellingBase):

    """
    Airborne FDEM modelling including variable bird height
    """

    def __init__(self, frequencies, coilspacing, nlay=2, verbose=False):
        """ Set up class by frequencies and geometries """
        pg.ModellingBase.__init__(self, verbose)
        self.nlay_ = nlay  # real layers (actually one more!)
        self.FOP_ = pg.FDEM1dModelling(nlay + 1, frequencies, coilspacing, 0.0)
        self.mesh_ = pg.createMesh1D(nlay, 2)  # thicknesses and resistivities
        self.mesh_.cell(0).setMarker(2)
        self.setMesh(self.mesh_)

    def response(self, model):
        thk = model(0, self.nlay)  # all thicknesses including bird height
        res = model(self.nlay - 1, self.nlay * 2)
        res[0] = 10000.
        return self.FOP_.response(pg.cat(thk, res))


class FDEM():
    """ Class for managing Frequency Domain EM data and their inversions """
    def __init__(self, x=None, freqs=None,
                 coilSpacing=None, inphase=None, outphase=None,
                 filename=None, scaleFreeAir=False):
        """
            Initialize data class and load data. Provide filename or data.
            If filename is given, data is loaded, overwriting settings.

            Parameters
            ----------
            x: array
                Array of measurement positions

            freq: array
                Measured frequencies

            coilSpacing : float
                Distance between 2 two coils

            inphase : array
                real part of |amplitude| * \exp^{i phase}

            outphase : array
                imaginary part of |amplitude| * \exp^{i phase}

            filename : str
                Filename to read from. Supported: .xyz (MaxMin), *.txt (Emsys)

            scaleFreeAir : bool
                Scale inphase and outphase data by free air (primary) solution

        """
        if isinstance(x, str) and freqs is None:  # string/file init
            filename = x

        self.x = x
        self.frequencies = freqs
        self.coilSpacing = coilSpacing

        self.IP = inphase
        self.OP = outphase
        self.ERR = None

        self.height = 1.0

        if filename:
            # check if filename extension is TXT or CSV
            fl = filename.lower()
            if fl.rfind('.txt') > 0 or fl.rfind('.csv') > 0:
                self.importEmsysAsciiData(filename)
            else:
                self.importMaxminData(filename)
        if np.any(self.frequencies):
            self.isActiveFreq = self.frequencies > 0.0
            self.activeFreq = np.nonzero(self.isActiveFreq)[0]

        if scaleFreeAir:
            freeAirSolution = self.FOP().freeAirSolution()
            self.IP /= freeAirSolution
            self.OP /= freeAirSolution

    def __repr__(self):
        if self.x is None:
            part1 = "<FDEMdata: sounding with {:d} frequencies".format(
                    len(self.frequencies))
        else:
            part1 = "<FDEMdata: {:d} soundings with {:d} frequencies".format(
                    len(self.x), len(self.frequencies))
        if self.coilSpacing is not None:
            cs = self.coilSpacing
            if hasattr(cs, '__iter__'):
                if len(cs) > 1:
                    part2 = "coil spacing is {:f}-{:f} m>".format(min(cs),
                                                                  max(cs))
            else:
                part2 = "coil spacing is {:f} m".format(cs)
            return part1 + ' , ' + part2
        else:
            return part1

    def importEmsysAsciiData(self, filename, verbose=False):
        """Import data from emsys text export:
            yields: positions, data, frequencies, error and geometry
        """
        cols = (1, 4, 6, 8, 9, 12, 15, 16)
        xx, sep, f, pf, ip, op, hmod, q = np.loadtxt(filename, skiprows=1,
                                                     usecols=cols, unpack=True)
        err = q / pf * 100.  # percentage of primary field

        if len(np.unique(sep)) > 1:
            print("Warning! Several coil spacings present in file!")

        self.coilSpacing = np.median(sep)
        f = np.round_(f)
        self.frequencies, mf, nf = np.unique(f, True, True)
        x, mx, nx = np.unique(xx, True, True)
        self.IP = np.ones((len(x), len(f))) * np.nan
        self.OP = np.ones((len(x), len(f))) * np.nan
        self.ERR = np.ones((len(x), len(f))) * np.nan

        for i in range(len(f)):
            # print(i, nx[i], nf[i])
            self.IP[nx[i], nf[i]] = ip[i]
            self.OP[nx[i], nf[i]] = op[i]
            self.ERR[nx[i], nf[i]] = err[i]

    def readHEMData(self, filename, takeevery=1, choosevcp=True):
        """ read RESOLVE type airborne EM data from .XYZ file """
        self.header = {}
        keyword = ''
        with open(filename) as f:
            for i, line in enumerate(f):
                if line[0] == '/':
                    line = line[1:].strip('\n').replace(',', '').replace('AND',
                                                                         '')
                    try:
                        result = [float(co) for co in line.split()]
                    except:
                        result = line.split()
                    if len(result) == 1:
                        result = result[0]
                    if keyword:
                        if isinstance(keyword, list):
                            for kw, res in zip(keyword, result):
                                self.header[kw] = res
                        else:
                            self.header[keyword] = result
                        keyword = ''
                    else:
                        keyword = result
                else:
                    break
            line = f.readline()
            print(line)
#            tmp = np.genfromtxt(fname=f, autostrip=True, comments='/',
#                skip_header=0, dtype=float, names=1, case_sensitive='lower',
#                missing_values='*', filling_values=-9999, skip_footer=1)
        tmp = np.genfromtxt(
            fname=filename, autostrip=True, comments='/',
            skip_header=i+1, dtype=float, names=True, case_sensitive='lower',
            missing_values='*', filling_values=-9999, skip_footer=1)
        # read properties from header
        if choosevcp:
            ivcp = np.nonzero(np.array(self.header['COILGEOMETRY']) == 1)[0]
        else:
            ivcp = range(len(self.header['FREQUENCY']))
        self.frequencies = np.array(self.header['FREQUENCY'])[ivcp]
        self.coilSpacing = np.array(self.header['COILSEPERATION'])[ivcp]
        # read properties from data block
        names = tmp.dtype.names
        pyproj = opt_import('pyproj', 'coordinate transformations EM lon,lat values')
        if 'lon' in names and 'lat' in names and pyproj is not None:
            utm = pyproj.Proj(proj='utm', zone=32, ellps='WGS84')  # projection
            x, y = utm(tmp['lon'], tmp['lat'])
        else:
            x, y = tmp['x'], tmp['y']
        self.pos = np.column_stack((x, y))[::takeevery]
        dx = np.sqrt(np.diff(self.pos[:, 0])**2 + np.diff(self.pos[:, 1])**2)
        self.x = np.hstack((0., np.cumsum(dx)))
        self.z = tmp['h_laser'][::takeevery]
        self.topo = tmp['topo'][::takeevery]
        IP = np.column_stack([tmp['real_'+str(i+1)] for i in ivcp])
        OP = np.column_stack([tmp['quad_'+str(i+1)] for i in ivcp])
        # better do a decimation or running average here
        self.IP = IP[::takeevery, :]
        self.OP = OP[::takeevery, :]
        self.isActiveFreq = self.frequencies > 0.0
        self.activeFreq = np.nonzero(self.isActiveFreq)[0]

    def importMaxminData(self, filename, verbose=False):
        """ Import MaxMin IPX format with pos, data, frequencies & geometry """
        delim = None
        fid = open(filename)

        for i, aline in enumerate(fid):
            if aline.split()[0][0].isdigit():  # number found
                break
            elif aline.find('COIL') > 0:  # [:6] == '/ COIL':
                self.coilSpacing = float(aline.replace(':', ': ').split()[-2])
#                pass
            elif aline.find('FREQ') > 0:  # [:6] == '/ FREQ':
                mya = aline[aline.find(':') + 1:].replace(',', ' ').split()
                myf = [float(aa) for aa in mya if aa[0].isdigit()]
                self.frequencies = np.array(myf)

        fid.close()

        if verbose:
            print("CS=", self.coilSpacing, "F=", self.frequencies)
        if aline.find(',') > 0:
            delim = ','

        nf = len(self.frequencies)
        if verbose:
            print("delim=", delim, "nf=", nf)
        A = np.loadtxt(filename, skiprows=i, delimiter=delim, comments='/').T
        x, y, self.IP, self.OP = A[0], A[1], A[
            2:nf * 2 + 2:2].T, A[3:nf * 2 + 2:2].T
        if max(x) == min(x):
            self.x = y
        else:
            self.x = x

    def deactivate(self, fr):
        """ Deactivate a single frequency """
        fi = np.nonzero(np.absolute(self.frequencies / fr - 1.) < 0.1)
        self.isActiveFreq[fi] = False
        self.activeFreq = np.nonzero(self.isActiveFreq)[0]

    def freq(self):
        """ Return active (i.e., non-deactivated) frequencies """
        return self.frequencies[self.activeFreq]

    def FOP(self, nlay=2):
        """
            Retrieve forward modelling operator using a block discretization

            Parameters
            ----------
            nlay : int
                Number of blocks

        """
        return pg.FDEM1dModelling(nlay, self.freq(), self.coilSpacing,
                                  -self.height)

    def FOPsmooth(self, zvec):
        """
            Retrieve forward modelling operator using fixed layers
            (smooth inversion)

            Parameters
            ----------
            zvec : array
        """
        return pg.FDEM1dRhoModelling(zvec, self.freq(), self.coilSpacing,
                                     -self.height)

    def selectData(self, xpos=0):
        """
            Retrieve inphase, outphase and error(if exist) vector from index
            or near given position

            Return: array, array, array|None
        """

        # check for index
        if isinstance(xpos, int) and (xpos < len(self.x)) and (xpos >= 0):
            n = xpos
        else:
            n = np.argmin(np.absolute(self.x - xpos))

        if self.ERR is not None:
            return (self.IP[n, self.activeFreq], self.OP[n, self.activeFreq],
                    self.ERR[n, self.activeFreq])
        else:
            return (self.IP[n, self.activeFreq], self.OP[n, self.activeFreq],
                    None)

    def error(self, xpos=0):
        """
            Return error vector
        """
        ip, op, err = self.selectData(xpos)
        return err

    def datavec(self, xpos=0):
        """
            Extract data vector (stack in and out phase) for given pos/no
        """

        ip, op, err = self.selectData(xpos)
        return np.hstack((ip, op))

    def errorvec(self, xpos=0, minvalue=0.0):
        """
            Extract error vector for a give position or sounding number
        """
        ip, op, err = self.selectData(xpos)
        return np.tile(np.maximum(err * 0.7071, minvalue), 2)
#        return pg.asvector(np.tile(np.maximum(err * 0.7071, minvalue), 2))

    def invBlock(self, xpos=0, nlay=2, noise=1.0,
                 stmod=30., lam=100., lBound=0., uBound=0., verbose=False):
        """
            Yield gimli inversion instance for block inversion
            inv(xpos,nlay) where nlay can be a FOP or a number of layers

            Parameters
            ----------
            xpos : array

            nLay : int
                Number of layers of the model to be determined

            noise : float
                Absolute data err in percent

            stmod : float or pg.RVector
                Starting model

            lam : float
                Global regularization parameter lambda.

            lBound : float
                Lower boundary for the model

            uBound : float
                Upper boundary for the model. 0 means no upper booundary

            verbose : bool
                Be verbose
        """

        self.transThk = pg.RTransLog()
        self.transRes = pg.RTransLogLU(lBound, uBound)
        self.transData = pg.RTrans()

        # EM forward operator
        if isinstance(nlay, pg.FDEM1dModelling):
            self.fop = nlay
        else:
            self.fop = self.FOP(nlay)

        data = self.datavec(xpos)

        self.fop.region(0).setTransModel(self.transThk)
        self.fop.region(1).setTransModel(self.transRes)

        if isinstance(noise, float):
            noiseVec = pg.RVector(len(data), noise)
        else:
            noiseVec = pg.asvector(noise)

        # independent EM inversion
        self.inv = pg.RInversion(data, self.fop, self.transData, verbose)
        if isinstance(stmod, float):  # real model given
            model = pg.RVector(nlay * 2 - 1, stmod)
            model[0] = 2.
        else:
            if len(stmod) == nlay * 2 - 1:
                model = pg.asvector(stmod)
            else:
                model = pg.RVector(nlay * 2 - 1, 30.)

        self.inv.setAbsoluteError(noiseVec)
        self.inv.setLambda(lam)
        self.inv.setMarquardtScheme(0.8)
        self.inv.setDeltaPhiAbortPercent(0.5)
        self.inv.setModel(model)
        self.inv.setReferenceModel(model)
        return self.inv

    def plotData(self, xpos=0, response=None, error=None, ax=None,
                 marker='bo-', rmarker='rx-', clf=True, addlabel='', nv=2):
        """
            Plot data as curves at given position
        """
        ip, op, err = self.selectData(xpos)

        if error is not None and err is not None:
            error = err

        fr = self.freq()

        if ax is None:
            fig, ax = plt.subplots(nrows=1, ncols=nv)
            ipax = ax[-2]
        else:
            ipax = ax[0]

        markersize = 4

        if error is not None:
            markersize = 2

        ipax.semilogy(ip, fr, marker, label='obs' + addlabel,
                      markersize=markersize)

        if error is not None and len(error) == len(ip):
            ipax.errorbar(ip, fr, xerr=error)

        # ipax.set_axis('tight')

        if error is not None:
            ipax.ylim((min(fr) * .98, max(fr) * 1.02))

        ipax.grid(True)
        ipax.set_xlabel('inphase [ppm]')
        ipax.set_ylabel('f [Hz]')

        if response is not None:
            rip = np.asarray(response)[:len(ip)]
            ipax.semilogy(rip, fr, rmarker, label='syn' + addlabel)

        ipax.legend(loc='best')

        opax = None

        if ax is None:
            opax = plt.subplot(1, nv, nv)
        else:
            opax = ax[-1]

        opax.semilogy(op, fr, marker, label='obs' + addlabel,
                      markersize=markersize)

        if error is not None and len(error) == len(ip):
            opax.errorbar(op, fr, xerr=error)

        if response is not None:
            rop = np.asarray(response)[len(ip):]
            opax.semilogy(rop, fr, rmarker, label='syn' + addlabel)

#        opax.set_axis('tight')

        if error is not None:
            opax.ylim((min(fr) * .98, max(fr) * 1.02))

        opax.grid(True)
        opax.set_xlabel('outphase [ppm]')
        opax.set_ylabel('f [Hz]')
        opax.legend(loc='best')
        # plt.subplot(1, nv, 1)
        return ax

    def plotDataOld(self, xpos=0, response=None,
                    marker='bo-', rmarker='rx-', clf=True):
        """
            plot data as curves at given position
        """
        ip, op = self.selectData(xpos)
        fr = self.freq()

        if clf:
            plt.clf()

        plt.subplot(121)
        plt.semilogy(ip, fr, marker, label='obs')
        plt.axis('tight')
        plt.grid(True)
        plt.xlabel('inphase [%]')
        plt.ylabel('f [Hz]')

        if response is not None:
            rip = np.asarray(response)[:len(ip)]
            plt.semilogy(rip, fr, rmarker, label='syn')

        plt.legend(loc='best')

        plt.subplot(122)
        plt.semilogy(op, fr, marker, label='obs')

        if response is not None:
            rop = np.asarray(response)[len(ip):]
            plt.semilogy(rop, fr, rmarker, label='syn')

        plt.axis('tight')
        plt.grid(True)
        plt.xlabel('outphase [%]')
        plt.ylabel('f [Hz]')
        plt.legend(loc='best')
        plt.show()

        return

    def showModelAndData(self, model, xpos=0, response=None, figsize=(8, 6)):
        """
        """
        fig, ax = plt.subplots(1, 3, figsize=figsize)

        model = np.asarray(model)
        nlay = (len(model) + 1) / 2

        thk = model[:nlay - 1]
        res = model[nlay - 1: 2 * nlay - 1]

        drawModel1D(ax[0], thk, res, plotfunction='semilogx')

        self.plotData(xpos, response, ax=ax[1:3], clf=False)
        return fig, ax

    def plotAllData(self, allF=True, orientation='horizontal', aspect=1000,
                    outname=None, show=False, figsize=(11, 6), everyx=1):
        """
            Plot data along a profile as image plots for IP and OP
        """

        if self.x is None:
            raise Exception("No measurement position array x given")

        freq = self.freq()
        np = 2

        if self.ERR is not None:
            np = 3

        if everyx is None:
            everyx = len(self.x) / 50

        fig, ax = plt.subplots(ncols=1, nrows=np, figsize=figsize)
        xfplot(ax[0], self.IP[:, self.activeFreq], self.x, freq,
               orientation=orientation, aspect=aspect, everyx=everyx)
        ax[0].set_title('inphase percent')

        xfplot(ax[1], self.OP[:, self.activeFreq], self.x, freq,
               orientation=orientation, aspect=aspect, everyx=everyx)
        ax[1].set_title('outphase percent')

        if self.ERR is not None:
            xfplot(ax[2], self.ERR[:, self.activeFreq], self.x, freq,
                   orientation=orientation)
            ax[2].set_title('error percent')

        if outname is not None:
            plt.savefig(outname)

        if show:
            plt.show()

        return

    def plotModelAndData(self, model, xpos, response,
                         modelL=None, modelU=None):
        self.plotData(xpos, response, nv=3)
        show1dmodel(model, color='blue')
        if modelL is not None and modelU is not None:
            pass
            # draw1dmodelErr(model, modelL, modelU)  # !!!!

        return

    def FOP2d(self, nlay):
        """ 2d forward modelling operator """
        return FDEM2dFOP(self, nlay)

    def inv2D(self, nlay, lam=100., resL=1., resU=1000., thkL=1.,
              thkU=100., minErr=1.0):
        """
            2d LCI inversion class
        """

        if isinstance(nlay, int):
            modVec = pg.RVector(nlay * 2 - 1, 30.)
            cType = 0  # no reference model
        else:
            modVec = nlay
            cType = 10  # use this as referencemodel
            nlay = (len(modVec) + 1) / 2

        # init forward operator
        self.f2d = self.FOP2d(nlay)

        # transformations
        self.tD = pg.RTrans()
        self.tThk = pg.RTransLogLU(thkL, thkU)
        self.tRes = pg.RTransLogLU(resL, resU)

        for i in range(nlay - 1):
            self.f2d.region(i).setTransModel(self.tThk)

        for i in range(nlay - 1, nlay * 2 - 1):
            self.f2d.region(i).setTransModel(self.tRes)

        # set constraints
        self.f2d.region(0).setConstraintType(cType)
        self.f2d.region(1).setConstraintType(cType)

        # collect data vector
        datvec = pg.RVector(0)

        for i in range(len(self.x)):
            datvec = pg.cat(datvec, self.datavec(i))

        # collect error vector
        if self.ERR is None:
            error = 1.0
        else:
            error = []
            for i in range(len(self.x)):
                err = np.maximum(self.ERR[i][self.activeFreq] * 0.701, minErr)
                error.extend(err)

        # generate starting model by repetition
        model = pg.asvector(np.repeat(modVec, len(self.x)))
        INV = pg.RInversion(datvec, self.f2d, self.tD)
        INV.setAbsoluteError(error)
        INV.setLambda(lam)
        INV.setModel(model)
        INV.setReferenceModel(model)

        return INV

if __name__ == "__main__":
    from optparse import OptionParser

    parser = OptionParser("usage: %prog [options] fdem",
                          version="%prog: " + pg.__version__)
    parser.add_option("-v", "--verbose", dest="verbose", action="store_true",
                      help="be verbose", default=False)
    parser.add_option("-n", "--nLayers", dest="nlay", help="number of layers",
                      type="int", default="4")
    parser.add_option("-x", "--xPos", dest="xpos", help="position/no to use",
                      type="float", default="0")
    parser.add_option("-l", "--lambda", dest="lam", help="init regularization",
                      type="float", default="30")
    parser.add_option("-e", "--error", dest="err", help="error estimate",
                      type="float", default="1")

    (options, args) = parser.parse_args()

    if options.verbose:
        __verbose__ = True

    A = NDMatrix(13, 6, 5)

    if len(args) == 0:
        parser.print_help()
        print("Please add a data file name.")
        sys.exit(2)
    else:
        datafile = args[0]

    fdem = FDEM(datafile)
    print(fdem)

    nlay = options.nlay
    xpos = options.xpos
    name = datafile.lower().rstrip('.xyz')
    fdem = FDEM(datafile)
    print(fdem)
    fdem.deactivate(56320.)  # do not use highest frequency
    fdem.plotAllData(outname=name + '-alldata.pdf')
    INV = fdem.invBlock(xpos=xpos, lam=options.lam, nlay=options.nlay,
                        noise=options.err, verbose=False)
    model = np.asarray(INV.run())
    INV.echoStatus()
    print("thk = ", model[:nlay - 1])
    print("res = ", model[nlay - 1:])
    fig, ax = fdem.showModelAndData(model, xpos, INV.response())
    INV = fdem.inv2D(options.nlay)
    INV.run()
#    fig.savefig(name+str(xpos)+'-result.pdf', bbox_inches='tight')
    plt.show()
