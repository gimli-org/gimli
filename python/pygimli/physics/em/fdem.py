#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
    Was macht das Ding
'''

import pygimli as g

from pygimli import FDEM1dModelling, RVector, asvector, RTrans, RTransLog, RTransLogLU, RInversion
from pygimli.viewer import show1dmodel

import matplotlib.pyplot as plt
import numpy as np

def importEmsysAsciiData(filename, verbose=False):
    """ pure import function reading in positions, data, frequencies, error and geometry """
    
    xx, sep, f, pf, ip, op, hmod, q = np.loadtxt(filename, skiprows=1, usecols=(1,4,6,8,9,12,15,16), unpack=True)

    err = q / pf  * 100. # percentage of primary field
    
    if len(np.unique(sep)) > 1:
        print("Warning! Several coil spacings present in file!")
    
    coilspacing = np.median(sep)
    
    f = np.round_(f)
    freq, mf, nf = np.unique(f, True, True)
    x, mx, nx = np.unique(xx, True, True)
    IP = np.ones((len(x), len(freq))) * np.nan
    OP = np.ones((len(x), len(freq))) * np.nan
    ERR = np.ones((len(x), len(freq))) * np.nan
    
    for i in range(len(f)):
        #print i, nx[i], nf[i]
        IP[nx[i],nf[i]] = ip[i]
        OP[nx[i],nf[i]] = op[i]
        ERR[nx[i],nf[i]] = err[i]
    
    return x, freq, coilspacing, IP, OP, ERR    

def importMaxminData(filename, verbose = False):
    """ 
        Pure import function reading in positions, data, frequencies and geometry 
    """
    
    delim = None
    fid = open(filename)
    coilspacing = 0.
    freq = []
    
    for i, aline in enumerate(fid):
        if aline.split()[0][0].isdigit(): #number found
            break
        elif aline.find('COIL') > 0:     #[:6] == '/ COIL':
            coilspacing = float(aline.replace(':',': ').split()[-2])
        elif aline.find('FREQ') > 0:   #[:6] == '/ FREQ':
            freq = np.array([float(aa) for aa in aline[aline.find(':')+1:].replace(',',' ').split() if aa[0].isdigit()])
    
    fid.close()
    
    if verbose: print("CS=", coilspacing, "F=", freq)
    if aline.find(',')>0: delim=','
    
    nf = len(freq)
    if verbose: print("delim=", delim, "nf=", nf)
    A = np.loadtxt(filename, skiprows=i, delimiter=delim).T
    x, y, IP, OP = A[0], A[1], A[2:nf*2+2:2].T, A[3:nf*2+2:2].T
    if max(x)==min(x):
        x = y

    return x, freq, coilspacing, IP, OP

def xfplot(ax, DATA, x, freq, everyx=5, orientation='vertical', aspect=40):
    """
        Plots a matrix according to x and frequencies
        
    """
    
    nt = list(range(0, len(x), everyx))
    plt.imshow(DATA.T, interpolation='nearest')
    plt.ylim(plt.ylim()[::-1])
    ax.set_xticks(nt)
    ax.set_xticklabels(["%g" % xi for xi in x[nt]])
    ax.set_yticks(list(range(0,len(freq)+1, 2)))
    ax.set_yticklabels(["%g" % freq[i] for i in range(0, len(freq), 2)])
    plt.colorbar(orientation=orientation, aspect=aspect)
    plt.xlabel('x [m]')
    plt.ylabel('f [Hz]')

class FDEM2dFOP(g.ModellingBase):
    """
    """
    def __init__(self, data, nlay=2, verbose=False):
        """
        """
        g.ModellingBase.__init__(self, verbose)
        self.nlay_ = nlay
        self.FOP_  = data.FOP(nlay)
        self.nx_   = len(data.x)
        self.nf    = len(data.freq())
        self.mesh_ = g.createMesh1D(self.nx_, 2*nlay-1)
        self.setMesh(self.mesh_)

    def response(self, model):
        """
        """
        modA = np.asarray(model).reshape((self.nlay_*2-1,self.nx_)).T
        resp = g.RVector(0)
        for modi in modA:
            resp = g.cat(resp, self.FOP_.response(modi))
        
        return resp
    
class FDEMData():
    """
        Managing fdem data
    """
    def __init__(self, x=None, freqs=None, 
                 coilSpacing=None, inphase=None, outphase=None,
                 filename=None, scaleFreeAir=False):
        """
            Initialize data class and load data. Either provide filename or data.
            If a filename is given we try to load the data and overwrite data settings
            
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
                Filename to read from. Supported so far: ..??
                
            scaleFreeAir : bool
                Scale inphase and outphase data by free air solution (primary field)
            
        """
        self.x = x
        self.frequencies = freqs
        self.coilSpacing = coilSpacing
        
        self.IP = inphase
        self.OP = outphase
        self.ERR = None
        
        self.height = 1.0
        
        if filename:
            # check if filename extension is TXT or CSV
            if filename.lower().rfind('.txt') > 0 or filename.lower().rfind('.csv') > 0:
                self.x, self.frequencies, self.coilSpacing, self.IP, self.OP, self.ERR = importEmsysAsciiData(filename)
            else:
                self.x, self.frequencies, self.coilSpacing, self.IP, self.OP = importMaxminData(filename)
            
        self.isActiveFreq = self.frequencies > 0.0
        self.activeFreq = np.nonzero(self.isActiveFreq)[0]
        
        if scaleFreeAir:
            freeAirSolution = self.FOP().freeAirSolution();
            self.IP /= freeAirSolution
            self.OP /= freeAirSolution
        

    def __repr__(self):
        if self.x is None:
            return "<FDEMdata: sounding with %d frequencies, " \
                   "coilspacing is %.1f>" % (len(self.frequencies),
                                             self.coilSpacing)
        else:
            return "<FDEMdata: %d soundings with each %d frequencies, " \
                   "coilspacing is %.1f>" % (len(self.x),
                                             len(self.frequencies),
                                             self.coilSpacing)
        
    def showInfos(self): 
        """
            Only for old scripts using it
        """
        print(__repr__(self))

    def deactivate(self, fr):
        """
            Deactivate a single frequency
        """
        fi = np.nonzero(np.absolute(self.frequencies / fr - 1.) < 0.1)
        
        self.isActiveFreq[fi] = False
        self.activeFreq = np.nonzero(self.isActiveFreq)[0]
        
    def freq(self):
        """
            Return active (i.e., nonzero frequencies)
        """
        return self.frequencies[self.activeFreq]
    
    def FOP(self, nlay=2):
        """
            Retrieve forward modelling operator using a block discretization
            
            Parameters
            ----------
            nlay : int
                Number of blocks
                
        """
        return FDEM1dModelling(nlay, self.freq(), self.coilSpacing, -self.height)
    
    def FOPsmooth(self, zvec):
        """
            Retrieve forward modelling operator using fixed layers 
            (smooth inversion)
            
            Parameters
            ----------
            zvec : array
                ???
        """
        return FDEM1dRhoModelling(asvector(zvec), asvector(self.freq()), self.coilSpacing, -self.height)
    
    def selectData(self, xpos=0):
        """ 
            Retrieve inphase, outphase and error(if exist) vector from index 
            or near given position 
            
            Return: array, array, array|None
        """
        
        if isinstance(xpos, int) and (xpos < len(self.x)) and (xpos >= 0): # index
            n = xpos
        else:
            n = np.argmin(np.absolute(self.x - xpos))
        
        if self.ERR is not None:
            return self.IP[n, self.activeFreq], self.OP[n, self.activeFreq], self.ERR[n, self.activeFreq]
        else:
            return self.IP[n, self.activeFreq], self.OP[n, self.activeFreq], None

    def error(self, xpos=0):
        """ 
            Return error vector 
        """
        ip, op, err = selectData(xpos)
        return err

    def datavec(self, xpos=0):
        """ 
            Extract data vector (stacking inphase and outphase
        """
            
        ip, op, err = self.selectData(xpos)
        return asvector(np.hstack((ip, op)))
    
    def errorvec(self, xpos=0, minvalue=0.0):
        """
            Extract error vector 
        """
        ip, op, err = self.selectData(xpos)
        return asvector(np.tile(np.maximum(err * 0.7071, minvalue), 2))
    
    def invBlock(self, xpos=0, nlay=2, noise=1.0,
                 stmod=30., lam=100., lBound=1., uBound=0., verbose=False):
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
            
            stmod : float
                Constant Starting model
                
            lam : float
                Global regularization parameter lambda.
                
            lBound : float
                Lower boundary for the model
                
            uBound : float
                Upper boundary for the model. 0 means no upper booundary
            
            verbose : bool
                Be verbose
        """
            
        self.transThk = RTransLog()
        self.transRes = RTransLogLU(lBound, uBound)
        self.transData = RTrans()
        
        # EM forward operator
        if isinstance(nlay, FDEM1dModelling):
            self.fop = nlay
        else:
            self.fop = self.FOP(nlay)
        
        data = self.datavec(xpos)
        
        self.fop.region(0).setTransModel(self.transThk)
        self.fop.region(1).setTransModel(self.transRes)
        
        if isinstance(noise, float):
            noiseVec = RVector(len(data), noise)
        else:
            noiseVec = asvector(noise)
        
        # independent EM inversion
        self.inv = RInversion(data, self.fop, self.transData, verbose)
        if isinstance(stmod, float): # real model given
            model = RVector(nlay * 2 - 1, stmod)
            model[0] = 2.
        else:
            if len(stmod) == nlay*2-1:
                model = asvector(stmod)
            else:
                model = RVector(nlay*2-1, 30.)
        
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
        
        ipax = None
        
        if ax is None:
            if clf: plt.clf()
            plt.subplot(1,nv,nv-1)
        else:
            ipax = ax[0]
        
        markersize = 4

        if error is not None:
            markersize = 2
        
        ipax.semilogy(ip, fr, marker, label='obs'+addlabel, markersize=markersize)

        if error is not None and len(error) == len(ip):
            ipax.errorbar(ip, fr, xerr=error)

        #ipax.set_axis('tight')

        if error is not None:
            ipax.ylim((min(fr)*.98,max(fr)*1.02))

        ipax.grid(True)
        ipax.set_xlabel('inphase [%]')
        ipax.set_ylabel('f [Hz]')
        
        if response is not None:
            rip = np.asarray(response)[:len(ip)]
            ipax.semilogy(rip, fr, rmarker, label='syn' + addlabel)
        
        ipax.legend(loc='best')
        
        opax = None
        
        if ax is None:
            opax = plt.subplot(1, nv, nv)
        else:
            opax = ax[1]
        
        opax.semilogy(op, fr, marker, label='obs'+addlabel,
                     markersize=markersize)
        
        if error is not None and len(error) == len(ip):
            opax.errorbar(op, fr, xerr=error)

        if response is not None:
            rop = np.asarray(response)[len(ip):]
            opax.semilogy(rop, fr, rmarker, label='syn'+addlabel)
        
        #opax.set_axis('tight')
        
        if error is not None:
            opax.ylim((min(fr) * .98, max(fr) * 1.02))
        
        opax.grid(True)
        opax.set_xlabel('outphase [%]')
        opax.set_ylabel('f [Hz]')
        opax.legend(loc='best')
        #plt.subplot(1, nv, 1)
        return 

        
    def plotDataOld(self, xpos=0, response = None,
                    marker='bo-', rmarker='rx-', clf=True):
        """
            plot data as curves at given position
        """
        ip, op = self.selectData(xpos)
        fr = self.freq()
        
        if clf: plt.clf()
        
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
    
    def showModelAndData(self, model, xpos=0, response=None):
        """
        """
        plt.clf()
        
        model = np.asarray(model)
        nlay = (len(model) + 1) / 2
        
        thk = model[:nlay-1]
        res = model[nlay-1: 2*nlay-1]
        
        ax1 = plt.subplot(131)
        show1dmodel(res, thk)
        
        ax2 = plt.subplot(132)
        ax3 = plt.subplot(133)
        
        self.plotData(xpos, response, ax=(ax2, ax3), clf=False)
         
    def plotAllData(self, allF=True, orientation='vertical', outname=None):
        """
            Plot data along a profile as image plots for IP and OP
        """
        
        if self.x is None:
            raise Exception("No measurement position array x given")
        
        freq = self.freq()
        nf = len(freq)
        np = 2
        
        if self.ERR is not None:
            np = 3
        
        plt.clf()
        ax1 = plt.subplot(np, 1, 1)
        xfplot(ax1, self.IP[:,self.activeFreq], self.x, freq,
               orientation=orientation)
        ax1.set_title('inphase percent')
        
        ax2 = plt.subplot(np,1,2)
        xfplot(ax2, self.OP[:,self.activeFreq], self.x, freq,
               orientation=orientation)
        ax2.set_title('outphase percent')
        
        if self.ERR is not None:
            ax3 = plt.subplot(np,1,3)
            xfplot(ax3, self.ERR[:,self.activeFreq], self.x, freq,
                   orientation=orientation)
            ax3.set_title('error percent')

        if outname is not None:
            plt.savefig(outname)
        
        plt.show()
        return

    def plotModelAndData(self, model, xpos, response, modelL=None, modelU=None):
        plotData(xpos, response, nv=3)
        show1dmodel(model, color='blue')
        if modelL is not None and modelU is not None:
            draw1dmodelErr(model, modelL, modelU)

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
            modVec = g.RVector(nlay*2 -1, 30.)
            cType = 0 # no reference model
        else:
            modVec = nlay
            cType = 10 # use this as referencemodel
            nlay = (len(modVec) + 1) / 2


        # init forward operator
        self.f2d = self.FOP2d(nlay)

        # transformations
        self.tD = g.RTrans()
        self.tThk = g.RTransLogLU(thkL, thkU)
        self.tRes = g.RTransLogLU(resL, resU)
        
        for i in range(nlay-1): self.f2d.region(i).setTransModel(self.tThk)
        
        for i in range(nlay-1, nlay*2-1): self.f2d.region(i).setTransModel(self.tRes)
        
        # set constraints
        self.f2d.region(0).setConstraintType(10)
        self.f2d.region(1).setConstraintType(10)

        # collect data vector
        datvec = g.RVector(0)
        
        for i in range(len(self.x)):
            datvec = g.cat(datvec, self.datavec(i))

        # collect error vector
        errVec = []
        
        for i in range(len(self.x)):
            errVec.extend(np.maximum(self.ERR[i][self.activeFreq] * 0.701, minErr))
            errVec.extend(np.maximum(self.ERR[i][self.activeFreq] * 0.701, minErr))

        # generate starting model by repetition
        model = g.asvector(np.repeat(modVec, len(self.x)))
        INV = g.RInversion(datvec, self.f2d, self.tD) 
        INV.setAbsoluteError(g.asvector(errVec))
        INV.setLambda(lam)
        INV.setModel(model)
        INV.setReferenceModel(model)
        
        return INV

if __name__ == "__main__":
    import sys
    from optparse import OptionParser

    parser = OptionParser("usage: %prog [options] fdem", version="%prog: " + g.versionStr() )
    parser.add_option("-v", "--verbose", dest="verbose", action="store_true"
                            , help="be verbose", default=False)
    
    (options, args) = parser.parse_args()

    if options.verbose:
        __verbose__ = True 
        
    if len(args) == 0:
        parser.print_help()
        print("Please add a mesh or model name.")
        sys.exit(2)
    else:
        if len(args) == 0:
            parser.print_help()
            print("Please add a mesh or model name.")
            sys.exit(2)
        else:
            fdem = FDEMData(args[0])


