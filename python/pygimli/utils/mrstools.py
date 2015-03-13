# -*- coding: utf-8 -*-
"""
mrstools - tools for magnetic resonance sounding (MRS)
"""

import pygimli as pg
import numpy as N
import pylab as P
from .base import draw1dmodel, rndig

class MRS1dBlockQTModelling(pg.ModellingBase):
    """
    MRS1dBlockQTModelling - pygimli modelling class for block-mono QT inversion
    f=MRS1dBlockQTModelling(lay, KR, KI, zvec, t, verbose = False )
    """
    def __init__(self, nlay, KR, KI, zvec, t, verbose=False ):
        """constructor."""
        mesh = pg.createMesh1DBlock(nlay, 2) # thk, wc, T2*
        pg.ModellingBase.__init__(self, mesh, verbose)
        self.KR_ = KR
        self.KI_ = KI
        self.zv_ = N.array(zvec)
        self.nl_ = nlay
        self.nq_ = len(KR)
        self.t_  = N.array(t)
        self.nt_ = len(t)

    def response(self, par):
        """compute response function."""
        nl = self.nl_
        thk = par(0, nl-1)
        wc  = par(nl-1, 2*nl-1)
        t2  = par(2*nl-1, 3*nl-1)
        zthk = N.cumsum(thk)
        zv = self.zv_
        lzv = len(zv)
        izvec = N.zeros(nl + 1, N.int32)
        rzvec = N.zeros(nl + 1)
        for i in range(nl-1):
            ii = (zv<zthk[i]).argmin()
            izvec[ i + 1 ] = ii
            if ii <= len(zv):
                rzvec[i+1] = (zthk[i] - zv[ii-1]) / (zv[ii] - zv[ii-1])
        izvec[-1] = lzv - 1
        A = N.zeros((self.nq_, self.nt_), dtype=complex)
        for i in range(nl):
            wcvec = N.zeros(lzv-1)
            wcvec[ izvec[i]:izvec[i+1] ] = wc[i]
            if izvec[i+1] < lzv:
                wcvec[ izvec[i+1] - 1 ] = wc[i] * rzvec[ i + 1 ]
            if izvec[i] > 0:
                wcvec[ izvec[i] - 1 ] = wc[i] * (1 - rzvec[i])
            amps = N.dot(self.KR_, wcvec) + N.dot(self.KI_, wcvec) * 1j

            for ii in range(len(A)):
                A += N.exp(-self.t_ / t2[i]) * amps[ii]

        return pg.asvector(N.abs(A).ravel())

def loadmrsproject(mydir):
    """
    load mrs project from given directory (zkernel.ve)

    (datacube.dat, KR/KI.bmat, zkernel.vec)
    """
    if mydir is None:
        mydir = '.'
    if mydir[-1] != '/':
        mydir = mydir + '/'
    # load files from directory
    zvec = pg.RVector(mydir + 'zkernel.vec')
    KR = pg.RMatrix(mydir + 'KR.bmat')
    KI = pg.RMatrix(mydir + 'KI.bmat')
    A = pg.RMatrix()
    pg.loadMatrixCol(A, mydir + 'datacube.dat')
    t = N.array(A[0])
    #data
    datvec = pg.RVector()
    for i in range(1, len(A)):
        datvec = pg.cat(datvec, A[i])
    #print len(A), len(t)
    return KR, KI, zvec, t, datvec

def qtblockmodelling(mydir, nlay,
                     startvec=None, lowerbound=None, upperbound=None):
    """Loads data from dir, creates forward operator

    Parameters
    ----------
    mydir :

    nlay :

    startvec :

    lowerbound :

    upperbound :

    Examples
    --------

    t,datvec,f=qtblockmodelling(mydir,nlay,startvec=(10,0.3,0.2),
                                lowerbound=(0.1,0,0.02),
                                upperbound(100.,0.45,,0.5))
    """
    KR, KI, zvec, t, datvec = loadmrsproject(mydir)
    if startvec is None:
        startvec = [10., 0.3, 0.2]
    if lowerbound is None:
        lowerbound = [0.1, 0.0, 0.02]
    if upperbound is None:
        upperbound = [100., 0.45, 0.5]
    # modelling operator
    f = MRS1dBlockQTModelling(nlay, KR, KI, zvec, t) #A[0])
    f.region(0).setStartValue(startvec[0])
    f.region(1).setStartValue(startvec[1])
    f.region(2).setStartValue(startvec[2])
    # Model transformations
    f.transTH = pg.RTransLogLU(lowerbound[0], upperbound[0])
    f.transWC = pg.RTransLogLU(lowerbound[1], upperbound[1])
    f.transT2 = pg.RTransLogLU(lowerbound[2], upperbound[2])
    f.region(0).setTransModel(f.transTH)
    f.region(1).setTransModel(f.transWC)
    f.region(2).setTransModel(f.transT2)
    return t, datvec, f

def showqtresultfit(thk, wc, t2, datvec, resp, t,
                    islog=True, clim=None, nu=3, nv=2):
    """ show mrs qt result and data fit
    showqtresultfit(thk,wc,t2,datvec,resp,t,islog=True,clim=None,nu=3,nv=2)
    """
    if clim is None:
        cma = max(datvec)
        cmi = min(datvec)
        if islog:
            cma = N.log10(cma)
            cmi = cma - 1.5
        clim = (cmi, cma)

    nt = len(t)
    nq = len(datvec) / nt
    si = (nq, nt)

#    P.clf()
#    P.subplot(nu,nv,1)

    fig = P.figure(1)
    ax1 = fig.add_subplot(nu, nv, 1)

    draw1dmodel(wc, thk, islog=False, xlab=r'$\theta$')
#    P.subplot(nu,nv,3)
    ax3 = fig.add_subplot(nu, nv, 3)
    draw1dmodel(t2, thk, xlab='T2* in ms')
    ax3.set_xticks([0.02, 0.05, 0.1, 0.2, 0.5])
    ax3.set_xticklabels(('0.02', '0.05', '0.1', '0.2', '0.5'))
#    P.subplot(nu,nv,2)
    ax2 = fig.add_subplot(nu, nv, 2)
    if islog:
        P.imshow(N.log10(N.array(datvec).reshape(si)),
                 interpolation='nearest', aspect='auto')
    else:
        P.imshow(N.array(datvec).reshape(si),
                 interpolation='nearest', aspect='auto')
    P.clim(clim)
#    P.subplot(nu,nv,4)
    ax4 = fig.add_subplot(nu, nv, 4)
    if islog:
        P.imshow(N.log10(resp.reshape(si)),
                 interpolation='nearest',aspect='auto')
    else:
        P.imshow(resp.reshape(si),
                 interpolation='nearest',aspect='auto')
    misfit = N.array(datvec - resp)
    P.clim(clim)
#    P.subplot(nu,nv,5)
    ax5 = fig.add_subplot(nu, nv, 5)
    P.hist(misfit, bins=30)
    P.axis('tight')
    P.grid(which='both')
    P.text(P.xlim()[0], N.mean(P.ylim()),
           ' std=%g nV' % rndig(N.std(misfit), 3))
#    P.subplot(nu,nv,6)
    ax6 = fig.add_subplot(nu, nv, 6)
    P.imshow(misfit.reshape(si), interpolation='nearest', aspect='auto')
    ax = [ ax1, ax2, ax3, ax4, ax5, ax6 ]
    return ax
