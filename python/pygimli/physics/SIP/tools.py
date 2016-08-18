#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    Spectral induced polarisation (SIP) tools
"""

from math import pi
import numpy as np
import pygimli as pg
from .models import ColeColeComplex, ColeColeComplexSigma, PeltonPhiEM
from .models import ColeColeAbs, ColeColePhi, DoubleColeColePhi


def fitCCEMPhi(f, phi,  ePhi=0.001, lam=1000., verbose=True,
               mpar=(0.2, 0, 1), taupar=(1e-2, 1e-5, 100),
               cpar=(0.25, 0, 1), empar=(1e-7, 1e-9, 1e-5)):
    """fit a Cole-Cole term with additional EM term to phase"""
    fCCEM = PeltonPhiEM(f)
    fCCEM.region(0).setParameters(*mpar)    # m (start,lower,upper)
    fCCEM.region(1).setParameters(*taupar)  # tau
    fCCEM.region(2).setParameters(*cpar)   # c
    fCCEM.region(3).setParameters(*empar)   # tau-EM
    ICC = pg.RInversion(phi, fCCEM, False)  # set up inversion class
    ICC.setAbsoluteError(ePhi)  # 1 mrad
    ICC.setLambda(lam)  # start with large damping and cool later
    ICC.setMarquardtScheme(0.8)  # lower lambda by 20%/it., no stop chi=1
    model = ICC.run()  # run inversion
    if verbose:
        ICC.echoStatus()
    return model, np.asarray(ICC.response())


def fitCCPhi(f, phi,  ePhi=0.001, lam=1000., verbose=True, robust=False,
             mpar=(0.2, 0, 1), taupar=(1e-2, 1e-5, 100), cpar=(0.25, 0, 1)):
    """fit a Cole-Cole term with additional EM term to phase"""
    fCCEM = ColeColePhi(f)
    fCCEM.region(0).setParameters(*mpar)    # m (start,lower,upper)
    fCCEM.region(1).setParameters(*taupar)  # tau
    fCCEM.region(2).setParameters(*cpar)   # c
    ICC = pg.RInversion(phi, fCCEM, False)  # set up inversion class
    ICC.setAbsoluteError(ePhi)  # 1 mrad
    ICC.setLambda(lam)  # start with large damping and cool later
    ICC.setMarquardtScheme(0.8)  # lower lambda by 20%/it., no stop chi=1
    ICC.setRobustData(robust)
#    ICC.setMaxIter(0)
    model = ICC.run()  # run inversion
    if verbose:
        ICC.echoStatus()
    return model, np.asarray(ICC.response())
#    return np.append(model, ICC.chi2()), np.asarray(ICC.response())


def fit2CCPhi(f, phi,  ePhi=0.001, lam=1000., verbose=True, robust=False,
              mpar1=(0.2, 0, 1), taupar1=(1e-2, 1e-5, 100), cpar1=(0.2, 0, 1),
              mpar2=(0.2, 0, 1), taupar2=(1e-2, 1e-5, 100), cpar2=(0.2, 0, 1)):
    """fit a Cole-Cole term with additional EM term to phase"""
    f2CC = DoubleColeColePhi(f)
    f2CC.region(0).setParameters(*mpar1)    # m (start,lower,upper)
    f2CC.region(1).setParameters(*taupar1)  # tau
    f2CC.region(2).setParameters(*cpar1)   # c
    f2CC.region(3).setParameters(*mpar2)    # m (start,lower,upper)
    f2CC.region(4).setParameters(*taupar2)  # tau
    f2CC.region(5).setParameters(*cpar2)   # c
    ICC = pg.RInversion(phi, f2CC, False)  # set up inversion class
    ICC.setAbsoluteError(ePhi)  # 1 mrad
    ICC.setLambda(lam)  # start with large damping and cool later
    ICC.setMarquardtScheme(0.8)  # lower lambda by 20%/it., no stop chi=1
    ICC.setRobustData(robust)
    ICC.setDeltaPhiAbortPercent(1)
#    ICC.setMaxIter(0)
    model = ICC.run()  # run inversion
    if verbose:
        ICC.echoStatus()
    return np.asarray(model), np.asarray(ICC.response())


def fitCCAbs(f, amp, error=0.01, lam=1000., mstart=None,
             taupar=(1e-2, 1e-5, 100), cpar=(0.5, 0, 1)):
    """fit amplitude spectrum by Cole-Cole model"""
    fCC = ColeColeAbs(f)
    tLog = pg.RTransLog()
    fCC.region(0).setStartValue(max(amp))
    if mstart is None:  # compute from amplitude decay
        mstart = 1. - min(amp) / max(amp)
    fCC.region(1).setParameters(mstart, 0, 1)    # m (start,lower,upper)
    fCC.region(2).setParameters(*taupar)  # tau
    fCC.region(3).setParameters(*cpar)   # c
    ICC = pg.RInversion(amp, fCC, tLog, tLog, False)  # set up inversion class
    ICC.setRelativeError(error)  # perr + ePhi/data)
    ICC.setLambda(lam)  # start with large damping and cool later
    ICC.setMarquardtScheme(0.8)  # lower lambda by 20%/it., no stop chi=1
    model = np.asarray(ICC.run())  # run inversion
    ICC.echoStatus()
    response = np.asarray(ICC.response())
    return model, response


def fitCCC(f, amp, phi, eRho=0.01, ePhi=0.001, lam=1000., mstart=None,
           taupar=(1e-2, 1e-5, 100), cpar=(0.5, 0, 1)):
    """fit complex spectrum by Cole-Cole model"""
    fCC = ColeColeComplex(f)
    tLog = pg.RTransLog()
    fCC.region(0).setStartValue(max(amp))
    if mstart is None:  # compute from amplitude decay
        mstart = 1. - min(amp) / max(amp)
    fCC.region(1).setParameters(mstart, 0, 1)    # m (start,lower,upper)
    fCC.region(2).setParameters(*taupar)  # tau
    fCC.region(3).setParameters(*cpar)   # c
    data = pg.cat(amp, phi)
    ICC = pg.RInversion(data, fCC, False)  # set up inversion class
    ICC.setTransModel(tLog)
    error = pg.cat(eRho*amp, pg.RVector(len(f), ePhi))
    ICC.setAbsoluteError(error)  # perr + ePhi/data)
    ICC.setLambda(lam)  # start with large damping and cool later
    ICC.setMarquardtScheme(0.8)  # lower lambda by 20%/it., no stop chi=1
    model = np.asarray(ICC.run())  # run inversion
    ICC.echoStatus()
    response = np.asarray(ICC.response())
    return model, response[:len(f)], response[len(f):]


def fitCCCC(f, amp, phi, error=0.01, lam=10., taupar=(1e-2, 1e-5, 100),
            cpar=(0.25, 0, 1), mpar=[0, 0, 1]):
    """fit complex spectrum by Cole-Cole model based on sigma"""
    fCC = ColeColeComplexSigma(f)
    tLog = pg.RTransLog()
    fCC.region(0).setStartValue(1./max(amp))
    if mpar[0] == 0:
        mpar[0] = 1. - min(amp)/max(amp)
    fCC.region(1).setParameters(*mpar)    # m (start,lower,upper)
    fCC.region(2).setParameters(*taupar)  # tau
    fCC.region(3).setParameters(*cpar)   # c
    data = pg.cat(1./amp * np.cos(phi), 1./amp * np.sin(phi))
    ICC = pg.RInversion(data, fCC, False)  # set up inversion class
    ICC.setTransModel(tLog)
    ICC.setAbsoluteError(data*error+max(data)*0.0001)  # perr + ePhi/data)
    ICC.setLambda(lam)  # start with large damping and cool later
    ICC.setMarquardtScheme(0.8)  # lower lambda by 20%/it., no stop chi=1
    model = np.asarray(ICC.run())  # run inversion
    ICC.echoStatus()
    response = np.asarray(ICC.response())
    rRe, rIm = response[:len(f)], response[len(f):]
    rAmp = 1./np.sqrt(rRe**2+rIm**2)
    return model, rAmp, np.arctan(rIm/rRe)


def KramersKronig(f, re, im, usezero=False):
    """Return real/imaginary parts retrieved by Kramers-Kronig relations

       formulas including singularity removal according to Boukamp (1993)
    """
    from scipy.integrate import simps

    x = f * 2. * pi
    im2 = np.zeros(im.shape)
    re2 = np.zeros(im.shape)
    re3 = np.zeros(im.shape)
    drdx = np.diff(re) / np.diff(x)
    dredx = np.hstack((drdx[0], (drdx[:-1] + drdx[1:]) / 2, drdx[-1]))
    didx = np.diff(im) / np.diff(x)
    dimdx = np.hstack((didx[0], (didx[:-1] + didx[1:]) / 2, didx[-1]))
    for num in range(len(x)):
        w = x[num]
        x2w2 = x**2 - w**2
        x2w2[num] = 1e-12
        fun1 = (re - re[num]) / x2w2
        fun1[num] = dredx[num] / 2 / w
        im2[num] = -simps(fun1, x) * 2. * w / pi
        fun2 = (im * w / x - im[num]) / x2w2
        re2[num] = simps(fun2, x) * 2. * w / pi + re[0]
        fun3 = (im * x - im[num] * w) / x2w2
        fun3[num] = (im[num] / w + dimdx[num]) / 2
        re3[num] = simps(fun3, x) * 2. / pi + re[-1]

    if usezero:
        return re2, im2
    else:
        return re3, im2


if __name__ == "__main__":
    pass
