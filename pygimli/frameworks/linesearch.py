# -*- coding: utf-8 -*-
"""pyGIMLi - Linesearch.

Linesearch procedures used by various inversion frameworks.
"""
import numpy as np
import pygimli as pg


def tauVector(taumin=0.01, taumax=1, logScale=False, n=21):
    """Generate a vector of tau values."""
    if logScale:
        return np.logspace(np.log10(taumin), np.log10(taumax), n)
    else:
        return np.linspace(taumin, taumax, n)


def lineSearchExact(inv, dM, taus=None, show=False, **kwargs):
    """Line search by exact forward response.

    Parameters
    ----------
    inv : pg.Inversion
        pygimli Inversion (or any derived class) instance
    taus : array
        array containing the tau values to test, alternatively:
    taumin : float
        minimum value
    taumax : float
        maximum value
    logScale : bool
        use logarithmic scaling, otherwise linear
    show : bool
        show curve
    """
    if taus is None:
        taus = tauVector(**kwargs)

    phis = np.zeros_like(taus)
    for i, tau in enumerate(taus):
        newModel = inv.modelTrans.update(inv.model, dM*tau)
        newResponse = inv.fop.response(newModel)
        phis[i] = inv.inv.getPhi(newModel, newResponse)

    if show:
        pg.plt.semilogx(taus, phis)

    return taus[np.argmin(phis)], newResponse

def lineSearchInter(inv, dM, taus=None, show=False, **kwargs):
    """Optimizes line search parameter by linear response interpolation.

    Parameters
    ----------
    inv : pg.Inversion
        pygimli Inversion (or any derived class) instance
    taus : array
        array containing the tau values to test, alternatively:
    taumin : float
        minimum value
    taumax : float
        maximum value
    logScale : bool
        use logarithmic scaling, otherwise linear
    show : bool
        show curve
    """
    if taus is None:
        taus = tauVector(**kwargs)

    phis = np.zeros_like(taus)
    dT = inv.dataTrans
    oldResponse = dT(inv.response)
    fullModel = inv.modelTrans.update(inv.model, dM)
    fullResponse = dT(inv.fop.response(fullModel))
    for i, tau in enumerate(taus):
        newModel = inv.modelTrans.update(inv.model, dM*tau)
        newResponse = dT.inv(oldResponse + (fullResponse - oldResponse) * tau)
        phis[i] = inv.inv.getPhi(newModel, newResponse)

    if show:
        pg.plt.semilogx(taus, phis)

    return taus[np.argmin(phis)], newResponse


def lineSearchInterOld(inv, dM, nTau=100, maxTau=1.0):
    """Optimizes line search parameter by linear response interpolation."""
    tD = inv.dataTrans
    tM = inv.modelTrans
    model = inv.model
    response = inv.response
    modelLS = tM.update(model, dM * maxTau)
    responseLS = inv.fop.response(modelLS)
    taus = np.arange(1, nTau+1) / nTau * maxTau
    phi = np.ones_like(taus) * inv.phi()
    phi[-1] = inv.phi(modelLS, responseLS)
    t0 = tD.fwd(response)
    t1 = tD.fwd(responseLS)
    for i, tau in enumerate(taus):
        modelI = tM.update(model, dM*tau)
        responseI = tD.inv(t1*tau+t0*(1.0-tau))
        phi[i] = inv.phi(modelI, responseI)

    return taus[np.argmin(phi)], responseLS

def lineSearchQuad(inv, dM, tau05=0.3, tau1=1):
    """Optimize line search by fitting parabola by Phi(tau) curve."""
    return 0.1, None

def lineSearch(inv, dm, method=None, **kwargs):
    """Carry out line search."""
    if method.lower().startswith("exact"):
        return lineSearchExact(inv, dm, **kwargs)
    elif method.lower().startswith("int"):
        return lineSearchInter(inv, dm, **kwargs)
    elif method.lower().startswith("quad"):
        return lineSearchQuad(inv, dm, **kwargs)
    else:
        tau, response = lineSearchInter(inv, dm, **kwargs)
        if tau > 0.01 and tau <= 1:
            return tau, response

        tau, response = lineSearchQuad(inv, dm, **kwargs)
        if tau > 0.01 and tau <= 1:
            return tau, response
        else:
            return 0.1, None

if __name__ == "__main__":
    pass