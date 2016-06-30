#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""For the manager look at BERT https://gitlab.com/resistivity-net/bert
"""

import numpy as np
import matplotlib.pyplot as plt

import pygimli as pg


def slownessWillie(phi, sat=1, vm=4000, vw=1600, va=330,
                   mesh=None, meshI=None, fill=None):
    r"""
    Return slowness :math:`s` after Wyllie time-average equation

    .. math::
        s = (1-\phi) \cdot\frac{1}{v_m} + \phi \cdot S \cdot\frac{1}{v_w} +
        \phi \cdot(1 - S) \cdot\frac{1}{v_a}

    * :math:`\phi` - porosity 0.0 --1.0
    * :math:`S`    - fluid saturation 0.0 --1.0 [sat]
    * :math:`v_m`  - velocity of matrix
    * :math:`v_w`  - velocity of water
    * :math:`v_a`  - velocity of air

    If mesh is not None the resulting values are calculated for each cell of
    the mesh.
    All parameter can be scalar, array of length mesh.cellCount()
    or callable(pg.cell).
    If meshI is not None the result is interpolated to meshI.cellCenters()
    and prolonged (if fill ==1).

    Examples
    --------

    WRITEME


    """
    # look at resistivity.py
    raise BaseException('TODO')


# Wyllie (time-average) equation
def wyllie(phi, sat=1, vm=4000, vw=1600, va=330):
    """Return slowness after Wyllie time-average equation"""
    return 1./vm * (1-phi) + phi * sat * 1./vw + phi * (1 - sat) * 1./va


def transFwdWylliePhi(vm=4000, vw=1600, va=330):
    """Wyllie transformation function porosity(slowness) """
    return pg.RTransLin(1./vw-1./vm, 1./vm)


def transInvWylliePhi(vm=4000, vw=1600, va=330):
    """Inverse Wyllie transformation function porosity(slowness) """
    a1, b1 = 1./vm, 1./vw-1./vm
    return pg.RTransLin(1./b1, -a1/b1)


def transFwdWyllieS(phi, vm=4000, vw=1600, va=330):
    """Wyllie transformation function slowness(saturation) """
    return pg.RTransLin((1/vw-1./va)*phi, (1-phi)/vm+phi/va)


def transInvWyllieS(phi, vm=4000, vw=1600, va=330):
    """Inverse Wyllie transformation function slowness(saturation) """
    a2, b2 = 1./vm * (1 - phi) + phi * 1./va, phi * (1./vw - 1./va)
    return pg.RTransLin(1./b2, -a2/b2)


def test_Wyllie():
    """Test Wyllie and Transformations"""
    import unittest

    phivec = np.arange(0, 0.5, 0.01)
    swvec = np.arange(0, 1, 0.01)

    phi0 = 0.4

    tIWPhi = transInvWylliePhi()
    tIWS = transInvWyllieS(phi0)
    tFWPhi = transFwdWylliePhi()
    tFWS = transFwdWyllieS(phi0)

    ax = plt.subplots()[1]
    # direct function
    rA = wyllie(phivec)
    rS = wyllie(phi0, swvec)

    ax.semilogy(phivec, rA, 'b-')
    ax.semilogy(swvec, rS, 'r-')
    # forward transformation
    fA = tFWPhi.trans(phivec)
    fS = tFWS.trans(swvec)
    np.testing.assert_allclose(rA, fA, rtol=1e-12)
    np.testing.assert_allclose(rS, fS, rtol=1e-12)

    ax.semilogy(phivec, fA, 'bx', markersize=10)
    ax.semilogy(swvec, fS, 'rx', markersize=10)

    # inverse transformation
    iA = tIWPhi.invTrans(phivec)
    iS = tIWS.invTrans(swvec)
    np.testing.assert_allclose(rA, iA, rtol=1e-12)
    np.testing.assert_allclose(rS, iS, rtol=1e-12)

    ax.semilogy(phivec, iA, 'bo')
    ax.semilogy(swvec, iS, 'ro')

    plt.show()


if __name__ == "__main__":
    pass
