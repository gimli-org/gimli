#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""    
    For the manager look at BERT https://gitlab.com/resistivity-net/bert
    
"""

import pygimli as pg
import numpy as np

def slownessWillie(phi, S=1, vm=4000, vw=1600, va=330,
                  mesh=None, meshI=None, fill=None):
    r"""
    Return slowness :math:`s` after Wyllie time-average equation
    
    .. math:: 
        s = (1-\phi) \cdot\frac{1}{v_m} + \phi \cdot S \cdot\frac{1}{v_w} + \phi \cdot(1 - S) \cdot\frac{1}{v_a}
        
    * :math:`\phi` - porosity 0.0 --1.0
    * :math:`S`    - fluid saturation 0.0 --1.0
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
    #look at resistivity.py
    IMPLEMENTME 



# Wyllie (time-average) equation
def Wyllie(phi, Sw=1, vm=4000, vw=1600, va=330):
    """ return slowness after Wyllie time-average equation """
    return 1./vm * (1-phi) + phi * Sw * 1./vw + phi * (1 - Sw) * 1./va


def transInvWylliePhi(vm=4000, vw=1600, va=330):
    """ inverse Wyllie transformation function porosity(slowness) """
    a1, b1 = 1./vm, 1./vw-1./vm
    return pg.RTransLin(1./b1, -a1/b1)


def transInvWyllieS(phi, vm=4000, vw=1600, va=330):
    """ inverse Wyllie transformation function slowness(saturation) """
    a2, b2 = 1./vm * (1 - phi) + phi * 1./va, phi * (1./vw - 1./va)
    return pg.RTransLin(1./b2, -a2/b2)


def transFwdWylliePhi(vm=4000, vw=1600, va=330):
    """ inverse Wyllie transformation function porosity(slowness) """
    return pg.RTransLin(1./vw-1./vm, 1./vm)


def transFwdWyllieS(phi, vm=4000, vw=1600, va=330):
    """ inverse Wyllie transformation function slowness(saturation) """
    return pg.RTransLin((1/vw-1./va)*phi, (1-phi)/vm+phi/va)

def test_Wyllie():
    phivec = np.arange(0, 0.5, 0.01)
    swvec = np.arange(0, 1, 0.01)

    phi0 = 0.4
    phi1 = np.ones_like(phivec)
    s1 = np.ones_like(swvec)
    tIWPhi = transInvWylliePhi()
    tIWS = transInvWyllieS(phi0)
    tFWPhi = transFwdWylliePhi()
    tFWS = transFwdWyllieS(phi0)

    fig, ax = plt.subplots()
    # direct function
    ax.semilogy(phivec, Wyllie(phivec), 'b-')
    ax.semilogy(swvec, Wyllie(phi0, swvec), 'r-')
    # forward transformation
    ax.semilogy(phivec, tFWPhi.trans(phivec), 'bx')
    ax.semilogy(swvec, tFWS.trans(swvec), 'rx')
    # inverse transformation
    ax.semilogy(phivec, tIWPhi.invTrans(phivec), 'b+')
    ax.semilogy(swvec, tIWS.invTrans(swvec), 'r+')

if __name__ == "__main__":
    pass


