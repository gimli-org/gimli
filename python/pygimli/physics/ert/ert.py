#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    Currently only petrophysical models.
    
    For the manager look at BERT https://gitlab.com/resistivity-net/bert
    
"""

import pygimli as pg
import numpy as np

def resistivityArchie(rFluid, porosity, a=1.0, m=2.0, S=1.0, n=2.0,
                      mesh=None, meshI=None, fill=None):
    """
    Return resistivity of rock for the petrophysical model Archies law.
    
    .. math:: 
        \rho = a\rho_{\text{fl}}\phi^{-m}\S_w^{-n}
        
    * :math:`\rho` - the electrical conductivity of the fluid saturated rock
    * :math:`\rho_{\text{fl}}` - electrical conductivity of the fluid
    * :math:`\phi` - porosity 0.0 --1.0
    * :math:`a` - Tortuosity factor. (common 1)
    * :math:`m` - Cementation exponent of the rock
            (usually in the range 1.3 -- 2.5 for sandstones)
    * :math:`n` - is the saturation exponent (usually close to 2)
     
    If mesh is not None the resulting values are calculated for each cell of 
    the mesh. 
    All parameter can be scalar, array of length mesh.cellCount() 
    or callable(pg.cell). If rFluid is non-steady n-step distribution 
    than rFluid can be a matrix of size(n, mesh.cellCount())
    If meshI is not None the result is interpolated to meshI.cellCenters() 
    and prolonged (if fill ==1).
    
    Note
    ----
    We experience some unstable nonlinear behavior.
    Until this is clarified all results are rounded to the precision 1e-6.
     
    Examples
    --------
    
    WRITEME
     
     
    """
    
    if mesh == None:
        return rFluid * a * porosity**(-m) * S**(-n)
    
    rB = None
    
    if rFluid.ndim == 1:
        rB = pg.RMatrix(1, len(rFluid))
        rB[0] = pg.solver.parseArgToArray(rFluid, mesh.cellCount(), mesh)
    elif rFluid.ndim == 2:
        rB = pg.RMatrix(len(rFluid), len(rFluid[0]))
        for i in range(len(rFluid)):
            rB[i] = rFluid[i]
     
    porosity = pg.solver.parseArgToArray(porosity, mesh.cellCount(), mesh)
    a = pg.solver.parseArgToArray(a, mesh.cellCount(), mesh)
    m = pg.solver.parseArgToArray(m, mesh.cellCount(), mesh)
    S = pg.solver.parseArgToArray(S, mesh.cellCount(), mesh)
    n = pg.solver.parseArgToArray(n, mesh.cellCount(), mesh)
    
    r = pg.RMatrix(len(rB), len(rB[0]))
    for i in range(len(r)):
        r[i] = rB[i] * a * porosity**(-m) * S**(-n)
        
    r.round(1e-6)
    if meshI == None:
        if len(r) == 1:
            return r[0].copy()
        return r
    
    rI = pg.RMatrix(len(r), meshI.cellCount())
    if meshI:
        pg.interpolate(mesh, r, meshI.cellCenters(), rI) 
        
    if fill:
        for i in range(len(rI)):
            rI[i] = pg.solver.fillEmptyToCellArray(meshI, rI[i])
        
    rI.round(1e-6)
    #print(rI)
    if len(rI) == 1:
        return rI[0]
    return rI


if __name__ == "__main__":
    pass

