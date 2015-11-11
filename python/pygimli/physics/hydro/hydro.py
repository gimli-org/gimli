#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    WRITEME
"""

import pygimli as pg
import numpy as np

def permeabiltyEngelhardtPitter(poro, q=3.5, s=5e-3,
                                mesh=None, meshI=None):
    r"""
    Empirical model for porosity to hydraulic permeability.
    
    Postulated for sand and sandstones. [EngelhardtPit1955]
    
    .. math:: 
        k & = 2\cdot 10^7 \frac{\phi^2}{(1-\phi)^2}* \frac{1}{S^2} \\
        S & = q\cdot s \\
        s & = \sum_{i=1}(\frac{P_i}{r_i})
        
    * :math:`\phi` - poro 0.0 --1.0
    * :math:`S` - in cm^-1  specific surface in cm^2/cm^3
    * :math:`q` - (3 for spheres, > 3 shape differ from sphere)
        3.5 sand
    * :math:`s` - in cm^-1 (s = 1/r for particles with homogeneous radii r)
    * :math:`P_i` - Particle ration with radii :math:`r_i` on 1cm^3 Sample
    
    Parameters
    ----------
    
    Returns
    -------
    k :
        in Darcy
    """
    poro = parseArgToArray(poro, mesh.cellCount(), mesh)
    q = parseArgToArray(q, mesh.cellCount(), mesh)
    s = parseArgToArray(s, mesh.cellCount(), mesh)
    
    S = q * s
    k = 2e7 * (poro**2 / (1.0-poro)**2) * 1.0/S**2 * physics.constants.Darcy
        
    if meshI:
        k = pg.interpolate(mesh, k, meshI.cellCenters()) 
        k = pg.solver.fillEmptyToCellArray(meshI, k)
    return k


if __name__ == "__main__":
    pass

