# -*- coding: utf-8 -*-

"""
Module containing submodules for various geophysical methods.
"""

from . em import FDEMData, TDEMData
from . sNMR import MRS
from . SIP import SIPSpectrum
from . traveltime import Refraction
#from . gravimetry import Gravimetry
# from . seismics import *

__all__ = ("FDEMData", "TDEMData", "MRS", "SIPSpectrum", "Refraction")

from math import pi

class constants:
    mu0 = 4.0 * pi * 1e-7
    
    G = 6.6742e-11  # [m^3/(kg s^2)]
    GmGal = G / 1e-5  # mGal
    
    Darcy = 9.86923e-13 #[m^2]
    
    g = 9.798 #[m/s^2]