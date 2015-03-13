# -*- coding: utf-8 -*-

"""
Module containing submodules for various geophysical methods (physics)
"""

from . em import FDEMData, TDEMData
from . sNMR import MRS
from . SIP import SIPSpectrum
from . traveltime import Refraction
#from . gravimetry import *
#from . seismics import *

__all__ = [FDEMData, TDEMData, MRS, SIPSpectrum, Refraction]
