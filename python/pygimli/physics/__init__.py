# -*- coding: utf-8 -*-

"""
Module containing submodules for various geophysical methods (physics)
"""

from . em import *
from . sNMR import *
from . SIP import *
from . traveltime import *
#from . gravimetry import *
#from . seismics import *

__all__ = ["FDEMData, TDEMData, MRS, SIPSpectrum, Refraction"]
