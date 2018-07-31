#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Module containing submodules for various geophysical methods.
"""

from .constants import Constants
from .em import FDEM, TDEM
from .ert import ERTManager, ERTModelling, VESManager
from .SIP import SIPSpectrum
from .sNMR import MRS
from .traveltime import Refraction

constants = Constants

# from . gravimetry import Gravimetry
# from . seismics import *

# __all__ = ["FDEM", "TDEM", "MRS", "SIPSpectrum", "Refraction"]
