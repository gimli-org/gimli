#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Module containing submodules for various geophysical methods.
"""

from math import pi

from .em import FDEM, TDEM
from .sNMR import MRS
from .SIP import SIPSpectrum

from .ert import ERTModelling

from .methodmanger import MethodManager
from .methodmanger import MeshMethodManager
from .traveltime import Refraction

# from . gravimetry import Gravimetry
# from . seismics import *

# __all__ = ["FDEM", "TDEM", "MRS", "SIPSpectrum", "Refraction"]

class Constants(object):
    """"TODO WRITEME"""
    # magnetic constant, vacuum permeability
    mu0 = 4.0 * pi * 1e-7  # [(kg * m) / (A^2 * s^2)]

    # electric constant, vacuum permittivity
    e0 = 8.85418781762e-12  # [(A^2 * s^4)/(kg m^3)]

    G = 6.6742e-11  # [m^3/(kg s^2)]
    GmGal = G / 1e-5  # mGal

    Darcy = 9.86923e-13  # [m^2]

    g = 9.798  # [m/s^2]
    # also a g function of latitude (and altitude)?

constants = Constants
