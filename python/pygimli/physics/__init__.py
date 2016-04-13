# -*- coding: utf-8 -*-

"""
Module containing submodules for various geophysical methods.
"""

from math import pi

from . em import FDEM, TDEM
from . sNMR import MRS
from . SIP import SIPSpectrum
# from . ert import resistivityArchie  # not in git yet (no ERT module)
# should rather be in pygimli/petrophysics as resistivity is not only ERT

# from . gravimetry import Gravimetry
# from . seismics import *

__all__ = ("FDEM", "TDEM", "MRS", "SIPSpectrum", "Refraction")


class constants:
    # magnetic constant, vacuum permeability
    mu0 = 4.0 * pi * 1e-7  # [(kg * m) / (A^2 * s^2)]

    # electric constant, vacuum permittivity
    e0 = 8.85418781762e-12  # [(A^2 * s^4)/(kg m^3)]

    G = 6.6742e-11  # [m^3/(kg s^2)]
    GmGal = G / 1e-5  # mGal

    Darcy = 9.86923e-13  # [m^2]

    g = 9.798  # [m/s^2]
    # also a g function of latitude (and altitude)?


class MethodManager(object):
    """
        General manager to maintenance a measurement method.

        The method manager holds one instance of a forward operator and a
        appropriate inversion method to handle simulation and reconstruction of
        common geophysical problems.

    """
    def __init__(self, verbose=True, debug=False, **kwargs):
        self.verbose = verbose
        self.debug = debug
        self.figs = {}

        self.fop = self.createFOP_(verbose)
        if self.fop is None:
            raise Exception("createFOP does not return valid forward operator")

        self.tD = None
        self.tM = None

        self.inv = self.createInv_(self.fop, verbose, debug)
        if self.inv is None:
            raise Exception("createINV does not return valid inversion")

        self.setVerbose(verbose)

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        """ String representation of the class """
        out = type(self).__name__ + " object"
        if hasattr(self, 'dataContainer'):
            out += "\n" + self.dataContainer.__str__()
        if hasattr(self, 'mesh'):
            out += "\n" + self.mesh.__str__()
        # some other stuff (data and model size)?
        return out
#        return "Method Manager: " + str(self.__class__)

    def setVerbose(self, verbose):
        """ make the class verbose (put output to the console) """
        self.verbose = verbose
        self.inv.setVerbose(verbose)
        self.fop.setVerbose(verbose)

#    @classmethod
    def createFOP_(self, verbose=False):
        """ Create forward operator working on refined mesh """
        return self.createFOP(verbose)

    def createInv_(self, fop, verbose=True, dosave=False):
        """ Create inversion instance, data- and model transformations. """
        return self.createInv(fop, verbose, dosave)

    # Data-related methods
    def createData(self, sensors, scheme):
        """ Create an empty data set """
        pass

    def setData(self, data):
        """ Set data """
        pass

    def checkData(self):
        """ Check data validity """
        pass

    def estimateError(self, absoluteError=0.001, relativeError=0.001):
        """ estimate error composed of an absolute and a relative part """
        pass

    def showData(self, ax=None, response=None):
        """ show data in form of travel time curves """
        pass

    # Mesh related methods
    def createMesh(self, ):
        """ Create a mesh aka the parametrization """
        pass

    def setMesh(self, ):
        """ Create a mesh aka the parametrization """
        pass

    def showMesh(self, ax=None):
        """ show mesh in given axes or in a new figure """
        pass

    # Work related methods
    def invert(self, data, values=None, verbose=0, **kwargs):
        """ Invert the data and fill the parametrization. """
        raise('implement me')
        pass

    def simulate(self, mesh, values, data=None):
        """ Run a simulation aka the forward task. """
        pass

    # Visualization stuff
    def show(self, data, values=None, axes=None,
             cMin=None, cMax=None, colorBar=1, **kwargs):
        """ Forward the visualization """
        pass

    def showResult(self, ax=None, cMin=None, cMax=None, logScale=False,
                   **kwargs):
        """ show resulting velocity vector """
        pass

    def saveResult(self, folder=None, size=(16, 10),
                   **kwargs):
        """
        Saves the results in the specified folder.
        """
        pass

    @staticmethod
    def createArgParser(dataSuffix='dat'):
        """
            Create default argument parser for the following options:

            -Q, --quiet

            -R, --robustData

            -B, --blockyModel

            -l, --lambda: options.lam

            -i, --maxIter: options.depth

            --depth: options.depth


        """
        import argparse

        parser = argparse.ArgumentParser(
            description="usage: %prog [options] *." + dataSuffix)
        parser.add_argument("-Q", "--quiet", dest="quiet",
                            action="store_true", default=False,
                            help="Be verbose.")
        parser.add_argument("-R", "--robustData", dest="robustData",
                            action="store_true", default=False,
                            help="Robust data (L1 norm) minimization.")
        parser.add_argument("-B", "--blockyModel", dest="blockyModel",
                            action="store_true", default=False,
                            help="Blocky model (L1 norm) regularization.")
        parser.add_argument('-l', "--lambda", dest="lam", type=float,
                            default=100,
                            help="Regularization strength.")
        parser.add_argument('-i', "--maxIter", dest="maxIter", type=int,
                            default=20,
                            help="Maximum iteration count.")
        parser.add_argument("--depth", dest="depth", type=float,
                            default=None,
                            help="Depth of inversion domain. [None=auto].")
        parser.add_argument('dataFileName')
        return parser

from . traveltime import Refraction
