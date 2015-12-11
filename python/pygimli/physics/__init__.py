# -*- coding: utf-8 -*-

"""
Module containing submodules for various geophysical methods.
"""

from math import pi
from . em import FDEM, TDEM
from . sNMR import MRS
from . SIP import SIPSpectrum
from . traveltime import Refraction
# from . gravimetry import Gravimetry
# from . seismics import *

__all__ = ("FDEM", "TDEM", "MRS", "SIPSpectrum", "Refraction")


class constants:
    """ Container for different physical constants and units. """
    # magnetic constant, vacuum permeability
    mu0 = 4.0 * pi * 1e-7  # [(kg * m) / (A^2 * s^2)]
    # electric constant, vacuum permittivity
    e0 = 8.85418781762e-12  # [(A^2 * s^4)/(kg m^3)]
    # gravitational constant
    G = 6.6742e-11  # gravitation constant [m^3/(kg s^2)]
    GmGal = G / 1e-5  # [mGal=1e-5 SI]
    g = 9.798  # norm gravitational acceleration [m/s^2]
    # flow and transport
    Darcy = 9.86923e-13  # [m^2]
    mDarcy = Darcy * 1e-3


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
        self.fop = self.createFOP(debug)
        if self.fop is None:
            raise Exception("createFOP does not return valid forward operator")
        self.tD = None
        self.tM = None
        self.inv = self.createInv(self.fop, verbose, debug)
        if self.inv is None:
            raise Exception("createINV does not return valid inversion")

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        """ String representation of the class """
        return "Method Manager: " + str(self.__class__)

    def createFOP(self, refine=True):
        """ Create forward operator working on refined mesh """
        raise Exception("Overload me!")

    def createInv(self, fop, verbose=True, dosave=False):
        """ Create inversion instance, data- and model transformations. """
        raise Exception("Overload me!")

    # Visualization stuff
    def show(self, data, values=None, axes=None,
             cMin=None, cMax=None, colorBar=1, **kwargs):
        """ Show the model """
        pass

    def showData(self, ax=None, response=None):
        """ show data in form of travel time curves """
        pass

    def showMesh(self, ax=None):
        """ show mesh in given axes or in a new figure """
        pass

    def showResult(self, ax=None, cMin=None, cMax=None, logScale=False,
                   **kwargs):
        """ show resulting velocity vector """
        pass

    def getDepth(self):
        """ get typical investigation depth """
        pass

    # Mesh-related methods
    def createMesh(self):
        """ Create a mesh aka the parametrization """
        pass
        # could take dataContainer.sensors and getDepth independent of method

    def setParaMesh(self, mesh):
        """ Set mesh for the inverse problem """
        pass

    # Data-related methods
    def createData(self, sensors):
        """ Create an empty data set """
        pass

    def setData(self, data):
        """ Set data """
        pass

    def importData(self, filename):
        """ Import data """
        pass

    def checkData(self):
        """ Check data validity """
        pass

    def estimateError(self, absoluteError=0.001, relativeError=0.001):
        """ estimate error composed of an absolute and a relative part """
        pass

    # Workflow-related methods
    def invert(self, data, values=None, verbose=0, **kwargs):
        """ Invert the data and fill the parametrization. """
        raise('implement me')
        pass

    def simulate(self, values):  # , mesh, data=None):
        """ Run a simulation aka the forward task. """
        return self.fop.reponse(values)

    def saveResult(self, folder=None, size=(16, 10),
                   **kwargs):
        """ Saves the results in the specified folder. """
        pass

    def showResultAndFit(self, figsize=(10, 15), **kwargs):
        """ """
        pass

    def saveFigures(self, name=None, ext='pdf'):
        """save all existing figures to files"""
        if name is None and hasattr(self, 'basename'):
            name = self.basename
        if name is None or not any(name):
            name = 'out'
        for key in self.figs:
            self.figs[key].savefig(name+'-'+key+'.'+ext, bbox_inches='tight')


class MethodManager1D(MethodManager):
    """ Base manager class for 1D soundings """
    def __init__(self, nLay=3, nProp=1, zVec=None, type='block',
                 verbose=True, debug=False, **kwargs):
        MethodManager.__init__(self,
                               verbose=verbose, debug=debug, **kwargs)
        self.nProp = nProp
        if zVec is not None or type == 'occam':  # Occam type inversion
            if zVec is None:
                pass
            self.zVec = zVec
            self.nlay = len(zVec)
            self.type = 'occam'
        else:
            self.mesh = None
            self.type = 'block'

    def showData(self):
        """ """
        pass

    def showResultAndFit(self, figsize=(10, 15), **kwargs):
        """ """
        pass
