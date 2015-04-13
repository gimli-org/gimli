# -*- coding: utf-8 -*-

"""
Module containing submodules for various geophysical methods.
"""

from . em import FDEM, TDEM
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


#class MethodManager(object):
    #"""
        #General manager to maintenance a measurement method.

        #The method manager holds one instance of a forward operator and a
        #appropriate inversion method to handle simulation and reconstruction of
        #common geophysical problems.

    #"""
    #def __init__(self, verbose=False):
        #self.fop = self.createFOP(verbose)
        #self.tD = None
        #self.tM = None
        #self.inv = self.createInv(verbose)

    #def __str__(self):
        #return self.__repr__()

    #def __repr__(self):
        #""" String representation of the class """
        #out = "Abstract MethodManager."
        #pass

    #def createFOP(self, refine=True):
        #""" Create forward operator working on refined mesh """
        #raise Exception("Overload me!")

    #def createInv(self, verbose=True, dosave=False):
        #""" Create inversion instance, data- and model transformations. """
        #raise Exception("Overload me!")

    ## Visualization stuff
    #def show(self, data, values=None, axes=None,
             #cMin=None, cMax=None, colorBar=1, **kwargs):
        #""" Forward the visualization """
        #pass

    #def showData(self, ax=None, response=None):
        #""" show data in form of travel time curves """
        #pass

    #def showMesh(self, ax=None):
        #""" show mesh in given axes or in a new figure """
        #pass

    #def showResult(self, ax=None, cMin=None, cMax=None, logScale=False,
                   #**kwargs):
        #""" show resulting velocity vector """
        #pass

    ## Mesh related methods
    #def createMesh(self, ):
        #""" Create a mesh aka the parametrization """
        #pass

    #def setParaMesh(self, mesh):
        #""" Set mesh for the inverse problem """
        #pass

    ## Data related methods
    #def createData(self, sensors, scheme):
        #""" Create an empty data set """
        #pass

    #def setData(self, data):
        #""" Set data """
        #pass

    #def importData(self, filename):
        #""" Import data """
        #pass

    #def checkData(self):
        #""" Check data validity """
        #pass

    #def estimateError(self, absoluteError=0.001, relativeError=0.001):
        #""" estimate error composed of an absolute and a relative part """
        #pass

    ## Work related methods
    #def invert(self, data, values=None, verbose=0, **kwargs):
        #""" Invert the data and fill the parametrization. """
        #raise('implement me')
        #pass

    #def simulate(self, mesh, values, data=None):
        #""" Run a simulation aka the forward task. """
        #pass

    #def saveResult(self, folder=None, size=(16, 10),
                   #**kwargs):
        #"""
        #Saves the results in the specified folder.
        #"""
        #pass
