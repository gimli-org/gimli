#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""TODO WRITEME"""


class MethodManager(object):
    """General manager to maintenance a measurement method.

        The method manager holds one instance of a forward operator and a
        appropriate inversion method to handle simulation and reconstruction of
        common geophysical problems.
    """
    def __init__(self, verbose=True, debug=False):
        """Constructor."""
        self.verbose = verbose
        self.debug = debug
        self.figs = {}
        self.dataToken_ = 'nan'
        self.errIsAbsolute = False

        self.mesh = None
        self.dataContainer = None

        self.fop = self.createFOP_(verbose)
        if self.fop is None:
            raise BaseException("createFOP does not return "
                                "valid forward operator")

        self.tD = None
        self.tM = None

        self.inv = self.createInv_(self.fop, verbose, debug)
        if self.inv is None:
            raise BaseException("createINV does not return valid inversion")

        self.setVerbose(verbose)

    def __str__(self):
        """TODO WRITEME."""
        return self.__repr__()

    def __repr__(self):
        """String representation of the class."""
        out = type(self).__name__ + " object"
        if hasattr(self, 'dataContainer'):
            out += "\n" + self.dataContainer.__str__()
        if hasattr(self, 'mesh'):
            out += "\n" + self.mesh.__str__()
        # some other stuff (data and model size)?
        return out
        # return "Method Manager: " + str(self.__class__)

    def setVerbose(self, verbose):
        """Make the class verbose (put output to the console)"""
        self.verbose = verbose
        self.inv.setVerbose(verbose)
        self.fop.setVerbose(verbose)

    def dataToken(self):
        """Token name for the data in a DataContainer."""
        if self.dataToken_ == 'nan':
            print("Warning! the Manager don't know the data token")
        return self.dataToken_

    def apparentData(self):
        """Convert data into apparent data."""
        raise BaseException("IMPLEMENTME in derived class")

#    @classmethod
    def createFOP_(self, verbose=False):
        """Create forward operator working on refined mesh."""
        return self.createFOP(verbose)

    def createInv_(self, fop, verbose=True, dosave=False):
        """Create inversion instance, data- and model transformations."""
        return self.createInv(fop, verbose, dosave)

    # Data-related methods
    def createData(self, sensors, scheme):
        """Create an empty data set."""
        pass

    def setData(self, data):
        """Set data."""
        pass

    def checkData(self):
        """Check data validity."""
        pass

    @staticmethod
    def estimateError(data, absoluteError=0.001, relativeError=0.001):
        """Estimate error composed of an absolute and a relative part."""
        pass

    def showData(self, axes=None, response=None, name='data'):
        """Show data."""
        pass

    # Work related methods
    def invert(self, **kwargs):
        """Invert the data and fill the parametrization."""
        raise BaseException('implement me in derived class' + str(**kwargs))

    def simulate(self, **kwargs):
        """Run a simulation aka the forward task."""
        raise BaseException('implement me in derived class' + str(**kwargs))

    # Visualization stuff
    def show(self, data, values=None, axes=None,
             cMin=None, cMax=None, colorBar=1, **kwargs):
        """Forward the visualization."""
        pass

    def showResult(self, ax=None, cMin=None, cMax=None, logScale=False,
                   **kwargs):
        """Show resulting vector."""
        pass

    def saveResult(self, folder=None, size=(16, 10),
                   **kwargs):
        """Save results in the specified folder."""
        pass

    @staticmethod
    def createArgParser(dataSuffix='dat'):
        """Create default argument parser.

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


class MeshMethodManager(MethodManager):
    """TODO WRITEME """
    def __init__(self, **kwargs):
        """Conszructor."""
        super(MeshMethodManager, self).__init__(**kwargs)

    # Mesh related methods
    def createMesh(self, ):
        """Create a mesh aka the parametrization."""
        pass

    def setMesh(self, ):
        """Create a mesh aka the parametrization."""
        pass

    def showMesh(self, ax=None):
        """Show mesh in given axes or in a new figure."""
        pass

    @staticmethod
    def createArgParser(dataSuffix='dat'):
        """Create argument parser for the manager."""
        parser = MethodManager.createArgParser(dataSuffix)

        parser.add_argument("--paraMaxCellSize", dest="maxCellArea",
                            type=float, default=None,
                            help="Maximum cell size for parametric cells "
                                 "in m² [None=auto].")
        parser.add_argument("--zWeight", dest="zWeight",
                            type=float, default=0.7,
                            help="Weight for vertical smoothness "
                                 "(1=isotrope). [0.7]")
        return parser
