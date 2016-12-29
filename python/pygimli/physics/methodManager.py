#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""TODO WRITEME"""

import pygimli as pg

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
        self.model = None

        self.mesh = None  # to be deleted if MethodManagerMesh is used # TODO
        self.dataContainer = None  # dto.

        self.fop = self.createFOP_(verbose)
        if self.fop is None:
            raise BaseException("createFOP does not return "
                                "valid forward operator")

        self.tD = None
        self.tM = None  # why not using a default RTransLog

        self.inv = self.createInv_(self.fop, verbose, debug)
        if self.inv is None:
            raise BaseException("createINV does not return valid inversion")

        self.setVerbose(verbose)

    def __str__(self):
        """String representation of the class."""
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

    def setDataToken(self, token):
        """Set the token name to identity the data in a DataContainer."""
        self.dataToken_ = token

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

    def model(self):
        """Return the actual model."""
        return self.model

    def createData(self, sensors, scheme):
        """Create an empty data set."""
        pass

    def setData(self, data):
        """Set data."""
        self.dataContainer = data

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

            -R, --robustData: options.robustData

            -B, --blockyModel: options.blockyModel

            -l, --lambda: options.lam

            -i, --maxIter: options.maxIter

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


class MethodManager1d(MethodManager):
    """Method Manager base class for managers on a 1d discretization."""
    def __init__(self, **kwargs):
        """Constructor."""
        super(MethodManager1d, self).__init__(**kwargs)
        self.nlay = kwargs.pop('nlay', 2)
        if 'nLayers' in kwargs:
            self.nlay = kwargs['nLayers']
        self.nProperties = kwargs.pop('nProperties', 1)
        self.Occam = kwargs.pop('Occam', False)  # member nameing!


class MethodManager1dProfile(MethodManager1d):  # does this make sense?
    """Method manager base class for 1D methods along a profile."""
    def __init__(self, **kwargs):
        """Constructor."""
        super(MethodManager1dProfile, self).__init__(**kwargs)


class MeshMethodManager(MethodManager):
    """Method Manager base class for managers using a (non-1D) mesh."""
    def __init__(self, **kwargs):
        """Constructor."""
        super(MeshMethodManager, self).__init__(**kwargs)
        self.mesh = None

    # Mesh related methods
    def createMesh(self, ):
        """Create a mesh aka the parametrization."""
        pass

    def setMesh(self, mesh, refine=True):
        """Set the internal mesh for this Manager.

        Inject the mesh in the internal fop und inv.

        Initialize RegionManager.
        For more than two regions the first is assumed to be background.

        Optional the forward mesh can be refined for higher numerical accuracy.

        Parameters
        ----------

        DOCUMENTME!!!

        """
        if isinstance(mesh, str):
            mesh = pg.load(mesh)

        if self.verbose:
            print(mesh)

        self.mesh = pg.Mesh(mesh)
        self.mesh.createNeighbourInfos()

        self.fop.setMesh(self.mesh)
        self.fop.regionManager().setConstraintType(1)
        if self.fop.regionManager().regionCount() > 1:
            self.fop.regionManager().region(1).setBackground(True)
#            self.fop.regionManager().regions().begin().second.setBackground(1)

        self.fop.createRefinedForwardMesh(refine)
        self.paraDomain = self.fop.regionManager().paraDomain()
        self.inv.setForwardOperator(self.fop)  # necessary? CR: check this

    def showMesh(self, ax=None):
        """Show mesh in given axes or in a new figure."""
        raise Exception('IMPLEMENTME')
        pass

    def setData(self, data):
        """Set data container from outside."""
        isinstance
        if not isinstance(data, pg.DataContainer):
            raise Exception('IMPLEMENTME')

            if isinstance(data,str):
                raise Exception('IMPLEMENTME')
        else:
            self.dataContainer = data

        if not self.dataContainer.allNonZero(self.dataToken()):
            raise BaseException("No or partial values for:", self.dataToken())

        if not self.dataContainer.allNonZero('err'):
            raise BaseException("No or partial values for err")

        self.fop.setData(self.dataContainer)

    @staticmethod
    def createArgParser(dataSuffix='dat'):
        """Create argument parser for the manager."""
        parser = MethodManager.createArgParser(dataSuffix)

        parser.add_argument("--paraMaxCellSize", dest="maxCellArea",
                            type=float, default=0.0,
                            help="Maximum cell size for parametric cells "
                                 "in m² [0.0] 0.0 means no restriction.")
        parser.add_argument("--zWeight", dest="zWeight",
                            type=float, default=0.7,
                            help="Weight for vertical smoothness "
                                 "(1=isotrope). [0.7]")
        return parser

    def invert(self, data=None, vals=None, err=None, mesh=None, **kwargs):
        """Run the full inversion.

        The data and error needed to be set before.
        The meshes will be created if necessary.

        Parameters
        ----------

        DOCUMENTME!!!


        **kwargs

            * lam : float [20]
                regularization parameter
            * zWeight : float [0.7]
                relative vertical weight
            * maxIter : int [20]
                maximum iteration number
            * robustdata : bool [False]
                robust data reweighting using an L1 scheme (IRLS reweighting)
            * blockymodel : bool [False]
                blocky model constraint using L1 reweighting roughness vector
            * startModelIsReference : bool [False]
                startmodel is the reference model for the inversion

            forwarded to createMesh

            * depth
            * quality
            * paraDX
            * maxCellArea

        """
        if 'verbose' in kwargs:
            self.setVerbose(kwargs.pop('verbose'))

        if data is not None:
            # setDataContainer would be better
            self.setData(data)

        if vals is not None:
            self.dataContainer.set(self.dataToken(), vals)

        if err is not None:
            self.dataContainer.set('err', vals)

        # check for data container here
        dataVals = self.dataContainer(self.dataToken())
        errVals = self.dataContainer('err')

        if mesh is not None:
            self.setMesh(mesh)

        if self.mesh is None:
            self.createMesh(depth=kwargs.pop('depth', None),
                            quality=kwargs.pop('quality', 34.0),
                            maxCellArea=kwargs.pop('maxCellArea', 0.0),
                            paraDX=kwargs.pop('paraDX', 0.3))


        self.inv.setData(dataVals)
        self.inv.setError(errVals, not self.errIsAbsolute)

        zWeight = kwargs.pop('zWeight', 0.7)
        if 'zweight' in kwargs:
            zWeight = kwargs.pop('zweight', 0.7)
            print("zweight option will be removed soon. Please use zWeight.")

        self.fop.regionManager().setZWeight(zWeight)

        self.inv.setLambda(kwargs.pop('lam', 20))
        self.inv.setMaxIter(kwargs.pop('maxIter', 20))
        self.inv.setRobustData(kwargs.pop('robustData', False))
        self.inv.setBlockyModel(kwargs.pop('blockyModel', False))
        self.inv.setRecalcJacobian(kwargs.pop('recalcJacobian', True))

        # TODO: ADD MORE KWARGS
        pc = self.fop.regionManager().parameterCount()

        startModel = kwargs.pop('startModel',
                                pg.RVector(pc, pg.median(dataVals)))

        self.inv.setModel(startModel)

        if kwargs.pop('startModelIsReference', False):
            self.inv.setReferenceModel(startModel)

        # Run the inversion
        if len(kwargs) > 0:
            print("Keyword arguments unknown:")
            print(kwargs)
            print("Warning! There are unknown kwargs arguments.")

        model = self.inv.run()
        self.model = model(self.paraDomain.cellMarkers())

        return self.model
