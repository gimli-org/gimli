#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Method Manager

Provide the end user interface for method (geophysical) dependent
modelling and inversion as well as data and model visualization.
"""

import numpy as np
import pygimli as pg

from pygimli.utils import prettyFloat as pf


# Discuss .. rename to Framework or InversionFramework since he only manages
# the union of Inversion/Modelling and RegionManager(later)
class MethodManager(object):
    """General manager to maintenance a measurement method.

    Method Manager are the interface to end-user interaction and can be seen as
    simple but complete application classes which manage all tasks of
    geophysical data processing.

    The method manager holds one instance of a forward operator and an
    appropriate inversion framework to handle modelling and data inversion.

    Method Manager also helps with data import and export,
    handle measurement data error estimation as well as model and data
    visualization.

    Attributes
    ----------
    verbose : bool
        Give verbose output.

    debug : bool
        Give debug output.

    fop : :py:mod:`pygimli.frameworks.Modelling`
        Forward Operator instance .. knows the physics.
        fop is initialized by
        :py:mod:`pygimli.manager.MethodManager.initForwardOperator`
        and calls a valid
        :py:mod:`pygimli.manager.MethodManager.createForwardOperator` method
        in any derived classes.

    inv : :py:mod:`pygimli.frameworks.Inversion`.
        Inversion framework instance .. knows the reconstruction approach.
        The attribute inv is initialized by default but can be changed
        overwriting
        :py:mod:`pygimli.manager.MethodManager.initInversionFramework`
    """
    def __init__(self, fop=None, fw=None, **kwargs):
        """Constructor."""
        self._fop = fop
        self._fw = fw
        self._verbose = kwargs.pop('verbose', False)
        self._debug = kwargs.pop('debug', False)

        # The inversion framework
        self._initInversionFramework(verbose=self._verbose,
                                     debug=self._debug)

        # The forward operator is stored in self._fw
        self._initForwardOperator(verbose=self._verbose, **kwargs)

        # maybe obsolete
        self.figs = {}
        self.errIsAbsolute = False

    def __hash__(self):
        """Create a hash for Method Manager"""
        return pg.utils.strHash(str(type(self))) ^ hash(self.fop)

    @property
    def verbose(self):
        return self._verbose

    @verbose.setter
    def verbose(self, v):
        self._verbose = v
        self.fw.verbose = self._verbose

    @property
    def debug(self):
        return self._debug

    @debug.setter
    def debug(self, v):
        self._debug = v
        self.fw.debug = self._debug

    @property
    def fw(self):
        return self._fw

    @property
    def fop(self):
        return self.fw.fop

    @property
    def inv(self):
        return self.fw

    @property
    def model(self):
        return self.fw.model

    def reinitForwardOperator(self, **kwargs):
        """Reinitialize the forward operator.

        Sometimes it can be useful to reinitialize the forward operator.
        Keyword arguments will be forwarded to 'self.createForwardOperator'.
        """
        self._initForwardOperator(**kwargs)

    def _initForwardOperator(self, **kwargs):
        """Initialize or re-initialize the forward operator.

        Called once in the constructor to force the manager to create the
        necessary forward operator member. Can be recalled if you need to
        changed the mangers own forward operator object. If you want an own
        instance of a valid FOP call createForwardOperator.
        """
        if self._fop is not None:
            fop = self._fop
        else:
            fop = self.createForwardOperator(**kwargs)

        if fop is None:
            pg.critical("It seems that createForwardOperator method "
                        "does not return a valid forward operator.")
        if self.fw is not None:
            self.fw.reset()
            self.fw.setForwardOperator(fop)
        else:
            pg.critical("No inversion framework defined.")

    def createForwardOperator(self, **kwargs):
        """Mandatory interface for derived classes.

        Here you need to specify which kind of forward operator FOP
        you want to use.
        This is called by any initForwardOperator() call.

        Parameters
        ----------
        **kwargs
            Any arguments that are necessary for your FOP creation.

        Returns
        -------
        Modelling
            Instance of any kind of :py:mod:`pygimli.framework.Modelling`.
        """
        # pg.critical("No forward operator defined, either give one or"
        #             "overwrite in derived class")
        useBert = kwargs.pop('useBert', False)
        verbose = kwargs.pop('verbose', False)
        if useBert:
            pg.verbose('Create BertModelling FOP')
            fop = BertModelling(sr=kwargs.pop('sr', True), verbose=verbose)
        else:
            pg.verbose('Create ERTModelling FOP')
            fop = ERTModelling(**kwargs)

        return fop

    def _initInversionFramework(self, **kwargs):
        """Initialize or re-initialize the inversion framework.

        Called once in the constructor to force the manager to create the
        necessary Framework instance.
        """
        self._fw = self.createInversionFramework(**kwargs)

        if self.fw is None:
            pg.critical("createInversionFramework does not return "
                        "valid inversion framework.")

    def createInversionFramework(self, **kwargs):
        """Create default Inversion framework.

        Derived classes may overwrite this method.

        Parameters
        ----------

        **kwargs
            Any arguments that are necessary for your creation.

        Returns
        -------
        Inversion
            Instance of any kind of :py:mod:`pygimli.framework.Inversion`.
        """
        if self._fw is None:
            return pg.frameworks.Inversion(**kwargs)
        else:
            return self._fw

    def estimateError(self, data, errLevel=0.01, absError=None):
        """Estimate data error.

        Create an error of estimated measurement error.
        On default it returns an array of constant relative errors.
        More sophisticated error estimation should be done
        in specialized derived classes.

        TODO check, rel or abs in return.

        Parameters
        ----------
        data : iterable
            Data values for which the errors should be estimated.

        errLevel : float (0.01)
            Error level in percent/100.

        absoluteError : float (None)
            TODO

        Returns
        -------
        err : array
            Returning array of size len(data)
        """
        if absError is not None:
            return absError + data * errLevel

        return np.ones(len(data)) * errLevel

    def simulate(self, model, **kwargs):
        """Run a simulation aka the forward task."""

        ra = self.fop.response(par=model)

        noiseLevel = kwargs.pop('noiseLevel', 0.0)
        if noiseLevel > 0:
            err = self.estimateError(ra, errLevel=noiseLevel)
            ra *= 1. + pg.math.randn(ra.size()) * err
            return ra, err

        return ra

    def dataCheck(self, data):
        """Overwrite for special checks to return data values"""
        # if self._dataToken == 'nan':
        #     pg.critical('self._dataToken nan, should be set in class', self)
        #     return data(self._dataToken)
        return data

    def _ensureData(self, data):
        """Check data validity"""
        if data is None:
            data = self.fw.dataVals

        vals = self.dataCheck(data)

        if vals is None:
            pg.critical("There are no data values.")

        if abs(min(vals)) < 1e-12:
            print(min(vals), max(vals))
            pg.critical("There are zero data values.")

        return vals

    def errorCheck(self, err, dataVals=None):
        """Return relative error. Default we assume 'err' are relative values.
        Overwrite is derived class if needed. """
        if isinstance(err, pg.DataContainer):
            if not err.haveData('err'):
                pg.error('Datacontainer have no "err" values. '
                         'Fallback set to 0.01')
            return err['err']

        return err

    def _ensureError(self, err, dataVals=None):
        """Check error validity"""
        if err is None:
            err = self.fw.errorVals

        vals = self.errorCheck(err, dataVals)

        if min(vals) <= 0:
            pg.critical("All error values need to be larger then 0."
                        " either give and err argument or fill dataContainer "
                        " with a valid 'err' ", min(vals), max(vals))

        return vals

    def preRun(self, *args, **kwargs):
        """Called just before the inversion run starts."""
        pass

    def postRun(self, *args, **kwargs):
        """Called just after the inversion run."""
        pass

    def invert(self, data=None, err=None, **kwargs):
        """Invert the data.

        Invert the data by calling self.inv.run() with mandatory data and
        error values.

        TODO
            *need dataVals mandatory? what about already loaded data

        Parameters
        ----------
        dataVals : iterable
            Data values to be inverted.

        errVals : iterable | float
            Error value for the given data.
            If errVals is float we assume this means to be a global relative
            error and force self.estimateError to be called.
        """
        dataVals = self._ensureData(data)
        errVals = self._ensureError(err, dataVals)

        self.preRun(**kwargs)
        self.fw.run(dataVals, errVals, **kwargs)
        self.postRun(**kwargs)

        return self.fw.model

    def showModel(self, model, ax=None, **kwargs):
        """Shows a model.

        Draw model into a given axes or show inversion result from last run.
        Forwards on default to the self.fop.drawModel function
        of the modelling operator.
        If there is no function given, you have to override this method.

        Parameters
        ----------
        ax : mpl axes
            Axes object to draw into. Create a new if its not given.

        model : iterable
            Model data to be draw.

        Returns
        -------
        ax, cbar
        """
        if ax is None:
            fig, ax = pg.plt.subplots()

        return self.fop.drawModel(ax, model, **kwargs)

    def showData(self, data, ax=None, **kwargs):
        """Shows the data.

        Draw data values into a given axes or show the data values from
        the last run.
        Forwards on default to the self.fop.drawData function
        of the modelling operator.
        If there is no given function given, you have to override this method.

        Parameters
        ----------
        ax : mpl axes
            Axes object to draw into. Create a new if its not given.

        data : iterable | pg.DataContainer
            Data values to be draw.

        """
        if ax is None:
            fig, ax = pg.plt.subplots()

        return self.fop.drawData(ax, data, **kwargs)

    def showResult(self, model=None, ax=None, **kwargs):
        """Show the last inversion result.

        TODO
        ----
         DRY: decide showModel or showResult

        Parameters
        ----------
        ax : mpl axes
            Axes object to draw into. Create a new if its not given.

        model : iterable [None]
            Model values to be draw. Default is self.model from the last run

        Returns
        -------
        ax, cbar
        """
        if model is None:
            model = self.model
        return self.showModel(model, ax=ax, **kwargs)

    def showFit(self, ax=None, **kwargs):
        """Show the last inversion data and response."""
        ax = self.showData(data=self.inv.dataVals,
                           error=self.inv.errorVals,
                           label='Data',
                           ax=ax, **kwargs)
        ax = self.showData(data=self.inv.response,
                           label='Response',
                           ax=ax, **kwargs)

        if not kwargs.pop('hideFittingAnnotation', False):
            ax.text(0.99, 0.005, r"rrms: {0}, $\chi^2$: {1}".format(
                pf(self.fw.inv.relrms()), pf(self.fw.inv.chi2())),
                    transform=ax.transAxes,
                    horizontalalignment='right',
                    verticalalignment='bottom',
                    fontsize=8)

        if not kwargs.pop('hideLegend', False):
            ax.legend()

        return ax

    def showResultAndFit(self, **kwargs):
        """Calls showResults and showFit."""

        fig = pg.plt.figure()
        ax = fig.add_subplot(1, 2, 1)

        self.showResult(ax=ax, model=self.model, **kwargs)

        ax1 = fig.add_subplot(2, 2, 2)
        ax2 = fig.add_subplot(2, 2, 4)

        self.showFit(axs=[ax1, ax2], **kwargs)

        fig.tight_layout()
        return fig

    @staticmethod
    def createArgParser(dataSuffix='dat'):
        """Create default argument parser.

        TODO move this to some kind of app class

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
#        parser.add_argument("-R", "--robustData", dest="robustData",
#                            action="store_true", default=False,
#                            help="Robust data (L1 norm) minimization.")
#        parser.add_argument("-B", "--blockyModel", dest="blockyModel",
#                            action="store_true", default=False,
#                            help="Blocky model (L1 norm) regularization.")
        parser.add_argument('-l', "--lambda", dest="lam", type=float,
                            default=100,
                            help="Regularization strength.")
        parser.add_argument('-i', "--maxIter", dest="maxIter", type=int,
                            default=20,
                            help="Maximum iteration count.")
#        parser.add_argument("--depth", dest="depth", type=float,
#                            default=None,
#                            help="Depth of inversion domain. [None=auto].")
        parser.add_argument('dataFileName')
        return parser


class ParameterInversionManager(MethodManager):
    """Framework to invert unconstrained parameters."""
    def __init__(self, funct=None, fop=None, **kwargs):
        """Constructor."""
        if fop is not None:
            if not isinstance(fop, pg.frameworks.ParameterModelling):
                pg.critical("We need a fop if type ",
                            pg.frameworks.ParameterModelling)

#        fop = pg.frameworks.ParameterModelling(fop, petro)

        super(ParameterInversionManager, self).__init__(fop, **kwargs)

    def createInversionFramework(self, **kwargs):
        """
        """
        return pg.frameworks.MarquardtInversion(**kwargs)

    def invert(self, data=None, err=None, **kwargs):
        """
        Parameters
        ----------
        limits: {str: [min, max]}
            Set limits for parameter by parameter name.
        startModel: {str: startModel}
            Set the start value for parameter by parameter name.
        """
        limits = kwargs.pop('limits', {})

        for k, v in limits.items():
            self.fop.setRegionProperties(k, limits=v)

        startModel = kwargs.pop('startModel', {})

        if isinstance(startModel, dict):
            for k, v in startModel.items():
                self.fop.setRegionProperties(k, startModel=v)
        else:
            kwargs['startModel'] = startModel

        return super(ParameterInversionManager, self).invert(data=data,
                                                             err=err,
                                                             **kwargs)


class MethodManager1d(MethodManager):
    """Method Manager base class for managers on a 1d discretization."""
    def __init__(self, fop=None, **kwargs):
        """Constructor."""
        super(MethodManager1d, self).__init__(fop, **kwargs)

    def createInversionFramework(self, **kwargs):
        """
        """
        return pg.frameworks.Block1DInversion(**kwargs)

    def invert(self, data=None, err=None, **kwargs):
        """ """
        return super(MethodManager1d, self).invert(data=data, err=err,
                                                   **kwargs)


class MeshMethodManager(MethodManager):
    def __init__(self, **kwargs):
        """Constructor."""
        super(MeshMethodManager, self).__init__(**kwargs)

    def paraModel(self, model=None):
        """Give the model parameter regarding the parameter mesh."""
        if model is None:
            model = self.fw.model

        if len(model) == self.fw.parameterCount:
            return model
        else:
            self.fop.paraModel(model)

    def setMesh(self, mesh, ignoreRegionManager=False, **kwargs):
        """ """
        if ignoreRegionManager:
            mesh = self.fop.createRefinedFwdMesh(mesh, **kwargs)

        self.fop.setMesh(mesh, ignoreRegionManager=ignoreRegionManager)

    def setData(self, data):
        """ """
        if isinstance(data, pg.DataContainer):
            self.fop.data = data
        else:
            pg.critical("setting data array is not yet implemented.")

    def invert(self, data=None, mesh=None, zWeight=1.0, startModel=None,
               **kwargs):
        """Run the full inversion.

        Parameters
        ----------
        data : pg.DataContainer

        mesh : pg.Mesh [None]

        zWeight : float [1.0]

        startModel : float | iterable [None]

            If set to None fop.createDefaultStartModel(dataValues) is called.

        Other Parameters
        ----------------
        forwarded to Inversion.run

        Returns
        -------
        model : array
            Model mapped for match the paraDomain Cell markers.
            The calculated model is in self.fw.model.
        """
        if data is not None:
            self.setData(data)

        if mesh is not None:
            self.setMesh(mesh)

        self.fop._refineP2 = kwargs.pop('refineP2', False)

        dataVals = self._ensureData(self.fop.data)
        errorVals = self._ensureError(self.fop.data, dataVals)

        if self.fop.mesh() is None:
            pg.critical('Please provide a mesh')

        # inversion will call this itsself as default behaviour
        # if startModel is None:
        #     startModel = self.fop.createStartModel(dataVals)

        # pg._g('invert-dats', dataVals)
        # pg._g('invert-err', errVals)
        # pg._g('invert-sm', startModel)

        kwargs['startModel'] = startModel

        self.fop.setRegionProperties('*', zWeight=zWeight)

        # Limits is no mesh related argument here or base??
        limits = kwargs.pop('limits', None)
        if limits is not None:
            self.fop.setRegionProperties('*', limits=limits)

        # pg._y(pg.pf(self.fop._regionProperties))

        self.preRun(**kwargs)
        self.fw.run(dataVals, errorVals, **kwargs)
        self.postRun(**kwargs)
        return self.paraModel(self.fw.model)

    def showFit(self, axs=None, **kwargs):
        """Show data and the inversion result model response."""
        orientation = 'vertical'
        if axs is None:
            fig, axs = pg.plt.subplots(nrows=1, ncols=2)
            orientation = 'horizontal'

        self.showData(data=self.inv.dataVals,
                      orientation=orientation,
                      ax=axs[0], **kwargs)
        axs[0].text(0.0, 1.03, "Data",
                    transform=axs[0].transAxes,
                    horizontalalignment='left',
                    verticalalignment='center')

        resp = None
        data = None
        if 'model' in kwargs:
            resp = self.fop.response(kwargs['model'])
            data = self._ensureData(self.fop.data)
        else:
            resp = self.inv.response
            data = self.fw.dataVals

        self.showData(data=resp,
                      orientation=orientation,
                      ax=axs[1], **kwargs)
        axs[1].text(0.0, 1.03, "Response",
                    transform=axs[1].transAxes,
                    horizontalalignment='left',
                    verticalalignment='center')

        axs[1].text(1.0, 1.03, r"rrms: {0}%, $\chi^2$: {1}".format(
                pg.pf(pg.utils.rrms(data, resp)*100),
                pg.pf(self.fw.chi2History[-1])),
                    transform=axs[1].transAxes,
                    horizontalalignment='right',
                    verticalalignment='center')

        # if not kwargs.pop('hideFittingAnnotation', False):
        #     axs[0].text(0.01, 1.0025, "rrms: {0}, $\chi^2$: {1}"
        #             .format(pg.utils.prettyFloat(self.fw.inv.relrms()),
        #                     pg.utils.prettyFloat(self.fw.inv.chi2())),
        #             transform=axs[0].transAxes,
        #             horizontalalignment='left',
        #             verticalalignment='bottom')

        return axs

    def coverage(self):
        """Return coverage vector considering the logarithmic transformation.
        """
        covTrans = pg.core.coverageDCtrans(self.fop.jacobian(),
                                           1.0 / self.inv.response,
                                           1.0 / self.inv.model)
        nCells = self.fop.paraDomain.cellCount()
        return np.log10(covTrans[:nCells] / self.fop.paraDomain.cellSizes())

    def standardizedCoverage(self, threshhold=0.01):
        """Return standardized coverage vector (0|1) using thresholding.
        """
        return 1.0*(abs(self.coverage()) > threshhold)


class PetroInversionManager(MeshMethodManager):
    def __init__(self, petro, mgr=None, **kwargs):
        petrofop = kwargs.pop('petrofop', None)

        if petrofop is None:
            fop = kwargs.pop('fop', None)

            if fop is None and mgr is not None:
                fop = mgr.fw.fop
                self.dataCheck = mgr.dataCheck
                self.errorCheck = mgr.errorCheck

            if fop is not None:
                if not isinstance(fop, pg.frameworks.PetroModelling):
                    petrofop = pg.frameworks.PetroModelling(fop, petro)

        if petrofop is None:
            print(mgr)
            print(fop)
            pg.critical('implement me')

        super(PetroInversionManager, self).__init__(fop=petrofop, **kwargs)




































# class MethodManager1dProfile(MethodManager1d):  # does this make sense?
#     """Method manager base class for 1D methods along a profile."""
#     def __init__(self, **kwargs):
#         """Constructor."""
#         pg.critical('to be removed')
#         super(MethodManager1dProfile, self).__init__(**kwargs)


# class MethodManager0(object):
#     """General manager to maintenance a measurement method.

#         The method manager holds one instance of a forward operator and a
#         appropriate inversion method to handle simulation and reconstruction of
#         common geophysical problems.
#     """
#     def __init__(self, verbose=True, debug=False):
#         """Constructor."""
#         # pg.critical('to be removed')
#         self.verbose = verbose
#         self.debug = debug
#         self.figs = {}
#         self.dataToken_ = 'nan'
#         self.errIsAbsolute = False
#         self.model = None

#         self.mesh = None  # to be deleted if MethodManagerMesh is used # TODO
#         self.dataContainer = None  # dto.

#         self.fop = self.createFOP_(verbose)
#         if self.fop is None:
#             raise BaseException("createFOP does not return "
#                                 "valid forward operator")

#         self.tD = None
#         self.tM = None  # why not using a default RTransLog

#         self.inv = self.createInv_(self.fop, verbose, debug)
#         if self.inv is None:
#             raise BaseException("createINV does not return valid inversion")

#         self.setVerbose(verbose)

#     def __str__(self):
#         """String representation of the class."""
#         return self.__repr__()

#     def __repr__(self):
#         """String representation of the instance."""
#         out = type(self).__name__ + " object"
#         if hasattr(self, 'dataContainer'):
#             out += "\n" + self.dataContainer.__str__()
#         if hasattr(self, 'mesh'):
#             out += "\n" + self.mesh.__str__()
#         # some other stuff (data and model size)?
#         return out
#         # return "Method Manager: " + str(self.__class__)

#     def setVerbose(self, verbose):
#         """Make the class verbose (put output to the console)"""
#         self.verbose = verbose
#         self.inv.setVerbose(verbose)
#         self.fop.setVerbose(verbose)

#     def setDataToken(self, token):
#         """Set the token name to identity the data in a DataContainer."""
#         self.dataToken_ = token

#     def dataToken(self):
#         """Token name for the data in a DataContainer."""
#         if self.dataToken_ == 'nan':
#             print("Warning! the Manager don't know the data token")
#         return self.dataToken_

#     def apparentData(self):
#         """Convert data into apparent data."""
#         raise BaseException("IMPLEMENTME in derived class")

# #    @classmethod
#     def createFOP_(self, verbose=False):
#         """Create forward operator working on refined mesh."""
#         return self.createFOP(verbose)

#     def createInv_(self, fop, verbose=True, dosave=False):
#         """Create inversion instance, data- and model transformations."""
#         return self.createInv(fop, verbose, dosave)

#     def model(self):
#         """Return the actual model."""
#         return self.model

#     def createData(self, sensors, scheme):
#         """Create an empty data set."""
#         pass

#     def setData(self, data):
#         """Set data."""
#         self.dataContainer = data

#     def checkData(self):
#         """Check data validity."""
#         pass

#     @staticmethod
#     def estimateError(data, absoluteError=0.001, relativeError=0.001):
#         """Estimate error composed of an absolute and a relative part."""
#         pass

#     def showData(self, axes=None, response=None, name='data'):
#         """Show data."""
#         pass

#     # Work related methods
#     def invert(self, **kwargs):
#         """Invert the data and fill the parametrization."""
#         raise BaseException('implement me in derived class' + str(**kwargs))

#     def simulate(self, **kwargs):
#         """Run a simulation aka the forward task."""
#         raise BaseException('implement me in derived class' + str(**kwargs))

#     # Visualization stuff
#     def show(self, data, values=None, axes=None,
#              cMin=None, cMax=None, colorBar=1, **kwargs):
#         """Forward the visualization."""
#         pass

#     def showResult(self, ax=None, cMin=None, cMax=None, logScale=False,
#                    **kwargs):
#         """Show resulting vector."""
#         pass

#     def saveResult(self, folder=None, size=(16, 10),
#                    **kwargs):
#         """Save results in the specified folder."""
#         pass

#     @staticmethod
#     def createArgParser(dataSuffix='dat'):
#         """Create default argument parser.

#         Create default argument parser for the following options:

#             -Q, --quiet

#             -R, --robustData: options.robustData

#             -B, --blockyModel: options.blockyModel

#             -l, --lambda: options.lam

#             -i, --maxIter: options.maxIter

#             --depth: options.depth


#         """
#         import argparse

#         parser = argparse.ArgumentParser(
#             description="usage: %prog [options] *." + dataSuffix)
#         parser.add_argument("-Q", "--quiet", dest="quiet",
#                             action="store_true", default=False,
#                             help="Be verbose.")
#         parser.add_argument("-R", "--robustData", dest="robustData",
#                             action="store_true", default=False,
#                             help="Robust data (L1 norm) minimization.")
#         parser.add_argument("-B", "--blockyModel", dest="blockyModel",
#                             action="store_true", default=False,
#                             help="Blocky model (L1 norm) regularization.")
#         parser.add_argument('-l', "--lambda", dest="lam", type=float,
#                             default=100,
#                             help="Regularization strength.")
#         parser.add_argument('-i', "--maxIter", dest="maxIter", type=int,
#                             default=20,
#                             help="Maximum iteration count.")
#         parser.add_argument("--depth", dest="depth", type=float,
#                             default=None,
#                             help="Depth of inversion domain. [None=auto].")
#         parser.add_argument('dataFileName')
#         return parser


# class MeshMethodManager0(MethodManager0):
#     """Method Manager base class for managers using a (non-1D) mesh."""
#     def __init__(self, **kwargs):
#         """Constructor."""
#         # pg.critical('to be removed')
#         super(MeshMethodManager0, self).__init__(**kwargs)
#         self.mesh = None
#         self.data = None

#     # Mesh related methods
#     def createMesh(self, **kwargs):
#         """Create a mesh aka the parametrization."""
#         pass

#     def setMesh(self, mesh, refine=True):
#         """Set the internal mesh for this Manager.

#         Inject the mesh in the internal fop und inv.

#         Initialize RegionManager.
#         For more than two regions the first is assumed to be background.

#         Optional the forward mesh can be refined for higher numerical accuracy.

#         Parameters
#         ----------

#         DOCUMENTME!!!

#         """
#         if isinstance(mesh, str):
#             mesh = pg.load(mesh)

#         if self.verbose:
#             print(mesh)

#         self.mesh = pg.Mesh(mesh)
#         self.mesh.createNeighbourInfos()

#         self.fop.setMesh(self.mesh)
#         self.fop.regionManager().setConstraintType(1)
#         if self.fop.regionManager().regionCount() > 1:
#             self.fop.regionManager().region(1).setBackground(True)
# #            self.fop.regionManager().regions().begin().second.setBackground(1)

#         self.fop.createRefinedForwardMesh(refine)
#         self.paraDomain = self.fop.regionManager().paraDomain()
#         self.inv.setForwardOperator(self.fop)  # necessary? CR: check this

#     def showMesh(self, ax=None):
#         """Show mesh in given axes or in a new figure."""
#         raise Exception('IMPLEMENTME')
#         pass

#     def setData(self, data):
#         """Set data container from outside."""
#         if not isinstance(data, pg.DataContainer):
#             raise Exception('IMPLEMENTME')

#             if isinstance(data, str):
#                 raise Exception('IMPLEMENTME')
#         else:
#             self.dataContainer = data

#         if not self.dataContainer.allNonZero(self.dataToken()):
#             raise BaseException("No or partial values for:", self.dataToken())

#         if not self.dataContainer.allNonZero('err'):
#             raise BaseException("No or partial values for err")

#         self.fop.setData(self.dataContainer)

#     @staticmethod
#     def createArgParser(dataSuffix='dat'):
#         """Create argument parser for the manager."""
#         parser = MethodManager.createArgParser(dataSuffix)

#         parser.add_argument("--paraMaxCellSize", dest="maxCellArea",
#                             type=float, default=0.0,
#                             help="Maximum cell size for parametric cells "
#                                  "in m² [0.0] 0.0 means no restriction.")
#         parser.add_argument("--zWeight", dest="zWeight",
#                             type=float, default=0.7,
#                             help="Weight for vertical smoothness "
#                                  "(1=isotrope). [0.7]")
#         return parser

#     def invert(self, data=None, vals=None, err=None, mesh=None, **kwargs):
#         """Run the full inversion.

#         The data and error needed to be set before.
#         The meshes will be created if necessary.

#         DOCUMENTME!!!

#         Parameters
#         ----------

#         lam : float [20]
#                 regularization parameter
#         zWeight : float [0.7]
#                 relative vertical weight
#         maxIter : int [20]
#                 maximum iteration number
#         robustdata : bool [False]
#                 robust data reweighting using an L1 scheme (IRLS reweighting)
#         blockymodel : bool [False]
#                 blocky model constraint using L1 reweighting roughness vector
#         startModelIsReference : bool [False]
#                 startmodel is the reference model for the inversion

#         forwarded to createMesh

#         depth

#         quality

#         paraDX

#         maxCellArea

#         """
#         if 'verbose' in kwargs:
#             self.setVerbose(kwargs.pop('verbose'))

#         if data is not None:
#             # setDataContainer would be better
#             self.setData(data)

#         if vals is not None:
#             self.dataContainer.set(self.dataToken(), vals)

#         if err is not None:
#             self.dataContainer.set('err', vals)

#         # check for data container here
#         dataVals = self.dataContainer(self.dataToken())
#         errVals = self.dataContainer('err')

#         if mesh is not None:
#             self.setMesh(mesh)

#         if self.mesh is None:
#             self.createMesh(depth=kwargs.pop('depth', None),
#                             quality=kwargs.pop('quality', 34.0),
#                             maxCellArea=kwargs.pop('maxCellArea', 0.0),
#                             paraDX=kwargs.pop('paraDX', 0.3))

#         self.inv.setData(dataVals)
#         self.inv.setError(errVals, not self.errIsAbsolute)

#         zWeight = kwargs.pop('zWeight', 0.7)
#         if 'zweight' in kwargs:
#             zWeight = kwargs.pop('zweight', 0.7)
#             print("zweight option will be removed soon. Please use zWeight.")

#         self.fop.regionManager().setZWeight(zWeight)

#         self.inv.setLambda(kwargs.pop('lam', 20))
#         self.inv.setMaxIter(kwargs.pop('maxIter', 20))
#         self.inv.setRobustData(kwargs.pop('robustData', False))
#         self.inv.setBlockyModel(kwargs.pop('blockyModel', False))
#         self.inv.setRecalcJacobian(kwargs.pop('recalcJacobian', True))

#         # TODO: ADD MORE KWARGS
#         pc = self.fop.regionManager().parameterCount()

#         startModel = kwargs.pop('startModel',
#                                 pg.Vector(pc, pg.math.median(dataVals)))

#         self.inv.setModel(startModel)
#         # self.fop.setStartModel(startModel)

#         if kwargs.pop('startModelIsReference', False):
#             self.inv.setReferenceModel(startModel)

#         # Run the inversion
#         if len(kwargs) > 0:
#             print("Keyword arguments unknown:")
#             print(kwargs)
#             print("Warning! There are unknown kwargs arguments.")

#         model = self.inv.run()
#         self.model = model(self.paraDomain.cellMarkers())

#         return self.model
