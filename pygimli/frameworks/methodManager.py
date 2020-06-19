#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Method Manager

Provide the end user interface for method (geophysical) dependent
modelling and inversion as well as data and model visualization.
"""

import numpy as np
import pygimli as pg

from pygimli.utils import prettyFloat as pf


def fit(funct, data, err=None, **kwargs):
    """Generic function fitter.

    Fit data to a given function.

    TODO
    ----
        * will not work for harmonic function, maybe add harmonic flag to forward it to harmfit
        * Dictionary support for funct to submit user data..

    Parameters
    ----------
    funct: callable
        Function with the first argmument as data space, e.g., x, t, f, Nr. ..
        Any following arguments are the parameters to be fit.
        Except if a verbose flag if used.
    data: iterable (float)
        Data values
    err: iterable (float) [None]
        Data error values in %/100. Default is 1% if None are given.

    Other Parameters
    ----------------
    *dataSpace*: iterable
        Keyword argument of the data space of len(data).
        The name need to fit the first argument of funct.

    Returns
    -------
    model: array
        Fitted model parameter.
    response: array
        Model response.

    Example
    -------
    >>> import pygimli as pg
    >>>
    >>> func = lambda t, a, b: a*np.exp(b*t)
    >>> t = np.linspace(1, 2, 20)
    >>> data = func(t, 1.1, 2.2)
    >>> model, response = pg.frameworks.fit(func, data, t=t)
    >>> print(pg.core.round(model, 1e-5))
    2 [1.1, 2.2]
    >>> _ = pg.plt.plot(t, data, 'o', label='data')
    >>> _ = pg.plt.plot(t, response, label='response')
    >>> _ = pg.plt.legend()
    """
    mgr = ParameterInversionManager(funct, **kwargs)
    model = mgr.invert(data, err, **kwargs)
    return model, mgr.fw.response


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

    def __init__(self, fop=None, fw=None, data=None, **kwargs):
        """Constructor."""
        self._fop = fop
        self._fw = fw
        # we hold our own copy of the data
        self._verbose = kwargs.pop('verbose', False)
        self._debug = kwargs.pop('debug', False)

        self.data = None
        if data is not None:
            if isinstance(data, str):
                self.load(data)
            else:
                self.data = data

        # The inversion framework
        self._initInversionFramework(verbose=self._verbose,
                                     debug=self._debug)

        # The forward operator is stored in self._fw
        self._initForwardOperator(verbose=self._verbose, **kwargs)

        # maybe obsolete
        self.figs = {}
        self.errIsAbsolute = False

    def __hash__(self):
        """Create a hash for Method Manager."""
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
        pg.critical("No forward operator defined, either give one or "
                    "overwrite in derived class")

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

    def load(self, fileName):
        """API, overwrite in derived classes."""
        pg.critical('API, overwrite in derived classes', fileName)

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
            absolute error in the unit of the data

        Returns
        -------
        err : array
            Returning array of size len(data)
        """
        if absError is not None:
            return absError + data * errLevel

        return np.ones(len(data)) * errLevel

    def simulate(self, model, **kwargs):
        # """Run a simulation aka the forward task."""

        ra = self.fop.response(par=model)

        noiseLevel = kwargs.pop('noiseLevel', 0.0)
        if noiseLevel > 0:
            err = self.estimateError(ra, errLevel=noiseLevel)
            ra *= 1. + pg.randn(ra.size(), seed=kwargs.pop('seed', None)) * err
            return ra, err

        return ra

    def checkData(self, data):
        """Overwrite for special checks to return data values"""
        # if self._dataToken == 'nan':
        #     pg.critical('self._dataToken nan, should be set in class', self)
        #     return data(self._dataToken)
        return data

    def _ensureData(self, data):
        """Check data validity"""
        if data is None:
            data = self.fw.dataVals

        vals = self.checkData(data)

        if vals is None:
            pg.critical("There are no data values.")

        if abs(min(vals)) < 1e-12:
            print(min(vals), max(vals))
            pg.critical("There are zero data values.")

        return vals

    def checkError(self, err, dataVals=None):
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

        vals = self.checkError(err, dataVals)

        if vals is None:
            pg.warn('No data array given, set Fallback set to 1%')
            vals = np.ones(len(dataVals)) * 0.01

        try:
            if min(vals) <= 0:
                pg.critical("All error values need to be larger then 0."
                            " either give and err argument or fill dataContainer "
                            " with a valid 'err' ", min(vals), max(vals))
        except Exception as e:
            pg.critical("Can't estimate data error")


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
        if data is not None:
            self.data = data
        else:
            data = self.data

        dataVals = self._ensureData(data)

        errVals = self._ensureError(err, dataVals)

        self.preRun(**kwargs)
        self.fw.run(dataVals, errVals, **kwargs)
        self.postRun(**kwargs)

        return self.fw.model

    def showModel(self, model, ax=None, **kwargs):
        """Show a model.

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

        ax, cBar = self.fop.drawModel(ax, model, **kwargs)
        return ax, cBar

    def showData(self, data=None, ax=None, **kwargs):
        """Show the data.

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

        Returns
        -------
        ax, cbar

        """
        if ax is None:
            fig, ax = pg.plt.subplots()

        if data is None:
            data = self.data

        return self.fop.drawData(ax, data, **kwargs), None

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
        ax, cBar = self.showData(data=self.inv.dataVals,
                                 error=self.inv.errorVals,
                                 label='Data',
                                 ax=ax, **kwargs)
        ax, cBar = self.showData(data=self.inv.response,
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

        return ax, cBar

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
        elif funct is not None:
            fop = pg.frameworks.ParameterModelling(funct)
        else:
            pg.critical('you should either give a valid fop or a function so I can create the fop for you')

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
        dataSpace = kwargs.pop(self.fop.dataSpaceName, None)

        if dataSpace is not None:
            self.fop.dataSpace = dataSpace

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
        """Constructor.

        Attribute
        ---------
        mesh: pg.Mesh
            Copy of the main Mesh. Will be distributet to inversion and the fop.
            You can overwrite it with invert(mesh=mesh).
        """
        super(MeshMethodManager, self).__init__(**kwargs)
        self.mesh = None

    @property
    def paraDomain(self):
        return self.fop.paraDomain

    def paraModel(self, model=None):
        """Give the model parameter regarding the parameter mesh."""
        if model is None:
            model = self.fw.model

        if len(model) == self.fw.parameterCount:
            return model
        else:
            self.fop.paraModel(model)

    def createMesh(self, data=None, **kwargs):
        """API, implement in derived classes."""
        pg.critical('no default mesh generation defined .. implement in '
                    'derived class')

    def applyMesh(self, mesh, ignoreRegionManager=False, **kwargs):
        """ """
        if ignoreRegionManager:
            mesh = self.fop.createRefinedFwdMesh(mesh, **kwargs)
        self.fop.setMesh(mesh, ignoreRegionManager=ignoreRegionManager)

    def applyData(self, data):
        """ """
        self.fop.data = data

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

        Keyword Arguments
        -----------------
        forwarded to Inversion.run

        Returns
        -------
        model : array
            Model mapped for match the paraDomain Cell markers.
            The calculated model is in self.fw.model.
        """
        if data is None:
            data = self.data

        if data is None:
            pg.critical('No data given for inversion')

        if mesh is None:
            mesh = self.mesh

        if mesh is None:
            mesh = self.createMesh(data, **kwargs)

        self.mesh = mesh

        self.applyData(data)
        self.applyMesh(mesh)

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
                # Check! why I can't use mgr.fop
                #fop = mgr.fop
                fop = mgr.createForwardOperator()
                self.checkData = mgr.checkData
                self.checkError = mgr.checkError

            if fop is not None:
                if not isinstance(fop, pg.frameworks.PetroModelling):
                    petrofop = pg.frameworks.PetroModelling(fop, petro)

        if petrofop is None:
            print(mgr)
            print(fop)
            pg.critical('implement me')

        super().__init__(fop=petrofop, **kwargs)


class JointPetroInversionManager(MeshMethodManager):
    def __init__(self, petros, mgrs):
        """Initialize with lists of managers and transformations"""
        self.mgrs = mgrs

        self.fops = [pg.frameworks.PetroModelling(m.fop, p)
                     for p, m in zip(petros, mgrs)]

        super().__init__(fop=pg.frameworks.JointModelling(self.fops))

        ## just hold a local copy
        self.dataTrans = pg.trans.TransCumulative()

    def checkError(self, err, data=None):
        """Collect error values."""
        if len(err) != len(self.mgrs):
            pg.critical("Please provide data for all managers")

        vals = pg.Vector(0)
        for i, mgr in enumerate(self.mgrs):
            # we get the data values again or we have to split data
            dataVals = mgr.checkData(self.fop._data[i])
            vals = pg.cat(vals, mgr.checkError(err[i], dataVals))
        return vals

    def checkData(self, data):
        """Collect data values."""
        if len(data) != len(self.mgrs):
            pg.critical("Please provide data for all managers")

        self.dataTrans.clear()
        vals = pg.Vector(0)

        for i, mgr in enumerate(self.mgrs):
            self.dataTrans.add(mgr.inv.dataTrans, data[i].size())
            vals = pg.cat(vals, mgr.checkData(data[i]))

        self.inv.dataTrans = self.dataTrans
        return vals

    def invert(self, data, **kwargs):
        """Run inversion"""
        limits = kwargs.pop('limits', [0., 1.])
        self.fop.modelTrans.setLowerBound(limits[0])
        self.fop.modelTrans.setUpperBound(limits[1])

        kwargs['startModel'] = kwargs.pop('startModel',
                                          (limits[1]+limits[0])/2.)

        return super().invert(data, **kwargs)
