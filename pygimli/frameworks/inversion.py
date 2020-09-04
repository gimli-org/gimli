# -*- coding: utf-8 -*-
"""pyGIMLi - Inversion Frameworks.

These are basic inversion frameworks that usually needs a forward operator to run.
"""
import numpy as np
import pygimli as pg

from pygimli.utils import prettyFloat as pf


class Inversion(object):
    """Basic inversion framework.

    Changes to prior Versions (remove me)

        * holds the starting model itself, fop only provide a creator for SM
        fop.createStartModel(dataValues)

    Attributes
    ----------
    verbose : bool
        Give verbose output
    debug : bool
        Give debug output
    startModel : array
        Holds the current starting model
    model : array
        Holds the last active model
    maxIter : int [20]
        Maximal interation number.
    stopAtChi1 : bool [True]
        Stop iteration when chi² is one. If set to false the iteration stops
        after maxIter or convergence reached (self.inv.deltaPhiAbortPercent())
    """
    def __init__(self, fop=None, inv=None, **kwargs):
        self._verbose = kwargs.pop('verbose', False)
        self._debug = kwargs.pop('debug', False)

        # If this class or its derived is a Framework the _inv holds another
        # Inversion which allows us (remove me)........
        # this will be probably removed in the future
        self.isFrameWork = False # check if needed
        self._stopAtChi1 = True

        self._preStep = None
        self._postStep = None

        self._inv = None
        self._fop = None

        self.reset()

        if inv is not None:
            self._inv = inv
            self.isFrameWork = True
        else:
            self._inv = pg.core.RInversion(self._verbose, self._debug)

        self._dataTrans = pg.trans.TransLin()
        self.axs = None # for showProgress only
        self.maxIter = kwargs.pop('maxIter', 20)

        if fop is not None:
            self.setForwardOperator(fop)

    def reset(self):
        """"""
        self._model = None
        self._startModel = None
        self._dataVals = None
        self._errorVals = None

    @property
    def inv(self):
        return self._inv
    @property
    def fop(self):
        return self._fop
    @fop.setter
    def fop(self, f):
        self.setForwardOperator(f)

    def setForwardOperator(self, fop):
        self._fop = fop
        # we need to initialize the regionmanager by calling it once
        self._fop.regionManager()
        self._inv.setForwardOperator(fop)

    @property
    def verbose(self):
        return self._verbose
    @verbose.setter
    def verbose(self, v):
        self._verbose = v
        self.inv.setVerbose(self._verbose)
        self.fop.setVerbose(self._verbose)

    @property
    def debug(self):
        return self._debug
    @debug.setter
    def debug(self, v):
        self._debug = v
        self.inv.setDoSave(self._debug)

    @property
    def dataTrans(self):
        return self._dataTrans
    @dataTrans.setter
    def dataTrans(self, dt):
        self._dataTrans = dt
        self.inv.setTransData(self._dataTrans)

    @property
    def modelTrans(self):
        return self.fop.modelTrans

    @property
    def startModel(self):
        """ Gives the current default starting model.

        Returns the current default starting model or
        call fop.createStartmodel() if non is defined.
        """
        if self._startModel is None:
            sm = self.fop.regionManager().createStartModel()
            if len(sm) > 0 and max(abs(np.atleast_1d(sm))) > 0.0:
                self._startModel = sm
                pg.info("Created startmodel from region infos:", sm)
            else:
                pg.verbose("No region infos for startmodel")

        if self._startModel is None:
            sm = self.fop.createStartModel(self.dataVals)
            pg.info("Created startmodel from forward operator:", sm)
            self._startModel = sm
        return self._startModel

    @startModel.setter
    def startModel(self, model):
        """
        model: [float] | float
            Model used as starting model.
            Float value is used as constant model.
        """
        if model is None:
            self._startModel = None
        elif isinstance(model, float) or isinstance(model, int):
            self._startModel = np.ones(self.parameterCount) * float(model)
            pg.info("Startmodel set from given value.", float(model))
        elif hasattr(model, '__iter__'):
            if len(model) == self.parameterCount:
                pg.info("Startmodel set from given array.", model)
                self._startModel = model
            else:
                pg.error("Startmodel size invalid {0} != {0}.".
                         format(len(model), self.parameterCount))

    @property
    def model(self):
        """The last active model."""
        if self._model is None:
            if hasattr(self.inv, 'model'):
                ### inv is RInversion()
                if len(self.inv.model()) > 0:
                    return self.inv.model()
                else:
                    raise pg.critical("There was no inversion run so there is no last model")
            else:
                return self.inv.model
        return self._model

    @model.setter
    def model(self, m):
        self._model = m

    @property
    def response(self):
        if len(self.inv.response()) > 0:
            return self.inv.response()
        else:
            raise Exception("There was no inversion run so there is no response yet")

    # backward compatibility
    @property
    def dataErrs(self):
        pg.warn('do not use')
        return self._errorVals
    @dataErrs.setter
    def dataErrs(self, v):
        pg.warn('do not use')
        self._errorVals = v

    @property
    def dataVals(self):
        return self._dataVals
    @dataVals.setter
    def dataVals(self, d):
        """Set mandatory data values.

        Values == 0.0. Will be set to Tolerance
        """
        self._dataVals = d

        if self._dataVals is None:
            pg._y(d)
            pg.critical("Inversion framework needs data values to run")

        # zero can be a valid data value
        #
        # if min(abs(self._dataVals)) < 1e-12:
        #     print(self._dataVals)
        #     pg.warn("Found zero data values. Setting them to a TOLERANCE value of 1e-12")
        #     pg.fixZero(self._dataVals, 1e-12)

    @property
    def errorVals(self):
        return self._errorVals
    @errorVals.setter
    def errorVals(self, d):
        """Set mandatory error values.

        Values == 0.0. Will be set to Tolerance
        """
        self._errorVals = d

        if self._errorVals is None:
            pg._y(d)
            pg.critical("Inversion framework needs error values to run")

        if min(abs(self._errorVals)) < 1e-12:
            print(self._errorVals)
            pg.warn("Found zero error values. Setting them to a Fallback value of 1")
            pg.fixZero(self._errorVals, 1)

    @property
    def parameterCount(self):
        pC = self.fop.regionManager().parameterCount()
        if pC == 0:
            pg.warn("Parameter count is 0")
        return pC

    @property
    def robustData(self):
        return self.inv.robustData()
    @robustData.setter
    def robustData(self, v):
        if self.inv is not None:
            self.inv.setRobustData(v)

    @property
    def maxIter(self):
        return self.inv.maxIter()
    @maxIter.setter
    def maxIter(self, v):
        if self.inv is not None:
            self.inv.setMaxIter(v)

    @property
    def stopAtChi1(self):
        return self._stopAtChi1
    @stopAtChi1.setter
    def stopAtChi1(self, b):
        self._stopAtChi1 = b

    @property
    def minDPhi(self):
        return self.inv.deltaPhiAbortPercent()
    @minDPhi.setter
    def minDPhi(self, dPhi):
        return self.setDeltaChiStop(dPhi)
    def setDeltaChiStop(self, it):
        self.inv.setDeltaPhiAbortPercent(it)

    def echoStatus(self):
        self.inv.echoStatus()

    def setPostStep(self, p):
        self._postStep = p

    def setPreStep(self, p):
        self._preStep = p

    def setData(self, data):
        QUESTION_ISNEEDED
        if isinstance(data, pg.DataContainer):
            raise Exception("should not been here .. its Managers job")
            self.fop.setData(data)
        else:
            self.dataVals = data

    def chi2(self, response=None):
        return self.phiData(response) / len(self.dataVals)

    def phiData(self, response=None):
        """ """
        if response is None:
            response = self.response

        dT = self.dataTrans
        dData = (dT.trans(self.dataVals) - dT.trans(response)) / \
                 dT.error(self.dataVals, self.errorVals)

        return pg.math.dot(dData, dData)

    def phiModel(self, model=None):
        """ """
        if model is None:
            model = self.model

        rough = self.inv.roughness(model)
        return pg.math.dot(rough, rough)

    def phi(self, model=None, response=None):
        """ """
        phiD = self.phiData(response)
        if self.inv.localRegularization():
            return phiD
        else:
            return phiD + self.phiModel(model) * self.inv.getLambda()

    def relrms(self):
        """Relative root-mean-square misfit of the last run."""
        return self.inv.relrms()

    def absrms(self):
        """Absolute root-mean-square misfit of the last run."""
        return self.inv.absrms()

    def run(self, dataVals, errorVals, **kwargs):
        """Run inversion.

        The inversion will always start from the starting model taken from
        the forward operator.
        If you want to run the inversion from a specified prior model,
        e.g., from a other run, set this model as starting model to the FOP
        (fop.setStartModel).
        Any self.inv.setModel() settings will be overwritten.

        Parameters
        ----------
        dataVals : iterable
            Data values
        errorVals : iterable
            Relative error values. dv / v

        Keyword Arguments
        -----------------
        maxIter : int
            Overwrite class settings for maximal iterations number.
        dPhi : float [1]
            Overwrite class settings for delta data phi aborting criteria.
            Default is 1%
        """
        self.reset()
        if self.isFrameWork:
            pg.critical('in use?')
            return self._inv.run(dataVals, errorVals, **kwargs)

        if self.fop is None:
            raise Exception("Need a valid forward operator for the inversion run.")

        maxIter = kwargs.pop('maxIter', self.maxIter)
        minDPhi = kwargs.pop('dPhi', self.minDPhi)

        self.verbose = kwargs.pop('verbose', self.verbose)
        self.debug   = kwargs.pop('debug', self.debug)
        self.robustData = kwargs.pop('robustData', False)

        lam = kwargs.pop('lam', 20)

        progress = kwargs.pop('progress', None)
        showProgress = kwargs.pop('showProgress', False)

        self.inv.setTransModel(self.fop.modelTrans)

        self.dataVals = dataVals
        self.errorVals = errorVals

        sm = kwargs.pop('startModel', None)
        if sm is not None:
            self.startModel = sm

        self.inv.setData(self._dataVals)
        self.inv.setRelativeError(self._errorVals)
        self.inv.setLambda(lam)

        # temporary set max iter to one for the initial run call
        maxIterTmp = self.maxIter
        self.maxIter = 1

        if self.verbose:
            pg.info('Starting inversion.')
            print("fop:", self.inv.fop())
            if isinstance(self.dataTrans, pg.trans.TransCumulative):
                print("Data transformation (cumulative):")
                for i in range(self.dataTrans.size()):
                    print("\t", i, self.dataTrans.at(i))
            else:
                print("Data transformation:", self.dataTrans)
            if isinstance(self.modelTrans, pg.trans.TransCumulative):
                print("Model transformation (cumulative):")
                for i in range(self.modelTrans.size()):
                    if i < 10:
                        print("\t", i, self.modelTrans.at(i))
                    else:
                        print(".", end='')
            else:
                print("Model transformation:", self.modelTrans)

            print("min/max (data): {0}/{1}".format(pf(min(self._dataVals)),
                                                    pf(max(self._dataVals))))
            print("min/max (error): {0}%/{1}%".format(pf(100*min(self._errorVals)),
                                                      pf(100*max(self._errorVals))))
            print("min/max (start model): {0}/{1}".format(pf(min(self.startModel)),
                                                          pf(max(self.startModel))))

        ### To ensure reproduceability of the run() call inv.start() will
        ### reset self.inv.model() to fop.startModel().
        self.fop.setStartModel(self.startModel)
        self.inv.setReferenceModel(self.startModel)

        if self.verbose:
            print("-" * 80)
        if self._preStep and callable(self._preStep):
                self._preStep(0, self)

        self.inv.start()
        self.maxIter = maxIterTmp

        if self._postStep and callable(self._postStep):
            self._postStep(0, self)

        if showProgress:
            self.showProgress(showProgress)

        lastPhi = self.phi()
        self.chi2History = [self.chi2()]
        self.modelHistory = [self.startModel]

        for i in range(1, maxIter):

            if self._preStep and callable(self._preStep):
                self._preStep(i, self)

            if self.verbose:
                print("-" * 80)
                print("inv.iter", i + 1, "... ", end='')

            try:
                self.inv.oneStep()
            except RuntimeError as e:
                print(e)
                pg.error('One step failed. '
                         'Aborting and going back to last model')

            if np.isnan(self.model).any():
                print(model)
                pg.critical('invalid model')

            resp = self.inv.response()
            chi2 = self.inv.chi2()

            self.chi2History.append(chi2)
            self.modelHistory.append(self.model)

            if showProgress:
                self.showProgress(showProgress)

            if self._postStep and callable(self._postStep):
                self._postStep(i, self)

            ### we need to check  the following before oder after chi2 calc??
            self.inv.setLambda(self.inv.getLambda() * self.inv.lambdaFactor())

            if self.robustData:
                self.inv.robustWeighting()

            if self.inv.blockyModel():
                self.inv.constrainBlocky()

            phi = self.phi()
            dPhi = phi/lastPhi

            if self.verbose:
                print("chi² = {0} (dPhi = {1}%) lam: {2}".format(
                            round(chi2, 2), round((1-dPhi)*100, 2), self.inv.getLambda()))

            if chi2 <= 1 and self.stopAtChi1:
                print("\n")
                if self.verbose:
                    pg.boxprint("Abort criterion reached: chi² <= 1 (%.2f)" % chi2)
                break

            if (dPhi > (1.0 - minDPhi / 100.0)) and i > 2:
            # if dPhi < -minDPhi:
                if self.verbose:
                    pg.boxprint("Abort criteria reached: dPhi = {0} (< {1}%)".format(
                                round((1-dPhi)*100, 2), minDPhi))
                break

            lastPhi = phi

        ### will never work as expected until we unpack kwargs .. any idea for
        # better strategy?
        # if len(kwargs.keys()) > 0:
        #     print("Warning! unused keyword arguments", kwargs)

        self.model = self.inv.model()
        return self.model

    def showProgress(self, style='all'):
        r"""Called if showProgress=True is set for the inversion run.

        TODO
            *Discuss .. its a useful function but breaks a little
                the FrameWork work only concept.
        """
        if self.axs is None:
            axs = None
            if style == 'all' or style == True:
                fig, axs = pg.plt.subplots(1, 2)
            elif style == 'Model':
                fig, axs = pg.plt.subplots(1, 1)
            self.axs = axs
        ax = self.axs

        if style == 'Model':
            for other_ax in ax.figure.axes:
                # pg._y(type(other_ax).mro())
                if type(other_ax).mro()[0] == type(ax):
                    # only clear Axes not Colorbars
                    other_ax.clear()

            self.fop.drawModel(ax, self.inv.model())
        else:
            # for other_ax in ax[0].figure.axes:
            #     other_ax.clear()
            for _ax in self.axs:
                _ax.clear()
                try:
                    pg.viewer.mpl.twin(_ax).clear()
                except:
                    pass

            self.fop.drawModel(ax[0], self.inv.model(),
                               label='Model')
            self.fop.drawData(ax[1], self._dataVals, self._errorVals,
                              label='Data')
            self.fop.drawData(ax[1], self.inv.response(),
                              label='Response')

            ax[1].text(0.99, 0.005,
                    "Iter: {0}, rrms: {1}%, $\chi^2$: {2}"
                        .format(self.inv.iter(),
                                pf(self.inv.relrms()),
                                pf(self.inv.chi2())),
                        transform=ax[1].transAxes,
                        horizontalalignment='right',
                        verticalalignment='bottom',
                        fontsize=8)

            ax[1].figure.tight_layout()
        pg.plt.pause(0.05)


class MarquardtInversion(Inversion):
    """Marquardt scheme (local damping with decreasing regularization strength
    """
    def __init__(self, fop=None, **kwargs):
        super(MarquardtInversion, self).__init__(fop, **kwargs)
        self.stopAtChi1 = False
        self.inv.setLocalRegularization(True)
        self.inv.setLambdaFactor(0.8)

    def run(self, dataVals, errorVals, **kwargs):
        r"""Parameters
        ----------
        **kwargs:
            Forwarded to the parent class.
            See: :py:mod:`pygimli.modelling.Inversion`
        """
        self.fop.regionManager().setConstraintType(0)
        self.fop.setRegionProperties('*', cType=0)

        self.model = super(MarquardtInversion, self).run(dataVals, errorVals, **kwargs)
        return self.model

class Block1DInversion(MarquardtInversion):
    """
    Attributes
    ----------
    nLayers : int

    """
    def __init__(self, fop=None, **kwargs):
        #pg.warn("move this to the manager")
        super(Block1DInversion, self).__init__(fop=fop, **kwargs)

    def setForwardOperator(self, fop):
        if not isinstance(fop, pg.frameworks.Block1DModelling):
            pg.critical('Forward operator needs to be an instance of '
                        'pg.modelling.Block1DModelling but is of type:', fop)

        return super(Block1DInversion, self).setForwardOperator(fop)

    def fixLayers(self, fixLayers):
        """Fix layer thicknesses.

        Parameters
        ----------
        fixLayers : bool | [float]
            Fix all layers to the last value or set the fix layer
            thickness for all layers
        """
        if fixLayers is False:
            self.fop.setRegionProperties(0, modelControl=1.0)
        elif fixLayers is not None:
            # how do we fix values without modelControl?
            # maybe set the region to be fixed here
            self.fop.setRegionProperties(0, modelControl=1e6)
            if hasattr(fixLayers, '__iter__'):
                if len(fixLayers) != self.fop.nLayers:
                    print("fixLayers:", fixLayers)
                    pg.error("fixlayers needs to have a length of nLayers-1="
                             + str(self.fop.nLayers-1))
                self.fop.setRegionProperties(0, startModel=fixLayers)

    def setLayerLimits(self, limits):
        """Set min and max layer thickness.

        Parameters
        ----------
        limits : False | [min, max]
        """
        if limits is False:
            self.fop.setRegionProperties(0, limits=[0.0, 0.0], trans='log')
        else:
            self.fop.setRegionProperties(0, limits=limits, trans='log')

    def setParaLimits(self, limits):
        """Set the limits for each parameter region."""
        for i in range(1, 1 + self.fop.nPara):
            if self.fop.nPara == 1:
                self.fop.setRegionProperties(i, limits=limits, trans='log')
            else:
                self.fop.setRegionProperties(i, limits=limits[i-1], trans='log')

    def run(self, dataVals, errorVals,
            nLayers=None, fixLayers=None, layerLimits=None, paraLimits=None,
            **kwargs):
        r"""

        Parameters
        ----------
        nLayers : int [4]
            Number of layers.
        fixLayers : bool | [thicknesses]
            See: :py:mod:`pygimli.modelling.Block1DInversion.fixLayers`
            For fixLayers=None, preset or defaults are uses.
        layerLimits : [min, max]
            Limits the thickness off all layers.
            For layerLimits=None, preset or defaults are uses.
        paraLimits : [min, max] | [[min, max],...]
            Limits the range of the model parameter. If you have multiple
            parameters you can set them with a list of limits.

        **kwargs:
            Forwarded to the parent class.
            See: :py:mod:`pygimli.modelling.MarquardtInversion`
        """
        if nLayers is not None:
            self.fop.nLayers = nLayers

        if layerLimits is not None:
            self.setLayerLimits(layerLimits)

        if fixLayers is not None:
            self.fixLayers(fixLayers)

        if paraLimits is not None:
            self.setParaLimits(paraLimits)

        self.model = super(Block1DInversion, self).run(dataVals, errorVals, **kwargs)
        return self.model


class MeshInversion(Inversion):
    """
    ** UNUSED ** TO BE REMOVED
    Attributes
    ----------

    zWeight

    """
    def __init__(self, fop=None, **kwargs):
        pg.critical('Obsolete .. to be removed.')
        super(MeshInversion, self).__init__(fop=fop, **kwargs)
        self._zWeight = 1.0

    def setForwardOperator(self, fop):
        if not isinstance(fop, pg.frameworks.MeshModelling):
            pg.critical('Forward operator needs to be an instance of '
                        'pg.modelling.MeshModelling but is of type:', fop)

        return super(MeshInversion, self).setForwardOperator(fop)

    def run(self, dataVals, errorVals, mesh=None, zWeight=None, **kwargs):
        """
        """
        if mesh is not None:
            self.fop.setMesh(mesh)

        # maybe move this to the fop
        if zWeight is None:
            zWeight = self._zWeight
        self.fop.setRegionProperties('*', zWeight=zWeight)
        # maybe move this to the fop

        pg.debug('run with: ', self.fop.regionProperties())
        #### more mesh related inversion attributes to set?

        # ensure the mesh is generated
        self.fop.mesh()

        self.model = super(MeshInversion, self).run(dataVals, errorVals, **kwargs)

        return self.model


class PetroInversion(Inversion):
    def __init__(self, petro, fop=None, **kwargs):
        """
        Parameters
        ----------
        """
        pg.critical('Obsolete .. to be removed.')
        if fop is not None:
            if not isinstance(fop, pg.frameworks.PetroModelling):
                fop = pg.frameworks.PetroModelling(fop, petro)

        super(PetroInversion, self).__init__(fop=fop, **kwargs)

    def setForwardOperator(self, fop):
        if not isinstance(fop, pg.frameworks.PetroModelling):
            pg.critical('Forward operator needs to be an instance of '
                        'pg.modelling.PetroModelling but is of type:', fop)

        return super(PetroInversion, self).setForwardOperator(fop)

    def run(self, dataVals, errorVals, **kwargs):
        """
        """
        if 'limits' in kwargs:
            limits = kwargs.pop('limits')

            if len(self.fop.regionManager().regionIdxs()) > 1:
                pg.critical('implement')
            else:
                self.fop.setRegionProperties('*', limits=limits)

        #ensure the mesh
        self.fop.mesh()
        return super(PetroInversion, self).run(dataVals, errorVals, **kwargs)


class LCInversion(Inversion):
    """2D Laterally constrained inversion LCI framework.
    """
    def __init__(self, fop=None, **kwargs):

        if fop is not None:
            f = pg.frameworks.LCModelling(fop, **kwargs)

        super(LCInversion, self).__init__(f, **kwargs)
        self.dataTrans = pg.trans.TransLog()
        #self.setDeltaChiStop(0.1)

    def prepare(self, dataVals, errorVals, nLayers=4, **kwargs):
        dataVec = pg.RVector()
        for d in dataVals:
            dataVec = pg.cat(dataVec, d)

        errVec = pg.RVector()
        for e in errorVals:
            errVec = pg.cat(errVec, e)

        self.fop.initJacobian(dataVals=dataVals, nLayers=nLayers,
                              nPar=kwargs.pop('nPar', None))

        ### self.fop.initJacobian resets prior set startmodels
        if self._startModel is not None:
            self.fop.setStartModel(self._startModel)

        rC = self.fop.regionManager().regionCount()

        if kwargs.pop('disableLCI', False):
            self.inv.setMarquardtScheme(0.7)
            #self.inv.setLocalRegularization(True)
            for r in self.fop.regionManager().regionIdxs():
                self.fop.setRegionProperties(r, cType=0)
        else:
            #self.inv.stopAtChi1(False)
            cType = kwargs.pop('cType', None)
            if cType is None:
                cType = [1] * rC

            zWeight = kwargs.pop('zWeight', None)
            if zWeight is None:
                zWeight = [0.0] * rC

            self.fop.setRegionProperties('*',
                                         cType=cType,
                                         zWeight=zWeight,
                                         **kwargs)
            self.inv.setReferenceModel(self.fop.startModel())

        return dataVec, errVec

    def run(self, dataVals, errorVals, nLayers=4, **kwargs):
        lam = kwargs.pop('lam', 20)
        dataVec, errVec = self.prepare(dataVals, errorVals, nLayers, **kwargs)
        print('#'*50)
        print(kwargs)
        print('#'*50)
        return super(LCInversion, self).run(dataVec, errVec, lam=lam, **kwargs)
