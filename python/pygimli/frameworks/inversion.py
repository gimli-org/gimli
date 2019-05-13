# -*- coding: utf-8 -*-
"""pyGIMLi - Inversion Frameworks.

These are basic inversion frameworks that usually needs a forward operator to run.
"""
import numpy as np
import pygimli as pg


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
    """
    def __init__(self, fop=None, inv=None, **kwargs):
        self._verbose = kwargs.pop('verbose', False)
        self._debug = kwargs.pop('debug', False)

        # If this class or its derived is a Framework the _inv holds another
        # Inversion which allows us ........
        # this will be probably removed in the future
        self.isFrameWork = False

        self._preStep = None
        self._postStep = None

        self._inv = None
        self._fop = None

        self.reset()
        
        if inv is not None:
            self._inv = inv
            self.isFrameWork = True
        else:
            self._inv = pg.Inversion(self._verbose, self._debug)

        self._dataTrans = pg.RTransLin()
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
        if self.isFrameWork:
            return self._inv.inv
        else:
            return self._inv
    @property
    def fop(self):
        return self._fop
    @fop.setter
    def fop(self, f):
        self.setForwardOperator(f)

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
            self._startModel = self.fop.createStartModel()

        if self._startModel is None:
            self._startModel = self.fop.createDefaultStartModel(self.dataVals)

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
        elif hasattr(model, '__iter__'):
            self._startModel = model

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

        if min(abs(self._dataVals)) < 1e-12:
            print(self._dataVals)
            pg.warn("Found zero data values. Setting them to a TOLERANCE value of 1e-12")
            pg.fixZero(self._dataVals, 1e-12)

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
    def maxIter(self):
        return self.inv.maxIter()
    @maxIter.setter
    def maxIter(self, v):
        if self.inv is not None:
            self.inv.setMaxIter(v)

    def echoStatus(self):
        self.inv.echoStatus()

    def setDeltaChiStop(self, it):
        self.inv.setDeltaPhiAbortPercent(it)

    def setForwardOperator(self, fop):
        self._fop = fop
        # we need to initialize the regionmanager by calling it once
        self._fop.regionManager()
        self._inv.setForwardOperator(fop)

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

        return pg.dot(dData, dData)

    def phiModel(self, model=None):
        """ """
        if model is None:
            model = self.model

        rough = self.inv.roughness(model)
        return pg.dot(rough, rough)

    def phi(self, model=None, response=None):
        """ """
        phiD = self.phiData(response)
        if self.inv.localRegularization():
            return phiD
        else:
            return phiD + self.phiModel(model) * self.inv.getLambda()

    def run(self, dataVals, errorVals, **kwargs):
        """Run inversion.

        The inversion will always starts from the starting model taken from
        the forward operator.
        If you want to run the inversion from a specified prior model,
        e.g., from a other run, set this model as starting model to the FOP
        (fop.setStartModel).
        Any self.inv.setModel() settings will be overwritten.

        Parameters
        ----------
        errorVals : iterable
            Relative error values. dv / v

        Other Parameters
        ----------------
        maxIter : int
            Overwrite class settings for maximal iterations number.
        """
        if self.isFrameWork:
            return self._inv.run(dataVals, errorVals, **kwargs)

        if self.fop is None:
            raise Exception("Need a valid forward operator for the inversion run.")

        maxIter = kwargs.pop('maxIter', self.maxIter)

        self.verbose = kwargs.pop('verbose', self.verbose)
        self.debug   = kwargs.pop('debug', self.debug)

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
            print("Data transformation:", self.dataTrans)
            print("Model transformation:", self.modelTrans)
            print("min/max (data): {0:.2f}, {1:.2f}".format(min(self._dataVals), max(self._dataVals)))
            print("min/max (error): {0:.3f}%, {1:.3f}%".format(100*min(self._errorVals), 100*max(self._errorVals)))
            print("min/max (start model): {0:.2f}, {1:.2f}".format(min(self.startModel), max(self.startModel)))

        ### To ensure reproducability of the run() call inv.start() will
        ### reset self.inv.model() to fop.startModel().
        self.fop.setStartModel(self.startModel)
        self.inv.setReferenceModel(self.startModel)

        self.inv.start()
        self.maxIter = maxIterTmp

        if showProgress:
            self.showProgress(showProgress)

        lastPhi = self.phi()
        chi2History = [self.chi2()]

        for i in range(1, maxIter):

            if self._preStep and callable(self._preStep):
                self._preStep(i, self.inv)

            if self.verbose:
                print("-" * 80)
                print("inv.iter", i + 1, "... ", end='')

            self.inv.oneStep()
            resp = self.inv.response()
            chi2 = self.inv.chi2()

            if showProgress:
                self.showProgress(showProgress)

            self.inv.setLambda(self.inv.getLambda() * self.inv.lambdaFactor())

            if self.inv.robustData():
                self.inv.robustWeighting()

            if self.inv.blockyModel():
                self.inv.constrainBlocky()

            chi2History.append(chi2)

            if self._postStep and callable(self._postStep):
                self._postStep(i, self.inv)

            phi = self.phi()
            dPhi = (1-lastPhi / phi) * 100.

            if self.verbose:
                print("chi² = {0} (dPhi = {1}%) lam: {2}".format(
                            round(chi2, 2), round(dPhi, 2), self.inv.getLambda()))

            if chi2 <= 1:
                print("\n")
                if self.verbose:
                    pg.boxprint("Abort criteria reached: chi² <= 1")
                break

            if abs(dPhi) < self.inv.deltaPhiAbortPercent():
                if self.verbose:
                    pg.boxprint("Abort criteria reached: dPhi = {0} (< {1}%)".format(
                                round(dPhi, 2), self.inv.deltaPhiAbortPercent()))
                break

            lastPhi = phi

        if len(kwargs.keys()) > 0:
            print("Warning! unused keyword arguments", kwargs)

        self.model = self.inv.model()
        return self.model

    def showProgress(self, style='all'):
        r"""Called if showProgress=True is set for the inversion run.

        TODO
            * think .. its a useful function but breaks a little
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
                other_ax.clear()

            self.fop.drawModel(ax, self.inv.model())
        else:
            # for other_ax in ax[0].figure.axes:
            #     other_ax.clear()
            for _ax in self.axs:
                _ax.clear()
                try:
                    pg.mplviewer.twin(_ax).clear()
                except:
                    pass

            self.fop.drawModel(ax[0], self.inv.model(),
                               #label='Model'
                               )

            self.fop.drawData(ax[1], self._dataVals, self._errorVals,
                              #label='Data'
                              )
            self.fop.drawData(ax[1], self.inv.response(),
                              #label='Response'
                              )

            ax[1].text(0.01, 1.005,
                    "iter: %d, rrms: %.2g, $\chi^2$: %.2g" %
                        (self.inv.iter(), self.inv.relrms(), self.inv.chi2()),
                        transform=ax[1].transAxes,
                        horizontalalignment='left',
                        verticalalignment='bottom')

            ax[1].figure.tight_layout()
        pg.plt.pause(0.05)


class MarquardtInversion(Inversion):
    """Marquardt scheme (local damping with decreasing regularization strength
    """
    def __init__(self, fop=None, **kwargs):
        super(MarquardtInversion, self).__init__(fop, **kwargs)
        self.inv.setLocalRegularization(True)
        self.inv.stopAtChi1(False)
        self.inv.setLambdaFactor(0.8)

    def run(self, data, error, **kwargs):
        r"""Parameters
        ----------
        **kwargs:
            Forwarded to the parent class.
            See: :py:mod:`pygimli.modelling.Inversion`
        """
        self.fop.regionManager().setConstraintType(0)
        self.fop.setRegionProperties('*', cType=0)

        self.model = super(MarquardtInversion, self).run(data, error, **kwargs)
        return self.model

class Block1DInversion(MarquardtInversion):
    """
    Attributes
    ----------
    nLayers : int

    """
    def __init__(self, fop=None, **kwargs):
        super(Block1DInversion, self).__init__(fop=fop, **kwargs)
        # attributes:
        # nLayers, layerLimits, fixLayers

        self._nLayers = 4

    def setForwardOperator(self, fop):
        if not isinstance(fop, pg.frameworks.Block1DModelling):
            pg.critical('Forward operator needs to be an instance of '
                        'pg.modelling.Block1DModelling but is of type:', fop)

        return super(Block1DInversion, self).setForwardOperator(fop)

    def setLayers(self, nLayers):
        """Set amount of Layers"""
        self._nLayers = nLayers
        self.fop.initModelSpace(nLayers=nLayers)

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
            self.fop.setRegionProperties(0, modelControl=1e6)
            if hasattr(fixLayers, '__iter__'):
                if len(fixLayers) != self.fop.nLayers:
                    print("fixLayers:", fixLayers)
                    raise Exception("fixlayers needs to have a length of nLayers-1=" + str(self.fop.nLayers-1))
                self.fop.setRegionProperties(0, startModel=fixLayers)

            # TODO DRY to self.fop.createStartModel
            self.fop.setStartModel(self.fop.regionManager().createStartModel())

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

    def run(self, dataVals, errVals,
            nLayers=None, fixLayers=None, layerLimits=None,
            **kwargs):
        r"""

        Parameters
        ----------
        nLayers : int [4]
            Number of layers.
        fixLayers : bool | [thicknesses]
            For fixLayers=None, preset or defaults are uses.
            See: :py:mod:`pygimli.modelling.Block1DInversion.fixLayers`
        layerLimits : [min, max]
            For layerLimits=None, preset or defaults are uses.
            Set minimum or maximum layer thickness.

        **kwargs:
            Forwarded to the parent class.
            See: :py:mod:`pygimli.modelling.MarquardtInversion`
        """
        if nLayers is None:
            nLayers = self._nLayers

        ## initialize model space if needed
        self.setLayers(nLayers)

        if layerLimits is not None:
            self.setLayerLimits(layerLimits)

        if fixLayers is not None:
            self.fixLayers(fixLayers)

        self.model = super(Block1DInversion, self).run(dataVals, errVals, **kwargs)
        return self.model


class MeshInversion(Inversion):
    """
    ** UNUSED ** TO BE REMOVED
    Attributes
    ----------

    zWeight

    """
    def __init__(self, fop=None, **kwargs):
        TO_BE_REMOVED
        super(MeshInversion, self).__init__(fop=fop, **kwargs)
        self._zWeight = 1.0

    def setForwardOperator(self, fop):
        if not isinstance(fop, pg.frameworks.MeshModelling):
            pg.critical('Forward operator needs to be an instance of '
                        'pg.modelling.MeshModelling but is of type:', fop)

        return super(MeshInversion, self).setForwardOperator(fop)

    def run(self, dataVals, errVals, mesh=None, zWeight=None, **kwargs):
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

        self.model = super(MeshInversion, self).run(dataVals, errVals, **kwargs)

        return self.model


class PetroInversion(Inversion):
    def __init__(self, petro, mgr=None, fop=None, **kwargs):
        self.mgr = mgr

        if self.mgr is not None:
            fop = self.mgr.createForwardOperator(**kwargs)

        if fop is not None:
            fop = pg.frameworks.PetroModelling(fop, petro)

        super(PetroInversion, self).__init__(fop=fop, **kwargs)

    def setForwardOperator(self, fop):
        if not isinstance(fop, pg.frameworks.PetroModelling):
            pg.critical('Forward operator needs to be an instance of '
                        'pg.modelling.PetroModelling but is of type:', fop)

        return super(PetroInversion, self).setForwardOperator(fop)

    def run(self, dataVals, errVals, **kwargs):
        """
        """
        if 'limits' in kwargs:
            limits = kwargs.pop('limits', [0., 1.])

            if len(self.fop.regionManager().regionIdxs()) > 1:
                pg.critical('implement')
            else:
                self.fop.setRegionProperties('*', limits=limits)

        return super(PetroInversion, self).run(dataVals, errVals, **kwargs)


class LCInversion(Inversion):
    """2D Laterally constrained inversion LCI framework.
    """
    def __init__(self, fop=None, **kwargs):

        if fop is not None:
            f = pg.frameworks.LCModelling(fop, **kwargs)

        super(LCInversion, self).__init__(f, **kwargs)
        self.dataTrans = pg.RTransLog()
        #self.setDeltaChiStop(0.1)

    def prepare(self, dataVals, errVals, nLayers=4, **kwargs):
        dataVec = pg.RVector()
        for d in dataVals:
            dataVec = pg.cat(dataVec, d)

        errVec = pg.RVector()
        for e in errVals:
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

    def run(self, dataVals, errVals, nLayers=4, **kwargs):
        lam = kwargs.pop('lam', 20)
        dataVec, errVec = self.prepare(dataVals, errVals, nLayers, **kwargs)
        print('#'*50)
        print(kwargs)
        print('#'*50)
        return super(LCInversion, self).run(dataVec, errVec, lam=lam, **kwargs)
