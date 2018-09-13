# -*- coding: utf-8 -*-
"""
Vertical electrical sounding (VES) manager class.
"""
import numpy as np
import matplotlib.pyplot as plt

import pygimli as pg
from pygimli.mplviewer import drawModel1D

from pygimli.frameworks import Modelling, Block1DModelling

from pygimli.manager import MethodManager1d

def hasTwin(ax):
    for other_ax in ax.figure.axes:
        if other_ax is ax:
            continue
        if other_ax.bbox.bounds == ax.bbox.bounds:
            return other_ax
    return None


class VESModelling(Block1DModelling):
    """Vertical Electrical Sounding (VES) forward operator.

    Attributes
    ----------
    am :
        Part of data basis. Distances between A and M electrodes.
        A is first power, M is first potential electrode.
    bm :
        Part of data basis. Distances between B and M electrodes.
        B is second power, M is first potential electrode.
    an :
        Part of data basis. Distances between A and N electrodes.
        A is first power, M is second potential electrode.
    bn :
        Part of data basis. Distances between B and N electrodes.
        B is second power, M is second potential electrode.
    ab2 :
        Half distance between A and B.
        Only used for output and auto generated.
    """
    def __init__(self, ab2=None, mn2=None, **kwargs):

        self.am = None
        self.bm = None
        self.an = None
        self.bn = None
        self.ab2 = None
        self.mn2 = None

        super(VESModelling, self).__init__(**kwargs)

        self.setDataBasis(ab2=ab2, mn2=mn2)

    def createStartModel(self, rhoa, nLayers):
        self.setLayers(nLayers)

        #print(self.mn2)
        #print(self.ab2)
        startThicks = np.logspace(np.log10(min(self.mn2)/2),
                                  np.log10(max(self.ab2)/5), nLayers-1)
        startThicks = pg.utils.diff(pg.cat([0.0], startThicks))

        #print(startThicks)
        #exit()

        # layer thickness properties
        self.setRegionProperties(0, startModel=startThicks, trans='log')

        # resistivity properties
        self.setRegionProperties(1, startModel=np.median(rhoa), trans='log')

        # find a better way on forced update this seems  to be in user space here
        sm = self.regionManager().createStartModel()
        self.setStartModel(sm)
        return sm

    def setDataBasis(self, ab2=None, mn2=None,
                     am=None, bm=None, an=None, bn=None,
                     **kwargs):
        """Set data basis, i.e., arrays for all am, an, bm, bn distances.

        Parameters
        ----------
        """
        if ab2 is not None and mn2 is not None:

            if isinstance(mn2, float):
                mn2 = np.ones(len(ab2))*mn2

            if len(ab2) != len(mn2):
                print("ab2", ab2)
                print("mn2", mn2)
                raise Exception("length of ab2 is unequal length of nm2")

            self.am = ab2 - mn2
            self.an = ab2 + mn2
            self.bm = ab2 + mn2
            self.bn = ab2 - mn2

        elif am is not None \
            and bm is not None \
            and an is not None \
            and bn is not None:
            self.am = am
            self.bm = bm
            self.an = an
            self.bn = bn

        if self.am is not None and self.bm is not None:
            self.ab2 = (self.am + self.bm) / 2
            self.mn2 = pg.abs((self.am - self.an)) / 2

            self.k = (2.0 * np.pi) / (1.0/self.am - 1.0/self.an -
                                    1.0/self.bm + 1.0/self.bn)

    def response(self, par):
        return self.response_mt(par, 0)

    def response_mt(self, par, i=0):

        if self.am is not None and self.bm is not None:
            nLayers = (len(par)+1) // 2
            fop = pg.DC1dModelling(nLayers, self.am, self.bm, self.an, self.bn)
        else:
            raise Exception("I have no data basis .. "
                            "don't know what to calculate.")

        return fop.response(par)

    def drawModel(self, ax, model, **kwargs):
        pg.mplviewer.drawModel1D(ax=ax,
                                 model=model,
                                 plot=kwargs.pop('plot', 'loglog'),
                                 xlabel='Resistivity [$\Omega$m]', **kwargs)

    def drawData(self, ax, data, error=None, label=None, **kwargs):
        """
        """
        ra = data
        raE = error

        col = kwargs.pop('color', 'green')
        if label == 'Response':
            col = 'blue'

        ab2 = kwargs.pop('ab2', self.ab2)
        mn2 = kwargs.pop('mn2', self.mn2)
        plot = kwargs.pop('plot', 'loglog')

        plot = getattr(ax, plot)
        plot(ra, ab2, 'x-', color=col, **kwargs)

        if raE is not None:
            ax.errorbar(ra, ab2,
                        xerr=ra * raE, elinewidth=2, barsabove=True,
                        linewidth=0, color='red')

        ax.set_ylim(max(ab2), min(ab2))
        ax.set_xlabel('Apparent resistivity [$\Omega$m]')
        ax.set_ylabel('AB/2 in [m]')
        ax.grid(True)


class VESCModelling(VESModelling):
    def __init__(self, **kwargs):
        super(VESCModelling, self).__init__(nBlocks=2, **kwargs)
        self.phiAxe = None

    def phaseModel(self, model):
        """Return the current phase model values."""
        nLay = (len(model) + 1) // 3
        return pg.cat(model[0:nLay-1], 1000. * model[nLay*2-1::])

    def resModel(self, model):
        """Return the resistivity model values."""
        nLay = (len(model) + 1) // 3
        return model[0:nLay*2-1]

    def createStartModel(self, rhoa, nLayers):
        self.setLayers(nLayers)

        startThicks = np.logspace(np.log10(min(self.mn2)/2),
                                  np.log10(max(self.ab2)/5), nLayers-1)
        startThicks = pg.utils.diff(pg.cat([0.0], startThicks))

        # layer thickness properties
        self.setRegionProperties(0, startModel=startThicks, trans='log')

        # resistivity properties
        self.setRegionProperties(1, startModel=np.median(rhoa),
                                 trans='log')

        self.setRegionProperties(2, startModel=np.median(rhoa[len(rhoa)//2::]),
                                 trans='log')

        sm = self.regionManager().createStartModel()
        self.setStartModel(sm)
        return sm

    def response_mt(self, par, i=0):
        """ Multithread response for parametrization.

            Returns [|rhoa|, +phi(rad)] for [thicks, res, phi(rad)]
        """

        if self.am is not None and self.bm is not None:
            nLayers = (len(par) + 1) // 3
            fop = pg.DC1dModellingC(nLayers, self.am, self.bm, self.an, self.bn)
        else:
            raise Exception("I have no data basis .. "
                            "don't know what to calculate.")

        return fop.response(par)

    def drawModel(self, ax, model, **kwargs):
        """Draw 1D VESC Modell."""
        a1 = None
        a2 = None

        if err is not None:
            if isinstance(err, float):
                err = np.ones(len(data))*err

            super(VESCModelling, self).drawData(ax, data[0:len(data)//2], err[0:len(data)//2],
                             label=label)
            paE = err[len(data)//2::]
        else:
            a1 = ax
            a2 = hasTwin(ax)
            if a2 is None:
                a2 = ax.twiny()

        super(VESCModelling, self).drawModel(a1, self.resModel(model), **kwargs)

        plot = kwargs.pop('plot', 'semilogy')
        if plot is 'loglog':
            plot = 'semilogy'
        elif plot is 'semilogx':
            plot = 'plot'

        pg.mplviewer.drawModel1D(ax=a2,
                                 model=self.phaseModel(model),
                                 plot=plot,
                                 color='green',
                                 xlabel='Phase [mrad]', **kwargs)

    def drawData(self, ax, data, error=None, label=None, ab2=None, mn2=None,
                 **kwargs):
        """Draw 1D VESC Data."""
        a1 = None
        a2 = None

        if hasattr(ax, '__iter__'):
            if len(ax) == 2:
                a1 = ax[0]
                a2 = ax[1]
        else:
            a1 = ax
            a2 = hasTwin(ax)
            if a2 is None:
                a2 = ax.twiny()

        if ab2 is not None and mn2 is not None:
            self.setDataBasis(ab2=ab2, mn2=mn2)

        ra = data[0:len(data)//2]
        phi = data[len(data)//2::] * 1000. #mRad

        phiE = None # abs err
        raE = None # rel err
        if error is not None:
            if type(error) is float:
                raE = np.ones(len(data)//2) * error
                phiE = np.ones(len(data)//2) * error
            else:
                raE = error[0:len(data)//2]
                phiE = error[len(data)//2::]

        super(VESCModelling, self).drawData(a1, ra, error=raE,
                                            color='black',
                                            label=label, **kwargs)

        a2.semilogy(phi, self.ab2, 'x-', color='green', **kwargs)

        if phiE is not None:
            a2.errorbar(phi, self.ab2,
                        xerr=phiE, elinewidth=2, barsabove=True,
                        linewidth=0, color='red')


        a2.set_ylim(max(self.ab2), min(self.ab2))
        a2.set_xlabel('Apparent phase [mRad]', color='green')
        a2.set_ylabel('AB/2 in [m]')
        a2.grid(True)


class VESManager(MethodManager1d):
    """Vertical electrical sounding (VES) manager class.

    Examples
    --------
    >>> # no need to import matplotlib. pygimli's show does
    >>> import numpy as np
    >>> import pygimli as pg
    >>> from pygimli.physics import VESManager
    >>> ab2 = np.logspace(np.log10(1.5), np.log10(100), 32)
    >>> mn2 = 1.0
    >>> # 3 layer with 100, 500 and 20 Ohmm
    >>> # and layer thickness of 4, 6, 10 m
    >>> # over a Halfspace of 800 Ohmm
    >>> synthModel = pg.cat([4., 6., 10.], [100., 5., 20., 800.])
    >>> VES = VESManager()
    >>> ra, err = VES.simulate(synthModel, ab2=ab2, mn2=mn2, noiseLevel=0.01)
    >>> ax = VES.showData(ra, err)
    >>> # _= VES.invert(ra, err, nLayer=4, showProgress=0, verbose=0)
    >>> # ax = VES.showModel(synthModel)
    >>> # ax = VES.showResult(ax=ax)
    >>> pg.wait()
    """
    def __init__(self, **kwargs):
        """Constructor

        Parameters
        ----------

        complex : bool
            Accept complex resistivities.

        Attributes
        ----------
        complex : bool
            Accept complex resistivities.
        """
        self.__complex = kwargs.pop('complex', False)

        super(VESManager, self).__init__(**kwargs)

        self.transData = None
        self.transRho = pg.TransLog()
        self.transPhi = pg.TransLin()
        #self.transPhi = pg.TransLog()

    @property
    def complex(self):
        return self.__complex

    @complex.setter
    def complex(self, c):
        self.__complex = c
        self.initForwardOperator()

    def createForwardOperator(self, **kwargs):
        """Create Forward Operator.

        Create Forward Operator based on complex attribute.
        """
        if self.complex:
            return VESCModelling(**kwargs)
        else:
            return VESModelling(**kwargs)

    def simulate(self, model, ab2=None, mn2=None, **kwargs):
        """Simulate measurement data.
        """
        if ab2 is not None and mn2 is not None:
            self._fw.fop.setDataBasis(ab2=ab2, mn2=mn2)

        return super(VESManager, self).simulate(model, **kwargs)

    def invert(self, data=None, err=None, ab2=None, mn2=None, **kwargs):
        """Invert measured data.

        Parameters
        ----------

        Returns
        -------
        """
        if ab2 is not None and mn2 is not None:
            self.fop.setDataBasis(ab2=ab2, mn2=mn2)

        if isinstance(data, pg.DataContainer):
            raise StandardException('Implement me')

        nData = len(data)//2
        if self.complex:
            self.transData = pg.TransCumulative()
            self.transData.add(self.transRho, nData)
            self.transData.add(self.transPhi, nData)
        else:
            self.transData = pg.TransLog()

        self.inv.transData = self.transData

        self.fop.setRegionProperties(0, limits=[1., 1000.])
        self.fop.setRegionProperties(1, limits=[5., 10000.])
        self.fop.setRegionProperties(2, limits=[0.1/1000, 6./1000.])

        errVals = err
        if len(err) == 2 and self.complex:
            errVals = pg.cat(np.ones(nData)*err[0],
                             np.abs(err[1] / data[nData::]))
            #assume [relative resistivity Err, absolute phase Err]

        #ensure data and error sizes here

        return super(VESManager, self).invert(dataVals=data, errVals=errVals,
                                              **kwargs)

    def loadData(self, fileName, **kwargs):
        """ Load simple data matrix
        """
        mat = np.loadtxt(fileName)
        if len(mat[0]) == 4:
            self.fop.setDataBasis(ab2=mat[:,0], mn2=mat[:,1])
            return mat.T
        if len(mat[0]) == 6:
            self.complex = True
            self.fop.setDataBasis(ab2=mat[:,0], mn2=mat[:,1])
            return mat[:,0], mat[:,1], np.array(pg.cat(mat[:,2], mat[:,4])), np.array(pg.cat(mat[:,3], mat[:,5]))

    def exportData(self, fileName, data=None, error=None):
        """Export data into simple ascii matrix.

        Usefull?
        """
        mn2 = np.abs((self.fop.am - self.fop.an) / 2.)
        ab2 = (self.fop.am + self.fop.bm) / 2.
        mat = None
        if data is None:
            data = self.inv.dataVals

        if error is None:
            error = self.inv.errorVals

        if self.complex:
            nData = len(data)//2
            mat = np.array([ab2, mn2,
                            data[:nData], error[:nData],
                            data[nData:], error[nData:]
                            ]).T
            np.savetxt(fileName, mat, header='ab/2\tmn/2\trhoa\terr\tphia\terrphi')
        else:
            mat = np.array([ab2, mn2, data, error]).T
            np.savetxt(fileName, mat, header='ab/2\tmn/2\trhoa\terr')


def test_VESManager(showProgress=False):
    """
        run from console with: python -c 'import pygimli.physics.ert.ves as pg; pg.test_VESManager(1)'
    """
    thicks = [2., 10.]
    res = [100., 5., 30]
    phi = [0., 20., 0.]

    # model fails
    thicks = [2., 6., 10.]
    res = [100., 500., 20., 800.]
    phi = [0., 20., 50., 0]

    synthModel = pg.cat(thicks, res)
    ab2 = np.logspace(np.log10(1.5), np.log10(100.), 25)

    mgr = VESManager(verbose=True, debug=False)
    #mgr.fop.setRegionProperties(0, limits=[0.5, 200], trans='log')
    #mgr.fop.setRegionProperties(1, limits=[3, 2000], trans='log')
    ra, err = mgr.simulate(synthModel, ab2=ab2, mn2=1.0, noiseLevel=0.01)
    mgr.exportData('synth.ves', ra, err)

    mgr.invert(ra, err, nLayers=4, lam=100,
               showProgress=showProgress)

    ### Test -- reinit with new parameter count
    mgr.invert(ra, err, nLayers=3,
               showProgress=showProgress)

    #np.testing.assert_array_less(mgr.inv.inv.chi2(), 1)

    ### Test -- reinit with new data basis
    ab2 = np.logspace(np.log10(1.5), np.log10(50.), 10)
    ra, err = mgr.simulate(synthModel, ab2=ab2, mn2=1.0, noiseLevel=0.01)

    mgr2 = VESManager(verbose=False, debug=False)
    mgr2.invert(ra, err, nLayers=3, ab2=ab2, mn2=1.0,
                showProgress=showProgress)

    #np.testing.assert_array_less(mgr2.inv.inv.chi2(), 1)

    ### Test -- reinit with complex resistivies
    mgr.complex = True
    synthModel =  pg.cat(synthModel, phi)

    ra, err = mgr.simulate(synthModel, ab2=ab2, mn2=1.0, noiseLevel=0.01)
    mgr.exportData('synthc.ves', ra, err)
    mgr.invert(ra, err,
               showProgress=showProgress)

    np.testing.assert_array_less(mgr.inv.inv.chi2(), 1)

    if showProgress:
        print("test done");
        pg.wait()


def VESManagerApp():
    """Call VESManager as console app"""

    parser = VESManager.createArgParser(dataSuffix='ves')
    options = parser.parse_args()

    verbose = not options.quiet
    if verbose:
        print("VES Manager console application.")
        print(options._get_kwargs())

    mgr = VESManager(verbose=verbose, debug=pg.debug())

    ab2, mn2, ra, err = mgr.loadData(options.dataFileName)

    mgr.showData(ra, err)
    mgr.invert(ra, err, ab2, mn2,
               maxIter=options.maxIter,
               lam=options.lam,
               )
    mgr.showResultAndFit()
    pg.wait()
































class __VESManager():  # Should be derived from 1DManager
    """Vertical electrical sounding (VES) manager class."""
    def __init__(self,
                 ab2,
                 z,
                 mn2=None,
                 Type='smooth',
                 verbose=False):
        """
        Parameters
        ----------

        ab2: array_like
            Vector of distances between the point of the sounding and the
            current electrodes.

        mn2: array_like [ab2/3]
            Vector of distances between the point of the sounding and the
            potential electrodes.

        z: array_like, case specific
            smooth: z discretisation [m]\n
            block: number of layers
        """
        self.verbose = verbose
        self.type = Type

        self.ab2 = ab2

        self.mn2 = mn2
        if mn2 is None:
            self.mn2 = ab2/3

        self.Z = z  # z discretisation or nlay

        self.FOP = None
        self.INV = None
        self.startmodel = None

        self.createModelTrans()  # creates default as fallback

    def createFOP(self):
        """Creates the forward operator instance."""
        if self.type == 'block':
            self.FOP = pg.DC1dModelling(self.Z, self.ab2, self.mn2)
            mesh1D = pg.createMesh1DBlock(self.Z)

        if self.type == 'smooth':
            self.FOP = pg.DC1dRhoModelling(self.Z, self.ab2, self.mn2)
            mesh1D = pg.createMesh1D(nCells=len(self.Z) + 1)

        self.FOP.setMesh(mesh1D)
        self.applyModelTrans()

    @staticmethod
    def simulate(synmodel, ab2=None, mn2=None, errPerc=3.):
        """Forward calculation with optional noise

        Simulates a synthetic data set of a vertical electric sounding and
        appends gaussian distributed noise.
        Block only for now.

        Parameters
        ----------

        ab2: array_like
            Vector of distances between the point of the sounding and the
            current electrodes.

        mn2: array_like [ab2/3]
            Vector of distances between the point of the sounding and the
            potential electrodes.

        errPerc: float [3.]
            Percentage Value for the gaussian noise. Default are 3 %.

        """
        thk = synmodel[0]
        res = synmodel[1]
        if mn2 is None:
            mn2 = ab2/3
        FOP = pg.DC1dModelling(len(res), ab2, mn2)
        syndata = FOP.response(thk + res)
        syndata = syndata * (pg.randn(len(syndata)) * errPerc / 100. + 1.)
        return syndata

    def createINV(self, data, relErrorP=3., startmodel=None, **kwargs):
        """Create inversion instance

        Parameters
        ----------

        data : array_like
            Data array you would like to fit with this inverse modelling
            approach.

        relErrorP : float [3.]
            Percentage value of the relative error you assume. Default 3. means
            a 3 % error is assumed. Affects the chi2 criteria during the
            inversion process and therefore the inversion result (Inversion
            tries to fit the given data within the given errors).

        startmodel : array_like [None]
            Optional possibility to define the starting model for the inversion
            routine. The default will be the mean of the given data.

        **kwargs : keyword arguments
        ----------------------------

        Keyword arguments are redirected to the block inversion instance only!

        lambdaFactor: float < 1.0 [0.8]
            Inversion in Marquardt scheme reduces the lambda from initial high
            values down to a certain minimum. The reduction per step is
            represented by this value. Default is a reduction to 80% of the
            previous step (By the way, the default start lambda is 1000).

        robust: boolean [False]
            Recalculation of the errors to reduce the weight of spikes in the
            data. Not necessary for synthetic or "good" field data.
        """
        if self.FOP is None:
            self.createFOP()

        self.INV = pg.RInversion(data,
                                 self.FOP,
                                 False)

        self.applyDataTrans()

        if self.type == 'block':
            print(kwargs.pop('lambdaFactor'))
            self.INV.setMarquardtScheme(kwargs.pop('lambdaFactor', 0.8))
            self.INV.setRobustData(kwargs.pop('robust', False))

        if self.type == 'smooth':
            pass

        self.INV.setRelativeError(relErrorP / 100.0)

        if startmodel is not None:
            self.startmodel = startmodel

        if self.startmodel is None:
            self.createStartModel(data)

        if self.type == 'smooth':
            self.INV.setModel(self.startmodel)

        else:
            self.FOP.region(0).setStartValue(self.startmodel[0])
            self.FOP.region(1).setStartValue(self.startmodel[1])

    def createModelTrans(self,
                         thkBounds=(10., 0.1, 30.),
                         modelBounds=(10, 1.0, 10000.),
                         trans=('log', 'log')):
        """Define model transformations for inversion."""
        self.thkBounds = thkBounds  # thk(start, min, max)
        self.rhoBounds = modelBounds  # rho (model) (start, min, max)
        self.trans = trans
        if self.FOP is not None:
            self.applyModelTrans()

    def applyModelTrans(self):
        """Pass previously given bounds for the model transformation to FOP."""
        if self.FOP is None:
            raise Exception('initialize forward operator before appending \
                             proper boundaries for the model transformation')

        if self.type == 'smooth':
            self.FOP.region(0).setParameters(self.rhoBounds[0],
                                             self.rhoBounds[1],
                                             self.rhoBounds[2], self.trans[1])
        elif self.type == 'block':

            self.FOP.region(0).setParameters(self.thkBounds[0],
                                             self.thkBounds[1],
                                             self.thkBounds[2], self.trans[0])

            self.FOP.region(1).setParameters(self.rhoBounds[0],
                                             self.rhoBounds[1],
                                             self.rhoBounds[2], self.trans[1])

    def applyDataTrans(self):
        """Apply a logarithmic transformation to the data."""
        self.dataTrans = pg.RTransLog()
        self.INV.setTransData(self.dataTrans)

    def createStartModel(self, data):
        """Creates a default start model for the inversion."""
        if self.type == 'smooth':
            self.startmodel = pg.RVector(len(self.Z) + 1, np.median(data))
        else:
            self.startmodel = [self.getDepth() / self.Z / 2, np.median(data)]

    def invert(self, data, startmodel=None, relErrorP=3., lam=None, **kwargs):
        """Run inversion

        Creates forward operator, initializes inversion scheme and run
        inversion based on the input parameters. kwargs are redirected to
        the createINV method.
        """
        self.dataToFit = data
        if self.INV is None:
            self.createINV(data,
                           relErrorP=relErrorP,
                           startmodel=startmodel,
                           **kwargs)

        if self.type == 'smooth' and lam is None:
            self.resDistribution = self.INV.runChi1()

        else:
            if lam is None:
                lam = 1000.0
            self.INV.setLambda(lam)
            self.resDistribution = self.INV.run()

        if self.verbose is False:
            self.INV.echoStatus()
        return self.resDistribution

    def splitBlockModel(self):
        """Returns the thickness and ressitivities of the model."""
        z = int(self.Z) - 1
        return self.resDistribution[:z], self.resDistribution[z:]

    def getDepth(self):
        """Rule-of-thumb for Wenner/Schlumberger."""
        return np.max(self.ab2) / 3.  # rule-of-thumb for Wenner/Schlumberger

    def showResults(self, ax=None, syn=None, color=None):
        """Shows the Results of the inversion."""
        if ax is None:
            fig, ax = plt.subplots(ncols=1, figsize=(8, 6))
        if syn is not None:
            drawModel1D(ax, syn[0], syn[1], color='b', label='synthetic',
                        plotfunction='semilogx')
        if color is None:
            color = 'g'
        if self.type == 'smooth':
            drawModel1D(ax, self.Z, self.resDistribution, color=color,
                        label=r'$\lambda$={:2.2f}'
                        .format(self.INV.getLambda()))
        else:
            thk, rho = self.splitBlockModel()
            drawModel1D(ax, thk, rho, color=color,
                        label=r'$\ chi^2$={:2.2f}'
                        .format(self.INV.getChi2()))
        ax.grid(True, which='both')
        ax.legend(loc='best')
        return ax

    def showFit(self, ax=None, color='g', marker='-', syn=True):
        """Visualizes the data fit."""
        if syn is True:
            ax.loglog(self.dataToFit, self.ab2, 'bx-',
                      label='measured/synthetic')
        ax.loglog(self.INV.response(), self.ab2, color, label='fitted')
        ax.set_ylim((max(self.ab2), min(self.ab2)))
        ax.grid(True, which='both')
        ax.set_xlabel(r'$\rho_a$ [$\Omega$m]')
        ax.set_ylabel('AB/2 [m]')
        ax.legend(loc='best')
        return ax

    def showResultsAndFit(self, syn=None):
        """Calls showResults and showFit."""
        fig, ax = plt.subplots(ncols=2, figsize=(8, 6))
        self.showResults(ax=ax[0], syn=syn)
        self.showFit(ax=ax[1])

if __name__ == '__main__':
    VESManagerApp()


# The End
