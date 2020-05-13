# -*- coding: utf-8 -*-
"""
Vertical electrical sounding (VES) manager class.
"""
import numpy as np

import pygimli as pg
# from pygimli.frameworks import Modelling, Block1DModelling
from pygimli.frameworks import Block1DModelling, MethodManager1d


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
        A is first power, N is second potential electrode.
    bn :
        Part of data basis. Distances between B and N electrodes.
        B is second power, N is second potential electrode.
    ab2 :
        Half distance between A and B.
    mn2 :
        Half distance between A and B.
        Only used for input (feeding am etc.).
    """
    def __init__(self, ab2=None, mn2=None, **kwargs):
        r"""Constructor
        """
        self.am = None
        self.bm = None
        self.an = None
        self.bn = None
        self.ab2 = None
        self.mn2 = None

        super(VESModelling, self).__init__(**kwargs)

        if 'dataContainerERT' in kwargs:
            data = kwargs['dataContainerERT']
            if isinstance(data, pg.DataContainerERT):
                kwargs['am'] = [data.sensorPosition(data('a')[i]).distance(
                        data('m')[i]) for i in range(data.size())]
                kwargs['an'] = [data.sensorPosition(data('a')[i]).distance(
                        data('n')[i]) for i in range(data.size())]
                kwargs['bm'] = [data.sensorPosition(data('b')[i]).distance(
                        data('m')[i]) for i in range(data.size())]
                kwargs['bn'] = [data.sensorPosition(data('b')[i]).distance(
                        data('n')[i]) for i in range(data.size())]

        self.setDataSpace(ab2=ab2, mn2=mn2, **kwargs)

    def createStartModel(self, rhoa):
        r"""
        """
        if self.nLayers == 0:
            pg.critical("Model space is not been initialized.")

        startThicks = np.logspace(np.log10(min(self.mn2)/2),
                                  np.log10(max(self.ab2)/5),
                                  self.nLayers - 1)
        startThicks = pg.utils.diff(pg.cat([0.0], startThicks))

        # layer thickness properties
        self.setRegionProperties(0, startModel=startThicks, trans='log')

        # resistivity properties
        self.setRegionProperties(1, startModel=np.median(rhoa), trans='log')

        return super(VESModelling, self).createStartModel()

    def setDataSpace(self, ab2=None, mn2=None,
                     am=None, bm=None, an=None, bn=None,
                     **kwargs):
        """Set data basis, i.e., arrays for all am, an, bm, bn distances.

        Parameters
        ----------
        """
        # Sometimes you don't have AB2/MN2 but provide am etc.
        self.am = am
        self.an = an
        self.bm = bm
        self.bn = bn
        if ab2 is not None and mn2 is not None:  # overrides am etc.
            if isinstance(mn2, float):
                mn2 = np.ones(len(ab2))*mn2

            if len(ab2) != len(mn2):
                print("ab2", ab2)
                print("mn2", mn2)
                raise Exception("length of ab2 is unequal length of mn2")

            self.am = ab2 - mn2
            self.an = ab2 + mn2
            self.bm = ab2 + mn2
            self.bn = ab2 - mn2

        elif (am is not None and bm is not None and an is not None
              and bn is not None):
            self.am = am
            self.bm = bm
            self.an = an
            self.bn = bn

        if self.am is not None and self.bm is not None:
            self.ab2 = (self.am + self.bm) / 2
            self.mn2 = abs(self.am - self.an) / 2

            self.k = (2.0 * np.pi) / (1.0 / self.am - 1.0 / self.an -
                                      1.0 / self.bm + 1.0 / self.bn)

    def response(self, par):
        return self.response_mt(par, 0)

    def response_mt(self, par, i=0):

        if self.am is not None and self.bm is not None:
            nLayers = (len(par)+1) // 2
            fop = pg.core.DC1dModelling(nLayers,
                                        self.am, self.bm, self.an, self.bn)
        else:
            pg.critical("No data space defined don't know what to calculate.")

        return fop.response(par)

    def drawModel(self, ax, model, **kwargs):
        pg.viewer.mpl.drawModel1D(ax=ax,
                                  model=model,
                                  plot=kwargs.pop('plot', 'loglog'),
                                  xlabel=r'Resistivity ($\Omega$m)', **kwargs)
        ax.set_ylabel('Depth in (m)')

    def drawData(self, ax, data, error=None, label=None, **kwargs):
        r"""Draw modeled apparent resistivity data.

        Parameters
        ----------
        ax: axes
            Matplotlib axes object to draw into.

        data: iterable
            Apparent resistivity values to draw.

        error: iterable [None]
            Adds an error bar if you have error values.

        label: str ['$\varrho_a$']
            Set legend label for the amplitude.

        Other parameters
        ----------------
        ab2: iterable
            Override ab2 that fits data size.
        mn2: iterable
            Override mn2 that fits data size.
        plot: function name
            Matplotlib plot function, e.g., plot, loglog, semilogx or semilogy
        """
        ab2 = kwargs.pop('ab2', self.ab2)
        # mn2 = kwargs.pop('mn2', self.mn2)
        plot = kwargs.pop('plot', 'loglog')

        ra = data
        raE = error

        style = dict(pg.frameworks.modelling.DEFAULT_STYLES.get(
            label, pg.frameworks.modelling.DEFAULT_STYLES['Default']))
        style.update(kwargs)
        a1 = ax
        plot = getattr(a1, plot)
        if label is None:
            label = r'$\varrho_a$'
        plot(ra, ab2, 'x-', label=label, **style)

        if raE is not None:
            a1.errorbar(ra, ab2,
                        xerr=ra * raE, barsabove=True,
                        **pg.frameworks.modelling.DEFAULT_STYLES.get('Error',
                            pg.frameworks.modelling.DEFAULT_STYLES['Default']),
                        label='_nolegend_')

        a1.set_ylim(max(ab2), min(ab2))
        a1.set_xlabel(r'Apparent resistivity ($\Omega$m)')
        a1.set_ylabel(r'AB/2 (m)')
        a1.grid(True)
        a1.legend()


class VESCModelling(VESModelling):
    """Vertical Electrical Sounding (VES) forward operator. (complex)

    Vertical Electrical Sounding (VES) forward operator for complex
    resistivity values. see: :py:mod:`pygimli.physics.ert.VESModelling`
    """
    def __init__(self, **kwargs):
        super(VESCModelling, self).__init__(nPara=2, **kwargs)
        self.phiAxe = None

    def phaseModel(self, model):
        """Return the current phase model values."""
        nLay = (len(model) + 1) // 3
        return pg.cat(model[0:nLay-1], 1000. * model[nLay*2-1::])

    def resModel(self, model):
        """Return the resistivity model values."""
        nLay = (len(model) + 1) // 3
        return model[0:nLay*2-1]

    def createStartModel(self, rhoa):
        startThicks = np.logspace(np.log10(min(self.mn2)/2),
                                  np.log10(max(self.ab2)/5),
                                  self._nLayers-1)
        startThicks = pg.utils.diff(pg.cat([0.0], startThicks))

        # layer thickness properties
        self.setRegionProperties(0, startModel=startThicks,
                                 trans='log')

        # resistivity properties
        self.setRegionProperties(1, startModel=np.median(rhoa),
                                 trans='log')

        self.setRegionProperties(2, startModel=np.median(rhoa[len(rhoa)//2::]),
                                 trans='log')

        sm = self.regionManager().createStartModel()
        return sm

    def response_mt(self, par, i=0):
        """ Multithread response for parametrization.

            Returns [|rhoa|, +phi(rad)] for [thicks, res, phi(rad)]
        """

        if self.am is not None and self.bm is not None:
            nLayers = (len(par) + 1) // 3
            fop = pg.core.DC1dModellingC(nLayers,
                                         self.am, self.bm, self.an, self.bn)
        else:
            pg.critical("No data basis known.")

        return fop.response(par)

    def drawModel(self, ax, model, **kwargs):
        """Draw 1D VESC Modell."""
        a1 = ax
        a2 = pg.viewer.mpl.createTwinY(ax)

        super(VESCModelling, self).drawModel(a1,
                                             model=self.resModel(model),
                                             **kwargs)

        plot = kwargs.pop('plot', 'semilogy')
        if plot == 'loglog':
            plot = 'semilogy'
        elif plot == 'semilogx':
            plot = 'plot'

        pg.viewer.mpl.drawModel1D(ax=a2,
                                 model=self.phaseModel(model),
                                 plot=plot,
                                 color='C2',
                                 xlabel='Phase (mrad)',
                                 **kwargs)

        a2.set_xlabel('neg. phase (mRad)', color='C2')

    def drawData(self, ax, data, error=None, labels=None, ab2=None, mn2=None,
                 **kwargs):
        r"""Draw modeled apparent resistivity and apparent phase data.

        Parameters
        ----------
        ax: axes
            Matplotlib axes object to draw into.

        data: iterable
            Apparent resistivity values to draw. [rhoa phia].

        error: iterable [None]
            Rhoa in Ohm m and phia in radiand.
            Adds an error bar if you have error values. [err_rhoas err_phia]
            The error of amplitudes are assumed to be relative and the error
            of the phases is assumed to be absolute in mrad.

        labels: str [r'$\varrho_a$', r'$\varphi_a$']
            Set legend labels for amplitude and phase.

        Other parameters:
        -----------------
        ab2: iterable
            Override ab2 that fits data size.
        mn2: iterable
            Override mn2 that fits data size.
        plot: function name
            Matplotlib plot function, e.g., plot, loglog, semilogx or semilogy
        """
        a1 = None
        a2 = None

        if hasattr(ax, '__iter__'):
            if len(ax) == 2:
                a1 = ax[0]
                a2 = ax[1]
        else:
            a1 = ax
            a2 = pg.viewer.mpl.createTwinY(ax)

        if ab2 is not None and mn2 is not None:
            self.setDataSpace(ab2=ab2, mn2=mn2)

        ra = data[0:len(data)//2]
        phi = data[len(data)//2::] * 1000.  # mRad

        phiE = None  # abs err
        raE = None  # rel err

        if error is not None:
            if type(error) is float:
                raE = np.ones(len(data)//2) * error
                phiE = np.ones(len(data)//2) * error
            else:
                raE = error[0:len(data)//2]
                phiE = error[len(data)//2::]

        if labels is None:
            labels = [r'$\varrho_a$', r'$\varphi_a$']

        label = kwargs.pop('label', 'Data')

        style = dict(pg.frameworks.modelling.DEFAULT_STYLES.get(
            label, pg.frameworks.modelling.DEFAULT_STYLES['Default']))
        style.update(kwargs)

        super(VESCModelling, self).drawData(a1, ra, error=raE,
                                            label=labels[0], **style)

        style['Color'] = 'C2'

        a2.semilogy(phi, self.ab2, label=labels[1], **style)

        if phiE is not None:
            a2.errorbar(phi, self.ab2,
                        xerr=phiE,
                        **pg.frameworks.modelling.DEFAULT_STYLES.get('Error',
                              pg.frameworks.modelling.DEFAULT_STYLES['Default']),
                        barsabove=True,
                        label='_nolegend_'
                        )

        a2.set_ylim(max(self.ab2), min(self.ab2))
        a2.set_xlabel('Apparent neg. phase (mRad)', color='C2')
        a2.set_ylabel('AB/2 in (m)')
        a2.legend()
        a2.grid(True)


class VESManager(MethodManager1d):
    r"""Vertical electrical sounding (VES) manager class.

    Examples
    --------

    >>> import numpy as np
    >>> import pygimli as pg
    >>> from pygimli.physics import VESManager
    >>> ab2 = np.logspace(np.log10(1.5), np.log10(100), 32)
    >>> mn2 = 1.0
    >>> # 3 layer with 100, 500 and 20 Ohmm
    >>> # and layer thickness of 4, 6, 10 m
    >>> # over a Halfspace of 800 Ohmm
    >>> synthModel = pg.cat([4., 6., 10.], [100., 5., 20., 800.])
    >>> ves = VESManager()
    >>> ra, err = ves.simulate(synthModel, ab2=ab2, mn2=mn2, noiseLevel=0.01)
    >>> ax = ves.showData(ra, error=err)
    >>> # _= ves.invert(ra, err, nLayer=4, showProgress=0, verbose=0)
    >>> # ax = ves.showModel(synthModel)
    >>> # ax = ves.showResult(ax=ax)
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
        self._complex = kwargs.pop('complex', False)

        super(VESManager, self).__init__(**kwargs)

        self.inv.setDeltaChiStop(1)

        self.dataTrans = None
        self.rhoaTrans = pg.trans.TransLog()
        self.phiaTrans = pg.trans.TransLin()

    @property
    def complex(self):
        return self._complex

    @complex.setter
    def complex(self, c):
        self._complex = c
        self.reinitForwardOperator()

    def createForwardOperator(self, **kwargs):
        """Create Forward Operator.

        Create Forward Operator based on complex attribute.
        """
        if self.complex:
            return VESCModelling(**kwargs)
        else:
            return VESModelling(**kwargs)

    def simulate(self, model, ab2=None, mn2=None, **kwargs):
        """Simulate measurement data."""
        if ab2 is not None and mn2 is not None:
            self._fw.fop.setDataSpace(ab2=ab2, mn2=mn2)

        return super(VESManager, self).simulate(model, **kwargs)

    def preErrorCheck(self, err, dataVals=None):
        """Called before the validity check of the error values."""
        err = np.atleast_1d(err)
        if self.complex:
            if len(err) == 2:
                nData = len(dataVals) // 2
                err = pg.cat(np.ones(nData)*err[0],
                             np.abs(err[1] / dataVals[nData:]))
        else:
            if len(err) == 1:
                err = np.ones(nData)*err[0]

        return err

    def invert(self, data=None, err=None, ab2=None, mn2=None, **kwargs):
        """Invert measured data.

        Parameters
        ----------

        Keyword Arguments
        ----------------
        **kwargs
            Additional kwargs inherited from %(MethodManager1d.invert) and
            %(Inversion.run)

        Returns
        -------
        model : pg.Vector
            inversion result
        """
        if ab2 is not None and mn2 is not None:
            self.fop.setDataSpace(ab2=ab2, mn2=mn2)

        if data is not None:
            if self.complex:
                nData = len(data)//2
                self.dataTrans = pg.trans.TransCumulative()
                self.dataTrans.add(self.rhoaTrans, nData)
                self.dataTrans.add(self.phiaTrans, nData)
            else:
                self.dataTrans = pg.trans.TransLog()

            self.inv.dataTrans = self.dataTrans

        if 'layerLimits' not in kwargs:
            kwargs['layerLimits'] = [min(self.fop.mn2)/5,
                                     max(self.fop.ab2)/2]

        if 'paraLimits' in kwargs and self.complex:
            pL = kwargs['paraLimits'][1]
            kwargs['paraLimits'][1] = [pL[0]/1000, pL[1]/1000]

        return super(VESManager, self).invert(data=data, err=err, **kwargs)

    def loadData(self, fileName, **kwargs):
        """ Load simple data matrix
        """
        mat = np.loadtxt(fileName)
        if len(mat[0]) == 4:
            self.fop.setDataSpace(ab2=mat[:, 0], mn2=mat[:, 1])
            return mat.T
        if len(mat[0]) == 6:
            self.complex = True
            self.fop.setDataSpace(ab2=mat[:, 0], mn2=mat[:, 1])
            return (mat[:, 0], mat[:, 1],
                    np.array(pg.cat(mat[:, 2], mat[:, 4])),
                    np.array(pg.cat(mat[:, 3], mat[:, 5])))

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
            np.savetxt(fileName, mat,
                       header=r'ab/2\tmn/2\trhoa\terr\tphia\terrphi')
        else:
            mat = np.array([ab2, mn2, data, error]).T
            np.savetxt(fileName, mat, header=r'ab/2\tmn/2\trhoa\terr')


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


if __name__ == '__main__':
    VESManagerApp()
