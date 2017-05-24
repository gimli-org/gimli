#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""pygimli functions for dc resistivity / SIP data."""

# TODO Please sort the content into SIP package!

import pylab as P
import numpy as N

import pygimli as pg
from pygimli.utils import rndig
import string


def astausgleich(ab2org, mn2org, rhoaorg):
    """shifts the branches of a dc sounding to generate a matching curve."""
    ab2 = P.asarray(ab2org)
    mn2 = P.asarray(mn2org)
    rhoa = P.asarray(rhoaorg)
    um = P.unique(mn2)
    for i in range(len(um) - 1):
        r0, r1 = [], []
        ac = P.intersect1d(ab2[mn2 == um[i]], ab2[mn2 == um[i + 1]])
        for a in ac:
            r0.append(rhoa[(ab2 == a) * (mn2 == um[i])][0])
            r1.append(rhoa[(ab2 == a) * (mn2 == um[i + 1])][0])

        if len(r0) > 0:
            fak = P.mean(P.array(r0) / P.array(r1))
            print(fak)
            if P.isfinite(fak) and fak > 0.:
                rhoa[mn2 == um[i + 1]] *= fak

    return rhoa  # formerly pg as vector


def loadSIPallData(filename, outnumpy=False):
    """load SIP data with the columns ab/2,mn/2,rhoa and PHI with the
    corresponding frequencies in the first row."""
    if outnumpy:
        A = N.loadtxt(filename)
        fr = A[0, 3:]
        ab2 = A[1:, 0]
        mn2 = A[1:, 1]
        rhoa = A[1:, 2]
        PHI = A[1:, 3:]
    else:
        A = pg.RMatrix()
        pg.loadMatrixCol(A, 'sch/dc.ves')
        ndata = A.cols()
        ab2 = A[0](1, ndata)
        mn2 = A[1](1, ndata)
        rhoa = A[2](1, ndata)
        PHI = pg.RMatrix()
        fr = []
        for i in range(3, A.rows()):
            fr.append(A[i][0])
            PHI.push_back(A[i](1, ndata))
    return ab2, mn2, rhoa, PHI, fr


def makeSlmData(ab2, mn2, rhoa=None, filename=None):
    """generate a pygimli data container from ab/2 and mn/2 array."""
    data = pg.DataContainer()
    data.resize(len(ab2))
    pos = N.unique(N.hstack((ab2, mn2)))

    for elx in N.hstack((-pos[::-1], pos)):
        data.createElectrode(elx, 0., 0.)

    if filename is not None:
        f = open(filename, 'w')
        f.write(str(len(pos) * 2) + '\n#x y z\n')
        for elx in N.hstack((-pos[::-1], pos)):
            f.write(str(elx) + '\t0\t0\n')
            f.write(str(len(ab2)) + '\n#a\tb\tm\tn\tk\trhoa\n')

    lpos = len(pos)
    iab = pos.searchsorted(ab2)
    imn = pos.searchsorted(mn2)

    if (filename is not None) & (rhoa is None):
        rhoa = N.ones(len(ab2))

    for i in range(len(iab)):
        # print -pos[iab[i]], -pos[imn[i]], pos[imn[i]], pos[iab[i]]
        k = (ab2[i]**2 - mn2[i]**2) * N.pi / mn2[i] / 2.0
        if filename is not None:
            f.write(str(lpos - iab[i]) + '\t' + str(lpos + iab[i] + 1) + '\t')
            f.write(str(lpos - imn[i]) + '\t' + str(lpos + imn[i] + 1) + '\t')
            f.write(str(rndig(k, 4)) + '\t' + str(rhoa[i]) + '\n')
        data.createFourPointData(i, int(lpos - iab[i]), int(lpos + iab[i] + 1),
                                 int(lpos - imn[i]), int(lpos + imn[i] + 1))

    if filename is not None:
        f.close()

    data.set('rhoa', pg.asvector(rhoa))
    return data


def showsounding(ab2, rhoa, resp=None, mn2=None, islog=True, xlab=None):
    """
        Display a sounding curve (rhoa over ab/2) and an additional response.
    """
    if xlab is None:
        xlab = r'$\rho_a$ in $\Omega$m'

    ab2a = N.asarray(ab2)
    rhoa = N.asarray(rhoa)
    if mn2 is None:
        if islog:
            l1 = P.loglog(rhoa, ab2, 'rx-', label='observed')
        else:
            l1 = P.semilogy(rhoa, ab2, 'rx-', label='observed')

        P.hold(True)
        if resp is not None:
            if islog:
                l2 = P.loglog(resp, ab2, 'bo-', label='simulated')
            else:
                l2 = P.semilogy(resp, ab2, 'bo-', label='simulated')

            P.legend((l1, l2), ('obs', 'sim'), loc=0)
    else:
        for unmi in N.unique(mn2):
            if islog:
                l1 = P.loglog(rhoa[mn2 == unmi], ab2a[mn2 == unmi],
                              'rx-', label='observed')
            else:
                l1 = P.semilogy(rhoa[mn2 == unmi], ab2a[mn2 == unmi],
                                'rx-', label='observed')

            P.hold(True)
            if resp is not None:
                l2 = P.loglog(resp[mn2 == unmi], ab2a[mn2 == unmi],
                              'bo-', label='simulated')
                P.legend((l1, l2), ('obs', 'sim'))

    P.axis('tight')
    P.ylim((max(ab2), min(ab2)))
    locs = P.yticks()[0]
    if len(locs) < 2:
        locs = N.hstack((min(ab2), locs, max(ab2)))
    else:
        locs[0] = max(locs[0], min(ab2))
        locs[-1] = min(locs[-1], max(ab2))

    a = []
    for l in locs:
        a.append('%g' % rndig(l))

    P.yticks(locs, a)

    locs = P.xticks()[0]

    a = []
    for l in locs:
        a.append('%g' % rndig(l))

    P.xticks(locs, a)

    P.grid(which='both')
    P.xlabel(xlab)
    P.ylabel('AB/2 in m')
    # P.legend()
    P.show()
    return


def showsip1ddata(PHI, fr, ab2, mn2=None, cmax=None, ylab=True, cbar=True):
    """display SIP phase data as image plot."""
    P.cla()
    ax = P.gca()
    pal = P.cm.get_cmap()
    pal.set_under('w')
    pal.set_bad('w')
    if isinstance(PHI, pg.RVector):
        PHI = N.asarray(PHI)

    im = P.imshow(PHI.reshape((len(ab2), len(fr))),
                  interpolation='nearest', cmap=pal)
    if cmax is None:
        cmax = N.max(PHI)

    im.set_clim((0., cmax))

    ax.xaxis.set_label_position('top')
    P.xlabel('f in Hz')

    a = []
    df = 1
    for f in fr[::df]:
        a.append("%g" % rndig(f))

    P.xticks(N.arange(0, len(fr), df), a)
    xtl = ax.get_xticklabels()
    for i, xtli in enumerate(xtl):
        xtli.set_rotation('vertical')

    if ylab:
        a = []
        yla = 'AB/2'
        if mn2 is None:
            for i in range(len(ab2)):
                a.append(str(ab2[i]))
        else:
            yla = yla + '-MN/2'
            for i in range(len(ab2)):
                a.append('%g%g' % (rndig(ab2[i]), rndig(mn2[i])))

        P.yticks(N.arange(len(ab2)), a)
        P.ylabel(yla + ' in m')

    if cbar:
        P.colorbar(aspect=40, shrink=0.6)

    P.ylim((len(ab2) - 0.5, -0.5))
    P.show()
    P.ylim((len(ab2) - 0.5, -0.5))
    return


def showsip1dmodel(M, tau, thk, res=None, z=None,
                   cmin=None, cmax=None, islog=True):
    """
        Display an SIP Debye block model as image.
    """
    if z is None:
        z = N.cumsum(N.hstack((0., thk)))

    P.cla()
    pal = P.cm.get_cmap()
    pal.set_under('w')
    pal.set_bad('w')
    if isinstance(M, pg.RVector):
        M = N.asarray(M)

    if islog:
        M = N.log10(M)

    M = M.reshape((len(z), len(tau)))
    im = P.imshow(M, interpolation='nearest', cmap=pal)
    if cmax is None:
        cmax = N.max(M)

    if cmax is None:
        cmax = N.max(M)

    im.set_clim((cmin, cmax))

    a = []
    for t in tau[::2]:
        a.append("%g" % rndig(t * 1000, 2))

    P.xticks(N.arange(0, len(tau), 2), a)

    a = []
    for zi in z:
        a.append(str(zi))

    P.yticks(N.arange(len(z)) - 0.5, a)
    P.xlabel(r'$\tau$ in ms')
    P.ylabel('z in m')
    P.ylim((len(z) - 0.5, -0.5))
    P.colorbar(orientation='horizontal', aspect=40, shrink=0.6)

    if res is not None:
        xl = P.xlim()[1]
        for i in range(len(res)):
            P.text(xl, i, ' %g $\Omega$m' % rndig(res[i], 2))

    lgm = N.zeros((len(z), 1))
    tch = N.zeros((len(z), 1))
    lgt = N.log(tau)
    if islog:
        M = 10**M

    for n in range(len(M)):
        m = N.abs(M[n])
        tch[n] = N.sum(m)
        lgm[n] = N.exp(N.sum(m * lgt) / N.sum(m))

    tpos = N.interp(N.log(lgm), N.log(tau), N.arange(len(tau)))
    P.plot(tpos, N.arange(len(z)), 'w*')

    P.title('logarithmized spectral chargeability')
    P.show()
    return lgm, tch


class DebyeModelling(pg.ModellingBase):

    """forward operator for Debye decomposition."""

    def __init__(self, fvec, tvec=None, zero=False, verbose=False):

        if tvec is None:
            tvec = N.logspace(-4, 0, 5)

        mesh = pg.createMesh1D(len(tvec))

        if zero:
            mesh.cell(0).setMarker(-1)
            mesh.cell(len(tvec) - 1).setMarker(1)

        pg.ModellingBase.__init__(self, mesh, verbose)
        self.f_ = pg.asvector(fvec)
        self.t_ = tvec
        self.zero_ = zero

    def response(self, par):
        """phase spectrum as function of spectral chargeabilities."""
        y = pg.RVector(len(self.f_), 0.0)
        for (t, p) in zip(self.t_, par):
            wt = self.f_ * 2.0 * P.pi * t
            y = y + wt / (wt * wt + 1.) * p

        return y


def DebyeDecomposition(fr, phi, maxfr=None, tv=None, verbose=False,
                       zero=False, err=0.25e-3, lam=10., blocky=False):
    """Debye decomposition of a phase spectrum."""
    if maxfr is not None:
        idx = (fr <= maxfr) & (phi >= 0.)
        phi1 = phi[idx]
        fr1 = fr[idx]
        print("using frequencies from ", N.min(fr), " to ", N.max(fr), "Hz")
    else:
        phi1 = phi
        fr1 = fr

    if tv is None:
        tmax = 1. / N.min(fr1) / 2. / N.pi * 4.
        tmin = 1. / N.max(fr1) / 2. / N.pi / 8.
        tvec = N.logspace(N.log10(tmin), N.log10(tmax), 30)
    else:
        tvec = tv

    f = DebyeModelling(fr1, tvec, zero=zero)
    tvec = f.t_
    tm = pg.RTransLog()
    start = pg.RVector(len(tvec), 1e-4)
    if zero:
        f.region(-1).setConstraintType(0)  # smoothness
        f.region(0).setConstraintType(1)  # smoothness
        f.region(1).setConstraintType(0)  # min length
        f.regionManager().setInterRegionConstraint(-1, 0, 1.)
        f.regionManager().setInterRegionConstraint(0, 1, 1.)
        f.region(-1).setTransModel(tm)
        f.region(0).setTransModel(tm)
        f.region(1).setTransModel(tm)
        f.region(-1).setModelControl(1000.)
        f.region(1).setModelControl(1000.)
    else:
        f.regionManager().setConstraintType(1)  # smoothness

    inv = pg.RInversion(pg.asvector(phi1 * 1e-3), f, verbose)
    inv.setAbsoluteError(pg.RVector(len(fr1), err))
    inv.setLambda(lam)
    inv.setModel(start)
    inv.setBlockyModel(blocky)
    if zero:
        inv.setReferenceModel(start)
    else:
        inv.setTransModel(tm)

    mvec = inv.run()
    resp = inv.response()

    return tvec, mvec, N.array(resp) * 1e3, idx


class DoubleColeColeModelling(pg.ModellingBase):

    """
        Modelling using two Cole-Cole terms
    """

    def __init__(self, mesh, fvec, si=1.0, verbose=False):
        pg.ModellingBase.__init__(self, mesh, verbose)
        self.f_ = fvec
        self.si_ = si

    def response(self, par):
        """yields phase response response of double Cole Cole model."""
        y = pg.RVector(self.f_.size(), 0.0)
        wti = self.f_ * par[1] * 2.0 * P.pi
        wte = self.f_ * par[4] * 2.0 * P.pi
        for i in range(0, y.size()):
            cpI = 1. / (N.power(wti[i] * 1j, par[2]) + 1.)
            cpE = 1. / (N.power(wte[i] * 1j, par[5]) + 1.)
            y[i] = - N.imag(cpI) * par[0] - N.imag(cpE) * par[3] * self.si_
#            y[i] = - par[0] - N.imag(cpE) * par[3] * self.si_

        return y


def read1resfile(filename, readsecond=False, dellast=True):
    """read Radic instrument res file containing a single spectrum."""
    f = open(filename, 'r')
    line = f.readline()
    fr = []
    rhoa = []
    phi = []
    drhoa = []
    dphi = []
    while True:
        line = f.readline()
        if string.rfind(line, 'Freq') > -1:
            break

    if readsecond:
        while True:
            if string.rfind(f.readline(), 'Freq') > -1:
                break

    while True:
        line = f.readline()
        b = line.split('\t')
        if len(b) < 5:
            break

        fr.append(string.atof(b[0]))
        rhoa.append(string.atof(b[1]))
        phi.append(-string.atof(b[2]) * P.pi / 180.)
        drhoa.append(string.atof(b[3]))
        dphi.append(string.atof(b[4]) * P.pi / 180.)

    f.close()
    if dellast:
        fr.pop(0)
        rhoa.pop(0)
        phi.pop(0)
        drhoa.pop(0)
        dphi.pop(0)

    return fr, rhoa, pg.asvector(phi), drhoa, pg.asvector(dphi)


def ReadAndRemoveEM(filename, readsecond=False, doplot=False,
                    dellast=True, ePhi=0.5, ePerc=1., lam=2000.):
    """
        Read res1file and remove EM effects using a double-Cole-Cole model
        fr,rhoa,phi,dphi = ReadAndRemoveEM(filename, readsecond/doplot bools)
    """
    fr, rhoa, phi, drhoa, dphi = read1resfile(filename,
                                              readsecond,
                                              dellast=dellast)
    # forward problem
    mesh = pg.createMesh1D(1, 6)  # 6 independent parameters
    f = DoubleColeColeModelling(mesh, pg.asvector(fr), phi[2] / abs(phi[2]))
    f.regionManager().loadMap("region.control")
    model = f.createStartVector()

    # inversion
    inv = pg.RInversion(phi, f, True, False)
    inv.setAbsoluteError(phi * ePerc * 0.01 + ePhi / 1000.)
    inv.setRobustData(True)

    # inv.setCWeight(pg.RVector(6, 1.0)) # wozu war das denn gut?
    inv.setMarquardtScheme(0.8)
    inv.setLambda(lam)
    inv.setModel(model)
    erg = inv.run()
    inv.echoStatus()
    chi2 = inv.chi2()
    mod0 = pg.RVector(erg)
    mod0[0] = 0.0  # set IP term to zero to obtain pure EM term
    emphi = f(mod0)
    resid = (phi - emphi) * 1000.

    if doplot:
        s = "IP: m= " + str(rndig(erg[0])) + " t=" + str(rndig(erg[1])) + \
            " c =" + str(rndig(erg[2]))
        s += "  EM: m= " + str(rndig(erg[3])) + " t=" + str(rndig(erg[4])) + \
            " c =" + str(rndig(erg[5]))
        fig = P.figure(1)
        fig.clf()
        ax = P.subplot(111)
        P.errorbar(
            fr,
            phi *
            1000.,
            yerr=dphi *
            1000.,
            fmt='x-',
            label='measured')
        ax.set_xscale('log')
        P.semilogx(fr, f(mod0) * 1000., label='EM term (CC)')
        P.errorbar(fr, resid, yerr=dphi * 1000., label='IP term')
        ax.set_yscale('log')
        P.xlim((min(fr), max(fr)))
        P.ylim((0.1, max(phi) * 1000.))
        P.xlabel('f in Hz')
        P.ylabel(r'-$\phi$ in mrad')
        P.grid(True)
        P.title(s)
        P.legend(loc=2)  # ('measured','2-cole-cole','residual'))
        fig.show()

    return N.array(fr), N.array(rhoa), N.array(resid), N.array(
        phi) * 1e3, dphi, chi2, N.array(emphi) * 1e3
