#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""View ERT data
"""

import numpy as np
from numpy import ma

import pygimli as pg
from pygimli.mplviewer.dataview import showValMapPatches


def showERTData(data, vals=None, **kwargs):
    """Plot ERT data as pseudosection matrix (position over separation).

    Creates figure, axis and draw a pseudosection.

    Parameters
    ----------

    data : :gimliapi:`BERT::DataContainerERT`

    **kwargs :

        * axes : matplotlib.axes
            Axes to plot into. Default is None and a new figure and
            axes are created.
        * vals : Array[nData]
            Values to be plotted. Default is data('rhoa').
    """
    var = kwargs.pop('var', 0)
    if var > 0:
        import pybert as pb
        pg._g(kwargs)
        return pb.showData(data, vals, var=var, **kwargs)

    # remove ax keyword global
    ax = kwargs.pop('ax', None)

    if ax is None:
        fig = pg.plt.figure()
        ax = None
        axTopo = None
        if 'showTopo' in kwargs:
            ax = fig.add_subplot(1, 1, 1)
#            axs = fig.subplots(2, 1, sharex=True)
#            # Remove horizontal space between axes
#            fig.subplots_adjust(hspace=0)
#            ax = axs[1]
#            axTopo = axs[0]
        else:
            ax = fig.add_subplot(1, 1, 1)

    pg.checkAndFixLocaleDecimal_point(verbose=False)

    if vals is None:
        vals = 'rhoa'

    if isinstance(vals, str):
        if data.haveData(vals):
            vals = data(vals)
        else:
            pg.critical('field not in data container: ', vals)

        kwargs['cMap'] = kwargs.pop('cMap', pg.utils.cMap('rhoa'))
    
    kwargs['logScale'] = kwargs.pop('logScale', min(vals) > 0.0)

    ax, cbar = drawERTData(ax, data, vals=vals, **kwargs)
    
    #TODO here cbar handling like pg.show

    if 'xlabel' in kwargs:
        ax.set_xlabel(kwargs['xlabel'])
    if 'ylabel' in kwargs:
        ax.set_ylabel(kwargs['ylabel'])

    if 'showTopo' in kwargs:
        # if axTopo is not None:
        print(ax.get_position())
        axTopo = pg.plt.axes([ax.get_position().x0,
                           ax.get_position().y0,
                           ax.get_position().x0+0.2,
                           ax.get_position().y0+0.2])

        x = pg.x(data)
        x *= (ax.get_xlim()[1] - ax.get_xlim()[0]) / (max(x)-min(x))
        x += ax.get_xlim()[0]
        axTopo.plot(x, pg.z(data), '-o', markersize=4)
        axTopo.set_ylim(min(pg.z(data)), max(pg.z(data)))
        axTopo.set_aspect(1)

    # ax.set_aspect('equal')
    # plt.pause(0.1)
    pg.mplviewer.updateAxes(ax)
    return ax, cbar
    

def drawERTData(ax, data, vals=None, **kwargs):
    """Plot ERT data as pseudosection matrix (position over separation).

    Parameters
    ----------
    data : pybert.DataContainerERT
        data container with sensorPositions and a/b/m/n fields
    vals : iterable of data.size() [data('rhoa')]
        vector containing the vals to show
    ax : mpl.axis
        axis to plot, if not given a new figure is created
    cMin/cMax : float
        minimum/maximum color vals
    logScale : bool
        logarithmic colour scale [min(A)>0]
    label : string
        colorbar label

    **kwargs:
        * dx : float
            x-width of individual rectangles
        * ind : integer iterable or IVector
            indices to limit display
        * circular : bool
            Plot in polar coordinates when plotting via patchValMap
    Returns
    -------
    ax:
        The used Axes
    cbar:
        The used Colorbar or None
    """
    if vals is None:
        vals = data('rhoa')
        
    valid = data.get("valid").array().astype("bool")
    vals = ma.array(vals, mask=~valid)

    ind = kwargs.pop('ind', None)

    if ind is not None:
        vals = vals[ind]
        mid, sep = midconfERT(data, ind)
    else:
        mid, sep = midconfERT(data, circular=kwargs.get('circular', False))

    var = kwargs.pop('var', 0)
    cbar = None

    dx = kwargs.pop('dx', np.median(np.diff(np.unique(mid))))*2
    ax, cbar, ymap = showValMapPatches(vals, xVec=mid, yVec=sep, 
                                       dx=dx, ax=ax, **kwargs)
        
    if kwargs.get('circular', False):
        sM = np.mean(data.sensors(), axis=0)
        a = np.array([np.arctan2(s[1]-sM[1], s[0]-sM[0]) for s in data.sensors()])
        p = list(range(len(a)))
        # p.append(0)
        ax.plot(np.cos(a)[p], np.sin(a)[p], 'o', color='black')

        for i in range(len(a)):
            ax.text(1.15 * np.cos(a[i]),
                    1.15 * np.sin(a[i]), str(i+1),
                    horizontalalignment='center',
                    verticalalignment='center')
        ax.set_axis_off()
        ax.set_aspect(1)

    else:
        ytl = generateConfStr(np.sort([int(k) for k in ymap]))
        if 'DD1' in ytl and 'WB2' in ytl and 'DD2'not in ytl:
            ytl[ytl.index('DD1')] = 'WB1'
        if 'WA1' in ytl and 'SL2' in ytl and 'WA2'not in ytl:
            ytl[ytl.index('WA1')] = 'SL1'

        yt = ax.get_yticks()
        yt = np.unique(yt.clip(0, len(ytl)-1))
#        if yt[0] == yt[1]:
#            yt = yt[1:]
        dyt = np.diff(yt)
        if dyt[-1] < dyt[-2]:
            yt = yt[:-1]
        ax.set_yticks(yt)
        ax.set_yticklabels([ytl[int(yti)] for yti in yt])
    return ax, cbar


def midconfERT(data, ind=None, rnum=1, circular=False):
    """Return the midpoint and configuration key for ERT data.

    Return the midpoint and configuration key for ERT data.

    Parameters
    ----------
    data : pybert.DataContainerERT
        data container with sensorPositions and a/b/m/n fields

    ind : []
        Documentme

    rnum : []
        Documentme

    circular : bool
        Return midpoint in degree (rad) instead if meter.

    Returns
    -------
    mid : np.array of float
        representative midpoint (middle of MN, AM depending on array)
    conf : np.array of float
        configuration/array key consisting of
        1) array type (Wenner-alpha/beta, Schlumberger, PP, PD, DD, MG)
            00000: pole-pole
            10000: pole-dipole or dipole-pole
            30000: Wenner-alpha
            40000: Schlumberger or Gradient
            50000: dipole-dipole or Wenner-beta
        2) potential dipole length (in electrode spacings)
        3) separation factor (current dipole length or (di)pole separation)
    """
#    xe = np.hstack((pg.x(data.sensorPositions()), np.nan))  # not used anymore
    x0 = data.sensorPosition(0).x()
    xe = pg.x(data.sensorPositions()) - x0
    ux = pg.unique(xe)

    if len(ux) * 2 > data.sensorCount():  # 2D with topography case
        dx = np.array(pg.utils.diff(pg.utils.cumDist(data.sensorPositions())))
        dxM = pg.mean(dx)
        if min(pg.y(data)) != max(pg.y(data)) or \
           min(pg.z(data)) != max(pg.z(data)):
            # Topography case
            if (pg.max(abs(dx-dxM)) < dxM*0.9):
                # if the maximum spacing < meanSpacing/2 we assume equidistant
                # spacing and no missing electrodes
                dx = np.ones(len(dx)) * dxM
            else:
                # topography with probably missing electrodes
                dx = np.floor(dx/np.round(dxM))*dxM
                pass
        if max(dx) < 0.5:
            print("Detecting small distances, using mm accuracy")
            rnum = 3
        xe = np.hstack((0., np.cumsum(np.round(dx, rnum)), np.nan))

        de = np.median(np.diff(xe[:-1])).round(rnum)
        ne = np.round(xe/de)
    else:  # 3D (without topo) case => take positions directly
        de = np.median(np.diff(ux)).round(1)
        ne = np.array(xe/de, dtype=int)

    # a, b, m, n = data('a'), data('b'), data('m'), data('n')
    # check if xe[a]/a is better suited (has similar size)
    if circular:
        # for circle geometry
        center = np.mean(data.sensorPositions(), axis=0)
        r = data.sensors()[0].distance(center)
        s0 = data.sensors()[0]-center
        s1 = data.sensors()[1]-center
        p0 = np.arctan2(s0[1], s0[0])
        p1 = np.arctan2(s1[1], s1[0])
        if p1 > p0:
            #rotate left
            x=np.cos(np.linspace(0, 2*np.pi, data.sensorCount()+1)+p0)[:-1] * r
            y=np.sin(np.linspace(0, 2*np.pi, data.sensorCount()+1)+p0)[:-1] * r
        else:
            x=np.cos(np.linspace(2*np.pi, 0, data.sensorCount()+1)+p0)[:-1] * r
            y=np.sin(np.linspace(2*np.pi, 0, data.sensorCount()+1)+p0)[:-1] * r

        a = np.array([np.arctan2(y[i], x[i]) for i in data['a']])
        b = np.array([np.arctan2(y[i], x[i]) for i in data['b']])
        m = np.array([np.arctan2(y[i], x[i]) for i in data['m']])
        n = np.array([np.arctan2(y[i], x[i]) for i in data['n']])

        a = np.unwrap(a) % (np.pi*2)
        b = np.unwrap(b) % (np.pi*2)
        m = np.unwrap(m) % (np.pi*2)
        n = np.unwrap(n) % (np.pi*2)

    else:
        a = np.array([ne[int(i)] for i in data('a')])
        b = np.array([ne[int(i)] for i in data('b')])
        m = np.array([ne[int(i)] for i in data('m')])
        n = np.array([ne[int(i)] for i in data('n')])

    if ind is not None:
        a = a[ind]
        b = b[ind]
        m = m[ind]
        n = n[ind]

    anan = np.isnan(a)
    a[anan] = b[anan]
    b[anan] = np.nan
    ab, am, an = np.abs(a-b), np.abs(a-m), np.abs(a-n)
    bm, bn, mn = np.abs(b-m), np.abs(b-n), np.abs(m-n)

    if circular:
        for v in [ab, mn, bm, an]:
            v[v > np.pi] = 2*np.pi - v[v > np.pi]

    # 2-point (default) 00000
    sep = np.abs(a-m)
    mid = (a+m) / 2

    # 3-point (PD, DP) (now only b==-1 or n==-<1, check also for a and m)
    imn = np.isfinite(n)*np.isnan(b)
    mid[imn] = (m[imn]+n[imn]) / 2
    sep[imn] = np.minimum(am[imn], an[imn]) + 10000 + 100 * (mn[imn]-1) + \
        (np.sign(a[imn]-m[imn])/2+0.5) * 10000
    iab = np.isfinite(b)*np.isnan(n)
    mid[iab] = (a[iab]+b[iab]) / 2  # better 20000 or -10000?
    sep[iab] = np.minimum(am[iab], bm[iab]) + 10000 + 100 * (ab[iab]-1) + \
        (np.sign(a[iab]-n[iab])/2+0.5) * 10000
    #  + 10000*(a-m)

    # 4-point alpha: 30000 (WE) or 4000 (SL)
    iabmn = np.isfinite(a) & np.isfinite(b) & np.isfinite(m) & np.isfinite(n)
    ialfa = np.copy(iabmn)
    ialfa[iabmn] = (ab[iabmn] >= mn[iabmn]+2)  # old
    mnmid = (m[iabmn] + n[iabmn]) / 2
    ialfa[iabmn] = np.sign((a[iabmn]-mnmid)*(b[iabmn]-mnmid)) < 0

    mid[ialfa] = (m[ialfa] + n[ialfa]) / 2
    spac = np.minimum(bn[ialfa], bm[ialfa])
    abmn3 = np.round((3*mn[ialfa]-ab[ialfa])*10000)/10000
    sep[ialfa] = spac + (mn[ialfa]-1)*100*(abmn3 != 0) + \
        30000 + (abmn3 < 0)*10000
    # gradient

    # %% 4-point beta
    ibeta = np.copy(iabmn)
    ibeta[iabmn] = (bm[iabmn] >= mn[iabmn]) & (~ialfa[iabmn])

    if circular:
        # print(ab[ibeta])
        ibeta = np.copy(iabmn)

        def _averageAngle(vs):
            sumsin = 0
            sumcos = 0

            for v in vs:
                sumsin += np.sin(v)
                sumcos += np.cos(v)

            return np.arctan2(sumsin, sumcos)

        abC = _averageAngle([a[ibeta], b[ibeta]])
        mnC = _averageAngle([m[ibeta], n[ibeta]])

        mid[ibeta] = _averageAngle([abC, mnC])

        # special case when dipoles are completely opposite
        iOpp = abs(abs((mnC - abC)) - np.pi) < 1e-3
        mid[iOpp] = _averageAngle([b[iOpp], m[iOpp]])

        minAb = min(ab[ibeta])
        sep[ibeta] = 50000 + (np.round(ab[ibeta]/minAb)) * 100 + \
            np.round(np.minimum(np.minimum(am[ibeta], an[ibeta]),
                                np.minimum(bm[ibeta], bn[ibeta])) / minAb)
    else:
        mid[ibeta] = (a[ibeta] + b[ibeta] + m[ibeta] + n[ibeta]) / 4

        sep[ibeta] = 50000 + (ab[ibeta]-1) * 100 + np.minimum(
            np.minimum(am[ibeta], an[ibeta]), np.minimum(bm[ibeta], bn[ibeta]))

    # %% 4-point gamma
    # multiply with electrode distance and add first position
    if not circular:
        mid *= de
        mid += x0
    return mid, sep


def generateConfStr(yy):
    """Generate configuration string to characterize array."""
    types = ['PP', 'PD', 'DP', 'WA', 'SL', 'DD']  # base types
    spac = yy % 100  # source-receiver distance
    dip = np.round(yy//100) % 100  # MN dipole length
    typ = np.round(yy//10000)
    # check if SL is actually GR (multi-gradient)

    # check if DD-n-n should be renamed
    rendd = (np.mean(spac / (dip+1)) < 2.1)
    keys = []
    for s, d, t in zip(spac, dip, typ):
        key = types[t]
        if d > 0:
            if rendd and d+1 == s and t == 5:
                key = 'WB'
            else:
                key = key + str(d+1) + '-'
        key = key + "{:2d}".format(s)  # str(s)
        keys.append(key)

    return keys
