#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    Spectral induced polarisation (SIP) plotting tools
"""

import matplotlib.pyplot as plt


def showAmplitudeSpectrum(ax, freq, amp, ylabel=r'$\rho$ in $\Omega$m',
                          grid=True, marker='+', ylog=True, **kwargs):
    """show amplitude spectrum"""
    lab = kwargs.pop('label', 'obs')
    ax.semilogx(freq, amp, marker=marker, label=lab, **kwargs)
    if ylog:
        ax.set_yscale('log')
    ax.set_ylim(min(amp) * .99, max(amp * 1.01))
    ax.set_xlabel('f in Hz')
    ax.set_ylabel(ylabel)
    ax.grid(grid)


def showPhaseSpectrum(ax, freq, phi, ylabel=r'$\phi$ in mrad',
                      grid=True, marker='+', ylog=False, **kwargs):
    """show phase spectrum"""
    if 'label' not in kwargs:
        kwargs['label'] = 'obs'
    ax.semilogx(freq, phi, marker=marker, **kwargs)
    if ylog:
        ax.set_yscale('log')
    ax.set_xlabel('f in Hz')
    ax.set_ylabel(ylabel)
    ax.grid(grid)


def plotSpectrum(ax, freq, vals, ylabel=r'$\phi$ in mrad',
                 grid=True, marker='+', ylog=True, **kwargs):
    """show phase spectrum"""
    if 'label' not in kwargs:
        kwargs['label'] = 'obs'
    ax.loglog(freq, vals, marker=marker, **kwargs)
    if ylog:
        ax.set_yscale('log')
    ax.set_xlabel('f in Hz')
    ax.set_ylabel(ylabel)
    ax.grid(grid)


def showSpectrum(freq, amp, phi, nrows=2, ylog=None, ax=None):
    """show amplitude and phase spectra in two subplots"""
    if ax is None:
        fig, ax = plt.subplots(nrows=nrows, sharex=(nrows == 2))
    else:
        fig = ax[0].figure
    showAmplitudeSpectrum(ax[0], freq, amp, ylog=ylog)
    showPhaseSpectrum(ax[1], freq, phi, ylog=ylog)
    return fig, ax


if __name__ == "__main__":
    pass
