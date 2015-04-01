#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    Spectral induced polarisation (SIP) plotting tools
"""

import matplotlib.pyplot as plt


def showAmplitudeSpectrum(ax, freq, amp, ylabel=r'$\rho_a$ in $\Omega$m',
                          grid=True, marker='+', ylog=True, **kwargs):
    """ show amplitude spectrum """
    lab = kwargs.pop('label', 'obs')
    ax.semilogx(freq, amp, marker=marker, label=lab, **kwargs)
    if ylog:
        ax.set_yscale('log')
    ax.set_ylim(min(amp) * .99, max(amp * 1.01))
    ax.set_xlabel('f in Hz')
    ax.set_ylabel(ylabel)
    ax.grid(grid)


def showPhaseSpectrum(ax, freq, phi, ylabel=r'$\phi_a$ in mrad',
                      grid=True, marker='+', ylog=False, **kwargs):
    """ show phase spectrum """
    ax.semilogx(freq, phi, marker=marker, label='obs', **kwargs)
    if ylog:
        ax.set_yscale('log')
    ax.set_xlabel('f in Hz')
    ax.set_ylabel(ylabel)
    ax.grid(grid)


def plotSpectrum(ax, freq, vals, ylabel=r'$\phi_a$ in mrad',
                 grid=True, marker='+', ylog=True, **kwargs):
    """ show phase spectrum """
    ax.loglog(freq, vals, marker=marker, label='obs', **kwargs)
    if ylog:
        ax.set_yscale('log')
    ax.set_xlabel('f in Hz')
    ax.set_ylabel(ylabel)
    ax.grid(grid)


def showSpectrum(freq, amp, phi, nrows=2, ylog=None):
    """ show amplitude and phase spectra in two subplots """
    fig, ax = plt.subplots(nrows=nrows, sharex=(nrows == 2))
    showAmplitudeSpectrum(ax[0], freq, amp, ylog=ylog)
    showPhaseSpectrum(ax[1], freq, phi, ylog=ylog)
    return fig, ax


if __name__ == "__main__":
    pass