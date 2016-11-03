#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    Some additional infos here?
"""
import os
from pygimli.physics.SIP import SIPSpectrum


if __name__ == "__main__":

    sip = SIPSpectrum(os.path.dirname(__file__) + '/sipexample.txt')
    sip.showData()  # show amplitude and phase spectra
    sip.showDataKK()  # show data along with Kramers-Kronig calculated (data check)
    sip.showData(znorm=True)  # show normalized imaginary/real parts instead

    if 0:  # Pelton approach: fit Cole-Cole model along with EM Cole-Cole term
        sip.fitCCEM()
    else:  # determine permittivity to correct values with
        sip.removeEpsilonEffect()
        sip.fitColeCole(useCond=False)
    # %%

    # this fails
    # sip.fitDebyeModel(showFit=True)
    # %% create titles and plot data, fit and model
    sip.showAll(save=True)
