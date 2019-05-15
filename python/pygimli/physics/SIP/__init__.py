#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Spectral induced polarization (SIP) measurements and fittings."""

from .tools import (isComplex, toComplex, toPolar, squeezeComplex)

from .importData import (load)
from .sipspectrum import (SpectrumManager, SIPSpectrum)

from .models import (ColeColeRho, ColeColeSigma, ColeColePhi, DoubleColeColePhi,
                    tauRhoToTauSigma)

from .plotting import showSpectrum, drawPhaseSpectrum, drawAmplitudeSpectrum

#__all__ = [name for name in dir() if '_' not in name]
#__all__ = [SIPSpectrum, ColeColeRho, ColeColePhi, DoubleColeColePhi]
