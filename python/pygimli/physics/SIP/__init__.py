#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    Spectral induced polarization (SIP) measurements and fittings
"""

from .models import ColeColeRho, ColeColePhi, DoubleColeColePhi

from .sipspectrum import SIPSpectrum

#__all__ = [name for name in dir() if '_' not in name]
#__all__ = [SIPSpectrum, ColeColeRho, ColeColePhi, DoubleColeColePhi]
