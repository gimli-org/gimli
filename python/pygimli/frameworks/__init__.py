# -*- coding: utf-8 -*-
"""Unified and method independent inversion frameworks."""


from .modelling import (Modelling, Block1DModelling, MeshModelling,
                        PetroModelling, LCModelling, ParameterModelling)

from .inversion import (Inversion, MarquardtInversion,
                        Block1DInversion,
                        LCInversion)

from .methodManager import (MethodManager, MethodManager1d,
                            ParameterInversionManager,
                            MeshMethodManager, PetroInversionManager)

# will be removed very soon
# from .methodManager import (MethodManager0,
#                             MeshMethodManager0)

from .resolution import computeR

from .harmfit import HarmFunctor, harmfit, harmfitNative

__all__ = ['HarmFunctor', 'harmfitNative', 'harmfit', 'computeR']
