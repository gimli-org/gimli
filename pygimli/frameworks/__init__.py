# -*- coding: utf-8 -*-
"""Unified and method independent inversion frameworks."""


from .modelling import (Modelling, Block1DModelling, MeshModelling,
                        JointModelling, PriorModelling,
                        PetroModelling, LCModelling, ParameterModelling)

from .inversion import (Inversion, MarquardtInversion,
                        Block1DInversion,
                        LCInversion)

from .methodManager import (fit, MethodManager, MethodManager1d,
                            ParameterInversionManager,
                            MeshMethodManager,
                            PetroInversionManager,
                            JointPetroInversionManager)

# will be removed very soon
# from .methodManager import (MethodManager0,
#                             MeshMethodManager0)

# from .resolution import computeR

from .harmfit import HarmFunctor, harmfit, harmfitNative

__all__ = ['HarmFunctor', 'harmfitNative', 'harmfit']
