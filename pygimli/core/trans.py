# -*- coding: utf-8 -*-
"""Some specialization to the trans functions."""

from .core import pgcore
import numpy as np

__TransCumulative_addForGC__ = pgcore.RTransCumulative.add


def __TransCumulative_addForGC_MP__(self, T, *args):
    """Don't use directly.

    Monkey patch to keep the GC happy until redesign.
    """
    if not hasattr(self, '__trans__'):
        self.__trans__ = []
    self.__trans__.append(T)

    if len(args) == 1:
        # be sure avoid auto conversion from int to IndexArray
        if isinstance(args[0], int):
            return __TransCumulative_addForGC__(self, T, size=args[0])
    return __TransCumulative_addForGC__(self, T, *args)


pgcore.RTransCumulative.add = __TransCumulative_addForGC_MP__

# Aliases
Trans = pgcore.RTrans
TransLinear = pgcore.RTransLinear
TransLin = pgcore.RTransLin
TransPower = pgcore.RTransPower
TransLog = pgcore.RTransLog
TransLogLU = pgcore.RTransLogLU
TransCotLU = pgcore.RTransCotLU
TransCumulative = pgcore.RTransCumulative


class TransSymLog(pgcore.RTrans):
    """Transformation using a bilogarithmic scaling."""

    def __init__(self, tol=1e-12):
        """Forward transformation."""
        super().__init__()
        self.tol = tol

    def trans(self, x):
        """Forward transformation."""
        return pgcore.log(1 + np.abs(x) / self.tol) * np.sign(x)

    def invTrans(self, y):
        """Inverse transformation."""
        return (pgcore.exp(np.abs(y)) - 1.) * self.tol * np.sign(y)

    def deriv(self, x):
        """Derivative of the transformation."""
        return 1. / (np.abs(x) / self.tol + 1) / self.tol
