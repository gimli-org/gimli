# -*- coding: utf-8 -*-
"""Some specialization to the trans functions."""

from .core import pgcore

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
