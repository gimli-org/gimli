# -*- coding: utf-8 -*-
"""Some specialization to the trans functions."""

from pygimli.core import _pygimli_ as pg

__TransCumulative_addForGC__ = pg.RTransCumulative.add


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


pg.RTransCumulative.add = __TransCumulative_addForGC_MP__

# Aliases
Trans = pg.RTrans
TransLinear = pg.RTransLinear
TransLin = pg.RTransLin
TransPower = pg.RTransPower
TransLog = pg.RTransLog
TransLogLU = pg.RTransLogLU
TransCotLU = pg.RTransCotLU
TransCumulative = pg.RTransCumulative
