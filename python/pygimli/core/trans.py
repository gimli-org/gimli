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
    return __TransCumulative_addForGC__(self, T, *args)

pg.RTransCumulative.add = __TransCumulative_addForGC_MP__
