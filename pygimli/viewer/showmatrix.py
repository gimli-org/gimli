# -*- coding: utf-8 -*-
"""Generic matrix visualization tools."""

import matplotlib.pyplot as plt
import numpy as np

import pygimli as pg

from .mpl import createColorBar
from .mpl.matrixview import drawBlockMatrix, drawSparseMatrix


def showMatrix(mat, ax=None, **kwargs):
    """Show various pyGIMLi matrices using matplotlib.

    Parameters
    ----------
    mat : pyGIMLi matrix
    ax : mpl.axes
    **kwargs : forwarded to mpl plotting commands
    """
    if ax is None:
        ax = plt.subplots()[1]

    if isinstance(mat, pg.core.RSparseMapMatrix):
        gci = drawSparseMatrix(ax, mat)
        cBar = None
    elif isinstance(mat, pg.core.BlockMatrix):
        gci = drawBlockMatrix(ax, mat)
        cBar = createColorBar(gci,
                              size=kwargs.pop('size', 0.2),
                              pad=kwargs.pop('pad', None),
                              label=kwargs.pop('label', None))
        ids = mat.entries_matrixID()
        cBar.set_ticks(np.array(ids))
        cBar.set_ticklabels(["%d" % lab for lab in ids])
        cBar.set_label("Matrix ID")
        ax.set_title(mat.__str__() + "\n")

    else:
        pg.error("Matrix type not supported yet.")
    return ax, cBar
