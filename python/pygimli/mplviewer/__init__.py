# -*- coding: utf-8 -*-

from .dataview import *  

from .meshview import *
from .colorbar import *
from .overlayimage import *

import matplotlib.pyplot as plt
import numpy as np

goldenMean = 1.618 #(1.0 + math.sqrt(5.0)) / 2.0

def setOutputStyle(dim='w', paperMargin=5, xScale=1.0, yScale=1.0,
                   fontsize=9, scale=1, usetex=True):
    """
    """

    if dim == 'w':
        dim = 0
    else:
        dim = 1

    a4 = [21.0, 29.7]

    inches_per_cm = 1. / 2.54
    inches_per_pt   = 1.0 / 72.27 # pt/inch (latex)
    golden_mean = (1.0 + math.sqrt(5.0)) / 2.0

    textwidth       = (a4[0] - paperMargin) * inches_per_cm

    fig_width   = textwidth * xScale # fig width in inches
    fig_height  = textwidth * yScale # fig height in inches

    fig_size    = [fig_width * scale, fig_height * scale]

    #print "figsize:", fig_size
#    fig.set_size_inches(fig_size)

    #from matplotlib import rc
    ##rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    ##rc('font',**{'family':'serif','serif':['Palatino']})

    params = { 'backend': 'ps',
                'font.size'         : fontsize * scale,
                #'font.weight'       : 'bold',
                'axes.labelsize'    : fontsize * scale,
                'text.fontsize'     : fontsize * scale,
                'legend.fontsize'   : fontsize * scale,
                'xtick.labelsize'   : fontsize * scale,
                'ytick.labelsize'   : fontsize * scale,
                #font.sans-serif     : Bitstream Vera Sans, Lucida Grande, Verdana, Geneva, Lucid, Arial, Helvetica
                #'font.cmb10'     : 'cmb10',
                #'font.family'         : 'cursive',
                'font.family'          : 'sans-serif',
                #'font.sans-serif'   : 'Helvetica',
                'text.usetex'          : usetex,
                'figure.figsize'       : fig_size,
                'xtick.major.pad'      : 4 * scale,
                'xtick.minor.pad'      : 4 * scale,
                'ytick.major.pad'      : 4 * scale,
                'ytick.minor.pad'      : 4 * scale,
                'xtick.major.size'     : 4 * scale,     # major tick size in points
                'xtick.minor.size'     : 2 * scale,     # minor tick size in points
                'ytick.major.size'     : 4 * scale,     # major tick size in points
                'ytick.minor.size'     : 2 * scale,     # minor tick size in points
                'lines.markersize'     : 6 * scale,
                'lines.linewidth'      : 0.6 * scale
            }

    plt.rcParams.update(params)

# def setOutPutStyle
