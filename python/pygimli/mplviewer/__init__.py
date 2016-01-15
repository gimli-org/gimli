# -*- coding: utf-8 -*-
"""
    Viewer interface .. depends on matplotlib
"""

holdAxes_ = 0

def updateAxes(ax, a=None):
    """
        for internal use
    """
    if not holdAxes_:
        try:
            mpl.pyplot.pause(0.01)
        except Exception as e:
            #print(e)
            pass


from .dataview import *
from .meshview import *
from .colorbar import *
from .overlayimage import *

import matplotlib.pyplot as plt
import numpy as np

def hold(val=1):
    pg.mplviewer.holdAxes_ = val

def showLater(val=1):
    raise('do not use .. use show(hold=1) to keep pics in background')
    import matplotlib.pyplot as plt
    if val == 1:
        plt.ion()
    else:
        plt.ioff()
        plt.show()

def wait():
    plt.pause(0.01)
    plt.show()

goldenMean = 1.618  # (1.0 + math.sqrt(5.0)) / 2.0


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
    inches_per_pt = 1.0 / 72.27  # pt/inch (latex)
    golden_mean = (1.0 + math.sqrt(5.0)) / 2.0

    textwidth = (a4[0] - paperMargin) * inches_per_cm

    fig_width = textwidth * xScale  # fig width in inches
    fig_height = textwidth * yScale  # fig height in inches

    fig_size = [fig_width * scale, fig_height * scale]

    # print "figsize:", fig_size
#    fig.set_size_inches(fig_size)

    #from matplotlib import rc
    # rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    # rc('font',**{'family':'serif','serif':['Palatino']})

    params = {'backend': 'ps',
              'font.size': fontsize * scale,
              #'font.weight'       : 'bold',
              'axes.labelsize': fontsize * scale,
              'font.size': fontsize * scale,
              'legend.fontsize': fontsize * scale,
              'xtick.labelsize': fontsize * scale,
              'ytick.labelsize': fontsize * scale,
              # font.sans-serif     : Bitstream Vera Sans, Lucida Grande, Verdana, Geneva, Lucid, Arial, Helvetica
              #'font.cmb10'     : 'cmb10',
              #'font.family'         : 'cursive',
              'font.family': 'sans-serif',
              #'font.sans-serif'   : 'Helvetica',
              'text.usetex': usetex,
              'figure.figsize': fig_size,
              'xtick.major.pad': 4 * scale,
              'xtick.minor.pad': 4 * scale,
              'ytick.major.pad': 4 * scale,
              'ytick.minor.pad': 4 * scale,
              'xtick.major.size': 4 * scale,     # major tick size in points
              'xtick.minor.size': 2 * scale,     # minor tick size in points
              'ytick.major.size': 4 * scale,     # major tick size in points
              'ytick.minor.size': 2 * scale,     # minor tick size in points
              'lines.markersize': 6 * scale,
              'lines.linewidth': 0.6 * scale
              }

    plt.rcParams.update(params)

# def setOutPutStyle

def createAnimation(fig, animate, nFrames, dpi, out):
    """
        Create animation for the content of a given matplotlib figure.
    
        Until I know a better place.
    """
    anim = animation.FuncAnimation(fig, animate,
                                   frames=nFrames,
                                   interval=0.001, repeat=False)
    anim.save(out + ".mp4", writer=None, fps=20, dpi=dpi, codec=None,
              bitrate=24*1024, extra_args=None, metadata=None,
              extra_anim=None, savefig_kwargs=None)
    try:
        print("Create frames ... ")
        os.system('mkdir -p anim-' + out)
        os.system('ffmpeg -i ' + out + '.mp4 anim-' + out + '/movie%d.jpg')
    except:
        pass


import matplotlib.animation as animation

def saveAnimation(mesh, data, out, vData=None, plc=None, label='',
                  cMin=None, cMax=None, logScale=False, cmap=None, **kwargs):
    """
        Create and save an animation for a given mesh with a set of field data.
        
        Until I know a better place.
    """
    dpi=92
    scale=1
    fig = plt.figure(facecolor='white',
                     figsize=(scale*800/dpi, scale*490/dpi), dpi=dpi)  
    ax = fig.add_subplot(1,1,1)
        
    gci = pg.mplviewer.drawModel(ax, mesh, data=data[0],
                                 cMin=cMin, cMax=cMax, cmap=cmap,
                                 logScale=logScale)
    
    cbar = pg.mplviewer.createColorbar(gci, label=label, pad=0.55)
    ax.set_ylabel('Depth [m]')
    ax.set_xlabel('$x$ [m]')
        
    ticks = ax.yaxis.get_majorticklocs()
    tickLabels = []
    for t in ticks:
        tickLabels.append(str(int(abs(t))))

    ax.set_yticklabels(tickLabels)
    
    if plc:
        pg.show(plc, axes=ax)
    
    plt.tight_layout()
    plt.pause(0.001)
    
    def animate(i):
        print(out + ": Frame:", i, "/", len(data))
        
        if not vData is None:
            ax.clear()
            pg.mplviewer.holdAxes_ = 1
            pg.mplviewer.drawModel(ax, mesh, data=data[i],
                                 cMin=cMin, cMax=cMax, cmap=cmap,
                                 logScale=logScale)
            pg.mplviewer.drawStreams(ax, mesh, vData[i], **kwargs)
        else:
            pg.mplviewer.setMappableData(gci, data[i], 
                                         cMin=cMin, cMax=cMax,
                                         logScale=logScale)
            
            
        plt.pause(0.1)
    createAnimation(fig, animate, int(len(data)), dpi, out)
    