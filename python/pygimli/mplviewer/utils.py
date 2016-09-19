# coding=utf-8
"""Plotting utilities used througout the mplviewer package."""

import os

import matplotlib.animation as animation
import matplotlib.pyplot as plt

import pygimli as pg

holdAxes__ = 0


def updateAxes(ax, a=None):
    """For internal use."""
    if not holdAxes__:
        try:
            plt.pause(0.1)
        except BaseException as _:
            print(ax, a)


def hold(val=1):
    """TODO WRITEME."""
    pg.mplviewer.holdAxes__ = val


def wait():
    """TODO WRITEME."""
    plt.pause(0.1)
    plt.show()


def adjustWorldAxes(ax):
    """Set some common default properties for an axe."""
    ax.set_ylabel('Depth [m]')
    ax.set_xlabel('$x$ [m]')

    ticks = ax.yaxis.get_majorticklocs()
    tickLabels = []
    for t in ticks:
        tickLabels.append(str(int(abs(t))))

    ax.set_yticklabels(tickLabels)
    plt.tight_layout()
    plt.pause(0.01)


def saveFigure(fig, filename, pdfTrim=False):
    """Save figure as pdf."""
    fig.savefig(filename + '.pdf', bbox_inches='tight')
    #pdfTrim=1
    if pdfTrim:
        try:
            print("trying pdf2pdfS ... ")
            os.system('pdf2pdfBB ' + filename + '.pdf')
            os.system('pdf2pdfS ' + filename + '.pdf')
        except BaseException as _:
            print("fail local convert. Should be no problem.")


def saveAxes(ax, filename, adjust=False):
    """Save axes as pdf."""
    if adjust:
        adjustWorldAxes(ax)

    plt.pause(0.01)
    saveFigure(ax.figure, filename)

def setOutputStyle(dim='w', paperMargin=5, xScale=1.0, yScale=1.0, fontsize=9,
                   scale=1, usetex=True):
    """Set preferred output style."""
    if dim == 'w':
        dim = 0
    else:
        dim = 1

    a4 = [21.0, 29.7]

    inches_per_cm = 1. / 2.54
    # inches_per_pt = 1.0 / 72.27  # pt/inch (latex)
    # goldenMean = (1.0 + np.sqrt(5.0)) / 2.0

    textwidth = (a4[0] - paperMargin) * inches_per_cm

    fig_width = textwidth * xScale  # fig width in inches
    fig_height = textwidth * yScale  # fig height in inches

    fig_size = [fig_width * scale, fig_height * scale]

    # print "figsize:", fig_size
    # fig.set_size_inches(fig_size)

    # from matplotlib import rc
    # rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    # rc('font',**{'family':'serif','serif':['Palatino']})

    params = {
        'backend': 'ps',
        # 'font.weight'       : 'bold',
        'ax.labelsize': fontsize * scale,
        'font.size': fontsize * scale,
        'legend.fontsize': fontsize * scale,
        'xtick.labelsize': fontsize * scale,
        'ytick.labelsize': fontsize * scale,
        # font.sans-serif     : Bitstream Vera Sans, ...
        # 'font.cmb10'     : 'cmb10',
        # 'font.family'         : 'cursive',
        'font.family': 'sans-serif',
        # 'font.sans-serif'   : 'Helvetica',
        'text.usetex': usetex,
        'figure.figsize': fig_size,
        'xtick.major.pad': 4 * scale,
        'xtick.minor.pad': 4 * scale,
        'ytick.major.pad': 4 * scale,
        'ytick.minor.pad': 4 * scale,
        'xtick.major.size': 4 * scale,  # major tick size in points
        'xtick.minor.size': 2 * scale,  # minor tick size in points
        'ytick.major.size': 4 * scale,  # major tick size in points
        'ytick.minor.size': 2 * scale,  # minor tick size in points
        'lines.markersize': 6 * scale,
        'lines.linewidth': 0.6 * scale
    }
    plt.rcParams.update(params)


def setPlotStuff(fontsize=7, dpi=None):
    """TODO merge with setOutputStyle.

    Change ugly name.
    """
    from matplotlib import rcParams

    # print(rcParams.keys())

    # rcParams['ax.labelsize'] = fontsize  # REMOVED IN MPL.1.5
    # rcParams['ax.titlesize'] = fontsize  # REMOVED IN MPL.1.5
    # rcParams['ax.linewidth'] = 0.3  # REMOVED IN MPL.1.5
    rcParams['font.size'] = fontsize
    rcParams['xtick.labelsize'] = fontsize
    rcParams['ytick.labelsize'] = fontsize
    rcParams['legend.fontsize'] = fontsize
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['Helvetica']  # ['Times New Roman']
    rcParams['text.usetex'] = False

    #    rcParams['figure.figsize'] = 7.3, 4.2
    rcParams['xtick.major.size'] = 3
    rcParams['xtick.major.width'] = 0.3
    rcParams['xtick.minor.size'] = 1.5
    rcParams['xtick.minor.width'] = 0.3
    rcParams['ytick.major.size'] = rcParams['xtick.major.size']
    rcParams['ytick.major.width'] = rcParams['xtick.major.width']
    rcParams['ytick.minor.size'] = rcParams['xtick.minor.size']
    rcParams['ytick.minor.width'] = rcParams['xtick.minor.width']

    if dpi is not None:
        rcParams['figure.dpi'] = dpi
        rcParams['savefig.dpi'] = dpi


def createAnimation(fig, animate, nFrames, dpi, out):
    """Create animation for the content of a given matplotlib figure.

    Until I know a better place.
    """
    anim = animation.FuncAnimation(fig, animate, frames=nFrames,
                                   interval=0.001, repeat=False)
    anim.save(out + ".mp4", writer=None, fps=20, dpi=dpi, codec=None,
              bitrate=24 * 1024, extra_args=None, metadata=None,
              extra_anim=None, savefig_kwargs=None)
    try:
        print("Create frames ... ")
        os.system('mkdir -p anim-' + out)
        os.system('ffmpeg -i ' + out + '.mp4 anim-' + out + '/movie%d.jpg')
    except BaseException as _:
        pass


def saveAnimation(mesh, data, out, vData=None, plc=None, label='', cMin=None,
                  cMax=None, logScale=False, cmap=None, **kwargs):
    """Create and save an animation for a given mesh with a set of field data.

    Until I know a better place.
    """
    dpi = 92
    scale = 1
    fig = plt.figure(facecolor='white',
                     figsize=(scale * 800 / dpi, scale * 490 / dpi), dpi=dpi)
    ax = fig.add_subplot(1, 1, 1)

    gci = pg.mplviewer.drawModel(ax, mesh, data=data[0], cMin=cMin, cMax=cMax,
                                 cmap=cmap, logScale=logScale)

    pg.mplviewer.createColorbar(gci, label=label, pad=0.55)

    if plc:
        pg.show(plc, ax=ax)

    adjustWorldAxes(ax)

    def animate(i):
        """TODO WRITEME."""
        print(out + ": Frame:", i, "/", len(data))

        if vData is not None:
            ax.clear()
            pg.mplviewer.holdAxes_ = 1
            pg.mplviewer.drawModel(ax, mesh, data=data[i], cMin=cMin,
                                   cMax=cMax, cmap=cmap, logScale=logScale)
            pg.mplviewer.drawStreams(ax, mesh, vData[i], **kwargs)
        else:
            print(min(data[i]), max(data[i]))
            pg.mplviewer.setMappableData(gci, data[i], cMin=cMin, cMax=cMax,
                                         logScale=logScale)

        plt.pause(0.001)

    createAnimation(fig, animate, int(len(data)), dpi, out)
