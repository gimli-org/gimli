# CODING=Utf-8
"""Plotting utilities used throughout the viewer.mpl package."""

import os
import atexit

import numpy as np

# TODO expensive import costs 75% of total time
# see: python -X importtime -c 'import pygimli'
# import matplotlib.animation as animation
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
    
import pygimli as pg
from pygimli.utils import prettyFloat, prettyTime

holdAxes__ = 0


def updateFig(fig, force=False, sleep=.0001):
    """For internal use."""
    if not holdAxes__:
        try:
            fig.canvas.draw_idle()
            if force:
                fig.canvas.flush_events()
                #fig.canvas.draw()
                #pg.plt.show(block=False)
                pg.plt.pause(sleep)
                #time.sleep(sleep)
        except BaseException as e:
            print(fig, e)
            pg.warn("Exception raised", e)


def updateAxes(ax, force=False):
    """For internal use."""
    updateFig(ax.figure, force=force)


def hold(val=1):
    """TODO WRITEME."""
    globals()[holdAxes__] = val


def waitOnExit():
    backend = matplotlib.get_backend()
    if not 'inline' in backend:
        if 'Qt' in backend or 'Wx' in backend:

            if len(plt.get_fignums()) > 0:
                pg.info('Showing pending widgets on exit. '
                        'Close all figures or Ctrl-C to quit the programm')
                pg.wait()

# this can't be changed after import
if pg.rc['waitOnExit'] is True:
    atexit.register(waitOnExit)


def wait(**kwargs):
    """TODO WRITEME."""
    # plt.pause seems to be broken in mpl:2.1
    # ax.canvas.draw_onIdle()
    updateAxes(plt.gca())
    kp = kwargs.pop('untilKeyPressed', False)
    if kp == True:
        plt.waitforbuttonpress(**kwargs)
    else:
        plt.show(**kwargs)


def saveFigure(fig, filename, pdfTrim=False):
    """Save figure as pdf."""
    if '.pdf' in filename:
        filename = filename[0:filename.find('.pdf')]
    fig.savefig(filename + '.pdf', bbox_inches='tight')
    # pdfTrim=1
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

    updateAxes(ax, force=True)
    saveFigure(ax.figure, filename)


def insertUnitAtNextLastTick(ax, unit, xlabel=True, position=-2):
    """Replace the last-but-one tick label by unit symbol."""
    if xlabel:
        labels = ax.get_xticklabels()
        labels[position] = unit
        ax.set_xticklabels(labels)
    else:
        labels = ax.get_yticklabels()
        labels[-position] = unit
        ax.set_yticklabels(labels)


def adjustWorldAxes(ax):
    """Set some common default properties for an axe."""
    ax.set_ylabel('Depth (m)')
    ax.set_xlabel('$x$ (m)')

    renameDepthTicks(ax)
    plt.tight_layout()
    updateAxes(ax)


def renameDepthTicks(ax):
    """Switch signs of depth ticks to be positive"""
    @ticker.FuncFormatter
    def major_formatter(x, pos):
        return prettyFloat(-x) % x

    ax.yaxis.set_major_formatter(major_formatter)
    updateAxes(ax)


def setPrettyTimeTicks(ax):
    """Set x axis ticks with pretty time."""
    @ticker.FuncFormatter
    def major_formatter(x, pos):
        return pg.utils.prettyTime(x) % x

    ax.xaxis.set_major_formatter(major_formatter)
    pg.viewer.mpl.updateAxes(ax)



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

    gci = pg.viewer.mpl.drawModel(ax, mesh, data=data[0], cMin=cMin, cMax=cMax,
                                  cMap=cmap, logScale=logScale)

    pg.viewer.mpl.createColorbar(gci, label=label, pad=0.55)

    if plc:
        pg.show(plc, ax=ax)

    adjustWorldAxes(ax)

    def animate(i):
        """TODO WRITEME."""
        print(out + ": Frame:", i, "/", len(data))

        if vData is not None:
            ax.clear()
            pg.viewer.mpl.holdAxes_ = 1
            pg.viewer.mpl.drawModel(ax, mesh, data=data[i], cMin=cMin,
                                   cMax=cMax, cMap=cmap, logScale=logScale)
            pg.viewer.mpl.drawStreams(ax, mesh, vData[i], **kwargs)
        else:
            print(min(data[i]), max(data[i]))
            pg.viewer.mpl.setMappableData(gci, data[i], cMin=cMin, cMax=cMax,
                                         logScale=logScale)

        plt.pause(0.001)

    createAnimation(fig, animate, int(len(data)), dpi, out)


def plotLines(ax, line_filename, linewidth=1.0, step=1):
    """Read lines from file and plot over model."""
    xz = np.loadtxt(line_filename)
    n_points = xz.shape[0]
    if step == 2:
        for i in range(0, n_points, step):
            x = xz[i:i + step, 0]
            z = xz[i:i + step, 1]
            ax.plot(x, z, 'k-', linewidth=linewidth)
    if step == 1:
        ax.plot(xz[:, 0], xz[:, 1], 'k-', linewidth=linewidth)


def twin(ax):
    """Return the twin of ax if exist."""
    for other_ax in ax.figure.axes:
        if other_ax is ax:
            continue
        if other_ax.bbox.bounds == ax.bbox.bounds:
            return other_ax
    return None


def createTwinX(ax):
    """Utility function to create (or return existing) twin x axes for ax."""
    return _createTwin(ax, 'twinx')


def createTwinY(ax):
    """Utility function to create (or return existing) twin x axes for ax."""
    return _createTwin(ax, 'twiny')


def _createTwin(ax, funct):
    """Utility function to create (or return existing) twin x axes for ax."""
    tax = None
    for other_ax in ax.figure.axes:
        if other_ax is ax:
            continue
        if other_ax.bbox.bounds == ax.bbox.bounds:
            tax = other_ax

    if tax is None:
        tax = getattr(ax, funct)()

    return tax
