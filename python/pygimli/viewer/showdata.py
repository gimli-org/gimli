# -*- coding: iso-8859-1 -*-

try:
    import pygimli as pg
    import pygimli.mplviewer
except ImportError:
    import sys
    sys.stderr.write('''ERROR: cannot import the library 'pygimli'. Ensure that pygimli is in your PYTHONPATH ''')
    sys.exit(1)

import matplotlib.pyplot as plt

def showData(data, **kwargs):
    """
    Syntactic sugar for drawData.

    Creates figure, axis and Show.
    """
    MOVETOBERT
    fig = plt.figure()
    a = fig.add_subplot(111)
    drawData(a, data, data('rhoa'), **kwargs)
    a.set_aspect('equal')
    plt.show()
    return a

def drawData(axes, data, vals, pseudotype='A_M', cMin=None, cMax=None,
             showCbar=True, linear=False, label=""):
    """"""
    MOVETOBERT
    try:
        sheme = getattr(pg.mplviewer.Pseudotype, pseudotype)
    except:
        print("no pseudotype ", pseudotype, " found. Falling back to A_M")
        print(dir(pg.mplviewer.Pseudotype))
        sheme = pg.mplviewer.Pseudotype.A_M

    gci = pg.mplviewer.drawDataAsMatrix(axes, data, vals, pseudotype = sheme
                                        , logScale = not linear)

    if showCbar and (cMin is not cMax):
        pg.mplviewer.createColorbar(gci, cMin=cMin, cMax=cMax,
                                   nLevs=5, label=label)

    return gci