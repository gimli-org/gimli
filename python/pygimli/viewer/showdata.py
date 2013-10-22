# -*- coding: iso-8859-1 -*-

try:
    import pygimli as g
    import pygimli.mplviewer
except ImportError:
    import sys
    sys.stderr.write('''ERROR: cannot import the library 'pygimli'. Ensure that pygimli is in your PYTHONPATH ''')
    sys.exit(1)

import pylab

def showData(data, **kwargs):
    """
        Syntactic sugar for drawData. Creates figure, axis and Show.
    """
    MOVETOBERT
    fig = pylab.figure()
    a = fig.add_subplot(111)
    drawData(a, data, data('rhoa'), **kwargs)
    a.set_aspect('equal')
    pylab.show()
    return a

def drawData(axes, data, vals, pseudotype='A_M', cMin=None, cMax=None,
             showCbar=True, linear=False, label=""):
    """
    """
    MOVETOBERT
    try:
        sheme = getattr(g.mplviewer.Pseudotype, pseudotype)
    except:
        print "no pseudotype ", pseudotype, " found. Falling back to A_M"
        print dir(g.mplviewer.Pseudotype)
        sheme = g.mplviewer.Pseudotype.A_M

    gci = g.mplviewer.drawDataAsMatrix(axes, data, vals, pseudotype = sheme
                                        , logScale = not linear)

    if showCbar and (cMin is not cMax):
        g.mplviewer.createColorbar(gci, cMin=cMin, cMax=cMax,
                                   nLevs=5, label=label)

    return gci