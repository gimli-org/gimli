# -*- coding: utf-8 -*-

try:
    import pygimli as g
    from pygimli.mplviewer import drawMesh, drawModel, drawField
except ImportError:
    sys.stderr.write('''ERROR: cannot import the library 'pygimli'. Ensure that pygimli is in your PYTHONPATH ''')
    sys.exit( 1 )

import pylab as P

def showMesh( mesh, data = None, showLater = False, *args, **kwargs):
    ''
    ' Syntactic sugar, short-cut to create axes and plot node or cell values '
    ''
    fig = P.figure()
    a = fig.add_subplot( 111 )
    
    if data is None:
        drawMesh( a, mesh )
    else:
        if len( data ) == mesh.cellCount():
            drawModel( a, mesh, data, *args, **kwargs )
        elif len( data ) == mesh.nodeCount():
            drawField( a, mesh, data, *args, **kwargs )

    a.set_aspect( 'equal')
        
    if not showLater:
        P.show()

        
    #fig.show()
    #fig.canvas.draw()
    return a
#def showMesh( ... )

