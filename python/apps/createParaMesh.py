#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    This program is part of pygimli
    Visit http://www.resistivity.net for further information or the latest version.
"""

import sys, traceback
import getopt
from os import system

try:
    import pygimli as g
except ImportError:
    sys.stderr.write('''ERROR: cannot import the library 'pygimli'. Ensure that pygimli is in your PYTHONPATH ''')
    sys.exit( 1 )

    
from pygimli.meshtools import createParaMesh2dGrid
from pygimli.viewer import showMesh
from pygimli.mplviewer import drawModel, drawSensors
    
def main( argv ):
    
    import argparse
    
    parser = argparse.ArgumentParser( description = "Create GIMLi parameter mesh for a given datafile" )
        
    parser.add_argument("-v", "--verbose", dest = "verbose", action = "store_true"
                            , help = "Be verbose.", default = False )
    parser.add_argument("-o", "--output", dest = "outFileName"
                            , help = "File name for the resulting mesh. Default is [datafile.bms]"
                            , metavar = "File" )
    parser.add_argument("-d", "--dimension", dest = "dimension"
                            , help = "Dimension of the mesh to be created"
                            , type = int
                            , default = "2" )
    parser.add_argument("--paraDX", dest = "paraDX"
                            , help = "Horizontal distance between sensors, relative regarding sensor distance"
                            , type = float
                            , default = "1" )
    parser.add_argument("--paraDepth", dest = "paraDepth"
                            , help = "Maximum depth (absolute meter) for parametric domain, 0[default] means 0.4 * max sensor range"
                            , type = float
                            , default = "0" )
    parser.add_argument("--grid", dest = "grid"
                            , help = "Use a regular grid for parameter domain"
                            , action = "store_true"
                            , default = False )
    parser.add_argument("--paraDZ", dest = "paraDZ"
                            , help = "Vertical distance to the first depth layer, relative regarding sensor distance. For --grid only"
                            , type = float
                            , default = "1" )
    parser.add_argument("--nLayers", dest = "nLayers"
                            , help = "Number of depth layers [default = 11] For --grid only"
                            , type = int
                            , default = "11" )
                            
    parser.add_argument("--paraBoundary", dest = "paraBoundary"
                            , help = "Boundary offset for parameter domain in absolute sensor distance 2[default]"
                            , type = float
                            , default = "2" )
                            
    parser.add_argument( "--show", dest = "show"
                            , help = "Open a viewer window and show the mesh"
                            , action = "store_true"
                            , default = False )

    parser.add_argument( 'dataFileName' )

    options = parser.parse_args()
    
    # create default output filename
    if options.outFileName == None:
        options.outFileName = options.dataFileName.replace('.dat','') 
    
    if options.verbose: print( options )

    sensorPositions = g.DataContainer( options.dataFileName ).sensorPositions()
    
    if options.verbose: print( "Sensorcount = ", len( sensorPositions ) )
    
    if options.dimension == 3:
        raise( Exception( "not yet implemented" ) )
    
    swatch = g.Stopwatch( True )
    if options.dimension == 2:
        
        mesh = createParaMesh2dGrid( sensors = sensorPositions, 
                                     paraDX = options.paraDX,
                                     paraDZ = options.paraDZ,
                                     paraDepth = options.paraDepth,
                                     nLayers = options.nLayers,
                                     boundary = -1,       
                                     paraBoundary = options.paraBoundary, 
                                     verbose = options.verbose, quality = 30.0, smooth = True)
    
    if options.verbose:
        print( "generation takes ", swatch.duration(), " s" )
        print( mesh )
        print( "writing " + options.outFileName + '.bms' )
    
    mesh.save( options.outFileName )
    
    if options.show:
        ax = showMesh( mesh, showLater = True  )
        drawModel( ax, mesh, mesh.cellMarker(), alpha = 0.2 )
        xl = ax.get_xlim()
        yl = ax.get_ylim()        
        dx = ( xl[1] - xl[0] ) * 0.01
        dy = ( yl[1] - yl[0] ) * 0.01

        ax.set_xlim( xl[ 0 ] - dx, xl[ 1 ] + dx )
        ax.set_ylim( yl[ 0 ] - dy, yl[ 1 ] + dy )
        
        drawSensors( ax, sensorPositions )
        
        import pylab as P
        P.show() 
    
if __name__ == "__main__":
    main( sys.argv[ 1: ] )