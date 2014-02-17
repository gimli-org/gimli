#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    This program is part of pygimli
    Visit http://www.resistivity.net for further information or the latest version.
"""

import sys
import getopt
from os import system

import pygimli as g
import pybert as b

def exportSens( meshfile, sensMatrix, dataID, tol = 1e-5, save = False):
    mesh = g.Mesh( meshfile )
    mesh.showInfos()

    S = g.RMatrix();
    g.load( S, sensMatrix )
    print S.rows()
    print S.cols()

    if dataID == -1:
        for i in range( 0, S.rows() ):
            d = b.prepExportSensitivityData( mesh , S[ i ], tol )
            name = "sens-" + str( i )
            if save:
                print name
                g.save( d, name + ".vec" )
            mesh.addExportData( name, d )
    else:
        d = b.prepExportSensitivityData( mesh , S[ dataID ], tol )
        name = "sens-" + str( dataID )
        if save:
            g.save( d, name + ".vec" )
        mesh.addExportData( name, d )

    mesh.exportVTK("sens");

def main( argv ):
    from optparse import OptionParser

    parser = OptionParser( "usage: %prog [options] mesh"
                            , version="%prog: " + g.__version__
                            )
    parser.add_option("-v", "--verbose", dest="verbose", action = "store_true", default = False
                            , help="Be verbose." )
    parser.add_option("-s", "--sensMatrix", dest = "sensMatrix", metavar = "File", default = 'sens.bmat'
                            , help = "Sensitivity matrix" )
    parser.add_option("-i", "--dataId", dest = "id", type = int, default = -1
                            , help = "Export the data for the given id. -1 for all data." )

    parser.add_option("-t", "--tolerance", dest = "tolerance", type = float, default = 1e-5
                            , help = "Tolerance for Data preparation" )
    parser.add_option("", "--save", dest = "save", action = "store_true", default = False
                            , help = "Save single sensitivity vector" )

    (options, args) = parser.parse_args()

    if len( args ) == 0:
        parser.print_help()
        print "Please add a mesh-file."
        sys.exit( 2 )
    else:
        meshfile = args[ 0 ];

    if ( options.verbose ):
        print "meshfile =", meshfile
        print "sensMatrix =", options.sensMatrix
        print "ith.", options.id
        print "tol=",  options.tolerance
        print "save=",  options.save


    exportSens( meshfile, options.sensMatrix, options.id, options.tolerance, options.save )

if __name__ == "__main__":
    main( sys.argv[ 1: ] )
