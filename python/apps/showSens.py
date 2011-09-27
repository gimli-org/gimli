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

def exportSens( meshfile, sensMatrix, dataID ):
    mesh = g.Mesh( meshfile )
    mesh.showInfos()

    S = g.RMatrix();
    g.load( S, sensMatrix )
    print S.rows()
    print S.cols()
    
    if dataID == -1:
        for i in range( 0, S.rows() ):
            mesh.addExportData( "sens-" + str(i), g.prepExportSensitivityData( mesh , S[ i ], 1e-6 ) )
    else:
        mesh.addExportData( "sens-" + str( dataID ), g.prepExportSensitivityData( mesh , S[ dataID ], 1e-6 ) )
        
    mesh.exportVTK("sens");

def main( argv ):
    from optparse import OptionParser

    parser = OptionParser( "usage: %prog [options] mesh"
                            , version="%prog: " + g.versionStr() 
                            )
    parser.add_option("-v", "--verbose", dest="verbose", action = "store_true", default = False
                            , help="Be verbose." )
    parser.add_option("-s", "--sensMatrix", dest = "sensMatrix", metavar = "File", default = 'sens.bmat'
                            , help = "Sensitivity matrix" )
    parser.add_option("-i", "--dataId", dest = "id", type = int, default = 0
                            , help = "Export the data for the given id. -1 for all data." )

    (options, args) = parser.parse_args()

    if len( args ) == 0:
        parser.print_help()
        print "Please add a mesh-file."
        sys.exit( 2 )
    else:
        meshfile = args[ 0 ];

    id = 0
    
    if ( options.verbose ):
        print "meshfile =", meshfile
        print "sensMatrix =", options.sensMatrix
        print "ith.", options.id
        
    exportSens( meshfile, options.sensMatrix, options.id )
            
if __name__ == "__main__":
    main( sys.argv[ 1: ] )