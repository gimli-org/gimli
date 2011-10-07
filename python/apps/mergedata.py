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

#from pybertlib.importer import importData

def merge2data( Data1, Data2, sensorTokens = ['a','b','m','n','s','g'] ):
    ''' merge two data containers '''
    # check if sensors are identical,
    # if not build unified sensor list and remap sensor indices
    senseq = ( Data1.sensorCount() == Data2.sensorCount() )
    #    if senseq:
        
    if not senseq:
        return None
    
    # now the sensors are the same and we just have to paste the data vectors
    Data = g.DataContainer( )
    for tok in sensorTokens:
        if Data1.dataMap().has_key( tok ):
            Data.registerSensorIndex( tok )
        
    Data.setSensorPositions( Data1.sensorPositions() )
    Data.resize( Data1.size() + Data2.size() )
    
    for key in Data1.dataMap().keys():
        if Data2.dataMap().has_key( key ):
            Data.set( key, g.cat( Data1.get( key ), Data2.get( key ) ) )
                    
    return Data

def merge( N1, N2, snap = 0.1 ):
    '''
        Merge two datacontainer into one, by copying the sensor positions and datapoints from N2 into N1.\n
        Double sensor positions will be reused and snapped to a grid with gridsize snap
    '''
    data = g.DataContainer( N1 );
    data.resize( N1.size() + N2.size() )

    newSensors = range( N2.sensorCount() )
    
    for e in N2.electrodes():
        newSensors[ e.id() ] = data.createSensor( e.pos(), snap ).id();
        #print e.id(), newElecs[ e.id() ]
            
    for i in range( 0, N2.size() ):
        if N2( i ).a() > -1:
            eaID = newSensors[ N2( i ).a() ]
        else:
            eaID = -1

        if N2( i ).b() > -1:
            ebID = newSensors[ N2( i ).b() ]
        else:
            ebID = -1;

        if N2( i ).m() > -1:
            emID = newSensors[ N2( i ).m() ]
        else:
            emID = -1;

        if N2( i ).n() > -1:
            enID = newSensors[ N2( i ).n() ]
        else:
            enID = -1

        newid = N1.size() + i

        data.createFourPointData( newid, eaID, ebID, emID, enID );

    for key, val in N1.dataMap():
        if key is not 'a' and key is not 'b' and key is not 'm' and key is not 'n':
            data.set( key,    g.cat( N1( key ),    N2( key ) ) )

    return data
# def merge( ... )
    
def loadProjectFile( projectfile, verbose = False ):
    '''
        A project file defines how multiple data files are imported and merged.
        The currently supported formats are:                                                
        A list of multiple row entries with the following formats:                                                  
                       
        fileName1                                                                                                               
        fileName2 startx starty endx endy                                                                   
                                                                                                                                                                                          
        You can comment out a row by adding '#'.
    '''
    dataList = []

    fi = open( projectfile, "r" )
    content = fi.readlines( )
    fi.close();

    for c in content:
        row = c.split( '\n' )[0].split()
        d = None
        if row[ 0 ] != '#' :
            if len( row ) == 1:
                ### filename only
                d = importData( row[0] )
            elif len( row ) == 5:
                ### filename xstart ystart xend yend
                d = importData( row[0] )
                start = g.RVector3( float( row[1] ), float( row[2] ) )
                end   = g.RVector3( float( row[3] ), float( row[4] ) ) 
                
                for i in range( d.electrodeCount() ):
                    d.setPosition( i, start + float(i)*(end-start)/(d.electrodeCount()-1.) )
            else:
                print "cannot interprete project format: len(row) = ", len( row )
                return dataList
            
            dataList.append( d )
            
            if verbose: 
                print "append: ", d
                print "from:" , d.electrodePositions()[ 0 ], "to:", d.electrodePositions()[ -1 ]
            
    return dataList

def usage( exitcode = 0, comment = '' ):
    print comment
    print __doc__
    exit( exitcode )

def main( argv ):
    
    # OptionParser is deprecated since python 2.7, use argparse 
    from optparse import OptionParser

    parser = OptionParser( "usage: %prog [options] project file"
                            , version="%prog: " + g.versionStr() 
                            , description = loadProjectFile.__doc__ + 
                            '                                            The import data function provide the following data formats:                  ' + importData.__doc__
                            )
    parser.add_option("-v", "--verbose", dest="verbose", action = "store_true", default = False
                            , help="Be verbose." )
    parser.add_option("-o", "--output", dest = "outFileName", metavar = "File", default = None
                            , help = "Filename for the resulting data file." )
    parser.add_option("-s", "--snap", dest = "snap", type = "float", default = 0.1 
                            , help = "Snap coordinates to gridsize" )
                            
    (options, args) = parser.parse_args()

    projectFileName = None

    if len( args ) == 0:
        parser.print_help()
        print "Please add a project-file."
        sys.exit( 2 )
    else:
        projectFileName = args[ 0 ];

    if options.outFileName is None:
        options.outFileName = projectFileName[0:projectFileName.find('.pro')] + '.dat'
          
    if options.verbose:
        print options, args
        print "verbose =", options.verbose
        print "project =", projectFileName
        print "output =", options.outFileName
        print "snap =", options.snap

    dataList = loadProjectFile( projectFileName, verbose = options.verbose )

    outdata  = g.DataContainer()

    for d in dataList:
        outdata = merge( outdata, d, options.snap )
        if options.verbose:
            print outdata

    if options.verbose:
        print "outdata:", options.outFileName
        
    outdata.save( options.outFileName, 'a b m n r rhoa err ip' )

if __name__ == "__main__":
    main( sys.argv[ 1: ] )