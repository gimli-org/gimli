#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    This program is part of pygimli
    See http://www.resistivity.net/gimli for further information.
"""

import sys

import numpy as np

try:
    import pygimli as pg
except ImportError:
    sys.stderr.write('ERROR: cannot import the library pygimli.' +
                     ' Ensure that pygimli is in your PYTHONPATH')
    sys.exit(1)


def merge(data1, data2, ContainerType=pg.DataContainer, snap=0.001):
    '''
        Merge two datacontainers into one by copying the sensor positions
        and datapoints from data2 into data1.\n
        Double sensor positions will be unified and snapped to a grid.
    '''
    data = ContainerType(data1)
    data.add(data2, snap)
    return data
# def merge( ... )


def loadProjectFile(projectfile, ContainerType=pg.DataContainer, verbose=False):
    """
        A project file defines how multiple data files are imported and merged.
        The currently supported formats are:\n\n

        A list of multiple row entries with the following formats:

        fileName
        fileName interpolationfile
        fileName startx starty endx endy
        fileName x1 y1 x2 y2 x3 y3 ...

        You can comment out a row by adding '#'.

        interpolationfile is a 3-column-ascii-file (dx x y)
    """
    dataList = []

    fi = open(projectfile, "r")
    content = fi.readlines()
    fi.close()

    for c in content:
        row = c.split('\n')[0].split()
        d = None

        if len(row) > 0 and row[0] != '#':
            if len(row) == 1:  # filename only
                d = ContainerType(row[0])
            elif len(row) == 2:
                # Thomas?? ist das von dir?? was macht das
                # kommt mir nur vage bekannt vor, benutzt aber extra-File
                d = ContainerType(row[0])

                xn = pg.x(d.sensorPositions())
                zn = pg.z(d.sensorPositions())

                tape, xt, yt = np.loadtxt(row[1], unpack=True)

                x3n = np.interp(xn, tape, xt)
                y3n = np.interp(xn, tape, yt)

                for i in range(d.sensorCount()):
                    d.setSensorPosition(i, pg.RVector3(x3n[i], y3n[i], zn[i]))

            elif len(row) == 3:  # filename xstart xend
                d = ContainerType(row[0])

                start = pg.RVector3(float(row[1]), 0.0)
                end = pg.RVector3(float(row[2]), 0.0)

                for i in range(d.sensorCount()):
                    pos = start + float(i)*(end-start)/(d.sensorCount() - 1.)
                    d.setSensorPosition(i, pos)

            elif len(row) == 5:  # filename xstart ystart xend yend
                d = ContainerType(row[0])

                start = pg.RVector3(float(row[1]), float(row[2]))
                end = pg.RVector3(float(row[3]), float(row[4]))

                for i in range(d.sensorCount()):
                    pos = start + float(i)*(end-start)/(d.sensorCount() - 1.)
                    d.setSensorPosition(i, pos)
            elif len(row) == 7:
                # filename xstart ystart zstart xend yend zend
                d = ContainerType(row[0])

                start = pg.RVector3(float(row[1]), float(row[2]), float(row[3]))
                end = pg.RVector3(float(row[4]), float(row[5]), float(row[6]))

                for i in range(d.sensorCount()):
                    pos = start + float(i)*(end-start)/(d.sensorCount() - 1.)
                    d.setSensorPosition(i, pos)
            else:
                print(("cannot interprete project format: len(row) = ",
                       len(row)))
                return dataList

            dataList.append(d)

            if verbose:
                print(("append: ", d))
                print(("from:" , d.sensorPositions()[ 0 ],
                       "to:", d.sensorPositions()[ -1 ]))

    return dataList

def usage( exitcode = 0, comment = '' ):
    print(comment)
    print(__doc__)
    exit( exitcode )

def main(argv):

    # OptionParser is deprecated since python 2.7, use argparse
    from optparse import OptionParser

    parser = OptionParser("usage: %prog [options] project-file",
                          version="%prog: " + g.__version__,
                          description = loadProjectFile.__doc__ + \
                          '\n The import data function provides the following data formats:\n'
                          ) #+ importData.__doc__)

    parser.add_option("-v", "--verbose", dest="verbose", action = "store_true", default = False
                            , help="Be verbose." )
    parser.add_option("-o", "--output", dest = "outFileName", metavar = "File", default = None
                            , help = "Filename for the resulting data file." )
    parser.add_option("-s", "--snap", dest = "snap", type = "float", default = 0.1
                            , help = "Snap coordinates to gridsize" )
    parser.add_option("-B", "--bert", dest = "bert", action = "store_true", default = False
                            , help = "Use BERT sensor indices (a b m n)" )

    (options, args) = parser.parse_args()

    projectFileName = None

    if len( args ) == 0:
        parser.print_help()
        print("Please add a project-file.")
        sys.exit( 2 )
    else:
        projectFileName = args[ 0 ];

    if options.outFileName is None:
        options.outFileName = projectFileName[0:projectFileName.find('.pro')] + '.dat'

    if options.verbose:
        print((options, args))
        print(("verbose =", options.verbose))
        print(("project =", projectFileName))
        print(("output =", options.outFileName))
        print(("snap =", options.snap))

    ContainerTyp = g.DataContainer

    if options.bert:
        import pybert as b
        ContainerTyp = b.DataContainerERT

    dataList = loadProjectFile( projectFileName, ContainerTyp, verbose = options.verbose )

    outdata  = dataList[ 0 ]

    if options.verbose:
        print("start merging ...")
        print(outdata)

    for d in dataList[1:]:
        outdata = merge( outdata, d, ContainerTyp, options.snap )

        if options.verbose:
            print(outdata)

    if options.verbose:
        print(("Write file to: ", options.outFileName))

    outdata.save( options.outFileName )

if __name__ == "__main__":
    main( sys.argv[ 1: ] )
