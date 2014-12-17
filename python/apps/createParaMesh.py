#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    This program is part of pygimli
    Visit http://www.resistivity.net for further information or the latest version.
"""

import sys  # , traceback
import getopt
from os import system

try:
    import pygimli as pg
except ImportError:
    sys.stderr.write('ERROR: cannot import the library pygimli. ' +
                     'Ensure that pygimli is in your PYTHONPATH.')
    sys.exit(1)

import matplotlib.pyplot as plt
from pygimli.meshtools import createParaMesh2dGrid
from pygimli.viewer import showMesh
from pygimli.mplviewer import drawModel, drawSensors


def main(argv):
    import argparse

    description = "Create GIMLi parameter mesh for a given datafile"
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument("-v", "--verbose", dest="verbose", action="store_true",
                        help="be verbose", default=False)
    parser.add_argument("-o", "--output", dest="outFileName", metavar="File",
                        help="File name for mesh. Default [datafile].bms")
    parser.add_argument("-d", "--dimension", dest="dimension", type=int,
                        help="Dimension of the mesh", default="2")
    parser.add_argument("--paraDX", dest="paraDX", type=float, default="1",
                        help="Horizontal distance between sensors, " +
                        "relative regarding sensor distance")
    parser.add_argument("--paraDepth", dest="paraDepth", type=float,
                        help="Maximum depth [m] for parametric domain, " +
                        "0[default] means 0.4 * max sensor range", default="0")
    parser.add_argument("--grid", dest="grid",
                        action="store_true", default=False,
                        help="Use a regular grid for parameter domain")
    parser.add_argument("--paraDZ", dest="paraDZ", type=float, default="1",
                        help="Vertical distance to the first depth layer, " +
                        "relative regarding sensor distance. For --grid only")
    parser.add_argument("--nLayers", dest="nLayers", type=int, default="11",
                        help="Number of depth layers [default=11] for --grid")
    parser.add_argument("--paraBoundary", dest="paraBoundary", type = float,
                        help="Parameter domain boundary in profile distance",
                        default="2")
    parser.add_argument("--show", dest="show", action="store_true",
                        default=False, help="Show the mesh in new window.")

    parser.add_argument('dataFileName')

    options = parser.parse_args()

    # create default output filename
    if options.outFileName is None:
        options.outFileName = options.dataFileName.replace('.dat', '')

    if options.verbose:
        print(options)

    sensorPositions = pg.DataContainer(options.dataFileName).sensorPositions()

    if options.verbose:
        print(("Sensorcount = ", len(sensorPositions)))

    if options.dimension == 3:
        raise Exception

    swatch = pg.Stopwatch(True)
    if options.dimension == 2:

        mesh = createParaMesh2dGrid(sensors=sensorPositions,
                                    paraDX=options.paraDX,
                                    paraDZ=options.paraDZ,
                                    paraDepth=options.paraDepth,
                                    nLayers=options.nLayers,
                                    boundary=-1,
                                    paraBoundary=options.paraBoundary,
                                    verbose=options.verbose, quality=30.0,
                                    smooth=True)

    if options.verbose:
        print(("generation takes ", swatch.duration(), " s"))
        print( mesh )
        print(("writing " + options.outFileName + ".bms"))

    mesh.save(options.outFileName)

    if options.show:
        ax = showMesh(mesh, showLater=True)
        drawModel(ax, mesh, mesh.cellMarker(), alpha=0.2)
        xl = ax.get_xlim()
        yl = ax.get_ylim()
        dx = (xl[1] - xl[0]) * 0.01
        dy = (yl[1] - yl[0]) * 0.01

        ax.set_xlim(xl[0] - dx, xl[1] + dx)
        ax.set_ylim(yl[0] - dy, yl[1] + dy)

        drawSensors(ax, sensorPositions)

        plt.show()

if __name__ == "__main__":
    main(sys.argv[1:])
