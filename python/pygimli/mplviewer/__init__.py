# -*- coding: utf-8 -*-
"""Viewer interface .. depends on matplotlib."""

import os

import matplotlib.animation as animation
import matplotlib.pyplot as plt

import pygimli as pg

# are the following is suitable for a drawing package?
from .boreholes import BoreHole, BoreHoles, create_legend
from .colorbar import (addCoverageAlpha, autolevel, cmapFromName,
                       createColorbar, findAndMaskBestClim, setCbarLevels,
                       setMappableData)
from .dataview import (drawSensorAsMarker, generateMatrix, patchMatrix,
                       patchValMap, plotDataContainerAsMatrix, plotMatrix,
                       plotVecMatrix)
from .meshview import (CellBrowser, createMeshPatches,
                       createParameterContraintsLines, createTriangles,
                       draw1DColumn, drawField, drawMesh, drawMeshBoundaries,
                       drawModel, drawMPLTri, drawParameterConstraints,
                       drawPLC, drawSelectedMeshBoundaries,
                       drawSelectedMeshBoundariesShadow, drawSensors,
                       drawStreamLines, drawStreams, insertUnitAtNextLastTick,
                       plotLines)
from .overlayimage import (cacheFileName, deg2MapTile, getMapTile, mapTile2deg,
                           underlayMap)
from .utils import (createAnimation, hold, saveAnimation, setOutputStyle,
                    setPlotStuff, updateAxes, wait)

__all__ = [
    "BoreHole", "BoreHoles", "create_legend", "addCoverageAlpha", "autolevel",
    "cmapFromName", "createColorbar", "findAndMaskBestClim", "setCbarLevels",
    "setMappableData", "drawSensorAsMarker", "generateMatrix", "patchMatrix",
    "patchValMap", "plotDataContainerAsMatrix", "plotMatrix", "plotVecMatrix",
    "CellBrowser", "createMeshPatches", "createParameterContraintsLines",
    "createTriangles", "draw1DColumn", "drawField", "drawMesh",
    "drawMeshBoundaries", "drawModel", "drawMPLTri", "hold", "wait",
    "setOutputStyle", "setPlotStuff", "createAnimation", "saveAnimation",
    "drawParameterConstraints", "drawPLC", "drawSelectedMeshBoundaries",
    "drawSelectedMeshBoundariesShadow", "drawSensors", "drawStreamLines",
    "drawStreams", "insertUnitAtNextLastTick", "plotLines", "cacheFileName",
    "deg2MapTile", "getMapTile", "mapTile2deg", "underlayMap", "updateAxes"
]
