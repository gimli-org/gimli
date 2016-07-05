# -*- coding: utf-8 -*-
"""Viewer interface .. depends on matplotlib."""

# are the following is suitable for a drawing package?
from .boreholes import BoreHole, BoreHoles, create_legend
from .colorbar import (addCoverageAlpha, autolevel, cmapFromName,
                       createColorbar, findAndMaskBestClim, setCbarLevels,
                       setMappableData)
# TODO example scripts for the following and refactor is needed
# maybe ploter should named show or draw
from .dataview import (drawSensorAsMarker,  # dups to meshview??
                       generateMatrix,
                       patchMatrix,
                       patchValMap,
                       plotDataContainerAsMatrix,
                       plotMatrix,
                       plotVecMatrix)

from .meshview import (CellBrowser,
                       createMeshPatches,
                       createParameterContraintsLines,
                       createTriangles,
                       draw1DColumn,
                       drawField,
                       drawMesh,
                       drawMeshBoundaries,
                       drawModel,
                       drawMPLTri,
                       drawParameterConstraints,
                       drawPLC,
                       drawSelectedMeshBoundaries,
                       drawSelectedMeshBoundariesShadow,
                       drawSensors,
                       drawStreamLines,
                       drawStreams,
                       insertUnitAtNextLastTick,
                       plotLines)

from .overlayimage import (cacheFileName, deg2MapTile, getMapTile, mapTile2deg,
                           underlayMap)
from .utils import (hold,
                    wait,
                    updateAxes,
                    createAnimation,
                    saveAnimation,
                    setOutputStyle,
                    setPlotStuff)

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
