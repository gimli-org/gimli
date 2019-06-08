# -*- coding: utf-8 -*-
"""Matplotlib drawing functions used by `pygimli.viewer`."""

# are the following suitable for a drawing package?
from .utils import (hold,
                    wait,
                    updateFig,
                    updateAxes,
                    renameDepthTicks,
                    insertUnitAtNextLastTick,
                    saveFigure,
                    saveAxes,
                    adjustWorldAxes,
                    createAnimation,
                    saveAnimation,
                    setOutputStyle,
                    setPlotStuff,
                    plotLines,
                    twin, createTwinX, createTwinY)

from .boreholes import BoreHole, BoreHoles, create_legend

from .colorbar import (createColorBar,
                       createColorBarOnly,
                       findColorBar,
                       updateColorBar,
                       addCoverageAlpha,
                       autolevel,
                       cmapFromName,
                       findAndMaskBestClim,
                       setCbarLevels,
                       setMappableValues,
                       setMappableData)

from .meshview import (CellBrowser,
                       createMeshPatches,
                       createTriangles,
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
                       drawStreams)

from .overlayimage import (cacheFileName,
                           deg2MapTile,
                           getMapTile,
                           mapTile2deg,
                           underlayMap,
                           underlayBKGMap)

# TODO example scripts for the following and refactor is needed
# maybe ploter should named show or draw
from .dataview import (drawSensorAsMarker,  # dups to meshview??
                       showVecMatrix,
                       drawVecMatrix,
                       showValMapPatches,
                       drawValMapPatches,
                       showMatrix,
                       drawMatrix,
                       generateMatrix,
                       patchMatrix, # deprectated (Naming)
                       patchValMap, # deprectated (Naming)
                       plotDataContainerAsMatrix, # need renaming
                       plotMatrix, # deprectated (Naming)
                       plotVecMatrix,# deprectated (Naming)
                       )

# which of these do we actually need?
from .modelview import (drawModel1D,
                        showmymatrix,  # needed ?
                        draw1DColumn, # needed or redundant ?
                        draw1dmodel,   # needed or redundant ?
                        show1dmodel,  # needed or redundant ?
                        draw1dmodelErr,  # needed or redundant ?
                        draw1dmodelLU,  # needed or redundant ?
                        showStitchedModels,
                        showfdemsounding)

__all__ = [
    "BoreHole", "BoreHoles", "create_legend", "addCoverageAlpha", "autolevel",
    "cmapFromName",
    "createColorBar", "createColorBarOnly", "findColorBar", "updateColorBar",
    "findAndMaskBestClim", "setCbarLevels", "saveFigure", "saveAxes",
    "setMappableData", "drawSensorAsMarker", "generateMatrix", "patchMatrix",
    "patchValMap", "plotDataContainerAsMatrix", "plotMatrix", "plotVecMatrix",
    "CellBrowser", "createMeshPatches",
    "createTriangles", "draw1DColumn", "drawField", "drawMesh",
    "drawMeshBoundaries", "drawModel", "drawMPLTri", "hold", "wait",
    "setOutputStyle", "setPlotStuff", "createAnimation", "saveAnimation",
    "drawParameterConstraints", "drawPLC", "drawSelectedMeshBoundaries",
    "drawSelectedMeshBoundariesShadow", "drawSensors", "drawStreamLines",
    "drawStreams", "insertUnitAtNextLastTick", "plotLines", "cacheFileName",
    "deg2MapTile", "getMapTile", "mapTile2deg", "underlayMap", "updateAxes"
]

# plt.subplots() resets locale setting to system default .. this went
# horrible wrong for german 'decimal_point': ','
# https://github.com/matplotlib/matplotlib/issues/6706
# Qt5Agg resets it after importing figure;
# TkAgg resets it after importing pyplot.
# until its fixed we should maybe silently initialize the qt5agg backend and
# refix the locale afterwards. If someone have a plan to do.
#checkAndFixLocaleDecimal_point(verbose=True)


def createColorbar(*args, **kwargs):
    print("createColorbar is DEPRECATED .. please use createColorBar instead.")
    return createColorBar(*args, **kwargs)
