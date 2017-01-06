# -*- coding: utf-8 -*-
"""Mesh based data transformation and mapping, interpolation, extrapolation."""

import numpy as np
import pygimli as pg


def nodeDataToCellData(mesh, data):
    """Convert node data to cell data.

    Convert node data to cell data via interpolation to cell centers.

    Parameters
    ----------
    mesh : :gimliapi:`GIMLI::Mesh`
        2D or 3D GIMLi mesh

    data : iterable [float]
        Data of len mesh.nodeCount().
        TODO complex, R3Vector, ndarray

    Examples
    --------
    """
    if len(data) != mesh.nodeCount():
        raise BaseException("Dimension mismatch, expecting nodeCount(): " +
                            str(mesh.nodeCount()) +
                            "got: " + str(len(data)), str(len(data[0])))

    return pg.interpolate(mesh, data, mesh.cellCenters())


def cellDataToNodeData(mesh, data, style='mean'):
    """Convert cell data to node data.

    Convert cell data to node data via non-weighted averaging (mean) of common
    cell data.

    Parameters
    ----------
    mesh : :gimliapi:`GIMLI::Mesh`
        2D or 3D GIMLi mesh

    data : iterable [float]
        Data of len mesh.nodeCount().
        TODO complex, R3Vector, ndarray

    style : str ['mean']
        Interpolation style.
        * 'mean' : non-weighted averaging
        TODO harmonic averaging
        TODO weighted averaging (mean, harmonic)
        TODO interpolation via cell centered mesh

    Examples
    --------
    """
    if len(data) != mesh.cellCount():
        raise BaseException("Dimension mismatch, expecting cellCount(): " +
                            str(mesh.cellCount()) +
                            "got: " + str(len(data)), str(len(data[0])))

    if style == 'mean':

        if mesh.dim() == 1:
            return pg.cellDataToPointData(mesh, data[0])
        elif mesh.dim() == 2:
            return np.array([pg.cellDataToPointData(mesh, data[:, 0]),
                            pg.cellDataToPointData(mesh, data[:, 1])])
        elif mesh.dim() == 3:
            return np.array([pg.cellDataToPointData(mesh, data[0]),
                            pg.cellDataToPointData(mesh, data[1]),
                            pg.cellDataToPointData(mesh, data[2])])
    else:
        raise BaseException("Style '" +style + "'not yet implemented."
                            "Currently styles available are: 'mean, '")


def nodeDataToBoundaryData(mesh, data):
    """
        Assuming [NodeCount, dim] data
        DOCUMENT_ME
    """

    if len(data) != mesh.nodeCount():
        raise BaseException("Dimension mismatch, expecting nodeCount(): " +
                            str(mesh.nodeCount()) +
                            " got: " + str(len(data)), str(len(data[0])))
    dim = len(data[0])
    ret = np.zeros((mesh.boundaryCount(), dim))
    if dim == 1:
        for b in mesh.boundaries():
            ret[b.id()] = b.pot(b.center(), data)  # / b.nodeCount()
    elif dim == 2:
        for b in mesh.boundaries():
            # v = b.data(b.center(), data)
            # interpolation is hell slow here .. check!!!!!!
            v2 = (data[b.node(0).id()] + data[b.node(1).id()])*0.5
            # print(v -v2)
            ret[b.id()] = [v2[0], v2[1]]
    else:
        for b in mesh.boundaries():
            ret[b.id()] = b.data(b.center(), data)  # / b.nodeCount()

    return ret


def cellDataToBoundaryData(mesh, data):
    """ TODO DOCUMENT_ME """
    if len(data) != mesh.cellCount():
        raise BaseException("Dimension mismatch, expecting cellCount(): " +
                            str(mesh.cellCount()) +
                            "got: " + str(len(data)), str(len(data[0])))

    CtB = mesh.cellToBoundaryInterpolation()

    if isinstance(data, pg.R3Vector()):
        return np.array([CtB * pg.x(data),
                         CtB * pg.y(data),
                         CtB * pg.z(data)]).T
    else:
        return CtB * data



# simplistics extrapolation
def fillEmptyToCellArray(mesh, vals, slope=True):
    """
    Prolongate empty cell values to complete cell attributes.

    It is possible that you have zero values that need to be filled with
    appropriate attributes. This function tries to fill the empty values
    successive prolongation of the non zeros.

    Parameters
    ----------
    mesh : :gimliapi:`GIMLI::Mesh`
        For each cell of mesh a value will be returned.

    vals : array
        Array of size cellCount().

    Returns
    -------
    atts : array
        Array of length mesh.cellCount()
    """
    atts = pg.RVector(mesh.cellCount(), 0.0)
    oldAtts = mesh.cellAttributes()
    mesh.setCellAttributes(vals)
    mesh.createNeighbourInfos()
    # std::vector< Cell * >
    # empties = []

    if slope:
        # search all cells with empty neighbours
        ids = pg.find(mesh.cellAttributes() != 0.0)

        for c in mesh.cells(ids):
            for i in range(c.neighbourCellCount()):
                nc = c.neighbourCell(i)

                if nc:
                    if nc.attribute() == 0.0:
                        # c.setAttribute(99999)

                        b = pg.findCommonBoundary(c, nc)
                        # search along a slope
                        pos = b.center() - b.norm()*1000.
                        sf = pg.RVector()
                        startCell = c

                        while startCell:

                            startCell.shape().isInside(pos, sf, False)
                            nextC = startCell.neighbourCell(sf)
                            if nextC:
                                if nextC.attribute() == 0.0:
                                    nextC.setAttribute(c.attribute())
                                else:
                                    break

                            startCell = nextC

    mesh.fillEmptyCells(mesh.findCellByAttribute(0.0), background=-1)
    atts = mesh.cellAttributes()
    mesh.setCellAttributes(oldAtts)
    return atts

if __name__ == "__main__":
    pass
