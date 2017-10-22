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

    return pg.interpolate(mesh, data, destPos=mesh.cellCenters())


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

        if isinstance(data, pg.RVector):
            print(pg.cellDataToPointData(mesh, data))
            return pg.cellDataToPointData(mesh, data)

        if mesh.dim() == 1:
            return pg.cellDataToPointData(mesh, data)
        elif mesh.dim() == 2:
            return np.array([pg.cellDataToPointData(mesh, data[:, 0]),
                             pg.cellDataToPointData(mesh, data[:, 1])]).T
        elif mesh.dim() == 3:
            return np.array([pg.cellDataToPointData(mesh, data[0]),
                            pg.cellDataToPointData(mesh, data[1]),
                            pg.cellDataToPointData(mesh, data[2])])
    else:
        raise BaseException("Style '" + style + "'not yet implemented."
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

    if isinstance(data, pg.R3Vector):
        ret = np.zeros((mesh.boundaryCount(), 3))
        for b in mesh.boundaries():
            ret[b.id()] = sum(data[b.ids()]) / b.nodeCount()
        return ret

    dim = len(data[0])
    ret = np.zeros((mesh.boundaryCount(), dim))
    if dim == 1:
        for b in mesh.boundaries():
            ret[b.id()] = b.pot(b.center(), data)  # / b.nodeCount()
    elif dim == 2:
        for b in mesh.boundaries():
            if b.nodeCount() > 2:
                print(b)
                raise BaseException("implement me")

            ret[b.id()] = sum(data[b.ids()]) / b.nodeCount()
            continue
            # v = b.data(b.center(), data)
            # interpolation is hell slow here .. check!!!!!!
            v2 = (data[b.node(0).id()] + data[b.node(1).id()])*0.5
            # print(v -v2)
            ret[b.id()] = [v2[0], v2[1]]
    else:
        for b in mesh.boundaries():
            ret[b.id()] = sum(data[b.ids()]) / b.nodeCount()
            continue
            raise BaseException("don't use this until further checking")
            ret[b.id()] = b.vec(b.center(), data)  # / b.nodeCount()

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


def tapeMeasureToCoordinates(tape, pos):
    """Interpolate 2D tape measured topography to 2D Cartesian coordinates.

    Tape and pos value are expected to be sorted along distance to the orign.

    TODO optional smooth curve with harmfit

    Parameters
    ----------
    tape : [[x,z]] | [RVector3] | R3Vector
        List of tape measured topography points with measured distance (x)
        from origin and height (z)

    pos : iterable
        Query positions along the tape measured profile

    Returns
    -------
    res : ndarray(N, 2)
        Same as pos but with interpolated height values.
        The Distance between pos points and res (along curve) points remains.

    Examples
    --------
    >>> # no need to import matplotlib. pygimli's show does
    >>> import numpy as np
    >>> import pygimli as pg
    >>> import pygimli.meshtools as mt
    >>> elec = np.arange(11.)
    >>> topo = np.array([[0., 0.], [3., 2.], [4., 2.], [6., 1.], [10., 1.]])
    >>> _= pg.plt.plot(topo[:,0], topo[:,1])
    >>> p = mt.tapeMeasureToCoordinates(topo, elec)
    >>> pg.plt.gca().plot(p[:,0], p[:,1], 'o') #doctest: +ELLIPSIS
    [...]
    >>> pg.plt.gca().set_aspect(1)
    >>> pg.plt.show()
    """
    if isinstance(tape, pg.R3Vector) or isinstance(tape, pg.stdVectorRVector3):
        xTape = pg.x(tape)
        zTape = pg.z(tape)
    else:
        xTape = tape[:, 0]
        zTape = tape[:, 1]

    t = pg.utils.cumDist(pos)
    # print(t)
    tTape = pg.utils.cumDist(tape)
    xt = np.interp(t, tTape, xTape)
    zt = np.interp(t, tTape, zTape)

    pg.plt.plot(xTape, zTape)
    pg.plt.plot(xt, zt, 'o')
    pg.wait()
    return np.vstack([xt, zt]).T


def interpolate(*args, **kwargs):
    r"""Interpolation convinience function.

    Convinience function to interpolate different kind of data.
    Currently supported interpolation schemes are:

    * Mesh based values to arbitrary points, based on finite element
      interpolation (pg.core)

      Parameters:
        args: :gimliapi:`GIMLI::Mesh`, ...
            Arguments forwarded to :gimliapi:`GIMLI::interpolate`
        kwargs:
            Arguments forwarded to :gimliapi:`GIMLI::interpolate`
      Returns:
        Interpolated values

    * 1D point set :math:`u(x)` for ascending :math:`x`.
      Find interpolation function :math:`I = u(x)` and
      returns :math:`u_{\text{i}} = I(x_{\text{i}})`
      (interpolation methods are [**linear** via matplotlib,
      cubic **spline** via scipy, fit with **harmonic** functions' via pygimli])

      Parameters:
        args: xi, x, u
            * :math:`x_{\text{i}}` - target sample points
            * :math:`x` - function sample points
            * :math:`u` - function values
        kwargs:
            * method : string
                Specify interpolation method 'linear, 'spline', 'harmonic'
      Returns:
        ui: array of length xi
            :math:`u_{\text{i}} = I(x_{\text{i}})`, with :math:`I = u(x)`

    To use the core functions :gimliapi:`GIMLI::interpolate` start with a
    mesh instance as first argument or use the appropriate keyword arguments.

    TODO

    * 2D parametric to points (method=['linear, 'spline', 'harmonic'])
    * 2D/3D point cloud to points/grids (Delauney, 'linear, 'spline', 'harmonic')
    * Mesh to points based on nearest neighbour values (pg.core)

    Examples
    --------
    >>> # no need to import matplotlib. pygimli's show does
    >>> import numpy as np
    >>> import pygimli as pg
    >>> fig, ax = pg.plt.subplots(1, 1, figsize=(10, 5))
    >>> u = np.array([1.0, 12.0, 3.0, -4.0, 5.0, 6.0, -1.0])
    >>> xu = np.array(range(len(u)))
    >>> xi = np.linspace(xu[0], xu[-1], 1000)
    >>> _= ax.plot(xu, u, 'o')
    >>> _= ax.plot(xi, pg.interpolate(xi, xu, u, method='linear'),
    ...         color='blue', label='linear')
    >>> _= ax.plot(xi, pg.interpolate(xi, xu, u, method='spline'),
    ...            color='red', label='spline')
    >>> _= ax.plot(xi, pg.interpolate(xi, xu, u, method='harmonic'),
    ...         color='green', label='harmonic')
    >>> _= ax.legend()
    """
    core = False

    if len(args) > 0:
        if isinstance(args[0], pg.Mesh):
            core = True
    if 'srcMesh' in kwargs:
        core = True

    if core:
        return pg.core._pygimli_.interpolate(*args, **kwargs)

    x = 0
    u = 0
    xi = 0

    if len(args) == 3:
        xi = args[0]
        x = args[1]
        u = args[2]

    method = kwargs.pop('method', 'linear')

    if 'linear' in method:
        return np.interp(xi, x, u)

    if 'harmonic' in method:
        coeff = kwargs.pop('nc', int(np.ceil(np.sqrt(len(x)))))
        from pygimli.frameworks import harmfitNative
        return harmfitNative(u, x=x, nc=coeff, xc=xi, err=None)[0]

    if 'spline' in method:
        if pg.io.opt_import("scipy", requiredFor="use interpolate splines."):
            from scipy import interpolate
            tck = interpolate.splrep(x, u, s=0)
            return interpolate.splev(xi, tck, der=0)
        else:
            return xi*0.


if __name__ == '__main__':
    # no need to import matplotlib. pygimli's show does
    # import pygimli as pg
    import pygimli.meshtools as mt
    elec = np.arange(11.)
    topo = np.array([[0., 0.], [3., 2.], [4., 2.], [6., 1.], [10., 1.]])
    pg.plt.plot(topo[:, 0], topo[:, 1])
    p = mt.tapeMeasureToCoordinates(topo, elec)
    pg.plt.plot(p[:, 0], p[:, 1], 'o')
    pg.plt.gca().set_aspect(1)
    pg.wait()
