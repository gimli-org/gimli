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
        Data of len mesh.cellCount().
        TODO complex, R3Vector, ndarray

    style : str ['mean']
        Interpolation style.
        * 'mean' : non-weighted averaging
        TODO harmonic averaging
        TODO weighted averaging (mean, harmonic)
        TODO interpolation via cell centered mesh

    Examples
    --------
    >>> import pygimli as pg
    >>> grid = pg.createGrid(x=(1,2,3),y=(1,2,3))
    >>> celldata = np.array([1, 2, 3, 4])
    >>> nodedata = pg.meshtools.cellDataToNodeData(grid, celldata)
    >>> print(nodedata.array())
    [ 1.   1.5  2.   2.   2.5  3.   3.   3.5  4. ]
    """
    if len(data) != mesh.cellCount():
        raise BaseException("Dimension mismatch, expecting cellCount(): " +
                            str(mesh.cellCount()) +
                            "got: " + str(len(data)), str(len(data[0])))

    if style == 'mean':

        if np.ndim(data) == 1:
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


def interpolateAlongCurve(curve, t, **kwargs):
    """Interpolate along curve.

    Return curve coordinates for a piecewise linear curve :math:`C(t) = {x_i,y_i,z_i}`
    at positions :math:`t`.
    Curve and :math:`t` values are expected to be sorted along distance from the
    origin of the curve.

    Parameters
    ----------
    curve : [[x,z]] | [[x,y,z]] | [:gimliapi:`GIMLI::RVector3`] | :gimliapi:`GIMLI::R3Vector`
        Discrete curve for 2D :math:`x,z` curve=[[x,z]], 3D :math:`x,y,z`

    t : 1D iterable
        Query positions along the curve in absolute distance

    kwargs :
        If kwargs are given an additional curve smoothing is applied using
        :py:mod:`pygimli.meshtools.interpolate`. The kwargs will be delegated.

    Returns
    -------

    p : np.array
        Curve positions at query points :math:`t`.
        Dimension of p match the size of curve the coordinates.

    Examples
    --------
    >>> # no need to import matplotlib. pygimli's show does
    >>> import numpy as np
    >>> import pygimli as pg
    >>> import pygimli.meshtools as mt
    >>> fig, axs = pg.plt.subplots(2,2)
    >>> topo = np.array([[-2., 0.], [-1., 0.], [0.5, 0.], [3., 2.], [4., 2.], [6., 1.], [10., 1.], [12., 1.]])
    >>> t = np.arange(15.0)
    >>> p = mt.interpolateAlongCurve(topo, t)
    >>> _= axs[0,0].plot(topo[:,0], topo[:,1], '-x', mew=2)
    >>> _= axs[0,1].plot(p[:,0], p[:,1], 'o', color='red') #doctest: +ELLIPSIS
    >>>
    >>> p = mt.interpolateAlongCurve(topo, t, method='spline')
    >>> _= axs[1,0].plot(p[:,0], p[:,1], '-o', color='black') #doctest: +ELLIPSIS
    >>>
    >>> p = mt.interpolateAlongCurve(topo, t, method='harmonic', nc=3)
    >>> _= axs[1,1].plot(p[:,0], p[:,1], '-o', color='green') #doctest: +ELLIPSIS
    >>>
    >>> pg.plt.show()
    >>> pg.wait()
    """
    xC = np.zeros(len(curve))
    yC = np.zeros(len(curve))
    zC = np.zeros(len(curve))

    tCurve = kwargs.pop('tCurve', None)
    if tCurve is None:
        tCurve = pg.utils.cumDist(curve)
    dim = 3

    ## extrapolate starting overlaps
    if min(t) < min(tCurve):
        d = pg.RVector3(curve[1]) - pg.RVector3(curve[0])
        #d[2] = 0.0
        d.normalise()
        curve = np.insert(curve, [0],
                          [curve[0] - np.array(d*(min(tCurve)-min(t)))[0:curve.shape[1]]],
                          axis=0)
        tCurve = np.insert(tCurve, 0, min(t), axis=0)

    ## extrapolate ending overlaps
    if max(t) > max(tCurve):
        d = pg.RVector3(curve[-2]) - pg.RVector3(curve[-1])
        #d[2] = 0.0
        d.normalise()
        curve = np.append(curve,
                          [curve[-1] - np.array(d*(max(t)-max(tCurve)))[0:curve.shape[1]]],
                          axis=0)
        tCurve = np.append(tCurve, max(t))

    if isinstance(curve, pg.R3Vector) or isinstance(curve, pg.stdVectorRVector3):
        xC = pg.x(curve)
        yC = pg.y(curve)
        zC = pg.z(curve)
    else:
        if curve.shape[1] == 2:
            xC = curve[:, 0]
            zC = curve[:, 1]
            dim = 2
        else:
            xC = curve[:, 0]
            yC = curve[:, 1]
            zC = curve[:, 2]

    if len(kwargs.keys()) > 0:
        #interpolate more curve points to get a smooth line
        dTi = min(pg.utils.dist(pg.utils.diff(curve))) / 10.
        ti = np.arange(min(tCurve), max(tCurve)+dTi, dTi)
        xC = pg.interpolate(ti, tCurve, xC, **kwargs)
        zC = pg.interpolate(ti, tCurve, zC, **kwargs)

        if dim == 3:
            yC = pg.interpolate(ti, tCurve, yC, **kwargs)
        tCurve = ti

    xt = interpolate(t, tCurve, xC)
    zt = interpolate(t, tCurve, zC)

    if dim == 2:
        return np.vstack([xt, zt]).T

    yt = interpolate(t, tCurve, yC)

    return np.vstack([xt, yt, zt]).T


def tapeMeasureToCoordinates(tape, pos):
    """Interpolate 2D tape measured topography to 2D Cartesian coordinates.

    Tape and pos value are expected to be sorted along distance to the origin.

    DEPRECATED will be removed use
    :py:mod:`pygimli.meshtools.interpolateAlongCurve` instead

    TODO optional smooth curve with harmfit
    TODO parametric
    TODO parametric + Topo: 3d

    Parameters
    ----------
    tape : [[x,z]] | [RVector3] | R3Vector
        List of tape measured topography points with measured distance (x)
        from origin and height (z)

    pos : iterable
        Array of query positions along the tape measured profile t[0 ..

    Returns
    -------
    res : ndarray(N, 2)
        Same as pos but with interpolated height values.
        The Distance between pos points and res (along curve) points remains.

    """
    pg.deprecated("tapeMeasureToCoordinates", "interpolateAlongCurve")
    return interpolateAlongCurve(tape, pos)


def interpolate(*args, **kwargs):
    r"""Interpolation convinience function.

    Convenience function to interpolate different kind of data.
    Currently supported interpolation schemes are:

    * Interpolate mesh based data from one mesh to another
     (syntactic sugar for the core based interpolate (see below))

      Parameters:
        args: :gimliapi:`GIMLI::Mesh`, :gimliapi:`GIMLI::Mesh`, iterable
            `outData = interpolate(outMesh, inMesh, vals)`
            Interpolate values based on inMesh to outMesh.
            Values can be of length inMesh.cellCount() interpolated to 
            outMesh.cellCenters() or inMesh.nodeCount() which are interpolated tp 
            outMesh.positions().            
            
      Returns:
        Interpolated values.

    * Mesh based values to arbitrary points, based on finite element
      interpolation (from gimli core).

      Parameters:
        args: :gimliapi:`GIMLI::Mesh`, ...
            Arguments forwarded to :gimliapi:`GIMLI::interpolate`
        kwargs:
            Arguments forwarded to :gimliapi:`GIMLI::interpolate`
            
        `interpolate(srcMesh, destMesh)`
            All data from inMesh are interpolated to outMesh

      Returns:
        Interpolated values

    * Interpolate along curve.
      Forwarded to :py:mod:`pygimli.meshtools.interpolateAlongCurve`

      Parameters:
        args: curve, t

        kwargs:
            Arguments forwarded to
            :py:mod:`pygimli.meshtools.interpolateAlongCurve`

    * 1D point set :math:`u(x)` for ascending :math:`x`.
      Find interpolation function :math:`I = u(x)` and
      returns :math:`u_{\text{i}} = I(x_{\text{i}})`
      (interpolation methods are [**linear** via matplotlib,
      cubic **spline** via scipy, fit with **harmonic** functions' via pygimli])
      Note, for 'linear' and 'spline' the interpolate contains all original
      coordinates while 'harmonic' returns an approximate best fit.
      The amount of harmonic coefficients can be specfied with the 'nc' keyword.

      Parameters:
        args: xi, x, u
            * :math:`x_{\text{i}}` - target sample points
            * :math:`x` - function sample points
            * :math:`u` - function values
        kwargs:
            * method : string
                Specify interpolation method 'linear, 'spline', 'harmonic'
            * nc : int
                Number of harmonic coefficients for the 'harmonic' method.

      Returns:
        ui: array of length xi
            :math:`u_{\text{i}} = I(x_{\text{i}})`, with :math:`I = u(x)`


    To use the core functions :gimliapi:`GIMLI::interpolate` start with a
    mesh instance as first argument or use the appropriate keyword arguments.

    TODO

    * 2D parametric to points (method=['linear, 'spline', 'harmonic'])
    * 2D/3D point cloud to points/grids ('Delauney', 'linear, 'spline', 'harmonic')
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
    pgcore = False
    if 'srcMesh' in kwargs:
        pgcore = True

    elif len(args) > 0:
        if isinstance(args[0], pg.Mesh):
            if len(args) == 2 and isinstance(args[1], pg.Mesh):
                return pg.core._pygimli_.interpolate(args[0], args[1],
                                                     **kwargs)

            if len(args) == 3 and isinstance(args[1], pg.Mesh):
                pgcore = False # (outMesh, inMesh, vals)
            else:
                pgcore = True

    if pgcore:
        if len(args) == 3: # args: outData = (inMesh, inData, outPos)

            if args[1].ndim == 2: # outData = (inMesh, mat, vR3)

                outMat = pg.Matrix()
                pg.core._pygimli_.interpolate(args[0],
                                              inMat=np.array(args[1]),
                                              destPos=args[2],
                                              outMat=outMat,
                                              **kwargs)
                return np.array(outMat)

        if len(args) == 4: # args: (inMesh, inData, outPos, outData)

            if args[1].ndim == 1 and args[2].ndim == 1 and args[3].ndim == 1:
                return pg.core._pygimli_.interpolate(args[0],
                                                     inVec=args[1],
                                                     x=args[2],
                                                     y=args[3],
                                                     **kwargs)

            if isinstance(args[1], pg.RMatrix) and \
               isinstance(args[3], pg.RMatrix):
                return pg.core._pygimli_.interpolate(args[0],
                                                     inMat=args[1],
                                                     destPos=args[2],
                                                     outMat=args[3],
                                                     **kwargs)
            if isinstance(args[1], pg.RVector) and \
               isinstance(args[3], pg.RVector):
                return pg.core._pygimli_.interpolate(args[0],
                                                     inVec=args[1],
                                                     destPos=args[2],
                                                     outVec=args[3],
                                                     **kwargs)

        if len(args) == 5:
            if args[1].ndim == 1 and args[2].ndim == 1 and \
               args[3].ndim == 1 and args[4].ndim == 1:
                return pg.core._pygimli_.interpolate(args[0],
                                                     inVec=args[1],
                                                     x=args[2],
                                                     y=args[3],
                                                     z=args[4],
                                                     **kwargs)

        return pg.core._pygimli_.interpolate(*args, **kwargs)
        # end if pg.core:

    if len(args) == 3:

        if isinstance(args[0], pg.Mesh): # args: (outMesh, inMesh, data)
            outMesh = args[0]
            inMesh = args[1]
            data = args[2]

            if isinstance(data, pg.R3Vector) or isinstance(data, pg.stdVectorRVector3):
                x = pg.interpolate(outMesh, inMesh, pg.x(data))
                y = pg.interpolate(outMesh, inMesh, pg.y(data))
                z = pg.interpolate(outMesh, inMesh, pg.z(data))
                return np.vstack([x, y, z]).T

            if isinstance(data, np.ndarray):
                if data.ndim == 2 and data.shape[1] == 3:
                    x = pg.interpolate(outMesh, inMesh, data[:,0])
                    y = pg.interpolate(outMesh, inMesh, data[:,1])
                    z = pg.interpolate(outMesh, inMesh, data[:,2])
                    return np.vstack([x, y, z]).T

            if len(data) == inMesh.cellCount():
                return pg.interpolate(srcMesh=inMesh, inVec=data,
                                      destPos=outMesh.cellCenters())
            elif len(data) == inMesh.nodeCount():
                return pg.interpolate(srcMesh=inMesh, inVec=data,
                                      destPos=outMesh.positions())
            else:
                print(inMesh)
                print(outMesh)
                raise Exception("Don't know how to interpolate data of size",
                                str(len(data)))

            print("data: ", data)
            raise Exception("Cannot interpret data: ", str(len(data)))

        else: #args: xi, x, u

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
                if pg.optImport("scipy", requiredFor="use interpolate splines."):
                    from scipy import interpolate
                    tck = interpolate.splrep(x, u, s=0)
                    return interpolate.splev(xi, tck, der=0)
                else:
                    return xi*0.

    if len(args) == 2: # args curve, t
        curve = args[0]
        t = args[1]
        return interpolateAlongCurve(curve, args, **kwargs)


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
