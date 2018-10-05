# coding=utf-8

from __future__ import print_function

import os
import sys

import pygimli as pg


def streamline(mesh, field, startCoord, dLengthSteps, dataMesh=None,
               maxSteps=1000, verbose=False, coords=(0, 1)):
    """
        Create a streamline from startCoord and following a vector field in up
        and down direction.
    """
    xd, yd, vd = streamlineDir(mesh, field, startCoord, dLengthSteps,
                           dataMesh=dataMesh, maxSteps=maxSteps, down=True,
                           verbose=verbose, coords=coords)

    c = mesh.findCell(startCoord)

    if c is not None:
        c.setValid(True)

    xu, yu, vu = streamlineDir(mesh, field, startCoord, dLengthSteps,
                           dataMesh=dataMesh, maxSteps=maxSteps, down=False,
                           verbose=verbose, coords=coords)

    return xd + xu[1:], yd + yu[1:], vd + vu[1:]


def streamlineDir(mesh, field, startCoord, dLengthSteps, dataMesh=None,
                  maxSteps=1000, down=True, verbose=False, coords=(0, 1)):
    """
        down = -1, up = 1, both = 0
    """
    xd = []
    yd = []
    vd = []
    counter = 0

    pot = None
    vx = None
    vy = None
    isVectorData = False

    if isinstance(field, pg.R3Vector):
        field = field.array()

    if hasattr(field[0], '__len__'):
        if min(field[:, 0]) == max(field[:, 0]) and \
           min(field[:, 1]) == max(field[:, 1]):
            raise BaseException("No data range streamline: min/max == ",
                                min(field[:, 0]))
        vx = pg.RVector(field[:, 0])
        vy = pg.RVector(field[:, 1])

        isVectorData = True
    else:
        if min(field) == max(field):
            raise BaseException("No data range for streamline: min/max == ",
                                min(field))

        if dataMesh is not None:
            if len(field) == dataMesh.nodeCount():
                pot = pg.RVector(field)
            elif len(field) == dataMesh.cellCount():
                pot = pg.cellDataToPointData(dataMesh, field)
            else:
                raise BaseException(
                    "Data length (%i) for streamline is "
                    "neighter nodeCount (%i) nor cellCount (%i)" %
                    (len(field), dataMesh.nodeCount(), dataMesh.nodeCount()))
        else:
            if len(field) == mesh.nodeCount():
                pot = pg.RVector(field)
            elif len(field) == mesh.cellCount():
                pot = pg.cellDataToPointData(mesh, field)
            else:
                raise BaseException(
                    "Data length (%i) for streamline is "
                    "neighter nodeCount (%i) nor cellCount (%i)" %
                    (len(field), mesh.nodeCount(), mesh.nodeCount()))

    direction = 1
    if down:
        direction = -1

    # search downward
    pos = pg.RVector3(startCoord)
    c = mesh.findCell(startCoord)
    dLength = c.center().dist(c.node(0).pos()) / dLengthSteps

    # stream line starting point
    if c is not None:
        xd.append(pos[coords[0]])
        yd.append(pos[coords[1]])
        vd.append(-1)

    lastC = c
    lastU = -direction * 1e99

    d = None
    while c is not None and len(xd) < maxSteps:

        # valid .. temporary check if there is already a stream within the cell
        if not c.valid():
            break

        if isVectorData:
            u = 0.
            if len(vx) == mesh.cellCount():
                d = pg.RVector3(vx[c.id()], vy[c.id()])
            elif len(vx) == mesh.nodeCount():
                d = pg.RVector3(c.pot(pos, vx), c.pot(pos, vy))
            elif dataMesh:
                cd = dataMesh.findCell(pos)
                if cd is None:
                    raise BaseException("Cannot find " + str(pos) +
                                        " dataMesh")
                if len(vx) == dataMesh.cellCount():
                    d = pg.RVector3(vx[cd.id()], vy[cd.id()])
                elif len(vx) == dataMesh.nodeCount() and \
                    len(vy) == dataMesh.nodeCount():
                    d = pg.RVector3(cd.pot(pos, vx), cd.pot(pos, vy))
                else:
                    print(dataMesh)
                    print(len(vx), len(vy))
                    raise BaseException("data size wrong")
            else:
                raise Exception
        else:
            if dataMesh:
                cd = dataMesh.findCell(pos)
                if not cd:
                    break

                d = cd.grad(pos, pot)
                u = cd.pot(pos, pot)
            else:
                d = c.grad(pos, pot)
                u = c.pot(pos, pot)
        # print "cell:", c.id(), u
        # always go u down
        dAbs = d.length()
        if dAbs == 0.0:
            print(d,
                  "check this in streamlineDir(",)
            break

        if down:
            if u > lastU:
                break
        else:
            if u < lastU:
                break

        # * min(1.0, ((startCoord - pos).length()))
        pos += direction * d / dAbs * dLength
        c = mesh.findCell(pos, False)

        # Change cell here .. set old cell to be processed
        if c is not None:

            xd.append(pos[coords[0]])
            yd.append(pos[coords[1]])
            # set the stating value here
            if vd[0] == -1:
                vd[0] = dAbs
            vd.append(dAbs)

            # If the new cell is different from the current we move into the
            # new cell and make the last to be invalid ..
            # the last active contains a stream element
            if c.id() != lastC.id():
                lastC.setValid(False)
                lastC = c
                dLength = c.center().dist(c.node(0).pos()) / dLengthSteps
        else:
            # There is no new cell .. the last active contains a stream element
            lastC.setValid(False)

        lastU = u
        if verbose:
            print(pos, u)

    # Stream line has stopped and the current cell (if there is one) ..
    # .. contains a stream element
    if c is not None:

        c.setValid(False)

    if down:
        xd.reverse(), yd.reverse(), vd.reverse()

    return xd, yd, vd


def boundaryPlaneIntersectionLines(boundaries, plane):
    """Create Lines from boundaries that intersect a plane."""
    lines = []

    for b in boundaries:
        ps = []
        for i, n in enumerate(b.shape().nodes()):
            line = pg.Line(n.pos(), b.shape().node(
                (i + 1) % b.shape().nodeCount()).pos())
            p = plane.intersect(line, 1e-8, True)
            if p.valid():
                ps.append(p)

        if len(ps) == 2:
            lines.append(list(zip([ps[0].x(), ps[1].x()],
                                  [ps[0].z(), ps[1].z()])))
    return lines
# def intersectionLines


def number_of_processors():
    """Return number of processors on multiple platoforms."""
    # Windows
    if os.name == 'nt':
        return int(os.getenv('NUMBER_OF_PROCESSORS'))
    # Linux
    elif sys.platform == 'linux2':
        retv = 0
        with open('/proc/cpuinfo', 'rt') as cpuinfo:
            for line in cpuinfo:
                if line[:9] == 'processor':
                    retv += 1
        return retv

    # Please add similar hacks for MacOSX, Solaris, Irix,
    # FreeBSD, HPUX, etc.
    else:
        raise RuntimeError('unknown platform')
