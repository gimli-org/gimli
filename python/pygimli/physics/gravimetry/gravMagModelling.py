#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import pygimli as pg
import time
# from geomagnetics import GeoMagT0  # , date

mu0 = pg.physics.constants.mu0
G = pg.physics.constants.GmGal # mGal

deltaACyl = lambda R__, rho__: 2. * np.pi * R__**2. * rho__
# [m^2 kg/m^3]=[kg/m]
deltaMSph = lambda R__, rho__: 4. / 3. * np.pi * R__**3. * rho__  # [kg]
rabs = lambda r__: np.asarray([np.sqrt(x__.dot(x__)) for x__ in r__])
gradR = lambda r__: (r__.T / rabs(r__))
adot = lambda M__, x__: np.asarray([(a__.dot(M__)) for a__ in x__])


def magnetization(lat, lon, suszept, dat=(2010, 1, 1)):
    """
    TODO
    """
    T0, I, D = GeoMagT0(lat, lon, 0, dat)
    # indizierte Magnetisierung
    Mi = 1. / mu0 * suszept * T0
    # remanente Magnetisierung
    Mr = 0.

    print(T0, I, D, "abs: ", np.sqrt(T0.dot(T0)))

    return Mr + Mi


def BZPoly(pnts, poly, M, openPoly=False):
    """
    TODO

    Parameters
    ----------

    pnts : list
        Measurement points [[p1x, p1z], [p2x, p2z],...]
    poly : list
        Polygon [[p1x, p1z], [p2x, p2z],...]
    M : [M_x, M_y, M_z]
        Magnetization = [M_x, M_y, M_z]

    """
    dg, dgz = calcPolyGz(pnts, poly, 1.0, openPoly)
    return poissonEoetvoes(adot(M, -dgz))


def BaZSphere(pnts, R, pos, M):
    """
    Magnetic anomaly for a sphere.

    Calculate the vertical component of the anomalous magnetic field Bz for a
    buried sphere at position pos with radius R for a given magnetization M at
    measurement points pnts.

    Parameters
    ----------
    pnts : [[x,y,z], ]
        measurement points -- array[x,y,z]
    R : float
        radius
    pos : [float, float, float]
        [x,y,z] -- sphere center
    M : [float, float, float]
        [Mx, My, Mz] -- magnetization

    """
    return poissonEoetvoes(adot(M,
                                gradGZSphere(pnts, R, rho=1.0, pos=pos)))


def BaZCylinderHoriz(pnts, R, pos, M):
    """
    Magnetic anomaly for a horizontal cylinder.
    
    Calculate the vertical component of the anomalous magnetic field Bz for a
    buried horizontal cylinder at position pos with radius R for a given magnetization M at
    measurement points pnts.

    Parameters
    ----------
    
    Parameters
    ----------
    pnts : [[x,y,z], ]
        measurement points -- array[x,y,z]
    R : float
        radius
    pos : [float, float, float]
        [x,y,z] -- sphere center
    M : [float, float, float]
        [Mx, My, Mz] -- magnetization

    """
    return poissonEoetvoes(adot(M, 
                                gradGZCylinderHoriz(pnts, R, rho=1.0, pos=pos)))


def poissonEoetvoes(dg):
    """
    TODO
    """
    return mu0 / (4.0 * np.pi * G) * dg


def uSphere(r, R, rho, pos=(0., 0., 0.)):
    """
    Gravitationspotential einer Sphere mit Radius R und Dichte rho an pos.

    .. math:: u = -G * dM * \frec{1}{r}

    Parameters
    ----------

    """
    return -G * deltaMSph(R, rho) * 1. / rabs(r - pos)


def gradUSphere(r, R, rho, pos=(0., 0., 0.)):
    """
    Gravitationsfeldstrke einer Sphere mit Radius R und Dichte rho an pos


    .. math:: g = -G[m^3/(kg s^2)] * dM[kg] * 1/r^2 1/m^2] * \
    grad(r)[1/1] = [m^3/(kg s^2)] * [kg] * 1/m^2 * [1/1] == m/s^2

    Parameters
    ----------

    Returns
    -------
    [gx, gy, gz] :

    """

    # gesucht eigentlich g_z aber nach unten als -z
    return [1., 1., -1.] * (gradR(r - pos) * -
                            G * deltaMSph(R, rho) * 1. / (rabs(r - pos)**2)).T
# def gSphere(...)


def gradGZSphere(r, R, rho, pos=(0., 0., 0.)):
    """
    TODO

    .. math:: g = -\\nabla u

    Parameters
    ----------

    Returns
    -------
        [\d g_z /\dx, \d g_z /\dy, \d g_z /\dz]
    """
    t = pos[2]

    gzxyz = np.asarray([-3.0 * t * r[:, 0],
                        -3.0 * t * r[:, 1],
                        +2.0 * t * t - r[:, 0]**2 - r[:, 1]**2])
    return (G * deltaMSph(R, rho) / rabs(r - pos)**5. * gzxyz).T


def uCylinderHoriz(pnts, R, rho, pos=(0., 0.)):
    """
    TODO

    Parameters
    ----------

    Returns
    -------

    """
    u = np.zeros(len(pnts))
    for i, r in enumerate(rabs(pnts - pos)):
        if r > R:
            u[i] = -2 * np.pi * G * R * R * rho * np.log(r / R)
        else:
            u[i] = -np.pi * G * rho(r * r - R * R)

    return u


def gradUCylinderHoriz(r, R, rho, pos=(0., 0.)):
    """
    Analytical solution

    Calculate horizontal component of gravimetric field in mGal

    .. math::
        g = -G[m^3/(kg s^2)] * dM[kg/m] * 1/r[1/m] * grad(r)[1/1] =
        [m^3/(kg s^2)] * [kg/m] * 1/m * [1/1] == m/s^2

    Parameters
    ----------
    R   :
        Cylinder radius in [meter]
    p   :
        Cylinder center (x, z)
    rho :
        Density in [kg/m^3]

    Returns
    -------

    """
    return [1., -1.0] * (gradR(r - pos) * -G *
                         deltaACyl(R, rho) * 1. / (rabs(r - pos))).T
# def gZylinderHoriz():


def gradGZCylinderHoriz(r, R, rho, pos=(0., 0.)):
    """
    TODO

    .. math:: g = -grad u(r), with r = [x,z], |r| = \sqrt(x^2+z^2)

    Parameters
    ----------
    r   :

    R   :

    rho :
        Density in [kg/m^3]

    Returns
    -------
    grad gz, [gz_x, gz_z]

    """
    t = pos[1]

    gz_xz = np.asarray([-2.0 * r[:, 0] * (t - r[:, 1]),
                        1.0 * (- r[:, 0]**2 + (t - r[:, 1])**2)])

    return (G * deltaACyl(R, rho) / rabs(r - pos)**4. * gz_xz).T
# def gZSphere(...)


def gzHalfPlateHoriz(pnts, t, rho, pos=(0.0, 0.0)):
    """
    Analytical solution

    g = -grad u,

    Parameters
    ----------
    pnts   :

    t   :

    rho :
        Density in [kg/m^3]

    Returns
    -------
    gz:
        z-component of g
        .. math:: \\nabla(\\partial u/\\partial \\vec{r})_z
    """
    gz = np.zeros(len(pnts))

    for i, q in enumerate(pnts):
        zz1 = q[1] - pos[1]
        xx1 = q[0] - pos[0]
#        gz[i] = G * rho * (np.pi - 2. * \
#         np.arctan((q[0] - pos[0]) / (q[2] - pos[2])))
#         -z direction
        gz[i] = -G * rho * t * (np.pi + 2.0 * np.arctan2(xx1, zz1)) * -1.
#        gz[i] = G * rho * np.pi

    return gz
# def gzPlatteHoriz(...)


def gradGZHalfPlateHoriz(pnts, t, rho, pos=(0.0, 0.0)):
    """
    TODO

    .. math:: g = -\\nabla u

    Parameters
    ----------

    pnts : array (:math:`n\\times 2`)
        n 2 dimensional measurement points
    t : float
        Plate thickness in :math:`[\\text{m}]`
    rho : float
        Density in :math:`[\\text{kg}/\\text{m}^3]`

    Returns
    -------

    gz : array
        Gradient of z-component of g
        :math:`\\nabla(\\frac{\\partial u}{\\partial \\vec{r}}_z)`

    """

    gz = np.zeros((len(pnts), 2))

    for i, q in enumerate(pnts):

        zz1 = q[1] - pos[1]
        xx1 = q[0] - pos[0]

        gz[i, 0] = -2.0 * G * rho * t * (zz1 / (xx1 * xx1 + zz1 * zz1))
        gz[i, 1] = +2.0 * G * rho * t * (xx1 / (xx1 * xx1 + zz1 * zz1))

    return gz
# def gzPlatteHoriz(...)


def lineIntegralZ_WonBevis(p1, p2):
    """
    WonBevis1987.

    Returns
    -------
        g = -grad u =(Fx, 0.0, Fz), dFz(Fzx, Fzy, Fzz)
    """
    dg = pg.RVector3(0.0, 0.0, 0.0)
    dgz = pg.RVector3(0.0, 0.0, 0.0)
    pg.lineIntegralZ_WonBevis(p1, p2, dg, dgz)
    return np.asarray((dg[0], dg[1], dg[2])), np.asarray(
        (dgz[0], dgz[1], dgz[2]))

    x1 = p1[0]
    z1 = p1[1]
    x2 = p2[0]
    z2 = p2[1]

    x21 = x2 - x1
    z21 = z2 - z1
    z21s = z21 * z21
    x21s = x21 * x21

    xz12 = x1 * z2 - x2 * z1

    if x1 == 0. and z1 == 0.:
        return np.asarray((0.0, 0.0, 0.0)), np.asarray((0.0, 0.0, 0.0))
    if x2 == 0. and z2 == 0.:
        return np.asarray((0.0, 0.0, 0.0)), np.asarray((0.0, 0.0, 0.0))

    theta1 = np.arctan2(z1, x1)
    theta2 = np.arctan2(z2, x2)

    r1s = x1 * x1 + z1 * z1
    r2s = x2 * x2 + z2 * z2
    r1 = np.sqrt(r1s)
    r2 = np.sqrt(r2s)

    r21s = x21s + z21s
    R2 = r21s

    rln = np.log(r2 / r1)

    p = (xz12 / r21s) * \
        ((x1 * x21 - z1 * z21) / r1s - (x2 * x21 - z2 * z21) / r2s)
    q = (xz12 / r21s) * \
        ((x1 * z21 + z1 * x21) / r1s - (x2 * z21 + z2 * x21) / r2s)

    Fz = 0.0
    Fx = 0.0
    Fzx = 0.0  # dFz/dx
    Fzz = 0.0  # dFz/dz

    if np.sign(z1) != np.sign(z2):
        if (x1 * z2 < x2 * z1) and z2 >= 0.0:
            theta1 = theta1 + 2. * np.pi

        if (x1 * z2 > x2 * z1) and z1 >= 0.0:
            theta2 = theta2 + 2. * np.pi

    if x1 * z2 == x2 * z1:
        return np.asarray((0., 0.0, 0.)), np.asarray((0., 0.0, 0.))

    th12 = (theta1 - theta2)

    if abs(x21) < 1e-4:
        # print("case 3")
        Fz = x1 * rln
        Fx = 0.0
        Fzz = -p
        Fzx = q - z21s / r21s * rln
        # print(Zz, Zx, R2, x1, z1, x2, z2)

    else:  # default
        B = z21 / x21
        A = (x21 * xz12) / R2

        Fz = A * (th12 + B * rln)
        Fx = A * (-th12 * B + rln)
        z21dx21 = z21 / x21
#        z21x21 = z21 * x21

        fz = (th12 + z21dx21 * rln) / r21s

        Fzz = -p + x21s * fz
        Fzx = q - x21 * z21 * fz

        # // check this!!!
        # fx = (th12 * z21dx21 - rln)/r21s

    # print(np.asarray((Fx, 0.0, Fz)), np.asarray((Fzx, 0.0, Fzz)))
    return np.asarray((Fx, 0.0, Fz)), np.asarray((Fzx, 0.0, Fzz))


def calcPolyGz(pnts, poly, density=1, openPoly=False, forceOpen=False):
    """
    Calculate 2D gravimetric response at given points for a polygon with
    relative density change.

    pnts must be numbered clockwise. Else change the sign of the result.
    Return values are in mGal.

    Bei der magnetischen Loesung fehlt vermutlich ein 1/4.Pi im won & Bevis
    (oetvoes beziehung gl (9) ..... !!check this!!
    """

    qpnts = pnts
    N = len(pnts)

    if np.size(pnts[0]) == 1:
        qpnts = list(zip(pnts, np.zeros(N)))

    if not forceOpen:
        if np.linalg.norm(poly[0] - poly[-1], 2) < 1e-8:
            openPoly = True

    gz = np.zeros((N, 3))
    gzz = np.zeros((N, 3))

    for i, p in enumerate(qpnts):
        for j in range(len(poly) - (openPoly)):
            a = poly[j]
            b = poly[(j + 1) % len(poly)]
        #    print "a, b", a, b
            gzi, gzzi = lineIntegralZ_WonBevis(a - p, b - p)

#            print(gzi, gzzi)
            gz[i, :] += -gzi * [1.0, 1.0, -1.0]
            gzz[i, :] += -gzzi

    return density * 2.0 * -G * gz, density * 2.0 * -G * gzz
# def calcPolydgdz()


def angle(p1, p2, p3, Un):
    """
    Finds the angle between planes O-p1-p2 and O-p2-p3, where p1,p2,p3
    are coordinates of three points, taken in ccw order as seen from origin O.
    This is used by gravMag for finding the solid angle subtended by a polygon
    at the origin. Un is the unit outward normal vector to the polygon.
    After :cite:`SinghGup2001`
    """

    # Check if face is seen from inside
    inout = np.sign(Un.dot(p1))

    x2 = p2[0]
    y2 = p2[1]
    z2 = p2[2]

    # seen from inside; interchange p1 and p3
    if inout > 0:
        x3 = p1[0]
        y3 = p1[1]
        z3 = p1[2]
        x1 = p3[0]
        y1 = p3[1]
        z1 = p3[2]
    elif inout < 0:
        x1 = p1[0]
        y1 = p1[1]
        z1 = p1[2]
        x3 = p3[0]
        y3 = p3[1]
        z3 = p3[2]
    else:
        ang = 0.0
        perp = 1.0
        return ang, perp

    # Normals
    n1 = np.asarray([y2 * z1 - y1 * z2, x1 * z2 - x2 * z1, x2 * y1 - x1 * y2])
    n2 = np.asarray(
        [y3 * z2 - y2 * z3, x2 * z3 - x3 * z2, x3 * y2 - x2 * y3]) * -1.0

    n1 = n1 / np.linalg.norm(n1)
    n2 = n2 / np.linalg.norm(n2)

    perp = np.sum([x3, y3, z3] * n1)

    # sign of perp is -ve if points p1 p2 p3 are in cw order
    perp = np.sign(perp)
    r = np.sum((n1 * n2))
    ang = np.arccos(r)

    if perp < 0:
        ang = 2.0 * np.pi - ang

    return ang, perp
# angle(...)


def gravMagBoundarySinghGup(boundary):
    """
    Calculate [Fx, Fy, FZ] and [Fzx, Fzy, Fzz] at Origin for a given boundary.
    After :cite:`SinghGup2001`

    """
    shape = boundary.shape()
    # print(shape)

    r = shape.center()
    u = shape.norm()
    di = r.dot(u)

    P = 0.
    Q = 0.
    R = 0.

    l = u[0]
    m = u[1]
    n = u[2]

    # Berechne Raumwinkel
    W = 0
    for i in range(shape.nodeCount()):
        p1 = shape.node(i).pos()
        p2 = shape.node((i + 1) % shape.nodeCount()).pos()
        p3 = shape.node((i + 2) % shape.nodeCount()).pos()

        a, p = angle(p1, p2, p3, u)
        W += a

    W -= (shape.nodeCount() - 2) * np.pi

    fsign = float(np.sign(u.dot(shape.node(0).pos())))
    Omega = -fsign * W

    for i in range(shape.nodeCount()):

        vr1 = shape.node(i).pos()
        vr2 = shape.node((i + 1) % shape.nodeCount()).pos()

        r1 = vr1.abs()

        L = (vr2 - vr1).abs()

        Lx = vr2[0] - vr1[0]
        Ly = vr2[1] - vr1[1]
        Lz = vr2[2] - vr1[2]

        b = 2. * (vr1[0] * Lx + vr1[1] * Ly + vr1[2] * Lz)

        b2 = b / L / 2.
        if r1 + b2 == 0:
            I = (1.0 / L) * np.log((L - r1) / r1)
        else:
            I = (1.0 / L) * \
                np.log((np.sqrt(L * L + b + r1 * r1) + L + b2) / (r1 + b2))

        # print(I, L, b, r1,  Lx, Ly, Lz)
        P += I * Lx
        Q += I * Ly
        R += I * Lz

#    print "norm:", u
#    print Omega, l, m, n, P, Q, R
#    exitd
#    print "r", r

    Fx = di * (l * Omega + n * Q - m * R)
    Fy = di * (m * Omega + l * R - n * P)
    Fz = di * (n * Omega + m * P - l * Q)

    # M = [25.575706359959149,    0.000000000000000,   30.479939937615899]
    M = [0, 0, -1.0]
    Pd = u.dot(M)

    Fzx = Pd * (l * Omega + n * Q - m * R)
    Fzy = Pd * (m * Omega + l * R - n * P)
    Fzz = Pd * (n * Omega + m * P - l * Q)

#    print Fx, Fy, Fz, Fzx, Fzy, Fzz
#    Fzx, Fzy, Fzz = Fz* u
#    exitd
    return np.asarray([Fx, Fy, Fz]), np.asarray([Fzx, Fzy, Fzz]),


def solveGravimetry(mesh, dDensity=None, pnts=None, complete=False):
    """
    Solve gravimetric response.

    TOWRITE

    Parameters
    ----------
    mesh : :gimliapi:`GIMLI::Mesh`

    dDensity : float | array
        Density difference.

        * float -- solve for positive boundary marker only.
            Assuming one inhomogeneity.
        * [[int, float]] -- solve for multiple positive boundaries TOIMPL
        * array -- solve for one delta density value per cell
        * None -- return per cell kernel matrix G TOIMPL

    complete : bool [False]
        If True return whole solution or matrix for [dgx, dgy, dgz] and ... TODO

    Returns
    -------
    """
    if pnts is None:
        pnts = [[0.0, 0.0]]

    mesh.createNeighbourInfos()

    Gdg = None
    Gdgz = None
    dg = None
    dgz = None

    if complete:
        Gdg = np.zeros((len(pnts), mesh.cellCount(), 3))
        Gdgz = np.zeros((len(pnts), mesh.cellCount(), 3))
        dg = np.zeros((len(pnts), 3))
        dgz = np.zeros((len(pnts), 3))
    else:
        dg = np.zeros(len(pnts))
        Gdg = np.zeros((len(pnts), mesh.cellCount()))

    #times = []

    dgi = None
    dgzi = None

    for i, p in enumerate(pnts):
        mesh.translate(-pg.RVector3(p))

        for b in mesh.boundaries():
            if b.marker() != 0 or \
                hasattr(dDensity, '__len__') or \
                dDensity is None:

                if mesh.dimension() == 2:
                    #tic = time.time()
                    if complete:
                        dgi, dgzi = lineIntegralZ_WonBevis(b.node(0).pos(),
                                                           b.node(1).pos())
                        #times.append(time.time() - tic)
                        dgi *= -2.0
                        dgzi *= -2.0
                    else:
                        dgi = pg.lineIntegralZ_WonBevis(b.node(0).pos(),
                                                        b.node(1).pos())
                        dgi *= -2.0 * G
                else:
                    if complete:
                        dgi, dgzi = gravMagBoundarySinghGup(b)
                    else:
                        raise Exception("TOIMPL")

                if complete:
                    dgi *= [1.0, 1.0, -1.0]
                    dgi *= -G
                    dgzi *= -G

                if hasattr(dDensity, '__len__') or dDensity is None:
                    cl = b.leftCell()
                    cr = b.rightCell()

                    if cl:
                        Gdg[i][cl.id()] += dgi
                        if complete:
                            Gdgz[i][cl.id()] += dgzi
                    if cr:
                        Gdg[i][cr.id()] -= dgi
                        if complete:
                            Gdgz[i][cr.id()] -= dgzi
                else:
                    dg[i] += dgi * dDensity
                    if complete:
                        dgz[i] += dgzi * dDensity

        mesh.translate(pg.RVector3(p))

    #import matplotlib.pyplot as plt
    #print("times:", sum(times), np.mean(times))
    #plt.plot(times)

    if dDensity is None:
        if complete:
            return Gdg.transpose([0, 2, 1]), Gdgz.transpose([0, 2, 1])
        return Gdg
    elif hasattr(dDensity, '__len__'):
        if complete:
            dg = Gdg.transpose([0, 2, 1]).dot(dDensity)
            dgz = Gdgz.transpose([0, 2, 1]).dot(dDensity)
            return dg, dgz
        else:
            return Gdg.dot(dDensity)

    if complete:
        return dg, dgz
    return dg

class GravimetryModelling(pg.ModellingBase):
    """Gravimetry modelling operator"""
    def __init__(self, verbose=True):
        super(GravimetryModelling, self).__init__(verbose)
        self._J = pg.RMatrix()
        #until reference counting we need to hold the reference here
        self.setJacobian(self._J)

    def createStartmodel(self):
        return pg.RVector(self.regionManger().parameterCount(), 0.0)

    def setSensorPositions(self, pnts):
        self.sensorPositions = pnts

    def response(self, dDensity):
        return solveGravimetry(self.regionManager().paraDomain(),
                               dDensity, pnts=self.sensorPositions,
                               complete=False)

    def createJacobian(self, model):
        Gdz = solveGravimetry(self.regionManager().paraDomain(),
                               dDensity=None,
                               pnts=self.sensorPositions,
                               complete=False)
        self._J.resize(len(Gdz), len(Gdz[0]))
        for i in range(len(Gdz)):
            self._J.setVal(Gdz[i], i)

if __name__ == "__main__":
    import sys
    print(sys.argv[1:])
    print("do some tests here")
#    print lineIntegralZ([-2,-2], [2,-2])
