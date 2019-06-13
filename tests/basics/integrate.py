#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
pygimli integration function
"""

from math import pi, sqrt  # , sin, cos

import pygimli as pg


class funct:
    def __init__(self, f):
        self.f = f

    def __call__(self, vec):
        ret = pg.Vector(len(vec))
        for i, arg in enumerate(vec):
            ret[i] = self.f(arg)
        return ret


def integrate(f, ent, order):
    """ integrate function """

    J = 0
    x = []
    w = []
    if type(ent) is list:

        a = ent[0]
        b = ent[1]
        xs = pg.IntegrationRules.instance().gauAbscissa(order)
        w = pg.IntegrationRules.instance().gauWeights(order)

        x = (b - a) / 2.0 * pg.x(xs) + (a + b) / 2.0
        J = (b - a) / 2.

    else:
        J = ent.shape().jacobianDeterminant()
        xs = pg.IntegrationRules.instance().abscissa(ent.shape(), order)
        w = pg.IntegrationRules.instance().weights(ent.shape(), order)

        x = [ent.shape().xyz(xsi) for xsi in xs]

    for xx in x:
        print xx
    print w
    print funct(f)(x)
    return J * sum(funct(f)(x) * w)
# def integrate( ... )


def check(name, ist, soll):
    """
    """
    result = 'fail'
    if abs(ist-soll) < 1e-8:
        result = 'ok'
    print(name + ": ", ist, result, ist - soll)
# def check(  )

check('1d.1', integrate(lambda x_: x_, [0, 1], 9), 0.5)

mesh = pg.Mesh(2)

# H. T. RATHOD1*, K. V. NAGARAJA2, B. VENKATESUDU3 AND N. L. RAMESH4.
# Gauss Legendre quadrature over a triangle.
# J. Indian Inst. Sci., Sept.-Oct. 2004, 84, 183-188

n0 = mesh.createNode(pg.RVector3(0.0, 0.0, 0.0))
n1 = mesh.createNode(pg.RVector3(1.0, 0.0, 0.0))
n2 = mesh.createNode(pg.RVector3(0.0, 1.0, 0.0))
n3 = mesh.createNode(pg.RVector3(1.0, 1.0, 0.0))
n4 = mesh.createNode(pg.RVector3(pi/2., 0.0, 0.0))
n5 = mesh.createNode(pg.RVector3(pi/2., pi/2., 0.0))

tri1 = mesh.createTriangle(n0, n1, n2)
tri2 = mesh.createTriangle(n0, n1, n3)
tri3 = mesh.createTriangle(n0, n4, n5)

order = 9
# table(2) iter 9 : 0.400001147 0.666332910 0.874071505 1.000000006 0.718126535
# check( 'tri1.f1', integrate( lambda arg__: sqrt( arg__[0] + arg__[1] ), tri1, order ), 0.400001147 ) # 0.4
# check( 'tri1.f2', integrate( lambda arg__: 1./sqrt( arg__[0] + arg__[1] ), tri1, order ), 0.666332910 ) # 0.666666667
# check( 'tri2.f3', integrate( lambda arg__: 1./sqrt( arg__[0]**2. + arg__[1]**2. ), tri2, order ), 0.874071505) # 0.881373587
# check( 'tri3.f4', integrate( lambda arg__: sin( arg__[0] ) + cos( arg__[1] ), tri3, order ), 1.000000006 ) # 1.0
# check( 'tri2.f5', integrate( lambda arg__: exp( abs( arg__[0] + arg__[1] -1. ) ), tri2, order ), 0.718126535) # 0.71828183

quad1 = mesh.createQuadrangle(n0, n1, n3, n2)
check('quad1.f1', integrate(lambda arg__: arg__[0] * arg__[0] + arg__[1],
                            quad1, order), 5./6.)
# check('quad1.f2', integrate( lambda arg__: 1.0, quad1, order ), 1.0 )

edge1 = mesh.createEdge(n0, n1)
edge2 = mesh.createEdge(n1, n3)
edge3 = mesh.createEdge(n3, n2)
edge4 = mesh.createEdge(n2, n0)

print edge1.norm()

i1 = integrate(lambda arg__: sqrt(arg__[0] * arg__[0] + arg__[1]), edge1, 3)
i2 = integrate(lambda arg__: sqrt(arg__[0] * arg__[0] + arg__[1]), edge2, 3)
i3 = integrate(lambda arg__: sqrt(arg__[0] * arg__[0] + arg__[1]), edge3, 3)
i4 = integrate(lambda arg__: sqrt(arg__[0] * arg__[0] + arg__[1]), edge4, 3)

print(i1, edge1.norm())
print(i2, edge2.norm())
print(i3, edge3.norm())
print(i4, edge4.norm())
print(i1+i2+i3+i4)
