#!/usr/bin/env python
import numpy as np
import pygimli as pg

# Dummy data container
data = pg.DataContainer()
data.createSensor([0.0, 0.0])
data.createSensor([1.0, 2.0])
data.resize(1)
data.set("s", pg.RVector(1, 1.0))
data.set("g", pg.RVector(1, 2.0))
data.registerSensorIndex("s")
data.registerSensorIndex("g")

# Without secondary nodes
mesh = pg.createGrid([0,1,2],[0,1,2])

# Slowness
slo = [1,2,1,4]

def test_withoutSecNodes():
    fop = pg.TravelTimeDijkstraModelling(mesh, data)
    t_normal = fop.response(slo)
    np.testing.assert_allclose(t_normal, 1 + np.sqrt(2))

def test_withSecNodes():
    mesh2 = mesh.createSecondaryNodes(n=3)
    fop = pg.TravelTimeDijkstraModelling(mesh2, data)
    t_refined = fop.response(slo)
    np.testing.assert_allclose(t_refined, np.sqrt(5)) # only works if n_secNodes is odd number

def test_Jacobian():
    mesh2 = mesh.createSecondaryNodes(n=5)
    fop = pg.TravelTimeDijkstraModelling(mesh2, data)
    fop.createJacobian(slo)
    J = fop.jacobian()
    np.testing.assert_allclose(J * slo, np.sqrt(5))
