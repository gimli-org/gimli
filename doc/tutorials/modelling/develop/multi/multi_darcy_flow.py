#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

Darcy flow 

.. math::

    \vec(v) - \grad p & = 0 \\
    \grad\cdot\vec(v) & = g \\
    u & = u_0\,\text{on}\,\Gamma_D \\
    \vec(v) \vec{n} & = g\,\text{on}\,\Gamma_N


On the unit square.
"""

import pygimli as pg
import numpy as np
from pygimli.solver import solve

import pygimli.solver.fenics
from pygimli.solver.fenics import dolfin as dlf


def darcyFlow(model, p0, verbose=0):
    if verbose:
        dlf.set_log_level(dlf.INFO)
    else:
        dlf.set_log_level(dlf.ERROR)


    # Create mesh and define function space
    #mesh = dlf.UnitSquareMesh(16, 16, 'crossed')
    mesh = None
    try:
        mesh = dlf.RectangleMesh(dlf.Point(0, -5), dlf.Point(10, 0), 40, 20, "crossed")
    except:
        # for older fenics versions 
        mesh = dlf.RectangleMesh(0, -5, 10, 0, 40, 20, "crossed")

    class PermeabiltyInv(dlf.Expression):
        
        def eval_cell(self, values, r, cell):
            if hasattr(model, '__iter__'):
                if r[1] < -3.5:
                    values[:] = 1./model[2]
                elif r[1] < -0.5 - dlf.DOLFIN_EPS:
                    values[:] = 1./model[1]
                else:
                    values[:] = 1./model[0]
            else:
                #print(cell)
                #print(dir(cell))
                cMid = dlf.Cell(mesh, cell.index).midpoint()
                
                #print(r, cMid[0], cMid[1])
                #print(cell.get_vertex_coordinates())
                
                #print(dlf.Cell(cell))
                #print(dir(cell.this.midpoint()))
                #print(cell, cell.midpoint())
                c = model.findCell((cMid[0], cMid[1]))
                values[:] = 1./c.attribute()
            #print(r, cell, values)
            #print(values[:])

    kinv11 = PermeabiltyInv()
    kinv12 = dlf.Constant(0.0)
    kinv21 = dlf.Constant(0.0)
    kinv22 = PermeabiltyInv()
    Kinv = dlf.as_matrix(((kinv11, kinv12), (kinv21, kinv22)))


    order = 1
    # Define function spaces
    BDM = dlf.FunctionSpace(mesh, "RT", order)
    DG = dlf.FunctionSpace(mesh, "DG", order-1)
    W = BDM * DG

    def topBoundary(r, on_boundary):
        return on_boundary and r[1] == 0.0# - dlf.DOLFIN_EPS
        
    def inflowBoundary(r, on_boundary):
        return on_boundary #and r[0] < 0.0 + dlf.DOLFIN_EPS

    class pBoundaryDir(dlf.Expression):
        def eval(self, values, r):
            values[:] = 0.0
            if r[0] < 1e-7 and r[1] < -0.49:
                values[:] = p0
            
            
    class vBoundaryDir(dlf.Expression):
        def eval(self, values, r):
            values[:]=0.0
            #values[:] = -(alpha/2 * x*y**2 + beta*x - alpha/6.0 * x**3.0)
        def shape(self):
            return (2,)
            
    vBoundaryDir()
    top_val = dlf.Constant((0.0, 0.0))
    # pressureBound
    bc = dlf.DirichletBC(W.sub(0), top_val, topBoundary)

    f = dlf.Constant("0")

    g = pBoundaryDir()
    #g = dlf.Expression("-(0.3/2 * x[0]*x[1]*x[1] + 1*x[0] - 0.3/6.0 * x[0]*x[0]*x[0])")

    # Define trial and test functions
    (v, p) = dlf.TrialFunctions(W)# the unknowns to be recover
    (w, q) = dlf.TestFunctions(W)
    n = dlf.FacetNormal(mesh)

    #A = dot(Kinv*v, u)*dx - div(v)*p*dx + q*div(u)*dx
    #q*f*dx - inner(v, pbar*n)*ds
    A = (dlf.dot(Kinv * v, w) - dlf.div(w)*p + q * dlf.div(v)) * dlf.dx
    L = f * q * dlf.dx - dlf.inner(w, g*n) * dlf.ds

    w = dlf.Function(W)
    dlf.solve(A == L, w, bc)
    (v, p) = w.split()

    #dlf.plot(mesh, Kinv)
    #pygimli part starts here
    gimlimesh = pg.solver.fenics.convertDolfinMesh(mesh)
    pre = p.compute_vertex_values()
    vel = v.compute_vertex_values().reshape(2, len(pre))
            
    Q = dlf.FunctionSpace(mesh,"DG",0) 
    kinv = dlf.Function(Q)                                                                 
    kinv.interpolate(PermeabiltyInv())

    #pg.show(model, model.cellMarker(), label='marker')
    #pg.wait()
    #gimlimesh.setCellMarker(pg.interpolate(model, model.cellMarker(), gimlimesh.cellCenters()))

    return gimlimesh, np.array(vel), np.array(pre), np.array(1./np.array(kinv.vector()))

if __name__ == '__main__':
    mesh, vel, p, k = darcyFlow(model=[1, 0.01, 0.001])

    fig = pg.plt.figure()
    ax1 = fig.add_subplot(1, 2, 1)
    pg.show(mesh, axes=ax1)
    ax, cbar = pg.show(mesh,  p, colorBar=1, axes=ax1)
    cbar.ax.set_xlabel('pressure in ??')
    pg.show(mesh, vel, axes=ax1)
    ax2 = fig.add_subplot(1, 2, 2)
    pg.show(mesh, axes=ax2)
    ax, cbar = pg.show(mesh, (np.sqrt(vel[0]**2+vel[1]**2)), logScale=1, colorBar=1, axes=ax2)
    cbar.ax.set_xlabel('velocity in m/s)')
    pg.show(mesh, vel, axes=ax2)

    ## show perm
    #fig = pg.plt.figure()
    #ax2 = fig.add_subplot(1, 1, 1)
    #pg.show(mesh, axes=ax2)
    #ax, cbar = pg.show(mesh, 1./np.array(K.vector()), logScale=1, colorBar=1, axes=ax2)
    #cbar.ax.set_xlabel('Permeabilty in ??')
    #pg.show(mesh, vel, axes=ax2)


pg.wait()

