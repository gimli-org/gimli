#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt

import pygimli as pg
# from solverFVM import (solveFiniteVolume,
#                        createFVPostProzessMesh, diffusionConvectionKernel()


def buildUpB(b, rho, dt, u, v, dx, dy):

    b[1:-1, 1:-1] = rho*(1/dt*((u[2:, 1:-1] - u[0:-2, 1:-1])/(2*dx) +
                         (v[1:-1, 2:]-v[1:-1, 0:-2])/(2*dy)) -
                         ((u[2:, 1:-1]-u[0:-2, 1:-1])/(2*dx))**2 -
                         2*((u[1:-1, 2:]-u[1:-1, 0:-2])/(2*dy) *
                            (v[2:, 1:-1]-v[0:-2, 1:-1])/(2*dx)) -
                         ((v[1:-1, 2:]-v[1:-1, 0:-2])/(2*dy))**2)

    return b


def presPoisson(p, dx, dy, b):
    pn = np.empty_like(p)
    pn = p.copy()

    for q in range(nit):
        pn = p.copy()
        p[1:-1, 1:-1] = ((pn[2:, 1:-1]+pn[0:-2, 1:-1])*dy**2 +
                         (pn[1:-1, 2:]+pn[1:-1, 0:-2])*dx**2) / \
                        (2*(dx**2+dy**2)) - \
            dx**2*dy**2/(2*(dx**2+dy**2))*b[1:-1, 1:-1]

        p[-1, :] = p[-2, :]   # dp/dy = 0 at y = 2
        p[0, :] = p[1, :]   # dp/dy = 0 at y = 0
        p[:, 0] = p[:, 1]   # dp/dx = 0 at x = 0
        p[:, -1] = 0        # p = 0 at x = 2

    return p


def cavityFlow(nt, u, v, dt, dx, dy, p, rho, nu):
    un = np.empty_like(u)
    vn = np.empty_like(v)
    b = np.zeros((ny, nx))

    for n in range(nt):
        un = u.copy()
        vn = v.copy()

        b = buildUpB(b, rho, dt, u, v, dx, dy)
        p = presPoisson(p, dx, dy, b)

        u[1:-1, 1:-1] = un[1:-1, 1:-1] - \
            un[1:-1, 1:-1]*dt/dx*(un[1:-1, 1:-1] - un[0:-2, 1:-1]) - \
            vn[1:-1, 1:-1]*dt/dy*(un[1:-1, 1:-1] - un[1:-1, 0:-2]) - \
            dt/(2*rho*dx)*(p[2:, 1:-1] - p[0:-2, 1:-1]) + \
            nu*(dt/dx**2*(un[2:, 1:-1]-2*un[1:-1, 1:-1] + un[0:-2, 1:-1]) +
                dt/dy**2*(un[1:-1, 2:]-2*un[1:-1, 1:-1] + un[1:-1, 0:-2]))

        v[1:-1, 1:-1] = vn[1:-1, 1:-1] - \
            un[1:-1, 1:-1]*dt/dx*(vn[1:-1, 1:-1] - vn[0:-2, 1:-1]) - \
            vn[1:-1, 1:-1]*dt/dy*(vn[1:-1, 1:-1] - vn[1:-1, 0:-2]) - \
            dt/(2*rho*dy)*(p[1:-1, 2:]-p[1:-1, 0:-2]) + \
            nu*(dt/dx**2*(vn[2:, 1:-1] - 2*vn[1:-1, 1:-1] + vn[0:-2, 1:-1]) +
                dt/dy**2*(vn[1:-1, 2:]-2*vn[1:-1, 1:-1] + vn[1:-1, 0:-2]))

        u[0, :] = 0
        u[:, 0] = 0
        u[:, -1] = 1
        u[-1, :] = 0
        v[0, :] = 0
        v[-1, :] = 0
        v[:, 0] = 0
        v[:, -1] = 0

    return u, v, p

nx = 41
ny = 41
nt = 500
nit = 50
c = 1
dx = 2.0 / (nx-1)
dy = 2.0 / (ny-1)
x = np.linspace(0, 2, nx)
y = np.linspace(0, 2, ny)
Y, X = np.meshgrid(y, x)

rho = 1
nu = .1
dt = .001

u = np.zeros((ny, nx))
v = np.zeros((ny, nx))
p = np.zeros((ny, nx))
b = np.zeros((ny, nx))

nt = 200
u, v, p = cavityFlow(nt, u, v, dt, dx, dy, p, rho, nu)
# fig = plt.figure(figsize=(11,7), dpi=100)
# ax1 = fig.add_subplot(1,3,1)
# ax2 = fig.add_subplot(1,3,2)
# ax3 = fig.add_subplot(1,3,3)

grid = pg.createGrid(x, y)
fig = plt.figure()
ax1 = fig.add_subplot(1, 3, 1)
ax2 = fig.add_subplot(1, 3, 2)
ax3 = fig.add_subplot(1, 3, 3)

pl = pg.logTransDropTol(np.array((p.T).flat), 1e-2)
ul = pg.logTransDropTol(np.array((u.T).flat), 1e-2)
vl = pg.logTransDropTol(np.array((v.T).flat), 1e-2)

pg.show(grid, pl, logScale=False, showLater=True, colorBar=True, axes=ax1,
        cmap='b2r')
pg.show(grid, ul, logScale=False, showLater=True, colorBar=True, axes=ax2)
pg.show(grid, vl, logScale=False, showLater=True, colorBar=True, axes=ax3)

vel = np.vstack([np.array((u.T).flat), np.array((v.T).flat)]).T

pg.viewer.mpl.drawStreams(ax1, grid, vel)

#im1 = ax1.contourf(X,Y,p,alpha=0.5)    ###plnttong the pressure field as a contour
#divider1 = make_axes_locatable(ax1)
#cax1 = divider1.append_axes("right", size="20%", pad=0.05)
#cbar1 = plt.colorbar(im1, cax=cax1)

#im2 = ax2.contourf(X,Y,u,alpha=0.5)    ###plnttong the pressure field as a contour
#divider2 = make_axes_locatable(ax2)
#cax2 = divider2.append_axes("right", size="20%", pad=0.05)
#cbar2 = plt.colorbar(im2, cax=cax2)

#im3 = ax3.contourf(X,Y,v,alpha=0.5)    ###plnttong the pressure field as a contour
#divider3 = make_axes_locatable(ax3)
#cax3 = divider3.append_axes("right", size="20%", pad=0.05)
#cbar3 = plt.colorbar(im3, cax=cax3)

#ax1.contour(X,Y,p)               ###plotting the pressure field outlines
#ax1.quiver(X[::2,::2],Y[::2,::2],u[::2,::2],v[::2,::2]) ##plotting velocity
#ax1.xlabel('X')
#ax1.ylabel('Y')

plt.show()
#drawMesh(ax, grid)