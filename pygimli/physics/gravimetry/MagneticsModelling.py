import numpy as np
import pygimli as pg
from pygimli.matrix import NumpyMatrix


def SolveGravMagHolstein(pgmesh, pnts, cmp, igrf, foot=np.inf):
    if pnts is None:
        pnts = [[0.0, 0.0]]

    kernel = np.zeros((len(pnts), pgmesh.cellCount(), len(cmp)))
    # org: this does not make sense as igrf is either 3 or 7 long
    # B_dir = np.array(igrf / np.linalg.norm(igrf))
    # fakt = igrf[6] / (4*np.pi)
    # rather do like this:
    if len(igrf) == 3:  # an X, Y, Z vector
        F = np.linalg.norm(igrf)
        fakt = F / (4*np.pi)
        B_dir = np.array(igrf) / F
    elif len(igrf) == 7:  # an IGRF vector (D, I, H, X, Y, Z, F)
        fakt = igrf[6] / (4*np.pi)
        igrf = np.array(igrf[3:6])
        B_dir = np.array(igrf) / np.linalg.norm(igrf)
    else:
        raise Exception("Could not use IGRF vector")

    b_list, n_list, c_list = [], [], []
    for bd in pgmesh.boundaries():
        b_list.append([bd.allNodes()[0].id(),
                       bd.allNodes()[1].id(),
                       bd.allNodes()[2].id()])
        c_list.append([bd.leftCell(), bd.rightCell()])

    b_list = np.array(b_list)
    lb = len(b_list)

    for nd in pgmesh.nodes():
        n_list.append(nd.pos())

    n_list = np.array(n_list)

    cl, cr = [], []
    for i, c in enumerate(c_list):
        if c[0]:
            cl.append([i, c[0].id()])
        if c[1]:
            cr.append([i, c[1].id()])

    cl = np.array(cl)
    cr = np.array(cr)

    for i, p in enumerate(pnts):
        temp = np.zeros((lb, len(cmp)))
        r1 = n_list[b_list] - p
        r2 = r1[:, [1, 2, 0], :]
        r0 = r2 - r1
        u = np.sum(np.cross(r1, r2), 1)
        u = u / np.expand_dims(np.linalg.norm(u, axis=1), axis=1)
        ut = np.tile(u, 3).reshape((lb, 3, 3))
        ll = np.linalg.norm(r0, axis=2)
        t = r0/np.expand_dims(ll, axis=2)
        lm = (np.sum(r1*t, 2) + np.sum(r2*t, 2)) / 2
        h = np.cross(t, ut)
        hn = np.sum(h*r1, 2)
        v = np.sum(ut*r1, 2)
        r1n = np.linalg.norm(r1, axis=2)
        r2n = np.linalg.norm(r2, axis=2)
        rm = (r1n+r2n)/2
        lumbda = ll/(2*rm)

        # gravitational field
        g = hn*np.arctanh(lumbda)-np.sign(v)*v*np.arctan2(
            hn*lumbda, (rm*(1-lumbda**2)+abs(v)))

        # magnetic field vector and gravity gradient tensor
        b = h*np.expand_dims(np.arctanh(lumbda), axis=2) - \
            ut*np.expand_dims(np.sign(v)*np.arctan2(
                hn*lumbda, (rm*(1-lumbda**2)+abs(v))), axis=2)

        # magnetic gradient tensor
        d = (-2*lumbda*hn) / (r1n*r2n*(1-lumbda**2))
        e = (-2*lumbda*lm) / (r1n*r2n)
        f = (-2*lumbda*v) / (r1n*r2n*(1-lumbda**2))

        h1 = np.expand_dims(h, axis=3)
        h2 = np.swapaxes(h1, 2, 3)
        t1 = np.expand_dims(t, axis=3)
        t2 = np.swapaxes(t1, 2, 3)
        u1 = np.expand_dims(ut, axis=3)
        u2 = np.swapaxes(u1, 2, 3)

        B = (h1*h2-u1*u2)*np.expand_dims(d, (2, 3)) + \
            (t1*h2+h1*t2)*np.expand_dims(e, (2, 3))/2 + \
            (h1*u2+u1*h2)*np.expand_dims(f, (2, 3))

        P = np.dot(u, B_dir)

        g_vec = u*np.expand_dims(np.sum(g, 1), axis=1)
        B_vec = np.expand_dims(P, 1) * np.sum(b, 1)
        B_tens = np.expand_dims(P, (1, 2)) * np.sum(B, 1)

        jj = 0
        if 'gx' in cmp:
            temp[:, jj] = g_vec[:, 0]
            jj += 1

        if 'gy' in cmp:
            temp[:, jj] = g_vec[:, 1]
            jj += 1

        if 'gz' in cmp:
            temp[:, jj] = g_vec[:, 2]
            jj += 1

        # if 'gxx' in cmp:
        #     temp[:, jj]=G_tens[:, 0, 0]
        # if 'gxy' in cmp:
        #     temp[:, jj]=G_tens[:, 0, 1]
        # if 'gxz' in cmp:
        #     temp[:, jj]=G_tens[:, 0, 2]
        # if 'gyy' in cmp:
        #     temp[:, jj]=G_tens[:, 1, 1]
        # if 'gyz' in cmp:
        #     temp[:, jj]=G_tens[:, 1, 2]
        # if 'gzz' in cmp:
        #     temp[:, jj]=G_tens[:, 2, 2]

        if 'TFA' in cmp:
            temp[:, jj] = fakt*B_vec.dot(B_dir)
            jj += 1

        if 'Bx' in cmp:
            temp[:, jj] = fakt*B_vec[:, 0]
            jj += 1

        if 'By' in cmp:
            temp[:, jj] = fakt*B_vec[:, 1]
            jj += 1

        if 'Bz' in cmp:
            temp[:, jj] = fakt*B_vec[:, 2]
            jj += 1

        if 'Bxx' in cmp:
            temp[:, jj] = fakt*B_tens[:, 0, 0]
            jj += 1

        if 'Bxy' in cmp:
            temp[:, jj] = fakt*B_tens[:, 0, 1]
            jj += 1

        if 'Bxz' in cmp:
            temp[:, jj] = fakt*B_tens[:, 0, 2]
            jj += 1

        if 'Byy' in cmp:
            temp[:, jj] = fakt*B_tens[:, 1, 1]
            jj += 1

        if 'Byz' in cmp:
            temp[:, jj] = fakt*B_tens[:, 1, 2]
            jj += 1

        if 'Bzz' in cmp:
            temp[:, jj] = fakt * B_tens[:, 2, 2]
            jj += 1

        kernel[i][cl[:, 1]] += temp[cl[:, 0]]
        kernel[i][cr[:, 1]] -= temp[cr[:, 0]]

    return kernel.transpose([0, 2, 1])


class GravMagModelling(pg.Modelling):

    def __init__(self, mesh, points, cmp, igrf, foot=None):
        super().__init__()
        self.mesh = mesh
        self._J = pg.Matrix()
        self.sensorPositions = points
        self.components = cmp
        self.igrf = igrf
        self.footprint = foot
        self.setJacobian(self._J)
        self.kernel = SolveGravMagHolstein(self.mesh,
                                           pnts=self.sensorPositions,
                                           cmp=self.components, igrf=self.igrf,
                                           foot=self.footprint)
        self._J = NumpyMatrix(self.kernel)

    def response(self, model):
        return self.kernel.dot(model)

    def createJacobian(self):
        pass
