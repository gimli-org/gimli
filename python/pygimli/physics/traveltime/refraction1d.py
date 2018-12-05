from math import sqrt
import numpy as np
import pygimli as pg


def refractionNlayer(offsets, thk, vel, muteDirect=0):
    """First arrival of n-layered medium."""
    s = 1./np.array(vel)
    tt = offsets * s[0]
    if muteDirect:
        tt[:] = max(tt)
    nlay = len(vel)
    for i in range(1, nlay):  # i-th refracted
        tn = offsets * s[i]  # slope
        for j in range(i):  # sum over intercepts
            dsq = np.maximum(s[j]**2 - s[i]**2, 0)
            tn += 2 * thk[j] * sqrt(dsq)

        tt = np.minimum(tt, tn)
    return tt


class RefractionNLayer(pg.ModellingBase):
    """Forward operator for 1D Refraction seismic with layered model."""
    def __init__(self, offset=0, nlay=3, vbase=1100, verbose=True):
        """Init forward operator. Model are velocity increases.

        Parameters
        ----------
        offset : iterable
            vector of offsets between shot and geophone
        nlay : int
            number of layers
        vbase : float
            base velocity (all values are above)
        """
        super().__init__(verbose=verbose)
        self.nlay = nlay
        self.offset = offset
        self.vbase = vbase
        mesh = pg.createMesh1DBlock(nlay)
        self.setMesh(mesh)

    def thkVel(self, model):
        """Return thickness and velocity vectors from model."""
        return model(0, self.nlay-1), np.cumprod(model(
                self.nlay-1, self.nlay*2-1)) * self.vbase

    def response(self, model):
        """Return forward response f(m)."""
        assert len(model) == self.nlay*2-1
        return refractionNlayer(self.offset, *self.thkVel(model))


class RefractionNLayerFix1stLayer(pg.ModellingBase):
    """FOP for 1D Refraction seismic with layered model (e.g. water layer)."""
    def __init__(self, offset=0, nlay=3, v0=1465, d0=200, muteDirect=0,
                 verbose=True):
        """Init forward operator for velocity increases with fixed 1st layer.

        Parameters
        ----------
        offset : iterable
            vector of offsets between shot and geophone
        nlay : int
            number of layers
        v0 : float
            first layer velocity (at the same time base velocity)
        d0 : float
            depth of first layer
        """
        super().__init__(verbose=verbose)
        self.nlay = nlay
        self.offset = offset
        self.v0 = v0
        self.mDirect = muteDirect
        mesh = pg.createMesh1DBlock(nlay)
        self.setMesh(mesh)

    def thkVel(self, model):
        """Return thickness and velocity vectors from model."""
        thk = pg.cat(pg.Vector(1, self.d0), model(0, self.nlay-1))
        relvel = model(self.nlay-1, self.nlay*2-1)
        vel = self.v0 * np.cumprod(pg.cat(pg.Vector(1, 1.0), relvel))
        return thk, vel

    def response(self, model):
        """Return forward response f(m)."""
        assert len(model) == self.nlay*2-1
        return refractionNlayer(self.offset, *self.thkVel(model), self.mDirect)


if __name__ == '__main__':
    pass
