# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import pygimli as pg


class FunctionModelling(pg.ModellingBase):
    """New modelling operator returning f(x) derived from modelling base."""

    def __init__(self, nc, xvec, verbose=False):
        """Constructor, nc: number of coefficients, xvec: abscissa."""
        pg.ModellingBase.__init__(self, verbose)
        self.x_ = xvec
        self.nc_ = nc
        self.regionManager().setParameterCount(nc)

    def response(self, par):
        """
           the main thing - the forward operator: return f(x)
        """
        y = pg.RVector(len(self.x_), par[0])
        for i in range(1, self.nc_):
            y += pg.pow(self.x_, i) * par[i]

        return y

    def startModel(self):
        """Define the starting model."""
        return pg.RVector(self.nc_, 0.5)

if __name__ == "__main__":
    nc = 1  # polynomial degree
    x = np.arange(10.0)
    y = 1.1 + 0.6 * x
    y += np.random.randn(len(y)) * 0.2

    # two coefficients and x-vector (first data column)
    f = FunctionModelling(nc + 1, x)

    # initialize inversion with data and forward operator and set options
    inv = pg.RInversion(y, f, True, True)

    # constant absolute error of 0.01 (not necessary, only for chi^2)
    inv.setAbsoluteError(0.01)

    # the problem is well-posed and does not need regularization
    inv.setLambda(0)

    # actual inversion run yielding coefficient model
    coeff = inv.run()

    # get actual response and write to file.
    pg.save(inv.response(), "resnp.out")

    # print result to screen and save coefficient vector to file
    s = "y = " + str(round(coeff[0] * 1000) / 1000)

    for i in range(1, nc+1):
        s = s + " + " + str(round(coeff[i] * 1000) / 1000) + " x^" + str(i)

    print(s)

    pg.save(coeff, "out.vec")

    plt.plot(x, y, 'rx', x, inv.response(), 'b-')
    plt.title(s)
    plt.xlabel("x")
    plt.ylabel("y")
    plt.legend(("measured", "fitted"), loc="upper left")
    plt.show()
