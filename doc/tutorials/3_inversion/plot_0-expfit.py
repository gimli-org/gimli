#!/usr/bin/env python
# encoding: utf-8

r"""
Simple fit
==========

This tutorial shows how to do the simplest inversion case, a curve fit, by
setting up an own forward operator. The function to be fitted is`

.. math::

    f(x) = A * e^{-x/X}

with the two unknown coefficients A (a signal  amplitude) and X (a decay rate).
Both A and X are assumed to be positive which is often the case for physical
properties. The easiest way to do this is via a logarithmic transformation of
the model vector (containing A and X) which is very easily done in pyGIMLi.

First we import the pygimli library under a short name pg and the numerics
library numpy. Additionally we load the python plotting module of the library
matplotlib. Both are contained in most python distributions and systems.
"""

import pygimli as pg
import numpy as np
import matplotlib.pyplot as plt


###############################################################################
# We set up the modelling operator, i.e. to return :math:`{\bf f}({\bf x})` for
# given model parameters A and X subsumed in a vector. In order to be able to
# use operator in inversion, we derive from the abstract modelling base class.
# The latter holds the main mimic of generating Jacobian and adminstrating the
# model, the regularization and so on. The only function to overwrite is
# **response()**. If no function **createJacobian** is provided, they are
# computed by brute force (forward calculations with altered parameters).


class ExpModelling(pg.Modelling):
    def __init__(self, xvec, verbose=False):
        super().__init__()
        self.x = xvec

    def response(self, model):
        return model[0] * pg.exp(-self.x / model[1])

    def createStartModel(self, dataVals):
        return pg.Vector([1.0, 3.0])

###############################################################################
# The init function saves the x vector and defines the parameterization, i.e.
# two independent parameters (a 1D mesh with 1 cell and 2 properties).
# The response function computes the function using A=model[0] and X=model[1]
# The function startModel defines a meaningful starting vector. There are other
# methods to set the starting model as inv.setModel() but this one is a default
# one for people who use the class and forget about a starting model.

# We first create an abscissa vector using numpy (note that pygimli also
# provides an exp function and generate synthetic data with two arbitrary A and
# X values.

x = np.arange(0, 1, 1e-2)
data = 10.5 * np.exp(- x / 550e-3)

###############################################################################
# We define an (absolute) error level and add Gaussian noise to the data.

error = 0.5
data += pg.randn(*data.shape)*error
relError = error / data

###############################################################################
# Next, an instance of the forward operator is created. We could use it for
# calculating the synthetic data using f.response([10.5, 0.55]) or just
# f([10.5, 0.55]). We create a real-valued (R) inversion passing the forward
# operator, the data. A verbose boolean flag could be added to provide some
# output the inversion, another one prints more and saves files for debugging.

f = ExpModelling(x)
inv = pg.Inversion(f)

###############################################################################
# We create a real-valued logarithmic transformation and apply it to the model.
# Similar could be done for the data which are by default treated linearly.
# We then set the error level that is used for data weighting. It can be a
# float number or a vector of data length. One can also set a relative error.
# Finally, we define the inversion style as Marquardt scheme (pure local damping
# with decreasing the regularization parameter subsequently) and start with a
# relatively large regularization strength to avoid overshoot.
# Finally run yields the coefficient vector and we plot some statistics.

tLog = pg.trans.TransLog()
f.modelTrans = tLog
inv._inv.setMarquardtScheme()
inv._inv.setLambda(100)
coeff = inv.run(data, relError, verbose=True)
print(inv.relrms(), inv.chi2())
print(coeff)

###############################################################################
# We see that after 5 iterations the absolute rms value equals the noise level
# corresponding to a chi-squared misfit value of 1 as it should be the case for
# synthetic data. The relative rms (in %) is less relevant here but can be for
# other applications. Additionally the ranges for model and model response are
# given and the objective function consisting of data misfit and model
# roughness times lambda. Note that due to the local regularization the second
# term does not contribute to Phi. Set verbose to True to see the whole course
# of inversion. The values of the coefficient vector (a GIMLi real vector) are
# as expected close (i.e. equivalent) to the synthetic model.

###############################################################################
# We finally create a plotting figure and plot both data and model response.

plt.figure()
plt.plot(x, data, 'rx', x, inv.response, 'b-')

###############################################################################
# The createMesh1D automatically attributed the markers 0 and 1 to the two
# model parameters A and X, respectively. Each marker leads to a region that
# can be individually treated, e.g. the starting value, lower or upper bounds,
# or all three at the same time (setParameters). This changes the model
# transformation which can of course be region-specific.

# f.region(0).setLowerBound(0.1)
# f.region(0).setStartModel(3)
# f.region(1).setParameters(0.3, 0.01, 1.0)

###############################################################################
# If these are set before the inversion is used, they are used automatically.
# We set the model by hand using the new starting model

# inv.setVerbose(True)
# inv.setModel(f.createStartModel())
# print(inv.run())
# inv.echoStatus()

###############################################################################
# The result is pretty much the same as before but for stronger equivalence or
# smoothness-constrained regularization prior information might help a lot.
