#!/usr/bin/env python
# encoding: utf-8

r"""
Polyfit
=======

This tutorial shows a flexible inversion with an own forward calculation that
includes an own jacobian. We start with fitting a polynomial of degree
:math:`P`

.. math::

    f(x) = p_0 + p_1 x + \ldots + p_P x^P = \sum\limits_{i=0}^{P} p_i x^i

to given data :math:`y`.
The unknown model is the coefficient vector :math:`{\bf m}=[p_0,\ldots,p_P]`.
The vectorized function for a vector
:math:`{\bf x}=[x_1,\ldots,x_N]^T`
can be written as matrix-vector product

.. _eq:yAx:
.. math::

  {\bf f} ({\bf x}) = {\bf A} {\bf x} \quad\mbox{with}\quad {\bf A}=
  \left[
    \begin{array}{cccc}
        1 & x_1    & \ldots & x_1^P \\
   \vdots & \vdots & \ddots & \vdots \\
        1 & x_N    & \ldots & x_N^P
  \end{array}
  \right] =
  [ {\bf 1}\quad {\bf x} \quad {\bf x}^2 \ldots {\bf x}^P ] \;.

We set up the modelling operator, i.e. to return :math:`{\bf f}({\bf x})` for
given :math:`p_i`, as a class derived from the modelling base class.
The latter holds the main mimic of generating Jacobian, gradients by brute
force. The only function to overwrite is \cw{response()}.

Python is a very flexible language for programming and scripting and has many
packages for numerical computing and graphical visualization.
For this reason, we built Python bindings and compiled the library pygimli.
As a main advantage, all classes can be used and derived.
This makes the use of GIMLi very easy for non-programmers.
All existing modelling classes can be used, but it is also easy to create new
modelling classes.

We exemplify this by the preceding example.

First, the library must be imported.

To avoid name clashes with other libraries we suggest to import `pygimli` and
alias it to an easy name (as usually done for numpy or matplotlib), e.g. by
"""

import pygimli as pg
import numpy as np
import matplotlib.pyplot as plt

###############################################################################
# The modelling class is derived from ModellingBase, a constructor is defined
# and the response function is defined. Due to the linearity of the problem we
# store the matrix :math:`{\bf A}`, which is also the Jacobian matrix and use
# it for the forward calculation. A second function is just added as reference.
# We overwrite the method createJacobian as we know it but do nothing in the
# actual computation. If :math:`{\bf J}` depends on :math:`{\bf m}` this
# function must be filled.


class FunctionModelling(pg.core.ModellingBase):
    def __init__(self, nc, xvec, verbose=False):
        pg.core.ModellingBase.__init__(self, verbose)
        self.x_ = xvec
        self.nc_ = nc
        nx = len(xvec)
        self.regionManager().setParameterCount(nc)
        self.jacobian().resize(nx, nc)
        for i in range(self.nc_):
            self.jacobian().setCol(i, pg.math.pow(self.x_, i))

    def response(self, model):
        return self.jacobian() * model

    def responseDirect(self, model):
        y = pg.Vector(len(self.x_), model[0])

        for i in range(1, self.nc_):
            y += pg.math.pow(self.x_, i) * model[i]

        return y

    def createJacobian(self, model):
        pass  # if J depends on the model you should work here

    def startModel(self):
        return pg.Vector(self.nc_, 0.5)


###############################################################################
# Let us create some synthetic data for some x values

x = np.arange(0., 10., 0.5)
y = 1.1 + 2.1 * x - 0.2 * x**2
noise = 0.1
y += np.random.randn(len(y)) * noise

###############################################################################
# We now start by setting up the modelling operator, and inversion and run it.

fop = FunctionModelling(3, x)

# initialize inversion with data and forward operator and set options
inv = pg.Inversion(y, fop)

# constant absolute error of 0.01 is 1% (not necessary, only for chi^2)
inv.setAbsoluteError(noise)

# the problem is well-posed and does not need any regularization
inv.setLambda(0)

# actual inversion run yielding coefficient model
coeff = inv.run()

inv.echoStatus()

print(coeff)

###############################################################################
# The result is easily plotted by

plt.plot(x, y, 'rx', x, inv.response(), 'b-')
plt.show()
