#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

Polyfit
-------

This tutorial shows how to treat the simplest inversion case with an own
forward calculation. As regression is the easiest inverse problem we start with
fitting a polynomial of degree :math:`P`

.. math::

    f(x) = p_0 + p_1 x + \ldots + p_P x^P = \sum\limits_{i=0}^{P} p_i x^i
    
to some existing data :math:`y`.
The unknown model is the coefficient vector :math:`\m=[p_0,\ldots,p_P]`.
The vectorized function for a vector
:math:`\arr{x}=\transpose{[x_1,\ldots,x_N]} \grad \u{1}`
can be written as matrix-vector product

.. _eq:yAx:
.. math::

  \f(\arr{x}) = \A \arr{x} \quad\mbox{with}\quad \A=\left[ \begin{array}{cccc}
  1 & x_1 & \ldots & x_1^P \\
  \vdots & \vdots & \ddots & \vdots \\ 1 & x_N & \ldots & x_N^P
  \end{array} \right] =
  [ {\bf 1}\quad \arr{x} \quad \arr{x}^2 \ldots \arr{x}^P ] \;.

We set up the modelling operator, i.e. to return :math:`\f(\arr{x})` for given
:math:`p_i`, as a class derived from the modelling base class.
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

print(pg.__version__)

"""
As a result, all :ref:`sec:api` objects (classes and functions) can be referred
to with a preceding `pg.`, e.g., printing the version string for gimli.

Next, the modelling class is derived from ModellingBase, a constructor is
defined and the response function is defined.

The pygimli library must once be imported (in this case under the name g) and
all classes (e.g. modelling operators) can be used by pg.classname, e.g.
pg.RVector is the already known vector of real (double) values.

The main program is very easy then and the code is very similar to C++.
Data are loaded, both forward operator and inversion are created.
Inversion options are set and it the result of run is save to a file.
That's it.

As a main advantage of Python, the actual computations can be easily combined
with post-processing or visualization, even building graphical user-interfaces.
In this code example we use matplotlib, a plotting library inside of pylab, a
compound of different routines for numerics and plotting, comparable to MatLab.
"""


class FunctionModelling(pg.ModellingBase):
    def __init__(self, nc, xvec, verbose=False):
        pg.ModellingBase.__init__(self, verbose)
        self.x_ = xvec
        self.nc_ = nc
        self.regionManager().setParameterCount(nc)

    def response(self, par):
        y = pg.RVector(len(self.x_), par[0])
        
        for i in range(1, self.nc_):
            y += pg.pow(self.x_, i) * par[i]
            
        return y

    def startModel(self):
        return pg.RVector(self.nc_, 0.5)


"""
Let us create some synthetic data for some x values
"""

x = np.arange(0., 10., 1)
y = 1.1 + 2.1 * x

"""
We now start by setting up the modelling operator, and inversion and run it.
"""

nP = 3

# two coefficients and x-vector (first data column)
fop = FunctionModelling(nP, x)

# initialize inversion with data and forward operator and set options
inv = pg.RInversion(y, fop)

# constant absolute error of 0.01 is 1% (not necessary, only for chi^2)
inv.setAbsoluteError(0.01)

# the problem is well-posed and does not need regularization
inv.setLambda(0)

# actual inversion run yielding coefficient model
coeff = inv.run()

print(coeff)

"""
The result is easily plotted by
"""

plt.plot(x, y, 'rx', x, inv.response(), 'b-')

"""

.. image:: PLOT2RST.current_figure

"""

plt.show()
