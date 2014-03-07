.. _sec:design:

Design
======

In applied geophysics, various physical processes and fields are used to gain 
information about subsurface parameters.
Fields and processes can be well studied and understood by simulation assuming 
a parameter distribution.
This so-called forward task can be done by using analytical solution and 
numerical integration, or solving partial differential equations with finite 
difference or finite element techniques.
In the recent years, very different studies have been presented that open up 
new methods to be used.

However, in almost all approaches the task is finally to derive subsurface 
parameters, i.e. the inverse problem has to be solved.
Very often this is ill-posed, i.e. a variety of solutions is fitting the data 
within error bounds.
Hence regularization methods have to be applied.
There exist numerous inversion and regularization schemes, which do not have to 
be reinvented.
Furthermore, resolution analysis is desirable in order to appraise the quality 
of the results.
The idea of :ref:`sec:GIMLi` is to present a very flexible framework for 
geophysical inversion and modelling such that it can be used in an abstract way 
for any forward operator.
All problems such as optimization of regularization parameters, line search are 
solved generally.
The GIMLi library is structured into four layers (Fig. :ref:`fig:gimliblock`) 
that are based on each other:

.. _fig:gimliblock:
.. figure:: tutorial/pics/gimliblock.*
    :align: center

    Scheme of the GIMLi library

.. describe:: The basic layer

    holds fundamental algebraic methods and mesh containers for model parameterisation

.. describe:: The modelling & region layer

    administrates the modelling classes that are based on a base class and the connection to constraints and transform functions

.. describe:: The inversion layer

     is a template class for minimisation with different methods, inverse solvers, line search, :math:`\lambda` optimisation and resolution analysis

.. describe:: Inversion frameworks

    sophisticated techniques are formulated abstractly, e.g. time-lapse strategies, laterally constrained or roll-along inversion, different kinds of joint inversion

External programs are, e.g., mesh generators and solvers for linear systems.
For generating quality constrained irregular 2d and 3d meshes, we usually use :term:`Triangle` :cite:`shewchuk96b` and :term:`TetGen` :cite:`tetgen`. 
However, there is a tutorial on how to incorporate :term:`Gmsh` meshes :cite:`GeuzaineRemacle2009`.
Regular meshes can be directly created.
For solving linear systems we use the open-source collection :term:`SuiteSparse` :cite:`ChenDavHag+2009`, which contains multi-frontal direct \& iterative solvers as well as reordering algorithms.

External forward operators can be easily linked against GIMLi.
As soon they meet the requirements, the inversion can be setup and run with 2 lines.
With some additional effort, applications can easily be created either using C++ or, more easily, using Python scripts.
So far, various forward calculations are already included in GIMLi:
    * different 1d electromagnetic methods: VES, FDEM, TDEM, MT, TDR
    * first-arrival traveltime (refraction)
    * gravimetry
    * various fitting functions for 
