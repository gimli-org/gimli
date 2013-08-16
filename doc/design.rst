.. _sec:design:

Design
======

In geophysics, various physical processes and fields are used to gain information about the subsurface parameters.
The fields and processes can be very well studied and understood by simulation assuming a parameter distribution.
This so-called forward task can be done by using analytical formulae and numerical integration, or numerical schemes deploying finite difference or finite element techniques.
In the recent years, very different studies have been presented that open up new methods to be used.

However, in almost all approaches the task is finally to derive subsurface parameters, i.e. the inverse problem has to be solved.
Very often this is ill-posed, i.e. a variety of solutions is fitting the data within error bounds.
Hence regularization methods have to be applied.
There exist numerous inversion and regularization schemes, which do not have to be reinvented.
Furthermore, often a resolution analysis is required in order to appraise the quality of the results.
The idea of :ref:`GIMLi` is to present a very flexible framework for geophysical inversion and modelling such that it can be used from any forward operator.
All problems such as optimization of regularization parameters, line search are solved generally.
The GIMLi library is structured into four layers (Fig. :ref:`fig:gimliblock`) that are based on each other:

.. describe:: The basic layer

    holds fundamental algebraic methods and mesh containers for model parameterisation

.. describe:: The modelling & region layer

    administrates the modelling classes that are based on a basis class and the connection to constraints and transform functions

.. describe:: The inversion layer

     is a template class for minimisation with different methods, inverse solvers, line search, :math:`\lambda` optimisation and resolution analysis

.. describe:: Inversion frameworks

    sophisticated techniques are formulated, e.g. time-lapse strategies, roll-along inversion or different kinds of joint inversion

.. _fig:gimliblock:
.. figure:: tutorial/pics/gimliblock.*
    :align: center

    Scheme of the GIMLi library

External programs are, e.g., mesh generators and solvers for linear systems.
For generating quality constrained irregular 2d and 3d meshes, we usually use \cw{Triangle}~\citep{triangle} and \cw{TetGen}~\citep{tetgen}.
For solving linear systems we recommend the open-source collection \cw{SuiteSparse} \citep{davis}, which contains multi-frontal direct\&iterative solvers as well as reordering algorithms.

External forward operators can be easily linked against GIMLi.
As soon they meet the requirements, the inversion can be setup and run with 2 lines.