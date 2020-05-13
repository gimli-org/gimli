About pyGIMLi
=============

Introduction
------------

pyGIMLi is an open-source library for
:ref:`modelling<sec:about_gimli_modelling>` and
:ref:`inversion<sec:about_gimli_inversion>` and in geophysics. The
object-oriented library provides management for structured and unstructured
meshes in 2D and 3D, finite-element and finite-volume solvers, various
geophysical forward operators, as well as Gauss-Newton based frameworks for
constrained, joint and fully-coupled inversions with flexible regularization.

What is pyGIMLi suited for?

* analyze, visualize and invert geophysical data in a reproducible manner
* forward modelling of (geo)physical problems on complex 2D and 3D geometries
* inversion with flexible controls on a-priori information and regularization
* combination of different methods in constrained, joint and fully-coupled inversions
* teaching applied geophysics (e.g. in combination with `Jupyter notebooks <http://jupyter-notebook.readthedocs.io/en/latest/notebook.html#notebook-documents>`_)

What is pyGIMLi **NOT** suited for?

* for people that expect a ready-made GUI for interpreting their data

.. _sec:authors:

Authors
-------

.. image:: https://img.shields.io/github/contributors/gimli-org/gimli.svg
   :target: https://github.com/gimli-org/gimli/graphs/contributors

We gratefully acknowledge all contributors to the pyGIMLi open-source project and look forward to `your contribution <https://pygimli.org/contrib.html>`_!

.. include:: ../AUTHORS.rst

.. _sec:about_gimli_inversion:

Inversion
---------

One main task of pyGIMli is to carry out inversion, i.e. error-weighted
minimization, for given forward routines and data. Various types of
regularization on meshes (1D, 2D, 3D) with regular or irregular arrangement are
available. There is flexible control of all inversion parameters. The default
inversion framework is based on the generalized Gauss-Newton method.

Please see :ref:`tut:inversion` for examples and more
details.

.. _sec:about_gimli_modelling:

Modelling
---------

pyGIMLi comes with various geophysical forward operators, which can directly be
used for a given problem. In addition, abstract finite-element and finite-volume
interfaces are available to solve custom PDEs on a given mesh. See
:mod:`pygimli.physics` for a collection of forward operators and
:mod:`pygimli.solver` for the solver interface.

The modelling capabilities of pyGIMLi inlcude:

* 1D, 2D, 3D discretizations
* linear and quadratic shape functions (automatic shape function generator for possible higher order)
* Triangle, Quads, Tetrahedron, Prism and Hexahedron, mixed meshes
* solver for elliptic problems (Helmholtz-type PDE)

Please see :ref:`tut:modelling` for examples and more details.

License
-------
pyGIMLi is distributed under the terms of the **Apache 2.0** license. Details on
the license agreement can be found `here
<https://www.pygimli.org/license.html>`_.
