.. _sec:GIMLI:

About GIMLi
===========

Introduction
------------

GIMLi is an open-source multi-method library for solving
:ref:`inverse<sec:about_gimli_inversion>` and
:ref:`forward<sec:about_gimli_modelling>` modelling tasks.

What GIMLi is good for?:

* creating inversion applications (C++) and scripts (Python) for existing modules
* add your own forward calculations and build a robust inversion quickly
* combining different geophysical methods in various ways
* doing modelling of different PDEs

What GIMLi is **NOT** good for?:

* for people that expect a ready-made GUI for interpreting their data

What do you need to use GIMLi?:

* either the source code and some time to build the system or binaries
* preferrably Python with numpy and matplotlib, or a C++ compiler

.. _sec:authors:

Authors
-------

* **Carsten Rücker**

  *Berlin University of Technology, Department of Applied Geophysics, Berlin, Germany*

  carsten@pygimli.org

* **Thomas Günther**

  *Leibniz Institute for Applied Geophysics, Hannover, Germany*

  thomas@pygimli.org

* **Florian Wagner**

  *Helmholtz Centre Potsdam, GFZ German Research Centre for Geosciences, Potsdam, Germany*

  florian@pygimli.org

.. _sec:about_gimli_inversion:

Inversion
---------

One main task of GIMli is to carry out inversion, i.e. error-weighted minimization, for given forward routines and data.
Various types of regularization on meshes (1D,2D,3D) with regular or irregular arrangement are available.
There is a flexible control of all inversion parameters.

Please see :ref:`inversion tutorial<tut:inversion>` for examples and more details.

.. _sec:about_gimli_modelling:

Modelling
---------

Currently there is a Finite Element Framework for solving partial differential equations (PDE)

* 1D, 2D, 3D discretizations
* linear and quadratic shape functions (automatic shape function generator for possible higher order)
* Triangle, Quads, Tetrahedron, Prism and Hexahedron, mixed meshes
* solver for elliptic problems (Helmholtz-type PDE)

Please see :ref:`modelling tutorial<tut:modelling>` for examples and more details.

.. References are not displayed correctly in tex, so we leave this out for the moment
.. only:: html

    References
    ----------

    .. bibliography:: doc/biblio.bib
        :style: mystyle
        :all:

.. only:: latex

    License
    -------
    GIMLi is distributed under the terms of the **GNU GENERAL PUBLIC LICENSE**
    of the Free Software Foundation, Inc. (http://fsf.org/)

    Details can be found on the `GIMLi website
    <http://www.pygimli.org/COPYING.html>`_.
