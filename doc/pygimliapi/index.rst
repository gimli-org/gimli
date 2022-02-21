.. _sec:pygimliapi:
=====================
pyGIMLi API Reference
=====================

.. toctree::
   :maxdepth: 1
   :hidden:

.. only:: latex

  .. note::

    In the following, all Python modules, functions and classes are documented.
    For a reference of the C++ core library of the code, please visit:
    http://pygimli.org/gimliapi/.

.. currentmodule:: pygimli

.. automodule:: pygimli

.. rubric:: Module overview

.. autosummary::
  :toctree: _generated
  :template: module.rst

  .. core
  frameworks
  .. math
  .. matrix
  meshtools
  physics
  solver
  testing
  .. trans
  utils
  viewer

.. only:: latex

  .. toctree::
    :glob:

    _generated/*

.. note::

  This documentation is valid for |version|. Check your installed version with

  .. code-block:: python

    import pygimli as pg
    print(pg.__version__)

  and consider updating (e.g., via the conda package manager). If you have a
  newer version, the current documentation of the *dev* branch before
  the next release is available at: https://dev.pygimli.org.
