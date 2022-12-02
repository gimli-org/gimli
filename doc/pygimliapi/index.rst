:html_theme.sidebar_secondary.remove: true

.. _sec:pygimliapi:

=====================
pyGIMLi API Reference
=====================

.. only:: latex

  .. note::

    In the following, all Python modules, functions and classes are documented.
    For a reference of the C++ core library of the code, please visit:
    http://pygimli.org/gimliapi/.


.. currentmodule:: pygimli

.. rubric:: Module overview

.. autosummary::
  :template: module.rst
  :toctree: _generated
  :recursive:

  frameworks
  math
  meshtools
  physics
  solver
  testing
  utils
  viewer

.. note::

  This documentation is valid for |version|. Check your installed version with

  .. code-block:: python

    import pygimli as pg
    print(pg.__version__)

  and consider updating (e.g., via the conda package manager). If you have a
  newer version, the current documentation of the *dev* branch before
  the next release is available at: https://dev.pygimli.org.

