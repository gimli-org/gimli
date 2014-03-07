.. _sec:doc_todo:

===================
Documentation TODOs
===================

* somewhere explain difference between gimli and pygimli
* somewhere explain difference python bind of gimli and pygimli scripts
* transfer gmsh tutorial

.. todolist::

.. admonition:: See also

    **Restructured Text (reST) and Sphinx CheatSheets**

    :ref:`rst-basics`
 
    http://docs.geoserver.org/trunk/en/docguide/sphinx.html

    http://openalea.gforge.inria.fr/doc/openalea/doc/_build/html/source/sphinx/rest_syntax.html
        
    http://packages.python.org/an_example_pypi_project/sphinx.html


    **Conventions**

    PEP 8:

    http://www.python.org/dev/peps/pep-0008/

    Numpy/Scipy Docstring convention:

    https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt#docstring-standard
   
    Matplotlib naming convention:

    http://matplotlib.org/1.3.0/faq/usage_faq.html#coding-styles

**Requirements**

Needed external extensions for building this docu

:: 

    pip install numpydoc pybtex sphinxcontrib-programoutput


Linking to Python API
---------------------

To link to a specific module or function in the Python API from the docs,
you can use

::

    See :py:func:`pygimli.meshtools.grid.appendTriangleBoundary` in :py:mod:`pygimli.meshtools` for more details.

Which results in

See :py:func:`pygimli.meshtools.grid.appendTriangleBoundary` in :py:mod:`pygimli.meshtools` for more details.

.. seealso:: http://sphinx-doc.org/domains.html#cross-referencing-python-objects

Docstring examples
------------------
Functions
"""""""""
.. autofunction:: pygimli.meshtools.mesh.readHydrus2dMesh
.. autofunction:: pygimli.meshtools.mesh.readGmsh

Classes
"""""""
.. warning:: Usage and source link do not show when displaying classes, why?
.. autoclass:: pygimli.mplviewer.CellBrowser
