.. _sec:doc_todo:

=====================
*Documentation TODOs*
=====================

.. admonition:: See also

    Restructured Text (reST) and Sphinx CheatSheets:
 
    http://docs.geoserver.org/trunk/en/docguide/sphinx.html

    http://openalea.gforge.inria.fr/doc/openalea/doc/_build/html/source/sphinx/rest_syntax.html
        
    http://packages.python.org/an_example_pypi_project/sphinx.html

    Numpy/Scipy Docstring convention:

    https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt#docstring-standard
   
Needed external extensions for building this docu

:: 

    pip install numpydoc pybtex sphinxcontrib-programoutput

* somewhere explain difference between gimli and pygimli
* somewhere explain difference python bind of gimli and pygimli scripts

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
.. autofunction:: pygimli.meshtools.mesh.readHydrus2dMesh
.. autofunction:: pygimli.meshtools.mesh.readGmsh

