.. _sec:glossary:

Glossary
========

.. glossary::
    GCC
        The GNU Compiler Collection (GCC). Default compiler for most linux
        systems. The Windows port is :term:`MinGW` See: http://gcc.gnu.org/

    GIMLi
        The eponymous software package. See :ref:`sec:GIMLi`

    pyGIMLi
        The eponymous software package. See :ref:`sec:pygimliapi`

    libGIMLi
        The eponymous software package. See :ref:`sec:GIMLi`

    Gmsh
        Gmsh: a three-dimensional finite element mesh generator with built-in
        pre- and post-processing facilities
        :cite:`GeuzaineRemacle2009`
        http://geuz.org/gmsh/.
        See: :py:func:`pygimli.meshtools.mesh.readGmsh`

    BERT
        Boundless electrical resistivity tomography http://www.resistivity.net/

    IPython
        An improved :term:`Python` shell that integrates nicely with
        :term:`Matplotlib`. See http://ipython.scipy.org/.

    Matplotlib
        Matplotlib :term:`Python` package displays publication quality results.
        It displays both 1D X-Y type plots and 2D contour plots for structured
        and unstructured data. It works on all common platforms and produces
        publication quality hard copies. http://matplotlib.org

    MinGW
        MinGW, a contraction of "Minimalist GNU for Windows", is a minimalist
        development environment for native Microsoft Windows applications. See:
        http://www.mingw.org/

    MSYS
        MSYS, a contraction of "Minimal SYStem", is a Bourne Shell command line
        interpreter system. Offered as an alternative to Microsoft's cmd.exe,
        this provides a general purpose command line environment, which is
        particularly suited to use with MinGW, for porting of many Open Source
        applications to the MS-Windows platform. See: http://www.mingw.org/

    NumPy
        The :mod:`numpy` :term:`Python` package provides array arithmetic
        facilities. See: http://docs.scipy.org/doc/numpy

    Paraview
        Is an open-source, multi-platform data analysis and visualization
        application. See: http://paraview.org/

    Python
        The programming language that :term:`pyGIMLi` (and your scripts) are
        written in. See: http://www.python.org/

    Pylab
        At the moment, the current combination of :term:`Python`,
        :term:`NumPy`, :term:`SciPy`, :term:`Matplotlib`, and IPython provide a
        compelling environment for numerical analysis and computation. See:
        http://wiki.scipy.org/PyLab

    SciPy
        Scientific Tools for Python is open-source software for mathematics,
        science, and engineering. See: http://wiki.scipy.org/SciPy

    Sphinx
        The tools used to generate the :term:`GIMLi` documentation. See:
        http://sphinx-doc.org

    SuiteSparse
        SuiteSparse is a single archive that contains packages for solving
        large sparse problems using Sparse Cholesky factorization.
        http://faculty.cse.tamu.edu/davis/suitesparse.html

    Triangle
        A Two-Dimensional Quality Mesh Generator and Delaunay Triangulator.
        :cite:`Shewchuk96b`
        http://www.cs.cmu.edu/~quake/triangle.html
        See: :py:func:`pygimli.meshtools.mesh.readTriangle`

    Tetgen
        A Quality Tetrahedral Mesh Generator and a 3D Delaunay Triangulator.
        :cite:`Si2004`
        http://tetgen.org/
        See: :py:func:`pygimli.meshtools.mesh.readTetgen`
