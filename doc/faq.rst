.. _sec:faq:

Frequently asked questions
==========================

.. note::

    This section is still under construction. If you have a have a question
    (about GIMLi) that is not covered here, please open a `GitHub issue
    <https://github.com/gimli-org/gimli/issues>`_.

.. contents::
    :local:
    :backlinks: top

General
-------

What is the difference between BERT and GIMLi?
..............................................

GIMLi is a general library for modeling and inversion in geophysics. BERT
builds upon this framework and provides numerous high-level and user-friendly
functions and applications specifically tailored for DC resistivity studies.

Installation
------------

Python 2 or Python 3?
.....................

Short answer: Python 3. Long answer: Currently, pygimli is functionable with
all major Python versions (2.7, 3.4, 3.5). When compiling from source, it is
important that *boost_python* is build against the same Python version you want
to use. However, many (scientific) Python projects recommend to use Python 3
for various reasons and we will also drop support for Python 2 eventually. In
most cases you can translate your existing Python 2 scripts to Python 3 by
running *2to3* (https://docs.python.org/2/library/2to3.html).

What do I have to do to use pygimli in Spyder?
..............................................

If you have successfully compiled pygimli, but can't use it in Spyder with the
command:

.. code:: bash

    python -c 'import pygimli as pg; print(pg.__version)'

    python: can't open file 'import pygimli as pg; print(pg.__version__)':
    [Errno 2] No such file or directory

This is probably related to a path issue. Make sure you set the the correct
ones:

.. code:: bash

    export PYTHONPATH=$PYTHONPATH:$HOME/src/gimli/gimli/python
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/src/gimli/build/lib
    export PATH=$PATH:$HOME/src/gimli/build/bin

Make sure to add these settings to your *.bashrc* to make them permanent.
