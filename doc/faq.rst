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

GIMLi is a general library for modelling and inversion in geophysics. BERT
builds upon this framework and provides numerous high-level and user-friendly
functions and applications specifically tailored for DC resistivity studies.

Installation
------------

Python 2 or Python 3?
.....................

Short answer: Python 3. Long answer: Currently, pygimli is functionable with all
major Python versions (2.7, 3.4, 3.5, 3.6). When compiling from source, it is
important that *boost_python* is build against the same Python version you want
to use. However, many (scientific) Python projects recommend to use Python 3 for
`various reasons <http://python-3-for-scientists.readthedocs.io/>`_ and we will
also drop support for Python 2 eventually. In most cases you can translate your
existing Python 2 scripts to Python 3 by running `2to3
<https://docs.python.org/2/library/2to3.html>`_.

What do I have to do to use pygimli in Spyder?
..............................................

For Linux users it should be sufficient to have proper environment settings:

.. code:: bash

  export PYTHONPATH=$PYTHONPATH:$HOME/src/gimli/gimli/python
  export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/src/gimli/build/lib
  export PATH=$PATH:$HOME/src/gimli/build/bin

Note, these setting will ony affect applications that have been started from the
current console.
To make these settings permanent, make sure to add them to your *.bashrc*.

Windows users probably need to tell spyder where your pygimli
installation can be found by an appropriate setting of
the `PYTHONPATH` in the spyder Menubar:/tools/PYTHONPATH manager.

If you start spyder from Windows startmenu you probably need to set the
correct `PATH` somewhere in your windows personal settings.

You can test you installation within the console via:

.. code:: bash

  python -c 'import pygimli as pg; print(pg.version())'

If there is something wrong with the environment settings
you get an error like this:

.. code:: bash

  python: can't open file 'import pygimli as pg; print(pg.version())':
  [Errno 2] No such file or directory


Weird findings
--------------

My script called sip.py and nothing works
.........................................

Rename your file to something different. One of the prerequisite library (pyQT)
import the file sip.py as module for they own and just stuck.

Segfault on import
..................

Try `python -s -c "import pygimli"`. The `-s` option ensures that only system
packages are used. This avoids conflicts with local (pip) packages.

CXXABI_1.3.9 not included in libstdc++.so.6
...........................................
When installing conda packages on older machines, the above error may occur. If so, 
check out the suggestion made 
`here <https://github.com/ContinuumIO/anaconda-issues/issues/5191#issuecomment-368243432>`_.
