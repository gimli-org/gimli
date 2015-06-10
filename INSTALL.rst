.. _sec:install:

Installation
============

.. note::

    Binary installers for different platforms will be available for download
    soon. See below to build GIMLi from source.

Quick Installation
------------------

The following should suffice to build GIMLi from source on most Linux Platforms:

    .. code-block:: bash

        mkdir -p ~/src/gimli && cd ~/src/gimli
        svn checkout https://svn.code.sf.net/p/libgimli/code/trunk

        mkdir -p build && cd build
        cmake ../trunk
        make pygimli

See below for more detailed compilation instructions.

Prerequisites
-------------

To build GIMLi from source, the following tools are required:

* subversion, git(gccxml), mercurial(pygccxml), wget, tar
* cmake >= 2.8.10
* gcc >= 4.4

    tested on:

    * gentoo x86_64: gcc-4.4.5, gcc-4.5.3, gcc-4.5.4
    * Debian 3.2.46-1 x86_64: gcc-4.7.2
    * Ubuntu
    * Arch Linux
    * CentOS
    * MinGW32: gcc-4.6.2/4, gcc-4.7.2
    * MinGW64: gcc-4.5.4 (without python bindings)

For pygimli you additionally need:

* Python (2 and 3 are supported)
* numpy
* matplotlib >=1.1.0

Optional Prerequisites
......................
These tools can be installed system-wide with your native package manager (i.e.
apt-get). If not found, the provided build scripts will automatically download
and compile the necessary components.

* libboost >=1.46 (thread, python)

  tested: 1.46, 1.48, 1.49, 1.51, 1.52, 1.53, 1.57

* blas and lapack for suitesparse (system or auto via cmake)
* SuiteSparse (http://www.cise.ufl.edu/research/sparse/SuiteSparse)
* cppunit
* procps
* triangle (http://www.cs.cmu.edu/~quake/triangle.html)

Prerequisites automatically installed by the installer
......................................................
These tools are required to create the Python bindings and are likely to be
outdated in your distribution and should therefore not be installed
system-wide. The build scripts will install them automatically.

* gccxml
* pygccxml
* pyplusplus


Example Installation on Vanilla Debian
--------------------------------------

First install some of the necessary stuff. For sure you will need subversion to
get the source files and some things for the tool-chain:

.. code-block:: bash

    sudo apt-get install subversion git cmake mercurial
    sudo apt-get install libboost-all-dev libblas-dev liblapack-dev

If you want to use the pyGIMLi (Python scripts, bindings and apps):

.. code-block:: bash

    sudo apt-get install python-numpy python-matplotlib

Create a directory for your installation, e.g., $HOME/src

.. code-block:: bash

    mkdir -p ~/src
    cd src
    mkdir -p gimli
    cd gimli

Checkout the current sources for libgimli:

.. code-block:: bash

    svn checkout https://svn.code.sf.net/p/libgimli/code/trunk

We use cmake (http://www.cmake.org/) for compilation. We recommend using a
build directory parallel to the trunk path:

.. code-block:: bash

    mkdir -p build

The main directory structure should looks like this:

.. code-block:: bash

    gimli/trunk
    gimli/build

Change to the build path

.. code-block:: bash

    cd build

and configure the build:

.. code-block:: bash

    cmake ../trunk

If the output complains some missing dependencies you should install, just
install these and repeat the the last step.

To build the library just run make

.. code-block:: bash

    make

To speed up the build process using more CPUs, use the -j flag, e.g.:

.. code-block:: bash

    make -j 8

The libraries will be installed in build/lib and some test applications are
installed in build/bin

If you want to build the python bindings call

.. code-block:: bash

    make pygimli

The _pygimli_.so library will be copied into the source path
../trunk/python/pygimli. To use the gimli installation there have to be set
some environment variables:

.. code-block:: bash

    export PYTHONPATH=$PYTHONPATH:$HOME/src/gimli/trunk/python
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/src/gimli/build/lib
    export PATH=$PATH:$HOME/src/gimli/build/bin

You can test the pygimli build with:

.. code-block:: bash

    python -c 'import pygimli as pg; print(pg.__version)'

You can test your libgimli build with:

.. code-block:: bash

    make check

Of course the test will be very silent if you don't have cppunit installed.

Installation on Windows
-----------------------

Windows using MinGW systems
...........................

First install :term:`MinGW` and :term:`MSYS` to get a proper :term:`gcc` and a nice console

    * mingw-4.5.0 & msys-1.0.15 automatic installer: http://sourceforge.net/projects/mingw/files/

        tested: mingw-get-inst-20100909.exe

There is a new graphical installation and maintenance tool for MinGW which you should check out
        http://sourceforge.net/p/mingw/news/2013/07/graphical-installer-interface----new-snapshot-available/


The installation is common to the linux way with some small differences.

Prepare the directory structure like described above:
If you don't have a proper boost installation you can install them yourself:

.. code-block:: bash

    sh glimli/trunk/python/buildScripts/buildBoostWin32.sh

If you don't have blas and lapack you can install it via script

.. code-block:: bash

    cd gimli/external
    make lapack

The build is performed via cmake. While calling cmake *Mingw* users should be preferable generate for msys makefiles:

.. code-block:: bash

    cmake -G 'MSYS Makefiles' ../trunk

cmake provide an interactive configuration and fine tuning, e.g., for adjusting the boost-include and boost-library paths.

.. code-block:: bash

    cmake-gui ../trunk

To build the library just run make

.. code-block:: bash

    make

just need to set the environment:

.. code-block:: bash

    export PYTHONPATH=$PYTHONPATH:$(HOME)/src/gimli/trunk/python
    export PATH=$PATH:$(HOME)/src/gimli/build/lib
    export PATH=$PATH:$(HOME)/src/gimli/build/bin

Example Installation on Ubuntu
..............................

.. code-block:: bash

    sudo apt-get install subversion git cmake mercurial
    sudo apt-get install libboost-all-dev libblas-dev liblapack-dev
    sudo apt-get install python-matplotlib python-numpy

    mkdir -p ~/src/gimli
    cd ~/src/gimli
    svn checkout https://svn.code.sf.net/p/libgimli/code/trunk

    mkdir -p build
    cd build
    cmake ../trunk
    make gimli
    make pygimli
