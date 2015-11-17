Installation on Linux
---------------------

The following should suffice to build GIMLi from source on most Linux Platforms:

.. code-block:: bash

    mkdir -p ~/src/gimli && cd ~/src/gimli
    git clone https://github.com/gimli-org/gimli.git trunk

    mkdir -p build && cd build
    cmake ../trunk
    make gimli pygimli apps

See below for more detailed compilation instructions.


Detailed Installation on Vanilla Debian
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

First install some of the necessary stuff. For sure you will need subversion to
get the source files and some things for the tool-chain:

.. code-block:: bash

    sudo apt-get install subversion git cmake mercurial
    sudo apt-get install libboost-all-dev libblas-dev liblapack-dev

If you want to use the pyGIMLi (Python scripts, bindings and apps):

.. code-block:: bash

    sudo apt-get install python-numpy python-matplotlib
    sudo apt-get install libedit-dev clang-3.6-dev llvm-3.6-dev python3-dev


Create a directory for your installation, e.g., $HOME/src

.. code-block:: bash

    mkdir -p ~/src
    cd src
    mkdir -p gimli
    cd gimli

Checkout the current sources for libgimli:

.. code-block:: bash

    git clone https://github.com/gimli-org/gimli.git trunk

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

If you want to use the C++ commandline applications, call

.. code-block:: bash

    make apps

Compiled binaries will be written to `build/bin`.

You can test the pygimli build with:

.. code-block:: bash

    python -c 'import pygimli as pg; print(pg.__version)'

You can test your libgimli build with:

.. code-block:: bash

    make check

Of course the test will be very silent if you don't have cppunit installed.


Example Installation on Ubuntu
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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

Troubleshooting
^^^^^^^^^^^^^^^

If you experience runtime problems on starting pygimli like:

.. code-block:: bash

    ImportError: /usr/lib/libboost_python.so: undefined symbol: PyClass_Type

It may happen that CMake estimates the wrong libboost_python version by choosing py2 version instead of py3.
You can force cmake to select the correct version with:

.. code-block:: bash

    cmake ../trunk -DBoost_PYTHON_LIBRARY=/usr/lib64/libboost_python3.so

If the build misses libedit:

.. code-block:: bash

    /usr/bin/ld: cannot find -ledit

Install *libedit*, e.g. 'apt-get install libedit' on Debian/Ubuntu.

Useful cmake settings
^^^^^^^^^^^^^^^^^^^^^

You can rebuild and update all local generated third party software by setting the CLEAN environment variable:

.. code-block:: bash

    CLEAN=1 cmake ../trunk

Use alternative c++ compiler.

.. code-block:: bash

    CC=clang CXX=clang++ cmake ../trunk

Define alternative python version.
On default the version of your active python version will be chosen.
You will need numpy and boost-python builds with your desired python version.

.. code-block:: bash

    cmake ../trunk -DPYVERSION=3.3

Build the library with debug and profiling flags

.. code-block:: bash

    cmake ../trunk -DCMAKE_BUILD_TYPE=Debug

Build the library with gcc build.in sanity check 

.. code-block:: bash

    cmake ../trunk -DCMAKE_BUILD_TYPE=Debug -DASAN=1


Usefull make commands
^^^^^^^^^^^^^^^^^^^^^

More verbose build output to view the complete command line:

.. code-block:: bash

    make VERBOSE=1
 


