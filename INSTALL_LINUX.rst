Installation on Linux
---------------------

The following should suffice to build GIMLi from source on most Linux Platforms:

.. code-block:: bash

    mkdir -p ~/src/gimli && cd ~/src/gimli
    svn checkout https://svn.code.sf.net/p/libgimli/code/trunk

    mkdir -p build && cd build
    cmake ../trunk
    make pygimli

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


Note for gcc version > 5.0
^^^^^^^^^^^^^^^^^^^^^^^^^^

gccxml will not support gcc versions above 5.0 so we have to switch the xml caster to castxml to build pygimli.
This should be done via cmake automatically.
As an additional prerequisite castxml needs > clang-3.6.0 and > llvm-3.6.0 installed.
castxml will be installed automatically.


Usefull cmake settings
^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

    CC=clang CXX=clang++ cmake ../trunk

Define alternative python version.
On default the version of your active python version will be chosen.
You will need numpy and boost-python builds with your desired python version.

.. code-block:: bash

    cmake ../trunk -DPYVERSION=3.3


If you experience runtime problems on starting pygimli it may happen that CMake
estimates the wrong libboost_python version by choosing py2 version instead of py3.
You can force cmake to select the correct version with:

.. code-block:: bash

    cmake ../trunk -DBoost_PYTHON_LIBRARY=/usr/lib64/libboost_python3.so

Chose a xml caster for the python bindings.
Either gccxml (default for linux and gcc < 5) or castxml else.

.. code-block:: bash

    cmake ../trunk -DCASTER='gccxml'


Usefull make commands
^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

    make rebuild_thirdparty


