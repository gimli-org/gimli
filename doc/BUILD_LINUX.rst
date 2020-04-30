.. _sec:install_lin:

Installation on Linux
---------------------

Pre-build binary install with conda
...................................

.. only:: html

  .. image:: https://anaconda.org/gimli/pygimli/badges/installer/conda.svg
      :target: https://conda.anaconda.org/gimli

  .. image:: https://anaconda.org/gimli/pygimli/badges/downloads.svg
      :target: https://anaconda.org/gimli/pygimli

.. raw:: html

    <br><br>

On Linux platforms, the most comfortable way to install pygimli is via the conda
package manager contained in the `Anaconda distribution
<https://www.continuum.io/downloads#linux>`_. Anaconda is scientific Python
distribution with more than 100 Python packages included (~400 Mb). You can also
use the `lightweight alternative Miniconda <https://conda.io/miniconda.html>`_
(~35 Mb) and only install the packages you like to use.

Install Miniconda (only once):

.. code-block:: bash

    wget http://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
    chmod +x miniconda.sh
    ./miniconda.sh -b -p $HOME/miniconda
    export PATH=$HOME/miniconda/bin:$PATH # Note: Add this to your .bashrc for permanent use

To avoid conflicts with other packages, we recommend to install pygimli in a
separate environment. Here we call this environment `mypgenv`, but you can give
it any name. Note that this environment has to be created only once.

.. code-block:: bash

    conda create -n mypgenv -c gimli -c conda-forge pygimli python=3.6 -y

If you want to use pygimli, you have to activate the environment. You can put
this line in your `~/.bashrc` file so that it is activated automatically if you
open a terminal.

.. code-block:: bash

    conda activate mypgenv

To test if everything works correctly you can do the following:

.. code-block:: bash

    python -c "import pygimli; pygimli.test(show=False, onlydoctests=True)"

After that you can use pygimli with your text editor of choice and a terminal.
Depending on your preferences, you can also install third-party software such as
MATLAB-like integrated development environment (https://www.spyder-ide.org):

.. code-block:: bash

    conda install -c conda-forge spyder

Or alternatively, the web-based IDE JupyterLab (https://jupyterlab.readthedocs.io).

.. code-block:: bash

    conda install -conda-forge jupyterlab

Update your pygimli installation frome time to time, if want to have the newest
functionality:

.. code-block:: bash

    conda update pygimli

The only drawback using conda is that you are bound to the rhythm we update the
binary packages. Conda also can be seen as a sandbox Linux inside your system
and it might be difficult to combine system python packages and conda pyGIMLi.
If you like to keep your pyGIMLi version more recent (including all possible
drawbacks of versions that are actively developed) you should compile pyGIMli
using your systems toolchain.

Compile your own with the curl installer
........................................

If you are not using Anaconda, you can build pyGIMLi from source in the current
directory via:

.. code-block:: bash

   curl -Ls install.pygimli.org | bash

This script accepts a few more options. For help see:

.. code-block:: bash

    curl -Ls install.pygimli.org | bash -s help

If something goes wrong please take a look at the error message. In most cases
there are missing or outdated packages. Please have a look at the prerequisites
tab.

If the installation fails you can try the following instructions for manual
installation.

Detailed Installation on Debian Stretch
.......................................

In order to build pygimli (and gimli) for Python 2.7 and Python 3.5, install
the required packages:

.. code-block:: bash

    sudo apt-get install wget subversion git cmake mercurial \
        libboost-all-dev libblas-dev liblapack-dev \
        python python-setuptools \
        python-numpy python-matplotlib \
        libedit-dev clang llvm-dev python3-dev \
        python3  python3-numpy python3-matplotlib \
        python3-setuptools

Create a directory for your installation, e.g., $HOME/src

.. code-block:: bash

    mkdir -p ~/src
    cd src
    mkdir -p gimli
    cd gimli

Checkout the current sources for libgimli:

.. code-block:: bash

    git clone https://github.com/gimli-org/gimli.git

We use `cmake <https://cmake.org>`_ for compilation. We recommend using a
build directory parallel to the gimli (trunk) path:

.. code-block:: bash

    mkdir -p build

The main directory structure should looks like this:

.. code-block:: bash

    gimli/gimli
    gimli/build

Change to the build path

.. code-block:: bash

    cd build

and configure the build for Python 2.7 with:

.. code-block:: bash

    cmake ../gimli

If you want to compile for Python 3.5, alternatively use:

.. code-block:: bash

    cmake ../gimli -DPYVERSION=3.5

If the output complains about missing dependencies, install these and repeat
the the last step. To build the library just run `make`.

.. code-block:: bash

    make

To speed up the build process using more CPUs, use the `-j` flag, e.g.:

.. code-block:: bash

    make -j 8

The libraries will be installed in **build/lib** and some test applications are
installed in build/bin. If you want to build the Python bindings, call:

.. code-block:: bash

    make pygimli

You might add J=8 (`make pygimli J=8`) for using 8 jobs in parallel to speed up
the build (adapt this to the number of real cores of the computer). The library
_pygimli_.so library will be copied into the source path
**../gimli/python/pygimli** in the subdirectory core.

To use the gimli installation you need to set some environment variables (this
example assumes that the **src** directory resides in your home directory):

.. code-block:: bash

    export PYTHONPATH=$PYTHONPATH:$HOME/src/gimli/gimli/python
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/src/gimli/build/lib
    export PATH=$PATH:$HOME/src/gimli/build/bin

If you want to use the C++ command line applications, call

.. code-block:: bash

    make apps

Compiled binaries will be written to `build/bin`.

You can do a quick test of the pygimli build and installation with:

.. code-block:: bash

    python -c 'import pygimli as pg; print(pg.__version__)'

You can test your gimli build with:

.. code-block:: bash

    make check

Note that the test will be very silent if you don't have *cppunit* installed.

If you install pytest with

.. code-block:: bash

    sudo apt-get install python-pytest python3-pytest

then you can run the internal test suite with

.. code-block:: bash

    python -c "import pygimli; pygimli.test()"

Using Docker to build in Debian stretch (for advanced users only!)
..................................................................

If you want to use a Docker container to build (and possibly use) pyGIMLi, you
can use the Dockerfile found in the **scripts/** subdirectory named
*Dockerfile_DebianStretch*. Please refer to the file for further instructions.

Example Installation on Ubuntu
..............................

.. code-block:: bash

    sudo apt-get install libc-dev subversion git cmake mercurial
    sudo apt-get install libboost-all-dev libblas-dev liblapack-dev libedit-dev
    sudo apt-get install python3-dev python3-matplotlib python3-numpy

    mkdir -p ~/src/gimli
    cd ~/src/gimli
    git clone https://github.com/gimli-org/gimli.git

    mkdir -p build
    cd build
    cmake ../gimli
    make -j 4 gimli
    make pygimli J=4

Troubleshooting
...............

If you experience runtime problems on starting pygimli like:

.. code-block:: bash

    ImportError: /usr/lib/libboost_python.so: undefined symbol: PyClass_Type

It may happen that CMake estimates the wrong libboost_python version by choosing py2 version instead of py3.
You can force cmake to select the correct version with:

.. code-block:: bash

    cmake ../gimli -DBoost_PYTHON_LIBRARY=/usr/lib64/libboost_python3.so

If the build misses libedit:

.. code-block:: bash

    /usr/bin/ld: cannot find -ledit

Install *libedit*, e.g. 'apt-get install libedit' on Debian/Ubuntu.


castXML
.......

castXML (https://github.com/CastXML/CastXML/) is needed to generate the code for the python bindings.
Some systems provide castxml binary so the build system should detect it if installed.
As fallback solution the build system tries to install castxml binaries or try to compile there own if the binaries don't work.
You can enforce the local binary installation with:

.. code-block:: bash

    cmake ../../src/castXML/ -DCASTXML_LOCAL=1
    make

or the local binary compilation with:

.. code-block:: bash

    cmake ../../src/castXML/ -DCASTXML_LOCALSRC=1
    make


If castXML build complains about missing clang or llvm command, please go into
$(GIMLISRC)/../thirdParty/build-XXX-XXX/castXML and try configure and build cmake manually

.. code-block:: bash

    CC=clang-3.6 CXX=clang++-3.6 cmake ../../src/castXML/
    make

If you build castXML manually you can provide this binary to cmake via

.. code-block:: bash

    cmake ../gimli -DCASTER_EXECUTABLE=$(PATH_TO_CASTXML)


Useful cmake settings
.....................

You can rebuild and update all local generated third party software by setting
the CLEAN environment variable:

.. code-block:: bash

    CLEAN=1 cmake ../gimli

Use alternative c++ compiler.

.. code-block:: bash

    CC=clang CXX=clang++ cmake ../gimli

Define alternative python version. On default the version of your active python
version will be chosen. You will need numpy and boost-python builds with your
desired python version.

.. code-block:: bash

    cmake ../gimli -DPYVERSION=3.3

Build the library with debug and profiling flags

.. code-block:: bash

    cmake ../gimli -DCMAKE_BUILD_TYPE=Debug

Build the library with gcc build.in sanity check

.. code-block:: bash

    cmake ../gimli -DCMAKE_BUILD_TYPE=Debug -DASAN=1


Useful make commands
....................

More verbose build output to view the complete command line:

.. code-block:: bash

    make VERBOSE=1