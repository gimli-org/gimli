.. _sec:install_mac:

Building on MAC OS
------------------

These installation instructions were proposed by Claudio Jordi (ETH Zurich) in a
`GitHub issue
<https://github.com/gimli-org/gimli/issues/46#issuecomment-357735129>`_. Since
we do not have access to Macs, this has not been tested thorougly. Please
comment in the issue and let us know if it works.
know it works in the corresponding GitHub issue.
For installation instructions for MacOS, please refer to:


In most cases, the following will suffice to compile pyGIMLi in the current
directory.

.. code:: bash

    #========================================
    # First install:
    # - Xcode (AppStore)
    # - XQuartz (https://www.xquartz.org)
    # - homebrew (https://brew.sh)
    #========================================

    # Install python3

    brew install python3

    brew install boost --with-python3

    brew install boost-python --with-python3

    # install some prerequisites that are not yet installed (might be more than what is here…)

    brew install mercurial
    brew install wget

    # install matplotlib, … using pip3 (for python3)

    pip3 install scipy
    pip3 install numpy
    pip3 install matplotlib

    # follow installation instructions from pygimli.org

    mkdir -p ~/src
    cd src
    mkdir -p gimli
    cd gimli
    git clone https://github.com/gimli-org/gimli.git
    mkdir -p build

    cd build

    # The version (here 3.6.4_2) needs to be set to the installed python3 version

    cmake ../gimli -DPYTHON_EXECUTABLE=/usr/local/bin/python3 -DPYTHON_INCLUDE_DIR=/usr/local/Cellar/python3/3.6.4_2/Frameworks/Python.framework/Versions/3.6/include/python3.6m -DPYTHON_LIBRARY=/usr/local/Cellar/python3/3.6.4_2/Frameworks/Python.framework/Versions/3.6/lib/libpython3.6.dylib -DPY_NUMPY=/usr/local/lib/python3.6/site-packages/numpy

    # This was needed for the compilation of some c++ stuff
    export CPLUS_INCLUDE_PATH=/usr/local/Cellar/python3/3.6.4_2/Frameworks/Python.framework/Versions/3.6/Headers

    make -j 8

    make pygimli

    curl -Ls install.pygimli.org | bash

.. note::

    Conda packages for Mac OS will follow soon.
