##GIMLi - Geophysical Inversion and Modelling Library##
        
An open-source multi-method library for solving inverse and forward tasks.

#### Prerequisites ####

We need some more or less important tools to build libgimli
 
    . subversion, git
    . cmake >= 2.8.10
    . gcc >= 4.4

    tested: gentoo x86_64: gcc-4.4.5, gcc-4.5.3, gcc-4.5.4
                Debian 3.2.46-1 x86_64: gcc-4.7.2
                Ubuntu 
                MinGW32: gcc-4.6.2/4, gcc-4.7.2
                MinGW64: gcc-4.5.4 (without python bindings)

    . libboost >=1.46 (thread, [python])

    tested: 1.46, 1.48, 1.49, 1.51, 1.53(mingw)
                
#### Optional Prerequisites -- some can be installed via provided scripts ####

    . blas and lapack for suitesparse (system or script)
    . SuiteSparse http://www.cise.ufl.edu/research/sparse/SuiteSparse/
      tested: SuiteSparse-4.2.1.tar.gz (via script)

    . triangle http://www.cs.cmu.edu/~quake/triangle.html (via script)
    . cppunit
    . procps
    . gccxml, pygccxml and pyplusplus (via script)
    . matplotlib >=1.1.0
    . doxygen        


#### Example Installation on Vanilla Debian ####

First install some of the necessary stuff. For sure you will need subversion to get the source files and some things for the tool-chain

    sudo apt-get install subversion git cmake

    sudo apt-get install libboost-all-dev libblas-dev liblapack-dev

Optional install useful stuff:

    sudo apt-get install libcppunit-dev libprocps0-dev

    sudo apt-get install python-matplotlib

    sudo apt-get install doxygen

Create a directory for your installation, e.g., $HOME/src

    mkdir -p ~/src
    cd src
    mkdir -p gimli
    cd gimli

Checkout the current sources for libgimli:
    
    svn checkout https://svn.code.sf.net/p/libgimli/code/trunk

Getting necessary and optional external libraries. 
We need suitesparse and triangle. If they are not installed system-wide we provide a small script to obtain them.

    cd trunk
    cd external
    make 

Check if there are now a lib and an include path. With some content in it (e.g. libtriangle.a and libcholmod.a)
Then you can go back to the main gimli path:
    
    cd ../..
    
#### Dependencies for python bindings ####

Python bindings are generally a good idea because some tools depend on it. 
If you just want the libgimli library you can skip this part

The python binding files are generated automatic by using gccxml, pygccxml and pyplusplus
If you cannot install them with your distribution we provide a script for it.

    cd ..
    sh gimli/trunk/python/buildScripts/buildToolChain.sh

if 'sh COMMAND.sh' complains about missing pushd or popd try 'bash COMMAND.sh':

    bash gimli/trunk/python/buildScripts/buildToolChain.sh


##### Building with cmake #####

We test a new build system using cmake http://www.cmake.org/ that hopefully avoid a lot of problems from the past.
In the first, cmake provide out of source build so we recommend using a build directory beside the trunk path:

    cd gimli
    mkdir -p build
    
the main directory structure should looks like this:

    gimli/trunk
    gimli/build

change to the build path:

    cd build

and configure the build:
    
    cmake ../trunk

If the output complains some missing dependencies you want to install .. just install these and repeat the the last step. 

To build the library just run make
    
    make

The libraries will be installed in build/lib and some test applications are installed in build/bin

If you want to build the python bindings call
    
    make pygimli

the _pygimli_.so library will be copied into the source path ../trunk/python/pygimli. 
To use the gimli installation there have to be set some environment variables:

    export PYTHONPATH=$PYTHONPATH:$HOME/src/gimli/trunk/python
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/src/gimli/build/lib
    export PATH=$PATH:$HOME/src/gimli/build/bin

You can test the pygimli build with:

    python -c 'import pygimli as g; print g.versionStr()'

You can test your libgimli build with:

    make check

Of course the test will be very silent if you don't have cppunit installed.

If you have doxygen installed, you can create the api documentation:

    make html

### Installation on Windows ###

##### Windows i.e. Mingw systems #####

First install mingw and msys to get a proper gcc an and the msys console

mingw-4.5.0 & msys-1.0.15 automatic installer: http://sourceforge.net/projects/mingw/files/
    tested: mingw-get-inst-20100909.exe

    http://sourceforge.net/p/mingw/news/2013/07/graphical-installer-interface----new-snapshot-available/

    
The installation is common to the linux way with some small differences.

Prepare the directory structure like described above:
If you don't have a proper boost installation you can install them yourself:

    sh glimli/trunk/python/buildScripts/buildBoostWin32.sh

If you don't have blas and lapack you can install it via script

    cd gimli/external
    make lapack

The build is performed via cmake. While calling cmake *Mingw* users should be preferable generate for msys makefiles:

    cmake -G 'MSYS Makefiles' ../trunk

cmake provide an interactive configuration and fine tuning, e.g., for adjusting the boost-include and boost-library paths.

    cmake-gui ../trunk 

To build the library just run make
    
    make

just need to set the environment:

    export PYTHONPATH=$PYTHONPATH:$(HOME)/src/gimli/trunk/python
    export PATH=$PATH:$(HOME)/src/gimli/build/lib
    export PATH=$PATH:$(HOME)/src/gimli/build/bin



##### Using cmake with CodeBlocks #####

First, for sure, you need codeblocks from: http://www.codeblocks.org/downloads/26
    
    tested: codeblocks-10.05-setup.exe

to come ....

#### Example Installation on Ubuntu ####

    sudo apt-get install subversion git cmake
    sudo apt-get install libboost-all-dev libblas-dev liblapack-dev
    sudo apt-get install libcppunit-dev
    sudo apt-get install python-matplotlib
    sudo apt-get install doxygen

    mkdir -p ~/src/gimli
    cd ~/src/gimli
    svn checkout https://svn.code.sf.net/p/libgimli/code/trunk
    cd trunk/external/
    make
    cd ../../../
    bash gimli/trunk/python/buildScripts/buildToolChain.sh    
    cd gimli
    mkdir build
    cd build
    cmake ../trunk
    make
    make pygimli
