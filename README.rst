========
Overview
========

    :ref:`GIMLi` - Geophysical Inversion and Modelling Library
    is an open-source multi-method library for solving inverse and forward tasks.

Tested with:
    gentoo x86_64 gcc-4.4.5
    gentoo x86_64 gcc-4.5.3
    mingw win32 gcc-4??
    mingw win64 gcc-4.5.4 (without python bindings)


BUILD (for Windows systems)

Please look into the directory mingw and the README file therein

BUILD (for posix systems such as Linux, Unix, MacOS):

Prerequisites:
You will need at least (incomplete for now):
>=automake-1.11
>=gcc-4.4.5 (just for some boost stuff that uses long long)


For the first time after initial svn checkout:

./autogen.sh
make

If you install from source-package
./configure
make

After further svn update:
make 


OPTIONAL EXTERNALS:

for the testsuite
cppunit (http://sourceforge.net/projects/cppunit)

support for light-weight linear solver:
amd (http://www.cise.ufl.edu/research/sparse/amd/)
ldl (http://www.cise.ufl.edu/research/sparse/ldl/)

support for advanced linear solver:
blas                    
lapack                  
colamd                  
camd                    
cholmod                 
(http://www.cise.ufl.edu/research/sparse/SuiteSparse/)

support for 2d mesh-generation:
triangle                

You can build triangle and suitesparse by calling the Makefile in the external directory


Windows toolchain:

mingw-4.5.0 & msys-1.0.15
    automatic installer: http://sourceforge.net/projects/mingw/files/
    testet: mingw-get-inst-20100909.exe
    
codeblocks
    from: http://www.codeblocks.org/downloads/26
    testet: codeblocks-10.05-setup.exe
    change codeblocks toolchain (setting/Compiler & Debugger../) 
    to c:/mingw or MinGW64 base directory

cmake >=2.8.7		

Add boost include path to the build:

CPPFLAGS='-I /c/home/boost/include' ./configure


Step by step: what I did on vanilla win7-64bin:
win7-64.prof, gcc-4.6.1, boost-1.48.0, python-2.7.2	

install:
	mingw (mingw-get-inst-20111118.exe) (toolchain)
	TortoiseSVN-1.7.4.22459-x64-svn-1.7.2 (libgimli sources)
	>=cmake-2.8.7 (blas, lapack, gccxml)
	python-2.7.2 (python binding, boost)
	codeblocks-10.05-setup
	http://www.sliksvn.com/en/download (cmd-line tool for pygccxml)

reboot

adjust ~/.profile:

> export PATH=/c/python27:$PATH

choose an installation path
e.g.: (note, path might be created first)
> cd /c/user/$USER/src

Getting libgimli sources:

> svn co https://libgimli.svn.sourceforge.net/svnroot/libgimli/trunk libgimli

Install third party dependencies: triangle, suitesparse

> cd libgimli/externals
> make

Install third party dependencies: blas, lapack

> cd libgimli/externals
> make lapack

copy blas and lapack into libgimli/mingw path

Install third party dependencies: boost
	
> wget -nc -nd http://sourceforge.net/projects/boost/files/boost/1.48.0/boost_1_48_0.tar.gz
> tar -xzvf boost_1_48_0.tar.gz
(maybee some errors using tar at command line, choose totalcommander instead)

fix boot sources that boost-python can be build

> sh libglimli/python/buildScripts/buildBoostWin32.sh

copy libboost-threads into libgimli/mingw path

build libgimli.dll by using codeblocks

> sh libglimli/python/buildScripts/buildToolChain.sh

build pygimli




