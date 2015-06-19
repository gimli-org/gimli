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
^^^^^^^^^^^^^^^^^^^^^^

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
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

These tools are required to create the Python bindings and are likely to be
outdated in your distribution and should therefore not be installed
system-wide. The build scripts will install them automatically.

* gccxml
* pygccxml
* pyplusplus
