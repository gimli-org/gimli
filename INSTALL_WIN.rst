Installation on Windows
-----------------------

First, install the msys2 environment to get a proper gcc an and the msys
console (should be 64bit on most systems): http://msys2.github.io/

You will also need a native python installation if you want to use pygimli. We
recommend WinPython: http://winpython.github.io/#releases

After you installed the msys2 64bit environment start the console once so it
builds your personal home directory where you find a .bashrc. Probably under:

.. code-block:: bash

    c:/msys64/home/YOUR_USERNAME

Edit .bashrc and add the winpython installation path to your default PATH.

.. code-block:: bash

    export PATH=$PATH:/c/PATH_TO_YOUR_WINPYTHON/WinPython-64bit-3.5.1.2/python-3.5.1/

This is necessary since the gimli installation needs to know your valid python
installation.

Ideally the following one-liner will suffice to compile pyGIMLi in the current directory.

**Note: The script will automatic take care any needs or updates of the msys environment. 
And also needs to modify some of the llvm system files.**

.. code:: bash

    curl -Ls install.pygimli.org | bash 

This script accept a few more options. See for help:

.. code:: bash

    curl -Ls install.pygimli.org | bash -s help

If everything runs fine, inclusive some tests, the script will tell you some 
additional PATH and PYTHONPATH settings for your .bashrc to use pygimli inside
the console or any IDEs like spyder.

If there goes something wrong ensure to take a look on the error message. 

You can alse try the following instructions for manual installation. 

Manual installation
...................

Be sure to have an updated msys environment. Run at least:

.. code-block:: bash

    pacman -Sy

to update your local package databases. See: https://sourceforge.net/p/msys2/wiki/MSYS2%20installation/
for further update instructions.

To get a complete working toolchain you need some packages installed.

.. code-block:: bash

    pacman -S make tar git subversion mercurial unzip wget patch

.. code-block:: bash

    pacman -S mingw-w64-i686-cmake mingw-w64-i686-gcc mingw-w64-i686-gcc-fortran
    pacman -S mingw-w64-i686-openblas mingw-w64-i686-doxygen

.. code-block:: bash

    pacman -S mingw-w64-x86_64-cmake mingw-w64-x86_64-gcc mingw-w64-x86_64-gcc-fortran
    pacman -S mingw-w64-x86_64-openblas mingw-w64-x86_64-doxygen
    pacman -S mingw-w64-x86_64-llvm mingw-w64-x86_64-clang


The rest of the installation is like the linux way with some small differences.

Prepare the directory structure like described above:

The build is performed via cmake. While calling cmake *Mingw* users should be
preferable generate for msys makefiles:

.. code-block:: bash

    cmake -G 'MSYS Makefiles' ../trunk

If cmake complains missing python stuff, ensure the python interpreter you want
to use, is in your execution path.

To build the library just run

.. code-block:: bash

    make

To build pygimli just run

.. code-block:: bash

    make pygimli

This will take some time and you can grab a coffee (or two).
If it finishs without any errors you just need to set the environment:

.. code-block:: bash

    export PYTHONPATH=$PYTHONPATH:$HOME/src/gimli/trunk/python
    export PATH=$PATH:$HOME/src/gimli/build/lib
    export PATH=$PATH:$HOME/src/gimli/build/bin

If you want to use the C++ commandline applications, call

.. code-block:: bash

    make apps

Compiled binaries will be written to `build/bin`.

You can test the pygimli build with:

.. code-block:: bash

    python -c 'import pygimli as pg; print(pg.__version__)'


Using cmake with CodeBlocks
...........................

First, for sure, you need codeblocks from: http://www.codeblocks.org/downloads/26

tested: 10.05, 12.11, 13.12, each without integrated mingw but a real mingw installation

To generate the codeblocks project files run

.. code-block:: bash

    cmake -G "CodeBlocks - MinGW Makefiles"

and open the libgimli.cbp with codeblocks. Set up your compiler and run Build All.

First install :term:`MinGW` and :term:`MSYS` to get a proper :term:`gcc` and a nice console

* mingw-4.5.0 & msys-1.0.15 automatic installer: http://sourceforge.net/projects/mingw/files/

tested: mingw-get-inst-20100909.exe

There is a new graphical installation and maintenance tool for MinGW which you should check out
http://sourceforge.net/p/mingw/news/2013/07/graphical-installer-interface----new-snapshot-available/
