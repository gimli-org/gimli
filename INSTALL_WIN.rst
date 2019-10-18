.. _sec:install_win:
Installation on Windows
-----------------------

You need a Python installation with the modules numpy and matplotlib.
There are minimalistic installers and pre-packaged distributions.
We recommend `Anaconda <http://www.continuum.io/>`_ 64bit but `WinPython <http://winpython.github.io/#releases>`_ will to the job too.

Binary installers
.................

To avoid the building the C++ core, we provide binary packages built against Anaconda,
but working as well with WinPython.
Make sure you have an appropriate Python version (3.5.x or 3.6.x, 64bit) installed.
There are two ways, exe installers and wheels. The latter can be easily installed by a
package manager like the WinPython Control Panel (located in the WinPython main directory).

..  image:: https://img.shields.io/badge/pyGIMLi1.0.11_win64-Python_3.7_Wheel-green.svg
   :target: https://github.com/gimli-org/gimli/releases/download/v1.0.11/pygimli-1.0.11-cp37-cp37m-win_amd64.whl
..  image:: https://img.shields.io/badge/pyGIMLi1.0.11_win64-Python_3.6_Wheel-green.svg
   :target: https://github.com/gimli-org/gimli/releases/download/v1.0.11/pygimli-1.0.11-cp36-cp36m-win_amd64.whl

..  image:: https://img.shields.io/badge/pyGIMLi1.0.10_win64-Python_3.7_Wheel-green.svg
   :target: https://github.com/gimli-org/gimli/releases/download/v1.0.10/pygimli-1.0.10-cp37-cp37m-win_amd64.whl
..  image:: https://img.shields.io/badge/pyGIMLi1.0.10_win64-Python_3.6_Wheel-green.svg
   :target: https://github.com/gimli-org/gimli/releases/download/v1.0.10/pygimli-1.0.10-cp36-cp36m-win_amd64.whl
..  image:: https://img.shields.io/badge/pyGIMLi1.0.9_win64-Python_3.6_Wheel-green.svg
   :target: https://github.com/gimli-org/gimli/releases/download/v1.0.9/pygimli-1.0.9-cp36-cp36m-win_amd64.whl

In the package manager just select the whl file and press Install to install or upgrade.
If there is no package manager you can install the wheel by calling pip on the command line
(e.g. on the Windows command prompt, in the WinPython control shell or under MSYS2/cygwin):

.. code-block:: bash

    pip install pygimli-1.0.10-cp37-cp37m-win_amd64.whl # if pip is not in path use: python -m pip

See also `Issue #76 <https://github.com/gimli-org/gimli/issues/76>`_ for screenshots.
Uninstall can be made also with the Control Panel or by pip uninstall.

Alternatively you can use the MSI installer where you can choose your Python installation

..  image:: https://img.shields.io/badge/pyGIMLi1.0.11_win64-MSI_for_Python_3.7-green.svg
   :target: https://github.com/gimli-org/gimli/releases/download/v1.0.11/pygimli-1.0.11.win-amd64-py37.msi
..  image:: https://img.shields.io/badge/pyGIMLi1.0.11_win64-MSI_for_Python_3.6-green.svg
   :target: https://github.com/gimli-org/gimli/releases/download/v1.0.11/pygimli-1.0.11.win-amd64-py36.msi

..  image:: https://img.shields.io/badge/pyGIMLi1.0.10_win64-MSI_for_Python_3.7-green.svg
   :target: https://github.com/gimli-org/gimli/releases/download/v1.0.10/pygimli-1.0.10.win-amd64-py37.msi
..  image:: https://img.shields.io/badge/pyGIMLi1.0.10_win64-MSI_for_Python_3.6-green.svg
   :target: https://github.com/gimli-org/gimli/releases/download/v1.0.10/pygimli-1.0.10.win-amd64-py36.msi
..  image:: https://img.shields.io/badge/pyGIMLi1.0.9_win64-MSI_for_Python_3.6-green.svg
   :target: https://github.com/gimli-org/gimli/releases/download/v1.0.9/pygimli-1.0.9.win-amd64-py36.msi

As a result you should be able to import pygimli in any script from any location.
Since the C++ core changes only rarely, you can also check out the code via git, copy the
pyd/dll files from the core directory of the binary distribution (lib/site-packages/pygimli)
to the core directory of the git clone, and set PYTHONPATH to the latter.

Building pyGIMLI from source
............................

First, you need a Linux-like command shell along with a gcc compiler.
Although there might be different solutions (Cygwin, Git Bash, MinGW/MSYS),
we only support the MSYS2 (Minimal System 2) hosted at http://www.msys2.org.
As computers and modern Windows (>=7) are 64bit we only test this.
Avoid installing into strange Windows folders, e.g., c:\ProgramData\mys64 is fine.

After installing MSYS, start the console once so it builds your personal home
directory where you find a .bashrc file, e.g. in

.. code-block:: bash

    c:\ProgramData\mys64\home\YOUR_USERNAME

Edit .bashrc so that the Anaconda or WinPython installation path is added to your default
PATH.

e.g.:

.. code-block:: bash

    export PATH=/c/ProgramData/Anaconda3

This is necessary since gimli needs to know valid python installation and
version. Ideally the following one-liner will suffice to compile pyGIMLi in the
current directory.

**Note: The script will automatically take care of requirements and updates of MSYS2.
And also needs to modify/patch some of the llvm system files.**

.. code:: bash

    curl -Ls install.pygimli.org | bash

This script accepts a few more options. For help see

.. code:: bash

    curl -Ls install.pygimli.org | bash -s help

If everything runs fine, including some tests, the script will tell you some
additional PATH and PYTHONPATH settings for your .bashrc to use pygimli inside
the console or any IDE like Spyder (coming along with WinPython).

If something goes wrong, please take a look on the error message.

You can alse try the following instructions for manual installation.

Manual installation
...................

Make sure to have an updated msys2 environment. Run at least:

.. code-block:: bash
    pacman -Su
    pacman -Sy

to update your local package databases. See https://sourceforge.net/p/msys2/wiki/MSYS2%20installation/
for further instructions.

To get a complete working toolchain you need some packages installed.

.. code-block:: bash

    pacman -S make tar git subversion mercurial unzip wget patch

.. code-block:: bash

    pacman -S mingw-w64-x86_64-cmake mingw-w64-x86_64-gcc mingw-w64-x86_64-gcc-fortran
    pacman -S mingw-w64-x86_64-openblas mingw-w64-x86_64-suitesparse
    pacman -S mingw-w64-x86_64-doxygen mingw-w64-x86_64-cppunit
    pacman -S mingw-w64-x86_64-llvm mingw-w64-x86_64-clang

The rest of the installation is like the linux way with some small differences.

Prepare the directory structure as described above:

The build is performed via cmake. While calling cmake *MSYS* users should tell
using the Unix makefile generator to find the correct gcc versions:

.. code-block:: bash

    cmake ../gimli -G 'Unix Makefiles'

If cmake complains about missing python stuff, make sure the Python interpreter
is in your execution path. 

To build the library, just run

.. code-block:: bash

    make -j2

You might add the option -jN to use a number of N CPUs in parallel.
To build pygimli, run

.. code-block:: bash

    make pygimli J=2

You might add J=N to use a number of N CPUs in parallel.
Building pygimli takes some time and you can grab a coffee (or two).
If it finishes without errors you just need to set the environment:
(note that pygimli is still built in-place, in pygimli/core)

.. code-block:: bash

    export PYTHONPATH=$PYTHONPATH:$HOME/src/gimli/gimli/python
    export PATH=$PATH:$HOME/src/gimli/build/lib
    export PATH=$PATH:$HOME/src/gimli/build/bin

Compiled binaries will be written to `build/bin`.

You can test the pygimli build with:

.. code-block:: bash

    python -c 'import pygimli as pg; print(pg.version())'


Using cmake with CodeBlocks
...........................

Codeblocks is a nice C++ IDE available on http://www.codeblocks.org/downloads/

Tested versions 13.12/16.01, each without integrated mingw but a real MinGW/MSYS.

To generate the codeblocks project files run

.. code-block:: bash

    cmake -G "CodeBlocks - MinGW Makefiles"

and open the libgimli.cbp with codeblocks. Set up your compiler and run Build All.
