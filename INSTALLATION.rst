Installation
============

On all platforms, we recommend to install pyGIMLi via the conda package manager
contained in the Anaconda distribution. For details on how to install Anaconda,
we refer to: https://docs.anaconda.com/anaconda/install/

To avoid conflicts with other packages, we recommend to install pygimli in a
separate environment. Here we call this environment `pg`, but you can give
it any name. Note that this environment has to be created only once.

Open a terminal (Linux & Mac) or the Anaconda Prompt (Windows) and type:

.. code-block:: bash

    conda create -n pg -c gimli -c conda-forge pygimli

If you want to use pygimli, you have to activate the environment. You can put
this line in your `~/.bashrc` file so that it is activated automatically if you
open a terminal.

.. code-block:: bash

    conda activate pg

To test if everything works correctly you can do the following:

.. code-block:: bash

    python -c "import pygimli; pygimli.test(show=False, onlydoctests=True)"

After that you can use pygimli with your text editor of choice and a terminal.
Depending on your preferences, you can also install third-party software such as
the MATLAB-like integrated development environment (https://www.spyder-ide.org):

.. code-block:: bash

    conda install -c conda-forge spyder

Or alternatively, the web-based IDE JupyterLab (https://jupyterlab.readthedocs.io).

.. code-block:: bash

    conda install -conda-forge jupyterlab

Update your pygimli installation frome time to time, if want to have the newest
functionality:

.. code-block:: bash

    conda update -c gimli -c conda-forge pygimli

The only drawback using conda is that you are bound to the rhythm we update the
binary packages. Conda also can be seen as a sandbox inside your system
and it might be difficult to combine system python packages and conda pyGIMLi.
If you like to keep your pyGIMLi version more recent (including all possible
drawbacks of versions that are actively developed) you should compile pyGIMli
using your systems toolchain.
