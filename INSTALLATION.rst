Installation
============

|conda| |downloads| |platforms| |condaversion| |latest|

On all platforms, we recommend to install pyGIMLi via the conda package manager
contained in the Anaconda distribution. For details on how to install Anaconda,
we refer to: https://docs.anaconda.com/anaconda/install/

A **step-by-step guide for Windows users** can be found `here
<https://www.pygimli.org/_downloads/pygimli_win_install_guide.pdf>`_.

To avoid conflicts with other packages, we recommend to install pygimli in a
separate environment. Here we call this environment `pg`, but you can give
it any name. Note that this environment has to be created only once.

Open a terminal (Linux & Mac) or the Anaconda Prompt (Windows) and type:

.. code-block:: bash

    conda create -n pg -c gimli -c conda-forge pygimli=1.2.4

If you are using Windows or Mac, a new environment named "pg" should be visible
in the Anaconda Navigator. If you want to use pyGIMLi from the command line, you
have to activate the environment. You can put this line in your `~/.bashrc` file
so that it is activated automatically if you open a terminal.

.. code-block:: bash

    conda activate pg

After that you can use pyGIMLi with your text editor of choice and a terminal.

Usage with Spyder or JupyterLab
-------------------------------

Depending on your preferences, you can also install third-party software such as
the MATLAB-like integrated development environment (https://www.spyder-ide.org):

.. code-block:: bash

    conda install -c conda-forge spyder

Or alternatively, the web-based IDE JupyterLab (https://jupyterlab.readthedocs.io):

.. code-block:: bash

    conda install -c conda-forge jupyterlab

Testing
-------

To test if everything works correctly you can do the following:

.. code-block:: bash

    python -c "import pygimli; pygimli.test(show=False, onlydoctests=True)"

Staying up-to-date
------------------

Update your pyGIMLi installation from time to time, if want to have the newest
functionality:

.. code-block:: bash

    conda update -c gimli -c conda-forge pygimli

The only drawback of using conda is that you are bound to the rhythm in which we
update the conda packages. In order to work with the latest Python codes you
should create an environment with the latest pyGIMLi C++ core only,

.. code-block:: bash

    conda create -n pgcore -c gimli -c conda-forge pgcore
    
retrieve the source code by git

.. code-block:: bash

    git clone https://github.com/gimli-org/gimli
    cd gimli

and install pygimli as a development package

.. code-block:: bash

    conda develop .

Alternatively you could set the PYTHONPATH variable but you would have to care
for dependencies by yourself.

Later you can just update the pygimli code by

.. code-block:: bash

    git pull
    
Only if you need recent changes to the C++ core, you have to compile
pyGIMLi using your systems toolchain as described in 
https://www.pygimli.org/compilation.html#sec-build

.. |conda| image:: https://anaconda.org/gimli/pygimli/badges/installer/conda.svg
   :target: https://anaconda.org/gimli/pygimli
.. |downloads| image:: https://anaconda.org/gimli/pygimli/badges/downloads.svg
   :target: https://anaconda.org/gimli/pygimli
.. |condaversion| image:: https://anaconda.org/gimli/pygimli/badges/version.svg
   :target: https://anaconda.org/gimli/pygimli
.. |latest| image:: https://anaconda.org/gimli/pygimli/badges/latest_release_date.svg
   :target: https://anaconda.org/gimli/pygimli
.. |platforms| image:: https://anaconda.org/gimli/pygimli/badges/platforms.svg
   :target: https://anaconda.org/gimli/pygimli
