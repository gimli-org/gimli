Installation
============

.. note::

  On all platforms, we recommend to install pyGIMLi via the conda package manager contained in the Anaconda distribution.

  To avoid conflicts with other packages, we recommend to install pygimli in a
  separate environment. Here we call this environment `pg`, but you can give
  it any name. Note that this environment has to be created only once.

  .. code-block:: bash

      conda create -n pg -c gimli/label/test -c conda-forge pygimli

  If you want to use pygimli, you have to activate the environment. You can put
  this line in your `~/.bashrc` file so that it is activated automatically if you
  open a terminal.

  .. code-block:: bash

      conda activate pg

  To test if everything works correctly you can do the following:

  .. code-block:: bash

      python -c "import pygimli; pygimli.test(show=False, onlydoctests=True)"
