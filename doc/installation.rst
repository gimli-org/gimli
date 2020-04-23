:tocdepth: 2

.. _sec:install:

Installation
============

.. note::

  On all platforms, we recommend to install pyGIMLi via the conda package manager contained in the Anaconda distribution.

  To avoid conflicts with other packages, we recommend to install pygimli in a
  separate environment. Here we call this environment `mypgenv`, but you can give
  it any name. Note that this environment has to be created only once.

  .. code-block:: bash

      conda create -n mypgenv -c gimli/label/test -c conda-forge pygimli

  If you want to use pygimli, you have to activate the environment. You can put
  this line in your `~/.bashrc` file so that it is activated automatically if you
  open a terminal.

  .. code-block:: bash

      conda activate mypgenv

  To test if everything works correctly you can do the following:

  .. code-block:: bash

      python -c "import pygimli; pygimli.test(show=False, onlydoctests=True)"


.. raw:: html

    Choose your operating system for more information, in particular on
    alternative ways of installation or compilation from source:<br><br>
    <ul class="nav nav-tabs">
        <li class="nav-item">
            <a class="nav-link active" href="#install_win" role="tab" data-toggle="tab">
                <i class="fab fa-windows"></i> Windows</a>
        </li>
        <li class="nav-item">
            <a class="nav-link" href="#install_lin" role="tab" data-toggle="tab">
                <i class="fab fa-linux"></i> Linux</a>
        </li>
        <li class="nav-item">
            <a class="nav-link" href="#install_mac" role="tab" data-toggle="tab">
                <i class="fab fa-apple"></i> Mac OS</a>
        </li>
        <li class="nav-item">
            <a class="nav-link" href="#install_pre" role="tab" data-toggle="tab">
                <i class="fas fa-info"></i> Prerequisites</a>
        </li>
    </ul>
    <div class="tab-content">

.. raw:: html

    <div role="tabpanel" class="tab-pane active" id="install_win">

.. include:: ../INSTALL_WIN.rst

.. raw:: html

    </div>
    </div>
    </div>
    <div role="tabpanel" class="tab-pane" id="install_lin">
    <div>
    <div>

.. include:: ../INSTALL_LINUX.rst

.. raw:: html

    </div>
    </div>
    </div>
    <div role="tabpanel" class="tab-pane" id="install_mac">
    <div>
    <div>

.. include:: ../INSTALL_MAC.rst

.. raw:: html

    </div>
    </div>
    <div role="tabpanel" class="tab-pane" id="install_pre">
    <div>

.. include:: ../dependencies.rst

.. raw:: html

    </div>
    </div>
