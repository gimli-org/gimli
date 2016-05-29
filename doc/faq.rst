.. _sec:faq:

Frequently asked questions
==========================

.. note::

    This section is still under construction. If you have a have a question (about GIMLi) that is not covered here
    write a email to mail@pygimli.org.

.. error::

    Can we create a nice Question environment?

Installation
------------

.. admonition:: Question: Cant't use pygimli in spyder?

    Hello
    I'm running Linux/Windows/WhatEver.. with version XXX

    I try to install pygimli using the following procedure :
 
    ::

        ...
        ...
        make pygimli
    
    but when i call the module in spyder or with the command:
    python -c 'import pygimli as pg; print(pg.__version)'
 
    python: can't open file 'import pygimli as pg; print(pg.__version__)':
    [Errno 2] No such file or directory

    What do you think is happening?


If

:: 

    make pygimli 

ended without any errors and you find a message something like this:

::

    [100%] Linking CXX shared module /home/carsten/src/gimli/gimli/python/pygimli/_pygimli_.so

Then the compilation procedure seems to be successful and you still need to set environment variables:

::

    export PYTHONPATH=$PYTHONPATH:$HOME/src/gimli/gimli/python
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/src/gimli/build/lib
    export PATH=$PATH:$HOME/src/gimli/build/bin

Please check if the installation works in general using the any terminal or console:

::

    python -c 'import pygimli as pg; print(pg.__version__)'


.. error::
    
    Carsten: Anyone with spyder can help out with more details on how to set environments there?

