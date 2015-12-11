GitHub notes for core developers
================================

Contributing to the code
------------------------

.. note::

    These steps are only relevant, if you have write privileges to the main
    GitHub repository. A more general contribution guide will follow soon.

1. Clone the GIMLi repository:

    .. code-block:: bash

        git clone https://github.com/gimli-org/gimli gimli-src && cd gimli-src

You should be on the master branch by default. Check by calling *git status*.
GIMLi encourages developers to keep the *master branch* clean. So please use
the default development branch called *dev*. Save `this script
<https://gist.github.com/florian-wagner/0be0b09bbabec6791ed3>` in .git/hooks
for automatic warnings in case you are about to push to the master branch.

2. Checkout the the dev branch:

    .. code-block:: bash

        git checkout dev

3. Make changes in your favorite editor.

4. Add new files to the staging area:

    .. code-block:: bash

       git add new_file1 new_file2

5. Make a commit with a *meaningful* message:

    .. code-block:: bash

       git commit -a -m "Added important new feature."

6. Pull the latest developments from GitHub:

    .. code-block:: bash

       git pull

7. Push to the origin development branch:

    .. code-block:: bash

       git push origin dev

Note that if you are on the `dev` branch, a `git push` should suffice, but the
command above is more explicit about the branch which should be pushed.

Version numbering
-----------------

Following the PEP conventions, pygimli holds a `__version__` string.

.. code-block:: python

    import pygimli as pg
    print(pg.__version__)
    v0.9-0-gf5a6772-modified

This string consists of:

    *Last git tag*-*Number of commits since last tag*-*last commit hash*

The string has a *-modified* suffix if pygimli has uncommitted changes.

To produce a new version, type:

.. code-block:: bash

   git tag -a -m "First official release" "v1.0" # tags last commit as v1.0
   git push --tags # pushes tags to GitHub


Read api documentation from shell:

.. code-block:: bash

    python -c "import pygimli as pg; help(pg.test)

Run api examples from shell:

.. code-block:: bash

    python -c "import pygimli as pg; pg.test(pg.meshtools.createCircle)"



