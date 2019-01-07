.. _sec:coding_guidelines:

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

2. Checkout the dev branch:

    .. code-block:: bash

        git checkout dev

3. Make changes in your favorite editor.

4. Add new files to the staging area:

    .. code-block:: bash

       git add new_file1 new_file2

5. Make a commit with a *meaningful* message:

    .. code-block:: bash

       git commit -a -m "Added important new feature."

6. Pull the latest developments from GitHub using automatic rebase:
  For more info see: http://kernowsoul.com/blog/2012/06/20/4-ways-to-avoid-merge-commits-in-git/
    
    .. code-block:: bash

       git pull --rebase

You can set the rebase behavior on default with:

    .. code-block:: bash

        git config --global branch.autosetuprebase always

7. Push to the origin development branch:

    .. code-block:: bash

       git push origin dev

Note that if you are on the `dev` branch, a `git push` should suffice, but the
command above is more explicit about the branch which should be pushed.

8. Show available tags and branches:

    .. code-block:: bash

       git tag
       git branch

9. If you plan bigger changes or developments create an own branch for it:

    .. code-block:: bash

        git branch my_new_test_branch
        git checkout my_new_test_branch

Work on your code in branch `my_new_test_branch` and test it.
If you are happy with your results merge them back
into the dev branch and delete the branch.

    .. code-block:: bash

        git checkout dev
        git merge my_new_test_branch
        git branch -d my_new_test_branch

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

    git tag -a -m "First official release" "v1.0.x" # tags last commit as v1.0.x
    git push --tags # pushes tags to GitHub

.. _sec:coding_guidelines:

Coding rules
------------

* General we try to use Pep 8 https://www.python.org/dev/peps/pep-0008/?

* All names should be literally and in CamelShape style.

* Classes starts with Upper Case Letter.

* Members and Methods always starts with Lower Case Letter.

* All class member (self.member) need to be initialized in the Constructor.

* (ugly²) Do not use multi initialize in one line, e.g., a, b, c = 0, 0, 0

* check for data types with 'if isinstance(var, type):' instead 'if type(var) == type:'

Behaviour by name for global functions:
.......................................

.. code-block:: python

    createFOO(...)
        """Always needs to return an instance of FOO.""

.. code-block:: python

    showFOO(Bar, ...)
        """Always open a window or optionally use a given Axes to show us Bar as Foo."""
        return ax, cbar

.. code-block:: python

    drawFOO(ax, Bar...)
        """Always need an Axes ax and draws Bar as Foo""
        return graphics_object

.. code-block:: python

    readFOO(fileName, *args):
        """Read object from disc."""
        return obj

.. code-block:: python

    importFOO(fileName, obj, *args):
        """Import object from disc into an existing object."""
        return obj

.. code-block:: python

    exportFOO(obj, fileName):
        """Export object to disc in foreign (FOOs) format."""
        return true

.. code-block:: python

    convertFOO(fooObj):
        """Convert Foo obj into gimli Obj"""
        return gimliObj

API Documentation and doctests
------------------------------

Use the following documentation syntax or see at:
https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt#sections

.. code-block:: python

  def foo(arg1, arg2):

  """Short description, i.e., one line to explain what foo does.

    Explain a little bit more verbose was foo does. Use references :cite:`Archie1942`

    Use links to pygimli api :py:mod:`pygimli.meshtools.appendTriangleBoundary`.

    Use math.

    .. math :
        a + \simga * \rho

    Explain all parameters.

    Parameters
    ----------

    arg1: type(arg1) use links to :gimliapi:`GIMLI::Mesh`
        Describe arg1.

    arg1: type(arg1)
        Describe arg1.

    Write examples that might be complete minimal script that will be executed and tested automatically.

    Attributes
    ----------
        For members

    Examples
    --------
    >>> import foo
    >>>
    >>>

    Yields
    ------
        int : result

    See Also
    --------
        average : Weighted average,
        e.g., Link to tutorials :ref:`tut:Modelling_BC` assuming there
        has been set a appropriate label in the tutorial.

    Raises
    ------
        LinAlgException
            If the matrix is not numerically invertable.


  """

.. _sec:coding_guidelines:

Coding Guidelines
-----------------

Use pylint or prospector to improve code quality.


We use: (You can find exceptions in .landscape.yml)

* pep8
* pep257
* pylint
* pyflakes

.. _sec:testing:

Testing
-------

Run specific API examples from shell:

.. code-block:: bash

  python -c "import pygimli as pg; pg.test(pg.meshtools.createCircle, show=True)"

Run a specific test from shell.

.. code-block:: bash

  python -c "import pygimli; from pygimli.physics.petro.resistivity import *; test_Archie()"

Run all tests

.. code-block:: bash

  python -c "import pygimli; pygimli.test(show=True)"

Run pylint from shell to check code:

.. code-block:: bash

    pylint --rcfile $GIMLIROOT/.pylintrc file.py

Run prospector to check code like landscape.io do:

.. code-block:: bash

    prospector --profile=$GIMLIROOT/.prospector.yml file.py

Read api documentation from shell:

.. code-block:: bash

  python -c "import pygimli as pg; help(pg.test)"

More information on pygimli's native testing function:

.. autofunction:: pygimli.test
