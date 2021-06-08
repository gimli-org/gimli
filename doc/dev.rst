.. _sec:coding_guidelines:

GitHub notes for pyGIMLi developers
===================================

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


Branches
--------

+--------------+----------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------------------------------------------------------------+
| Branch name  | Description                | Purpose                                                                                                                                                                                                                                | CI (Jenkins) Current behavior                                                                                                                                                               | CI (Jenkins) Wishlist                                                                                                         | Rules                                                                                                                  |
+==============+============================+========================================================================================================================================================================================================================================+=============================================================================================================================================================================================+===============================================================================================================================+========================================================================================================================+
| dev          | Main development branch    | This is where new versions are developed and pull requests are merged to.                                                                                                                                                              | Runs build, documentation and tests after each push to GitHub. Update http://dev.pygimli.org                                                                                                | Test pull requests, build nightly conda packages for test channel                                                             | Run pg.test() and make doc before pushing. Tag related issues in commit message.                                       |
+--------------+----------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------------------------------------------------------------+
| release      | Latest release             | This is where releases are "staged". This usually means creating a git tag and manually merging dev into release. Hot fixes and website typos can be directly commited here.                                                           | Starts in empty workspace, runs build, documentation and tests after each push to GitHub. If the build is successful, release is merged into master and http://www.pygimli.org is updated.  | Test "backward compatibility" (e.g., run example scripts from last release with this version). Test on Mac and Windows, too.  | Make sure the tag is annotated and the version string is following the principles described below.                     |
+--------------+----------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------------------------------------------------------------+
| master       | Latest and tested release  | Make sure that if people checkout the repository, they always have a working version.                                                                                                                                                  |                                                                                                                                                                                             | Build pgcore (if necessary) and pygimli conda packages for release.                                                           | Never push anything to master!                                                                                         |
+--------------+----------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------------------------------------------------------------+
| *            | Feature branches           | Larger endeavors and major restructuring should happen in dedicated feature branches (or forks), which are eventually merged to dev. This can also be useful if you want to give write access to others to jointly work on a feature.  | Automatic testing can be requested (florian@pygimli.org).                                                                                                                                   |                                                                                                                               | Start feature branch from dev. Inform other developers about your develpment (to avoid conflicts and redundant work).  |
+--------------+----------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------------------------------------------------------------+

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
https://www.sphinx-doc.org/en/master/usage/extensions/napoleon.html

.. code-block:: python

  def foo(arg1, arg2):

  """Short description, i.e., one line to explain what foo does.

    Explain a little bit more verbose was foo does. Use references :cite:`Archie1942`

    Use links to pygimli api :py:mod:pygimli.manager` for modules
    :py:func:pygimli.solver.solveFiniteElements for functions


    Use math.

    .. math :
        a + \simga * \rho

    Explain all parameters.

    Args
    ----
    arg1: type | use links to :gimliapi:`GIMLI::Mesh`
        Describe arg1.
    arg2: type
        Describe arg2.

    Keyword Args
    ------------
    args: type
        Description.

    Attributes
    ----------
        For members

    Returns
    -------
    type:

    Examples
    --------
    >>> import foo
    >>>
    >>>

    See Also
    --------
        average : Weighted average,
        e.g., Link to tutorials :ref:`tut:Modelling_BC` assuming there
        has been set a appropriate label in the tutorial.

    References
    ----------

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

Adding an example to the paper carousel
.......................................

1. Put a catchy figure in *doc/_static/usecases*. It should be 750X380 px. You can use this command:

.. code-block:: bash

  mkdir -p converted
  mogrify -path converted -format jpg -resize "750x380^" -gravity North -crop 750x380+0+0 +repage -quality 95 *.jpg

2. Add a new dictionary with additional information into *doc/paper_carousel.py*. The format should be like this:

.. code-block:: python

  dict(img="jordi2018.jpg",
       title="Improved geophysical images through geostatistical regularization",
       subtitle="Jordi et al. (2018), Geophysical Journal International",
       link="https://doi.org/10.1093/gji/ggy055")


Contributing to the Documentation
---------------------------------
.. note::

    sphinx is a required to add and build the documentation. To learn more about sphinx, 
    visit the following source: https://www.sphinx-doc.org/en/master/index.html

PyGIMLi's documentation is built on sphinx. If you would like to contribute to the documentation of functions that you have been using and/or 
have added. Please follow these steps: 

1. Clone the latest pygimli and checkout the dev branch as shown in "Contributing to Code" and navigate to the *doc* directory.

2. If you do not have pgcore installed, create an environment called *pgdev* that installs pgcore and activate it:

.. code-block:: bash

    conda create -n pgdev -c gimli -c conda-forge pgcore
    conda activate pgdev
    
3. In order to also install PyGIMLi requirements on the new "pydev" environment, execute the following line of code in the *doc* directory:

.. code-block:: bash

    while read requirement; do conda install --yes $requirement; done < requirements.txt

4. Please verify if you have the following sphinx components installed (these can be installed via conda-forge):

    * sphinxcontrib-programoutput 
    * sphinxcontrib-bibtex < 2.0.0 
    * sphinxcontrib-doxylink 
    * bibtexparse
    * sphinx-gallery

5. Make changes to documentation following the guidelines of API documentation. 

6. Build the final documentation 

.. code-block:: bash

    make html

7. The process should build all of the documentation so it may take some time to complete. When it finishes, the last line should state:

.. code-block:: bash

    Build finished. The HTML pages are in _build/html.

8. Open a Pull request in the dev branch of pygimli to merge your changes.  


ToDo-list
---------

.. todolist::
