.. _sec:coding_guidelines:

GitHub notes for pyGIMLi developers
===================================

Contributing to the code
------------------------

.. note::

    These steps are only relevant, if you have write privileges to the main
    GitHub repository. The general contribution guidlines can be found in _sec:contributing:.

1. Clone the GIMLi repository:

    .. code-block:: bash

        git clone https://github.com/gimli-org/gimli gimli-src && cd gimli-src

You should be on the master branch by default. Check by calling *git status*.
GIMLi encourages developers to keep the *master branch* clean. So please use
the default development branch called *dev*. Save `this script
<https://gist.github.com/florian-wagner/0be0b09bbabec6791ed3>`_ in .git/hooks
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