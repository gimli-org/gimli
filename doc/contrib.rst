.. _sec:contributing:

Contributing
============

As an open-source project, pyGIMLi always welcomes contributions from the community. Here, we offer guidance for 4 different ways of contributing with
increasing levels of required coding proficiency (A-D).

A. Submit a bug report
----------------------

If you experience issues with pyGIMLi or miss a certain feature, please `open a
new issue on GitHub <https://github.com/gimli-org/gimli/issues>`_. To do so,
you need to `create a GitHub account <https://github.com/join>`_.

B. Send us your example
-----------------------

We are constantly looking for interesting usage examples of pyGIMLi. If you have
used the package and would like to contribute your work to the :ref:`chapt:examples`, in your Python script are documented according to the `sphinx-gallery syntax
<https://sphinx-gallery.github.io/stable/tutorials/plot_parse.html#sphx-glr-tutorials-plot-parse-py>`_.

C. Quickly edit the documentation
---------------------------------

If you would like to edit a typo or improve the documentation of our website, you can submit fixes to the 
documentation page by clicking on "Improve this page" under the quick links menu on every page.
This will direct you to the Github page to make the edits directly on the ".rst" using Github (you will need a Github account)
Once you make your edits, fill out the fields under "Commit changes" (make sure to provide a good overview of the changes you have made). 
Note that you will make changes to a new branch under your fork since you will not have access to make changes in the existing branch. 
Once you click "Commit changes" this automatically creates a pull request which we will review and then merge them if everything is good to go!

If you would like to contribute to the API or make changes offline to the files under doc. Please make sure 
you have `sphinx <https://www.sphinx-doc.org/en/master/index.html>`_ installed and see the following *Contribute to Code* and *Documentation* for more information.


D. Contribute to code
---------------------

.. note::

    To avoid redundant work, please `contact us <mailto:mail@pygimli.org>`_  
    or check the `current issues <https://github.com/gimli-org/gimli/issues>`_ 
    before you start working on a non-trivial feature.

The preferred way to contribute to the pygimli code base is via a pull request
(PR) on GitHub. The general concept is explained `here
<https://guides.github.com/introduction/flow>`_ and involves the following steps:


1. Fork and clone the repository
++++++++++++++++++++++++++++++++

If you are a first-time contributor, you need `a GitHub account
<https://github.com/join>`_ and your own copy ("fork") of the code. To do so, go
to https://github.com/gimli-org/gimli and click the "Fork button" in the upper
right corner. This will create an identical copy of the complete code base under
your username on the GitHub server. You can navigate to this repository and clone it to your local disk:

.. code-block:: bash

  git clone https://github.com/YOUR_USERNAME/gimli

2. Prepare your environment
+++++++++++++++++++++++++++

To install the neccesary Python requirements in your computer, we recommend Anaconda and the conda package manager.
The cloned and forked repository contains an *environment.yml* file with the specification for all package requirements to build 
and test pyGIMLi. This file should exist in the directory where you have cloned the repository. In order to create a separate environment
where you will work. Run the following command on the same directory of the repository in your terminal: 

.. code-block:: bash

    conda env create
    conda activate pgdev
    
You will need to do this everytime you start a terminal. ``(pydev)`` should appear before your terminal location. 

3. Create a feature branch
++++++++++++++++++++++++++

Go to the source folder and create a feature branch to hold your changes. It is
advisable to give it a sensible name such as ``adaptive_meshes``.

.. code-block:: bash

    cd gimli
    git checkout -b adaptive_meshes

4. Start making your changes
++++++++++++++++++++++++++++

Go nuts! Add and modify files and regularly commit your changes with meaningful
commit messages. Remember that you are working in your own personal copy and in
case you break something, you can always go back. While coding, we encourage you
to follow a few :ref:`sec:coding_guidelines`.

.. code:: bash

    git add new_file1 new_file2 modified_file1
    git commit -m "Implemented adaptive meshes."

5. Test your code
+++++++++++++++++

Make sure that everything works as expected. New functions should always contain
a docstring with a test:

.. code:: python

  def sum(a, b):
      """Return the sum of `a` and `b`.

      Examples
      --------
      >>> a = 1
      >>> b = 2
      >>> sum(a,b)
      3
      """
      return a + b

When you run ``pg.test()`` the docstring test will be evaluated. See also the
section on :ref:`sec:testing`.

6. Documentation
++++++++++++++++

We use sphinx to build the web pages from these sources. To edit HTML file, navigate to the doc folder in the repository.
Once you have made edits. Run the following commands on your terminal. 

.. code:: bash

    make html

The process should build all of the documentation so it may take some time to complete. When it finishes, the last line should state:

.. code:: bash

    Build finished. The HTML pages are in _build/html.


7. Submit a pull request
++++++++++++++++++++++++

Once you implemented a functioning new feature, make sure your GitHub repository
contains all your commits:

.. code:: bash

  git push origin adaptive_meshes

After pushing, you can go to GitHub and you will see a green PR button. Describe
your changes in more detail. Once reviewed by the core developers, your PR will
be merged to the main repository.

.. _sec:Coding_rules:

Coding rules
------------

* General we try to use Pep 8 https://www.python.org/dev/peps/pep-0008/?

* All names should be literally and in CamelShape style.

* Classes starts with Upper Case Letter.

* Members and Methods always starts with Lower Case Letter.

* All class member (self.member) need to be initialized in the Constructor.

* (ugly²) Do not use multi initialize in one line, e.g., a, b, c = 0, 0, 0

* check for data types with 'if isinstance(var, type):' instead 'if type(var) == type:'

.. _sec:coding_guidelines:

Coding Guidelines
-----------------

Use pylint or prospector to improve code quality.

We use: (You can find exceptions in .landscape.yml)

* pep8
* pep257
* pylint
* pyflakes

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


