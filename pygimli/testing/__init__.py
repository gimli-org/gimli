# coding=utf-8
"""
Testing utilities

In Python you can call `pygimli.test()` to run all docstring
examples.

Writing tests for pygimli
-------------------------

Please check: https://docs.pytest.org/en/latest/
"""

import sys
from os.path import join, realpath

import matplotlib.pyplot as plt
import numpy as np

import pygimli as pg

import warnings

# __testingMode__ = False
#
# def setTestingMode(mode):
#     """Set pygimli testing mode.
#
#     Testing mode ensures a constant seed for the random generator if you use
#     pg.randn().
#     """
#     globals()[__testingMode__] = mode
#
# def testingMode():
#     """Determine pygimli testing mode.
#
#     Returns True if pygimli is in testing mode.
#     """
#     return globals()[__testingMode__]

def test(target=None, show=False, onlydoctests=False, coverage=False,
         htmlreport=False, abort=False, verbose=True):
    """Run docstring examples and additional tests.

    Examples
    --------
    >>> import pygimli as pg
    >>> # You can test everything with pg.test() or test a single function:
    >>> pg.test("pg.utils.boxprint", verbose=False)
    >>> # The target argument can also be the function directly
    >>> from pygimli.utils import boxprint
    >>> pg.test(boxprint, verbose=False)

    Parameters
    ----------
    target : function or string, optional
        Function or method to test. By default everything is tested.
    show : boolean, optional
        Show matplotlib windows during test run. They will be closed
        automatically.
    onlydoctests : boolean, optional
        Run test files in ../tests as well.
    coverage : boolean, optional
        Create a coverage report. Requires the pytest-cov plugin.
    htmlreport : str, optional
        Filename for HTML report such as www.pygimli.org/build_tests.html.
        Requires pytest-html plugin.
    abort : boolean, optional
        Return correct exit code, e.g. abort documentation build when a test
        fails.
    """
    # pg.setTestingMode(True)
    # Remove figure warnings
    np.random.seed(1337)
    plt.rcParams["figure.max_open_warning"] = 1000
    warnings.filterwarnings("ignore", category=UserWarning,
                            message='Matplotlib is currently using agg, which is a'
                                    ' non-GUI backend, so cannot show the figure.')

    printopt = np.get_printoptions()

    if verbose:
        pg.boxprint("Testing pygimli %s" % pg.__version__, sym="+", width=90)

    # Numpy compatibility (array string representation has changed)
    if np.__version__[:4] == "1.14":
        pg.warn("Some doctests will fail due to old numpy version.",
                "Consider upgrading to numpy >= 1.15")

    if target:
        if isinstance(target, str):
            # If target is a string, such as "pg.solver.solve"
            # the code below will overwrite target with the corresponding
            # imported function, so that doctest works.
            target = target.replace("pg.", "pygimli.")
            import importlib
            mod_name, func_name = target.rsplit('.', 1)
            mod = importlib.import_module(mod_name)
            target = getattr(mod, func_name)

        if show: # Keep figure openend if single function is tested
            plt.ioff()

        import doctest
        doctest.run_docstring_examples(target, globals(), verbose=verbose,
                                       optionflags=doctest.ELLIPSIS,
                                       name=target.__name__)
        return

    try:
        import pytest
    except ImportError:
        raise ImportError("pytest is required to run test suite. "
                          "Try 'sudo pip install pytest'.")

    old_backend = plt.get_backend()
    if not show:
        plt.switch_backend("Agg")
    else:
        plt.ion()

    cwd = join(realpath(__path__[0]), '..')

    excluded = [
        "gui", "physics/traveltime/example.py", "physics/em/fdemexample.py"
    ]

    if onlydoctests:
        excluded.append("testing")

    cmd = ([
        "-v", "-rsxX", "--color", "yes", "--doctest-modules", "--durations",
        "5", cwd
    ])
    for directory in excluded:
        cmd.extend(["--ignore", join(cwd, directory)])

    if coverage:
        pc = pg.optImport("pytest_cov", "create a code coverage report")
        if pc:
            cmd.extend(["--cov", "pygimli"])
            cmd.extend(["--cov-report", "term"])

    if htmlreport:
        ph = pg.optImport("pytest_html", "create a html report")
        if ph:
            cmd.extend(["--html", htmlreport])

    plt.close("all")
    exitcode = pytest.main(cmd)
    if abort:
        print("Exiting with exitcode", exitcode)
        sys.exit(exitcode)

    plt.switch_backend(old_backend)
    np.set_printoptions(**printopt)
