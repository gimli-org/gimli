# coding=utf-8
"""
Testing utilities
==================

In Python you can call pygimli.test(show=True) to run all docstring
examples.

Writing tests for pygimli
=========================

Please check: https://pytest.org/latest/example/index.html
"""

import os
import sys
from os.path import isdir, join, realpath

import matplotlib.pyplot as plt

from pygimli.io import opt_import


def test(target=None, show=False, onlydoctests=False, coverage=False,
         htmlreport=False, abort=False, verbose=True):
    """Run docstring examples and additional tests.

    Examples
    --------
    >>> import pygimli as pg
    >>> # You can test everything with pg.test() or test an individual function:
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

    old_backend = plt.get_backend()
    if not show:
        plt.switch_backend("Agg")

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

    cwd = join(realpath(__path__[0]), '..')

    excluded = [
        "gui",
        "physics/traveltime/example.py",
        # "physics/em/fdemexample.py"
    ]

    if onlydoctests:
        excluded.append("testing")

    cmd = (["-v", "-rsxX", "--color", "yes", "--doctest-modules", "--durations", 5, cwd])
    for directory in excluded:
        cmd.extend(["--ignore", join(cwd, directory)])

    if coverage:
        pc = opt_import("pytest_cov", "create a code coverage report")
        if pc:
            cmd.extend(["--cov", "pygimli"])
            cmd.extend(["--cov-report", "term"])

    if htmlreport:
        ph = opt_import("pytest_html", "create a html report")
        if ph:
            cmd.extend(["--html", htmlreport])

    exitcode = pytest.main(cmd)
    plt.switch_backend(old_backend)
    plt.close('all')
    if abort:
        sys.exit(exitcode)
