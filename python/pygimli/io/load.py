# -*- coding: utf-8 -*-
"""TODO Module docstring."""

import os.path
from importlib import import_module

import pygimli as pg
from pygimli.meshtools import readFenicsHDF5Mesh, readGmsh, readPLC


def optImport(module, requiredFor="use the full functionality"):
    """Import and return module only if it exists.

    If `module` cannot be imported, a warning is printed followed by the
    `requiredFor` string. Otherwise, the imported `module` will be returned.
    This function should be used to import optional dependencies in order to
    avoid repeated try/except statements.

    Parameters
    ----------
    module : str
        Name of the module to be imported.
    requiredFor : str, optional
        Info string for the purpose of the dependency.

    Examples
    --------
    >>> from pygimli.io import optImport
    >>> pg = optImport("pygimli")
    >>> pg.__name__
    'pygimli'
    >>> optImport("doesNotExist", requiredFor="do something special")
    No module named 'doesNotExist'.
    You need to install this optional dependency to do something special.
    """
    # set default message for common imports
    if not requiredFor and "matplotlib" in module:
        requiredFor = "visualize 2D content"

    if module.count(".") > 2:
        raise ImportError("Can only import modules and sub-packages.")

    try:
        mod = import_module(module)
    except ImportError:
        msg = ("No module named \'%s\'.\nYou need to install this optional "
               "dependency to %s.")
        print(msg % (module, requiredFor))
        mod = None

    return mod


opt_import = optImport


def load(fname, verbose=False):
    """General import function to load data and meshes from file.

    Parameters
    ----------
    fname : string
        Filename or folder of files to load.
    verbose : bool
        Be verbose.

    Examples
    --------
    >>> import os, tempfile
    >>> import pygimli as pg
    >>> from pygimli.io import load
    >>> fname = tempfile.mktemp(suffix=".bms")
    >>> pg.createGrid(range(3), range(3)).save(fname)
    1
    >>> mesh = load(fname)
    >>> os.remove(fname)
    >>> mesh.cellCount()
    4
    """
    import_routines = {
        # Data
        ".data": pg.DataContainer,
        ".ohm": pg.DataContainer,  # BERT compatibility
        ".shm": pg.DataContainer,  # BERT compatibility
        ".sgt": pg.DataContainer,
        # Vectors
        ".dat": pg.RVector,
        ".vector": pg.RVector,
        ".vec": pg.RVector,
        ".idx": pg.IVector,
        # Matrices
        ".bmat": pg.RMatrix,
        ".mat": pg.RMatrix,
        # Meshes
        ".poly": readPLC,
        ".bms": pg.Mesh,
        ".msh": readGmsh,
        ".vtk": pg.Mesh,
        ".h5": readFenicsHDF5Mesh  # fenics specs as default
    }

    if not os.path.exists(fname):
        raise Exception("File or directory named %s does not exist." % (fname))

    # recursive function call if fname is a folder of files
    if os.path.isdir(fname):
        files = os.listdir(fname)
        if verbose:
            print("Reading %s with %d files..." % (fname, len(files)))
        return [load(f) for f in files]

    suffix = os.path.splitext(fname)[1]

    if suffix in import_routines:
        try:
            return import_routines[suffix](fname)
        except BaseException as e:
            if verbose:
                import sys
                import traceback
                traceback.print_exc(file=sys.stdout)
                print(e)
                print("File extension %s seems to be not correct. "
                      "Trying auto-detect." % suffix)
    else:
        if verbose:
            print("File extension %s is unknown. Trying auto-detect." % suffix)

    for routine in import_routines.values():
        try:
            return routine(fname)
        except BaseException as _:
            # print(e)
            pass

    raise BaseException("File type of %s is unknown or file does not exist "
                        "and could not be imported." % fname)
