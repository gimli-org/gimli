# -*- coding: utf-8 -*-
"""TODO Module docstring."""

import os.path
from importlib import import_module
from urllib.request import urlretrieve

import numpy as np
import pygimli as pg
from pygimli.meshtools import readFenicsHDF5Mesh, readGmsh, readPLC, readSTL
from pygimli.utils import readGPX
from pygimli.utils import cache


# Example data repository
exampleDataRepository = ''.join((
    'https://raw.githubusercontent.com/',  # RAW files
    'gimli-org/example-data/',  # Organization and repository
    'master/'  # Branch
))


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
    >>> import pygimli as pg
    >>> pg = pg.optImport("pygimli")
    >>> pg.__name__
    'pygimli'
    >>> pg.optImport("doesNotExist", requiredFor="do something special")
    No module named 'doesNotExist'.
    You need to install this optional dependency to do something special.
    """
    # set default message for common imports
    if not requiredFor and "matplotlib" in module:
        requiredFor = "visualize 2D content"

    if module.count(".") > 2:
        raise ImportError("Can only import modules and sub-packages.")

    mod = None
    try:
        mod = import_module(module)
    except ImportError:
        msg = ("No module named \'%s\'.\nYou need to install this optional "
               "dependency to %s.")
        print(msg % (module, requiredFor))
        # exception would be better her but then the test fails
        # raise Exception(msg % (module, requiredFor))

    return mod


def opt_import(*args, **kwargs):
    pg.deprecated() # last vis: 20190903
    return optImport(*args, **kwargs)


def getConfigPath():
    r"""Return the user space path for configuration or cache files.

    For Windows: return 'C:\Documents and Settings\username\Appdata\gimli'

    For Linux return /home/username/.config/gimli

    """

    if 'APPDATA' in os.environ:
        return os.path.join(os.environ['APPDATA'], 'gimli')
    else:
        return os.path.join(os.environ['HOME'], '.config', 'gimli')

def load(fname, verbose=False, testAll=True, realName=None):
    """General import function to load data and meshes from file.

    Parameters
    ----------
    fname : string
        Filename or folder of files to load.
    testAll : bool [True]
        Test all filter when file suffix is unknown or loading fails.
    realName : str [None]
        Real file name.
        When fname is generic (i.e. without suffix) we can test the
        realName for type info.
    verbose : bool
        Be verbose.

    Examples
    --------
    >>> import os, tempfile
    >>> import pygimli as pg
    >>> fname = tempfile.mktemp(suffix=".bms")
    >>> pg.createGrid(range(3), range(3)).save(fname)
    1
    >>> mesh = pg.load(fname)
    >>> os.remove(fname)
    >>> mesh.cellCount()
    4
    """
    ImportFilter = {
        # Data
        ".dat": pg.DataContainerERT,
        ".data": pg.DataContainerERT,
        ".ohm": pg.DataContainerERT,  # BERT compatibility
        ".shm": pg.DataContainerERT,  # BERT compatibility
        ".sgt": pg.DataContainer,
        ".collect": pg.core.DataMap,
        # Vectors
        #".dat": pg.Vector,
        ".vector": pg.Vector,
        ".vec": pg.Vector,
        ".idx": pg.IVector,
        # Matrices
        ".bmat": pg.Matrix,
        ".mat": pg.Matrix,
        ".matrix": pg.matrix.SparseMapMatrix,
        # Meshes
        ".poly": readPLC,
        ".bms": pg.Mesh,
        ".msh": readGmsh,
        ".mod": pg.Mesh,
        ".vtk": pg.Mesh,
        ".stl": readSTL,
        ".h5": readFenicsHDF5Mesh,  # fenics specs as default
        # Misc
        ".gpx": readGPX,  # read gpx waypoints
        ".xy": np.loadtxt,  #
    }

    if not os.path.exists(fname):
        raise Exception("File or directory named %s does not exist." % (fname))

    # recursive function call if fname is a folder of files
    if os.path.isdir(fname):
        files = os.listdir(fname)
        if verbose:
            print("Reading %s with %d files..." % (fname, len(files)))
        return [load(f) for f in files]

    suffix = None
    if realName is not None:
        suffix = os.path.splitext(realName)[1]
    else:
        suffix = os.path.splitext(fname)[1]

    if suffix in ImportFilter:
        try:
            if verbose:
                print("Import {0} ({1})".format(fname, ImportFilter[suffix]))
            return ImportFilter[suffix](fname)
        except Exception as e:
            if verbose or pg.core.debug():
                import sys
                import traceback
                traceback.print_exc(file=sys.stdout)
                print(e)
                print("File extension %s seems to be not correct. "
                      "Trying auto-detect." % suffix)
    else:
        if verbose:
            print("File extension {0} is unknown. Trying auto-detect.".format(suffix))

    if testAll:
        for routine in ImportFilter.values():
            try:
                return routine(fname)
            except Exception as _:
                # print(e)
                pass

    raise Exception("File type of {0} is unknown or file does not exist "
                        "and could not be imported.".format(suffix))


def getExampleFile(path, cacheName=True):
    """Download and return temporary filename of file in example repository."""
    # TODO: Cache locally and check hash sums for potential file corruption
    url = exampleDataRepository + path
    fileName = path.split('/')[-1]
    if not os.path.exists(fileName):
        pg.info("Getting:", fileName)
        return urlretrieve(url)[0]
    else:
        pg.info("File allready exists:", fileName)
    return fileName

# maybe better move this in getExampleFile
# @cache
# def getExampleData(path, type=None):
#     """Get data from example repository."""
#     url = getExampleFile(path)
#     dat = load(url, realName=path)
#     return dat
