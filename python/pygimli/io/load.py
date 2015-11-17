# -*- coding: utf-8 -*-

import os.path

import pygimli as pg
from pygimli.meshtools import readGmsh, readPLC


def load(fname, verbose=False):
    """
    General import function to load data and meshes from file.

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
        ".idx": pg.IVector,
        # Matrices
        ".bmat": pg.RMatrix,
        ".mat": pg.RMatrix,
        # Meshes
        ".poly": readPLC,
        ".bms": pg.Mesh,
        ".msh": readGmsh,
        ".vtk": pg.Mesh
    }

    if not os.path.exists(fname):
        raise Exception("File or directory named %s does not exist." % fname)

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
        except Exception as e:
            if verbose:
                import sys
                import traceback
                traceback.print_exc(file=sys.stdout)
                print(e)
                print("File extension %s seems to be not correct. "
                    "Trying auto-detect." % suffix)
    else:
        if verbose:
            print(
                "File extension %s is not known. Trying auto-detect." %
                suffix)

    for routine in import_routines.values():
        try:
            return routine(fname)
        except Exception as e:
            # print(e)
            pass

    raise Exception(
        "File type of %s is unknown or file does not exist and could not be imported." %
        fname)
