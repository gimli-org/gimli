# -*- coding: utf-8 -*-

import os.path

import pygimli as pg
from pygimli.meshtools import readGmsh

def load(fname, verbose=False):
    """
    General import function to load data and meshes from file.

    Parameters
    ----------
    fname : string
        Filename or folder of files to load.
    verbose : bool
        Be verbose.
    """

    import_routines = {
            # Data
            ".data": pg.DataContainer,
            ".ohm": pg.DataContainer, # BERT compatibility
            # Vectors
            ".dat": pg.RVector,
            ".vector": pg.RVector,
            ".idx": pg.IVector,
            # Matrices
            ".bmat": pg.RMatrix,
            ".mat": pg.RMatrix,
            # Meshes
            ".bms": pg.Mesh,
            ".msh": readGmsh,
            ".poly": pg.Mesh,
            ".vtk": pg.Mesh
            }

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
        except:
            if verbose:
                print("File extension %s seems to be not correct. Trying auto-detect." % suffix)
    else:
        if verbose:
            print("File extension %s is not known. Trying auto-detect." % suffix)

    for routine in import_routines.values():
        try:
            return routine(fname)
        except:
            pass

    raise Exception("File type of %s is unknown and could not be imported." % fname)
