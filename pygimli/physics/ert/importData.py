#!/usr/bin/env python
# -*- coding: utf-8 -*-
import pygimli as pg


def load(fileName, verbose=False, **kwargs):
    """Shortcut to load ERT data.

    Import Data and try to assume the file format.
    Use pybert importer if installed.

    Parameters
    ----------
    fileName: str

    Returns
    -------
    data: pg.DataContainer

    """  
    
    data = pg.load(fileName)
    
    if isinstance(data, pg.DataContainerERT):   
        return data

    if verbose:
        pg.info("Try to import using pybert .. if available")

    pb = pg.optImport('pybert')
    data = pb.loadData(fileName)

    if isinstance(data, pg.DataContainerERT):   
        return data

    pg.critical("Can't import ERT data file.", fileName)

