#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""2.5D non-optimized totalfield forward operator for ERT.

Please use the BERT package for more advanced forward operator
https://gitlab.com/resistivity-net/bert
"""

import numpy as np

import pygimli as pg


def simulate(mesh, scheme, res, sr=True, useBert=True,
             verbose=False, **kwargs):
    """ERT forward calculation.

    Convenience function to use the ERT modelling operator
    if you like static functions.

    See :py:mod:`pygimli.ert.ERTManager.simulate` for description
    of the arguments.

    Parameters
    ----------
    mesh: :gimliapi:`GIMLI::Mesh` | str
        Modelling domain. Mesh can be a file name here.
    scheme: :gimliapi:`GIMLI::DataContainerERT` | str
        Data configuration. Scheme can be a file name here.
    res: see :py:mod:`pygimli.ert.ERTManager.simulate`
        Resistivity distribution.
    sr: bool [True]
        Use singularity removal technique.
    useBert: bool [True]
        Use Bert forward operator instead of the reference implementation.
    **kwargs:
        Forwarded to :py:mod:`pygimli.ert.ERTManager.simulate`
    """
    from .ertManager import ERTManager
    ert = ERTManager(useBert=useBert, sr=sr, verbose=verbose)

    if isinstance(mesh, str):
        mesh = pg.load(mesh)

    if isinstance(scheme, str):
        scheme = pg.physics.ert.load(scheme)

    return ert.simulate(mesh=mesh, res=res, scheme=scheme,
                        verbose=verbose, **kwargs)


@pg.cache
def createGeometricFactors(scheme, numerical=None, mesh=None,
                           h2=True, p2=True, verbose=False):
    """Create geometric factors for a given data scheme.

    Create geometric factors for a data scheme with and without topography.
    Calculation will be done analytical (only for half space geometry)
    or numerical.

    This function caches the result depending on scheme, mesh and pg.version()

    Parameters
    ----------
    scheme: :gimliapi:`GIMLI::DataContainerERT`
        Datacontainer of the scheme.
    numerical: bool | None [False]
        If numerical is None, False is assumed, we try to guess topography
        and warn if we think we found them.
        If set to True or False, numerical calculation will used respectively.
    mesh: :gimliapi:`GIMLI::Mesh` | str
        Mesh for numerical calculation. If not given, analytical geometric
        factors for halfspace earth are guessed or a default mesh will be
        created. The mesh will be h and p refined. If given topo is set to
        True. If the numerical effort is to high or the accuracy to low
        you should consider to calculate the factors manual.
    h2: bool [True]
        Default refinement to achieve high accuracy calculation
    p2: bool [True]
        Default refinement to achieve high accuracy calculation
    verbose: bool
        Give some output.
    """
    if numerical is None:
        numerical = False
        if min(pg.z(scheme)) != max(pg.z(scheme)):
            verbose = True
            pg.warn('Sensor z-coordinates not equal. Is there topography?')

    if numerical is False and mesh is None:
        if verbose:
            pg.info('Calculate analytical flat earth geometric factors.')

        return pg.core.geometricFactors(scheme, forceFlatEarth=True)

    if mesh is None:
        mesh = createInversionMesh(scheme)

    if verbose:
        pg.info('mesh', mesh)

    if h2 is True:
        m = mesh.createH2()
        if verbose:
            pg.info('h2 refine', m)

    if p2 is True:
        m = m.createP2()
        if verbose:
            pg.info('p2 refine', m)

    if verbose:
        pg.info('Calculate numerical geometric factors.')

    d = simulate(m, res=1.0, scheme=scheme, sr=False, useBert=True,
                 calcOnly=True, verbose=True)
    return 1./d['u']


def createInversionMesh(data, **kwargs):
    """Create default mesh for ERT inversion.

    Parameters
    ----------
    data: :gimliapi:`GIMLI::DataContainerERT`
        Data Container needs at least sensors to define the geometry of the
        mesh.

    Other Parameters
    ----------------
        Forwarded to :py:mod:`pygimli.meshtools.createParaMesh`

    Returns
    -------
    mesh: :gimliapi:`GIMLI::Mesh`
        Inversion mesh with default marker (1 for background,
        2 parametric domain)
    """
    mesh = pg.meshtools.createParaMesh(data.sensors(), **kwargs)
    return mesh


def createERTData(elecs, schemeName='none', **kwargs):
    """Create data scheme for compatibility (advanced version in BERT).

    Parameters
    ----------
    sounding : bool [False]
        Create a 1D VES Schlumberger configuration.
        elecs need to be an array with elecs[0] = mn/2 and elecs[1:] = ab/2.
    """
    if kwargs.pop('sounding', False):
        data = pg.DataContainerERT()
        data.setSensors(pg.cat(-elecs[::-1], elecs))

        nElecs = len(elecs)
        for i in range(nElecs-1):
            data.createFourPointData(i, i, 2*nElecs-i-1, nElecs-1, nElecs)

        return data

    if schemeName != "dd":
        import pybert as pb  # that's bad!!! TODO: remove pybert deps
        return pb.createData(elecs, schemeName, **kwargs)

    isClosed = kwargs.pop('closed', False)

    data = pg.DataContainerERT()
    data.setSensors(elecs)

    nElecs = len(elecs)
    a = []
    b = []
    m = []
    n = []
    eb = 0
    for i in range(nElecs):
        for j in range(eb + 2, nElecs):
            ea = i
            eb = ea + 1
            em = j
            en = em + 1

            if isClosed:
                en = en % nElecs

            if en < nElecs and en != ea:
                a.append(ea)
                b.append(eb)
                m.append(em)
                n.append(en)

    data.resize(len(a))
    data.add('a', a)
    data.add('b', b)
    data.add('m', m)
    data.add('n', n)
    data.set('valid', np.ones(len(a)))

    return data


def estimateError(data, absoluteError=0.001, relativeError=0.03,
                  absoluteUError=None, absoluteCurrent=0.1):
    """Estimate error composed of an absolute and a relative part.

    Parameters
    ----------
    absoluteError : float [0.001]
        Absolute data error in Ohm m. Need 'rhoa' values in data.

    relativeError : float [0.03]
        relative error level in %/100

    absoluteUError : float [0.001]
        Absolute potential error in V. Need 'u' values in data. Or
        calculate them from 'rhoa', 'k' and absoluteCurrent if no 'i'
        is given

    absoluteCurrent : float [0.1]
        Current level in A for reconstruction for absolute potential V

    Returns
    -------
    error : Array
    """
    if relativeError >= 0.5:
        print("relativeError set to a value > 0.5 .. assuming this "
              "is a percentage Error level dividing them by 100")
        relativeError /= 100.0

    if absoluteUError is None:
        if not data.allNonZero('rhoa'):
            pg.critical("We need apparent resistivity values "
                        "(rhoa) in the data to estimate a "
                        "data error.")
        error = relativeError + absoluteError / data('rhoa')
    else:
        u = None
        i = absoluteCurrent
        if data.haveData("i"):
            i = data('i')

        if data.haveData("u"):
            u = data('u')
        else:
            if data.haveData("r"):
                u = data('r') * i
            elif data.haveData("rhoa"):
                if data.haveData("k"):
                    u = data('rhoa') / data('k') * i
                else:
                    pg.critical("We need (rhoa) and (k) in the"
                                "data to estimate data error.")

            else:
                pg.critical("We need apparent resistivity values "
                            "(rhoa) or impedances (r) "
                            "in the data to estimate data error.")

        error = pg.abs(absoluteUError / u) + relativeError

    return error


if __name__ == "__main__":
    pass
