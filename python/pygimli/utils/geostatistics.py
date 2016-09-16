import time
from math import pi, sin, cos

import numpy as np
import scipy.linalg

import pygimli as pg


def covarianceMatrixVec(x, y, z=None, I=None, dip=0, strike=0, var=1):
    """Geostatistical covariance matrix.

    Parameters
    ----------

    """
    if I is None:
        I = [1, 1, 1]
    elif isinstance(I, (float, int)):
        I = [I, I, I]
    elif len(I) < 3:
        I.append(I[-1])

    if z is None:
        z = np.zeros_like(x)

    hx = x - x[:, np.newaxis]
    hy = y - y[:, np.newaxis]
    hz = z - z[:, np.newaxis]
    alpha = -dip * pi / 180  # rotation of operator
    beta = -strike * pi / 180
    Hx = (hx*cos(alpha)-hy*sin(alpha))*cos(beta) / I[0]
    Hy = (hx*sin(alpha)+hy*cos(alpha))*cos(beta) / I[1]
    Hz = hz * sin(beta) / I[2]
    CM = var*np.exp(-np.sqrt(Hx**2+Hy**2+Hz**2))  # Covariance matrix

    return CM


def covarianceMatrixPos(pos, **kwargs):
    """Position (R3Vector) based covariance matrix"""
    return covarianceMatrixVec(np.array(pg.x(pos)), np.array(pg.y(pos)),
                               np.array(pg.z(pos)), **kwargs)


def covarianceMatrix(mesh, nodes=False, **kwargs):
    """Geostatistical covariance matrix (cell or node) for given mesh.

    Parameters
    ----------
    mesh : gimliapi:`GIMLI::Mesh`
        Mesh
    nodes : bool [False]
        use node positions, otherwise (default) cell centers are used
    **kwargs

        I : float or list of floats
            correlation lengths (range) in individual directions
        dip : float
            dip angle of major axis (I[0])
        strike : float
            strike angle (for 3D)

    Returns
    -------
    Cm : np.array (square matrix of size cellCount/nodeCount)
        covariance matrix
    """
    if nodes:
        pos = mesh.positions()
    else:
        pos = mesh.cellCenters()
    return covarianceMatrixPos(pos, **kwargs)


def computeInverseRootMatrix(CM, thrsh=0.001, verbose=False):
    """Compute inverse square root (C^{-0.5} of matrix."""
    t = time.time()
    e_vals, e_vecs = scipy.linalg.eigh(CM)
    if verbose:
        print('(C) Calculation time for eigenvalue decomposition:\n%s sec' %
              (time.time()-t))

    t = time.time()
    A = scipy.linalg.inv(np.diag(np.sqrt(np.real(e_vals))))
    if verbose:
        print('(C) Calculation time for inv(diag(sqrt)):\n%s sec' %
              (time.time()-t))

    t = time.time()
    gemm = scipy.linalg.get_blas_funcs("gemm", [e_vecs, A])
    B = gemm(1, e_vecs, A)
    if verbose:
        print('(C) Calculation time for dot 1:\n%s sec' % (time.time()-t))

    t = time.time()
    gemm2 = scipy.linalg.get_blas_funcs("gemm", [B, np.transpose(e_vecs)])
    if verbose:
        print('gemm test:\n%s sec' % (time.time()-t))
    CM05 = gemm2(1, B, np.transpose(e_vecs))
    if verbose:
        print('(C) Calculation time for dot 2:\n%s sec' % (time.time()-t))

    if thrsh:
        nModel = len(CM)
        RCM05 = pg.RSparseMapMatrix(nModel, nModel)
        for i in range(nModel):
            for j in range(nModel):
                if np.abs(CM05[i][j]) > thrsh:
                    RCM05.setVal(i, j, CM05[i][j])
        return RCM05
    else:
        return CM05  # not making sense as constraints matrix
