"""
Methods to calculate the model resolution.
"""

import numpy as np

from pygimli.io import opt_import

lin = opt_import("scipy.linalg", "calculate model resolion matrices")
if lin:
    import scipy.sparse as sp
    from scipy.sparse import coo_matrix, csc_matrix
    from scipy.sparse.linalg import aslinearoperator, LinearOperator, lsqr



def computeR(J, C, alpha=0.5):
    r"""Return diagional of model resolution matrix.

    Calculates the formal model resolution matrix deterministically following:

    .. math::

        \mathbf{R_m} = (\mathbf{J}^T\mathbf{D}^T\mathbf{D}\mathbf{J} +
        \lambda \mathbf{C}^T\mathbf{C})^{-1} \mathbf{J}^T\mathbf{D}^T\mathbf{D}\mathbf{J}

    .. note::

        The current implementation assumes that :math:`\mathbf{D}` is the
        identitiy matrix, i.e. equal data weights.

    Parameters
    ----------
    J : array
        Jacobian matrix.
    C : array
        Constraint matrix.
    alpha : float
        Regularization strength :math:`\alpha`.
    """
    JTJ = J.T.dot(J)
    CM_inv = C.T.dot(C)
    Jsharp = lin.solve(JTJ + alpha * CM_inv, J.T)
    R = np.diag(Jsharp.dot(J))
    return R
