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

        \mathbf{R_m} = (\mathbf{J}^T\mathbf{D}^T\mathbf{D}\mathbf{J} + \alpha
        \mathbf{C}^T\mathbf{C})^{-1}
        \mathbf{J}^T\mathbf{D}^T\mathbf{D}\mathbf{J}

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

#
# # self-made imports
# import pygimli as g
# from bert_tools import *
#
# #general imports
# import numpy as np
# from scipy.io import loadmat
# import scipy.linalg as lin
# import scipy.sparse as sp
# from scipy.sparse import coo_matrix, csc_matrix
# from scipy.sparse.linalg import aslinearoperator, LinearOperator, lsqr
# import matplotlib.pyplot as plt
#
# def load_constraint_matrix(fname):
#     """ Load constraint matrix in sparse format """
#     global mesh
#     i, j, data = np.loadtxt(fname, unpack=True, dtype=int)
#     C = coo_matrix((data, (i,j)), shape=(mesh.boundaryCount(), mesh.cellCount()), dtype=int)
#     return C
#
# def Rmult(x):
#     """ Performs matrix-vector multiplication y = R * x """
#     global J
#     global C
#     global linop
#     r, s = C.shape
#     Jx = J.dot(x)
#     Jxpad = np.hstack((Jx, np.zeros(r)))
#     y = lsqr(linop, Jxpad, atol=1.0e-4, iter_lim=5000, show=False)
#     return y[0]
#
# def diagestrand(dim):
#     """
#         Diagonal estimation based on the algorithm of Bekas et al. (2007)
#     """
#     global Rmult
#
#     chi = 0.01  # stopping criterion (based on convergence tests)
#     d = np.zeros(dim)
#     t = d.copy()
#     q = d.copy()
#     norm = 999
#     k = 0
#
#     while norm > chi:
#         # draw random vector from normal distribution
#         v = np.random.randn(dim)
#         t += Rmult(v) * v
#         q += v * v
#         d_k = t / q
#
#         # evalute difference norm to check convergence
#         norm = np.linalg.norm(d - d_k)
#         d = d_k
#         k +=1
#
#     print "Stopped stochastic estimation after %d iterations with chi = %s." % (k, chi)
#     return d
#
# def compute_res(J, C, alpha=0.5):
#     """
#     Determinisctic computation of model resolution """
#
#     print "Starting deterministic estimation of resolution diagonal..."
#
#     JTJ = J.T.dot(J)
#     CM_inv = C.T.dot(C)
#     Jsharp = lin.solve(JTJ + alpha * CM_inv, J.T)
#     R = np.diag(Jsharp.dot(J))
#     return R
#
# def compare_res(J, C, alpha=0.5):
#     """ Compare determinstic to stochastic resolution diagonal estimate """
#
#     R_det = compute_res(J, C, alpha)
#     R_est = diagestrand(J.shape[1])
#
#     fig = plt.figure(figsize=(12,5))
#     ax1 = fig.add_subplot(121)
#     ax1.plot(R_det, R_est, 'k.', alpha=0.5)
#     ax1.plot((0,1), 'k--')
#     plt.xlabel('Deterministic estimate of resolution diagonal')
#     plt.ylabel('Stochastic estimate of resolution diagonal')
#     plt.axis('scaled')
#     plt.xlim(0, R_det.max()+0.1)
#     plt.ylim(0, R_det.max()+0.1)
#
#     ax2 = fig.add_subplot(122)
#     ax2.hist(R_det - R_est, 20, color='gray')
#     plt.xlabel('Deterministic - stochastic estimate of resolution diagonal')
#     plt.xlim(-0.15, 0.15)
#     plt.ylabel('Count')
#     plt.draw()
#     plt.show()
#
#     return fig
#
# if __name__ == '__main__':
#     alpha = 0.5
#     mesh = g.Mesh('meshParaDomain.bms')
#     J = loadsens('sens.bmat')
#     C = load_constraint_matrix('constraint.matrix')
#
#     p, q = J.shape
#     r, s = C.shape
#     linop = aslinearoperator(sp.vstack((J, C), 'csc'))
#
#     fig = compare_res(J,C)
#
#
#
# #diagR = diagestrand(q, 256)
# #
# #ntimes = 3
#
# #diags=np.zeros((ntimes, q))
#
# #for k in range(ntimes):
# ## use 256 random vectors for each estimate
# #numcols=256;
# #print 'iteration', k
#
# ## estimate the diagonal of r.
# #diagr=diagestrand(q, numcols);
#
# ## store the estimate.
# #diags[k,:] = diagr
#
# ## take median value.
# #diagrest = np.median(diags, 0)
