'''
Locality Preserving Projection (LPP)
Created on 2012/11/27
@author: du
'''

import numpy as np
from scipy.sparse import csr_matrix, lil_matrix
from ..base import BaseEstimator, TransformerMixin
from ..neighbors import kneighbors_graph
from ..utils.extmath import safe_sparse_dot
from ..utils import array2d, check_random_state, check_arrays
from scipy.linalg import eigh
from ..utils.arpack import eigsh

from .lpp_util import toSymmetrixMat

def adjacency_matrix(X, n_neighbors, mode='distance', use_ext=True):
    X = np.asanyarray(X)
    G = kneighbors_graph(X, n_neighbors, mode)
    G = G.tolil()
    nzx, nzy = G.nonzero()
    if use_ext:
        G = toSymmetrixMat(G, nzx, nzy)
    else:
        for xx, yy in zip(nzx, nzy):
            if G[yy, xx] == 0:
                G[yy, xx] = G[xx, yy]
    return G

def affinity_matrix(adj_mat, kernel_func="heat", kernel_param=1.0):
    W = csr_matrix(adj_mat)
    if kernel_func == "heat":
        np.exp(-W.data ** 2 / kernel_param, W.data)
    else:
        W.data = kernel_func(W.data)
    return W

def laplacian_matrix(afn_mat):
    col_sum = np.asarray(afn_mat.sum(0)).flatten()
    afn_mat = (-afn_mat).tolil()
    afn_mat.setdiag(col_sum)
    return afn_mat.tocsr(), col_sum

def auto_dsygv(M, N, k, k_skip=0, eigen_solver='autoack', tol=1E-6,
                  max_iter=100, random_state=None):
    """
    Find the null space of a matrix M.

    Parameters
    ----------
    M : {array, matrix, sparse matrix, LinearOperator}
        Input covariance matrix: should be symmetric positive semi-definite

    k : integer
        Number of eigenvalues/vectors to return

    k_skip : integer, optional
        Number of low eigenvalues to skip.

    eigen_solver : string, {'auto', 'arpack', 'dense'}
        auto : algorithm will attempt to choose the best method for input data
        arpack : use arnoldi iteration in shift-invert mode.
                    For this method, M may be a dense matrix, sparse matrix,
                    or general linear operator.
                    Warning: ARPACK can be unstable for some problems.  It is
                    best to try several random seeds in order to check results.
        dense  : use standard dense matrix operations for the eigenvalue
                    decomposition.  For this method, M must be an array
                    or matrix type.  This method should be avoided for
                    large problems.

    tol : float, optional
        Tolerance for 'arpack' method.
        Not used if eigen_solver=='dense'.

    max_iter : maximum number of iterations for 'arpack' method
        not used if eigen_solver=='dense'

    random_state: numpy.RandomState or int, optional
        The generator or seed used to determine the starting vector for arpack
        iterations.  Defaults to numpy.random.

    """
    if eigen_solver == 'auto':
        if M.shape[0] > 200 and k + k_skip < 10:
            eigen_solver = 'arpack'
        else:
            eigen_solver = 'dense'

    if eigen_solver == 'arpack':
        random_state = check_random_state(random_state)
        v0 = random_state.rand(M.shape[0])
        try:
            eigen_values, eigen_vectors = eigsh(M, k + k_skip, N, sigma=0.0,
                                                tol=tol, maxiter=max_iter,
                                                v0=v0)
        except RuntimeError as msg:
            raise ValueError("Error in determining null-space with ARPACK. "
                             "Error message: '%s'. "
                             "Note that method='arpack' can fail when the "
                             "weight matrix is singular or otherwise "
                             "ill-behaved.  method='dense' is recommended. "
                             "See online documentation for more information."
                             % msg)

        return eigen_vectors[:, k_skip:], eigen_values[k_skip:]
    elif eigen_solver == 'dense':
        if hasattr(M, 'toarray'):
            M = M.toarray()
        if hasattr(N, 'toarray'):
            N = N.toarray()
        eigen_values, eigen_vectors = eigh(
            M, N, eigvals=(0, k + k_skip - 1), overwrite_a=True)
#        index = np.argsort(np.abs(eigen_values))
#        return eigen_vectors[:, index], eigen_values
        return eigen_vectors, eigen_values
    else:
        raise ValueError("Unrecognized eigen_solver '%s'" % eigen_solver)

def lpp(X, n_neighbors, mode="distance", kernel_func="heat", kernel_param=1.0,
        k=2, eigen_solver='arpack', tol=1E-6, max_iter=100,
               random_state=None):
    print "making adjacency matrix"
    W = adjacency_matrix(X, n_neighbors, mode)
    print "making affinity matrix"
    W = affinity_matrix(W, kernel_func, kernel_param)
    print "making laplacian matrix"
    L, D = laplacian_matrix(W)
    print "eigen map"
    L = safe_sparse_dot(X.T, safe_sparse_dot(L, X))
    D = np.dot(X.T, D[:, np.newaxis] * X)
    return auto_dsygv(L, D, k, 0, eigen_solver, tol, max_iter, random_state)

class LPP(BaseEstimator, TransformerMixin):
    def __init__(self, n_neighbors=None, n_components=2, mode="distance",
                 kernel_func="heat", kernel_param=1.0, eigen_solver='arpack',
                 tol=1E-6, max_iter=100, random_state=None):
        self._n_neighbors = n_neighbors
        self._n_components = n_components
        self._mode = mode
        self._kernel_func = kernel_func
        self._kernel_param = kernel_param
        self._eigen_solver = eigen_solver
        self._tol = tol
        self._max_iter = max_iter
        self._random_state = random_state
        self._components = None

    def fit(self, X, y=None):
        if self._kernel_func == "heat" and self._kernel_param is None:
            self._kernel_param = 1.0 / X.shape[1]
        if self._n_neighbors is None:
            self._n_neighbors = max(int(X.shape[0] / 10), 1)
        self._components, _ = lpp(X, self._n_neighbors, self._mode,
                                  self._kernel_func, self._kernel_param,
                                  self._n_components, self._eigen_solver,
                                  self._tol, self._max_iter,
                                  self._random_state)
        return self

    def transform(self, X):
        return safe_sparse_dot(X, self._components)


def lem(X, n_neighbors, mode="distance", kernel_func="heat", kernel_param=1.0,
        k=2, eigen_solver='arpack', tol=1E-6, max_iter=100,
               random_state=None):
    W = adjacency_matrix(X, n_neighbors, mode)
    W = affinity_matrix(W, kernel_func, kernel_param)
    L, D = laplacian_matrix(W)
    D_temp = lil_matrix(L.shape)
    D_temp.setdiag(D)
    D = D_temp.tocsr()
    
    return auto_dsygv(L, D, k, 1, eigen_solver, tol, max_iter, random_state)

class LEM(BaseEstimator, TransformerMixin):
    def __init__(self, n_neighbors=None, n_components=2, mode="distance",
                 kernel_func="heat", kernel_param=1.0, eigen_solver='arpack',
                 tol=1E-6, max_iter=100, random_state=None):
        self._n_neighbors = n_neighbors
        self._n_components = n_components
        self._mode = mode
        self._kernel_func = kernel_func
        self._kernel_param = kernel_param
        self._eigen_solver = eigen_solver
        self._tol = tol
        self._max_iter = max_iter
        self._random_state = random_state
        self._components = None

    def fit(self, X, y=None):
        if self._kernel_func == "heat" and self._kernel_param is None:
            self._kernel_param = 1.0 / X.shape[1]
        if self._n_neighbors is None:
            self._n_neighbors = max(int(X.shape[0] / 10), 1)
        self.embedded, _ = lem(X, self._n_neighbors, self._mode,
                                  self._kernel_func, self._kernel_param,
                                  self._n_components, self._eigen_solver,
                                  self._tol, self._max_iter,
                                  self._random_state)
        return self

    def transform(self, X):
        return self.embedded

def test(n, use_ext):
    from ..datasets import make_s_curve, make_swiss_roll
    X, c = make_s_curve(2000)
    adjacency_matrix(X, n, use_ext=use_ext)
