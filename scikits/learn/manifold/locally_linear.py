"""Locally Linear Embedding"""

# Author: Fabian Pedregosa -- <fabian.pedregosa@inria.fr>
# License: BSD, (C) INRIA 2011

import numpy as np
from scipy.linalg import eigh, svd, qr
from scipy.sparse import linalg, eye, csr_matrix
from scipy.sparse.linalg import LinearOperator
from scipy_future import eigsh
from ..base import BaseEstimator
from ..utils import check_random_state
from ..neighbors import kneighbors_graph, BallTree, barycenter_weights

try:
    import pyamg
    pyamg_available = True
except ImportError:
    pyamg_available = False


def null_space(M, k, k_skip=1, eigen_solver='arpack', tol=1E-6, max_iter=100,
              random_state=None):
    """
    Find the null space of a matrix M.

    Parameters
    ----------
    M : array, matrix, sparse matrix, or LinearOperator
        Input covariance matrix: should be symmetric positive semi-definite

    k : number of eigenvalues/vectors to return

    k_skip : number of low eigenvalues to skip.

    eigen_solver : string ['arpack' | 'lobpcg' | 'dense']
        arpack : use arnoldi iteration in shift-invert mode.
                 For this method, M may be a dense matrix, sparse matrix,
                 or general linear operator.
        lobpcg : use locally optimized block-preconditioned conjugate gradient.
                 For this method, M may be a dense or sparse matrix.
                 A dense matrix M will be converted internally to a
                 csr sparse format.
        dense  : use standard dense matrix operations for the eigenvalue
                 decomposition.  For this method, M must be an array or matrix
                 type.  This method should be avoided for large problems.

    tol : tolerance for 'arpack' or 'lobpcg' methods.
          not used if eigen_solver=='dense'

    max_iter : maximum number of iterations for 'arpack' or 'lobpcg' methods
          not used if eigen_solver=='dense'
    """
    random_state = check_random_state(random_state)

    if eigen_solver == 'arpack':
        eigen_values, eigen_vectors = eigsh(M, k + k_skip, sigma=0.0,
                                            tol=tol, maxiter=max_iter)
        return eigen_vectors[:, k_skip:], np.sum(eigen_values[k_skip:])
    elif eigen_solver == 'lobpcg':
        # initial vectors for iteration
        X = np.random.rand(M.shape[0], k + k_skip)
        try:
            ml = pyamg.smoothed_aggregation_solver(M, symmetry='symmetric')
        except TypeError:
            ml = pyamg.smoothed_aggregation_solver(M, mat_flag='symmetric')
        prec = ml.aspreconditioner()

        # compute eigenvalues and eigenvectors with LOBPCG
        eigen_values, eigen_vectors = linalg.lobpcg(
            M, X, M=prec, largest=False, tol=tol, maxiter=max_iter)

        index = np.argsort(eigen_values)
        return (eigen_vectors[:, index[k_skip:]],
                np.sum(eigen_values[index[k_skip:]]))
    elif eigen_solver == 'dense':
        M = np.asarray(M)
        eigen_values, eigen_vectors = eigh(
            M, eigvals=(k_skip, k + k_skip), overwrite_a=True)
        index = np.argsort(np.abs(eigen_values))
        return eigen_vectors[:, index], np.sum(eigen_values)
    else:
        raise ValueError("Unrecognized eigen_solver '%s'" % eigen_solver)


def locally_linear_embedding(
    X, n_neighbors, out_dim, reg=1e-3, eigen_solver='arpack',
    tol=1e-6, max_iter=100, method='standard', H_tol=1E-4, M_tol=1E-12,
    random_state=None):
    """Perform a Locally Linear Embedding analysis on the data.

    Parameters
    ----------
    X : array-like or BallTree, shape [n_samples, n_features]
        Input data, in the form of a numpy array or a precomputed
        :class:`BallTree`.

    n_neighbors : integer
        number of neighbors to consider for each point.

    out_dim : integer
        number of coordinates for the manifold.

    reg : float
        regularization constant, multiplies the trace of the local covariance
        matrix of the distances.

    eigen_solver : {'arpack', 'lobpcg', 'dense'}
        arpack can handle both dense and sparse data efficiently
        lobpcg solver is usually faster than dense but depends on PyAMG.

    max_iter : integer
        maximum number of iterations for the lobpcg solver.

    method : string ['standard' | 'hessian' | 'modified']
        standard : use the standard locally linear embedding algorithm.
                   see reference [1]
        hessian  : use the hessian eigenmap method.  This method requires
                   n_neighbors > out_dim * (1 + (out_dim + 1) / 2.
                   see reference [2]
        modified : use the modified locally linear embedding algorithm.
                   see reference [3]
        ltsa     : use local tangent space alignment algorithm
                   see reference [4]

    hessian_tol : tolerance used for hessian eigenmapping method
                  only referenced if method == 'hessian'

    modified_tol : tolerance used for modified LLE method
                  only referenced if method == 'modified'

    Returns
    -------
    Y : array-like, shape [n_samples, out_dim]
        Embedding vectors.

    squared_error : float
       Reconstruction error for the embedding vectors. Equivalent to
       norm(Y - W Y, 'fro')**2, where W are the reconstruction weights.

    References
    ----------
      [1] Roweis, S. & Saul, L. Nonlinear dimensionality reduction by
          locally linear embedding.  Science 290:2323 (2000).
      [2] Donoho, D. & Grimes, C. Hessian eigenmaps: Locally linear embedding
          techniques for high-dimensional data. Proc Natl Acad Sci U S A.
          100:5591 (2003).
      [3] Zhang, Z. & Wang, J. MLLE: Modified Locally Linear Embedding
          Using Multiple Weights.
          http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.70.382
      [4] Zhang, Z. & Zha, H. Principal manifolds and nonlinear dimensionality
          reduction via tangent space alignment. Journal of Shanghai Univ.
          8:406 (2004)
    """

    if eigen_solver not in ('arpack', 'lobpcg', 'dense'):
        raise ValueError("unrecognized eigen_solver '%s'" % eigen_solver)

    if method not in ('standard', 'hessian', 'modified', 'ltsa'):
        raise ValueError("unrecognized method '%s'" % method)

    if hasattr(X, 'query'):
        # X is a ball tree
        balltree = X
        X = balltree.data
    else:
        balltree = BallTree(X)

    N, d_in = X.shape

    if out_dim > d_in:
        raise ValueError("output dimension must be less than or equal "
                         "to input dimension")
    if n_neighbors >= N:
        raise ValueError("n_neighbors must be less than number of points")

    M_sparse = (eigen_solver != 'dense')

    if method == 'standard':
        W = kneighbors_graph(
            balltree, n_neighbors=n_neighbors, mode='barycenter', reg=reg)

        # we'll compute M = (I-W)'(I-W)
        # depending on the solver, we'll do this differently
        if M_sparse:
            M = eye(*W.shape, format=W.format) - W
            M = np.dot(M.T, M).tocsr()
        else:
            M = (np.dot(W.T, W) - (W.T + W)).todense()
            M.flat[::M.shape[0] + 1] += 1  # W = W - I

    elif method == 'hessian':
        dp = out_dim * (out_dim + 1) / 2

        if n_neighbors <= out_dim + dp:
            raise ValueError("for method='hessian', n_neighbors must be "
                             "greater than out_dim*[1+(out_dim+1)/2]")

        neighbors = balltree.query(X, k=n_neighbors + 1, return_distance=False)
        neighbors = neighbors[:, 1:]

        X = np.asarray(X)

        Yi = np.empty((n_neighbors, 1 + out_dim + dp), dtype=np.float)
        Yi[:, 0] = 1

        M = np.zeros((N, N), dtype=np.float)

        for i in range(N):
            Gi = X[neighbors[i]]
            Gi -= Gi.mean(0)

            #build hessian estimator
            U, sig, VT = svd(Gi, full_matrices=0)
            Yi[:, 1:1 + out_dim] = U[:, :out_dim]

            j = 1 + out_dim
            for k in range(out_dim):
                Yi[:, j:j + out_dim - k] = U[:, k:k + 1] * U[:, k:out_dim]
                j += out_dim - k

            Q, R = qr(Yi)

            w = Q[:, out_dim + 1:]
            S = w.sum(0)

            S[np.where(abs(S) < H_tol)] = 1
            w /= S

            nbrs_x, nbrs_y = np.meshgrid(neighbors[i], neighbors[i])
            M[nbrs_x, nbrs_y] += np.dot(w, w.T)

        if M_sparse:
            M = csr_matrix(M)

    elif method == 'modified':
        if n_neighbors < out_dim:
            raise ValueError("modified LLE requires n_neighbors >= out_dim")

        neighbors = balltree.query(X, k=n_neighbors + 1, return_distance=False)
        neighbors = neighbors[:, 1:]

        #find the eigenvectors and eigenvalues of each local covariance
        # matrix. We want V[i] to be a [n_neighbors x n_neighbors] matrix,
        # where the columns are eigenvectors
        V = np.zeros((N, n_neighbors, n_neighbors))
        nev = min(d_in, n_neighbors)
        evals = np.zeros([N, nev])
        for i in range(N):
            V[i], evals[i], tmp = np.linalg.svd(X[neighbors[i]] - X[i])
        evals **= 2

        #find regularized weights: this is like normal LLE.
        # because we've already computed the SVD of each covariance matrix,
        # it's faster to use this rather than np.linalg.solve
        reg = 1E-3 * evals.sum(1)

        tmp = np.dot(V.transpose(0, 2, 1), np.ones(n_neighbors))
        tmp[:, :nev] /= evals + reg[:, None]
        tmp[:, nev:] /= reg[:, None]

        w_reg = np.zeros((N, n_neighbors))
        for i in range(N):
            w_reg[i] = np.dot(V[i], tmp[i])
        w_reg /= w_reg.sum(1)[:, None]

        #calculate eta: the median of the ratio of small to large eigenvalues
        # across the points.  This is used to determine s_i, below
        rho = evals[:, out_dim:].sum(1) / evals[:, :out_dim].sum(1)
        eta = np.median(rho)

        #find s_i, the size of the "almost null space" for each point:
        # this is the size of the largest set of eigenvalues
        # such that Sum[v; v in set]/Sum[v; v not in set] < eta
        s_range = np.zeros(N, dtype=int)
        evals_cumsum = np.cumsum(evals, 1)
        eta_range = evals_cumsum[:, -1:] / evals_cumsum[:, :-1] - 1
        for i in range(N):
            s_range[i] = np.searchsorted(eta_range[i, ::-1], eta)
        s_range += n_neighbors - nev  # number of zero eigenvalues

        #Now calculate M.
        # This is the [N x N] matrix whose null space is the desired embedding
        M = np.zeros((N, N), dtype=np.float)
        for i in range(N):
            s_i = s_range[i]

            #select bottom s_i eigenvectors and calculate alpha
            Vi = V[i, :, n_neighbors - s_i:]
            alpha_i = np.linalg.norm(Vi.sum(0)) / np.sqrt(s_i)

            #compute Householder matrix which satisfies
            #  Hi*Vi.T*ones(n_neighbors) = alpha_i*ones(s)
            # using prescription from paper
            h = alpha_i * np.ones(s_i) - np.dot(Vi.T, np.ones(n_neighbors))

            norm_h = np.linalg.norm(h)
            if norm_h < M_tol:
                h *= 0
            else:
                h /= norm_h

            #Householder matrix is
            #  >> Hi = np.identity(s_i) - 2*np.outer(h,h)
            #Then the weight matrix is
            #  >> Wi = np.dot(Vi,Hi) + (1-alpha_i) * w_reg[i,:,None]
            #We do this much more efficiently:
            Wi = (Vi - 2 * np.outer(np.dot(Vi, h), h)
                  + (1 - alpha_i) * w_reg[i, :, None])

            #Update M as follows:
            # >> W_hat = np.zeros( (N,s_i) )
            # >> W_hat[neighbors[i],:] = Wi
            # >> W_hat[i] -= 1
            # >> M += np.dot(W_hat,W_hat.T)
            #We can do this much more efficiently:
            nbrs_x, nbrs_y = np.meshgrid(neighbors[i], neighbors[i])
            M[nbrs_x, nbrs_y] += np.dot(Wi, Wi.T)
            Wi_sum1 = Wi.sum(1)
            M[i, neighbors[i]] -= Wi_sum1
            M[neighbors[i], i] -= Wi_sum1
            M[i, i] += s_i

        if M_sparse:
            M = csr_matrix(M)

    elif method == 'ltsa':
        neighbors = balltree.query(X, k=n_neighbors + 1, return_distance=False)
        neighbors = neighbors[:, 1:]

        M = np.zeros((N, N))

        for i in range(N):
            Xi = X[neighbors[i]]

            # compute out_dim largest eigenvalues of Xi * Xi^T
            v = np.linalg.svd(Xi - Xi.mean(0), full_matrices=True)[0]

            Gi = np.zeros((n_neighbors, out_dim + 1))
            Gi[:, 1:] = v[:, :out_dim]
            Gi[:, 0] = 1. / np.sqrt(n_neighbors)

            GiGiT = np.dot(Gi, Gi.T)

            nbrs_x, nbrs_y = np.meshgrid(neighbors[i], neighbors[i])
            M[nbrs_x, nbrs_y] -= GiGiT
            M[neighbors[i], neighbors[i]] += 1

    return null_space(M, out_dim, k_skip=1, eigen_solver=eigen_solver,
                      tol=tol, max_iter=max_iter,
                      random_state=random_state)


class LocallyLinearEmbedding(BaseEstimator):
    """Locally Linear Embedding

    Parameters
    ----------
    n_neighbors : integer
        number of neighbors to consider for each point.

    out_dim : integer
        number of coordinates for the manifold

    Attributes
    ----------
    `embedding_vectors_` : array-like, shape [out_dim, n_samples]
        Stores the embedding vectors

    `reconstruction_error_` : float
        Reconstruction error associated with `embedding_vectors_`
    """

    def __init__(self, n_neighbors=5, out_dim=2, random_state=None):
        self.n_neighbors = n_neighbors
        self.out_dim = out_dim
        self.random_state = random_state

    def fit(self, X, Y=None, reg=1e-3, eigen_solver='lobpcg', tol=1e-6,
            max_iter=100, **params):
        """Compute the embedding vectors for data X

        Parameters
        ----------
        X : array-like of shape [n_samples, n_features]
            training set.

        out_dim : integer
            number of coordinates for the manifold

        reg : float
            regularization constant, multiplies the trace of the local
            covariance matrix of the distances

        eigen_solver : {'lobpcg', 'dense'}
            use the lobpcg eigensolver or a dense eigensolver based on LAPACK
            routines. The lobpcg solver is usually faster but depends on PyAMG

        max_iter : integer
            maximum number of iterations for the lobpcg solver.

        Returns
        -------
        self : returns an instance of self.
        """
        self.random_state = check_random_state(self.random_state)
        self._set_params(**params)
        self.ball_tree = BallTree(X)
        self.embedding_, self.reconstruction_error_ = \
            locally_linear_embedding(
                self.ball_tree, self.n_neighbors, self.out_dim, reg=reg,
                eigen_solver=eigen_solver, tol=tol, max_iter=max_iter,
                random_state=self.random_state)
        return self

    def transform(self, X, reg=1e-3, **params):
        """
        Transform new points into embedding space.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]

        reg : float
            regularization constant, multiplies the trace of the local
            covariance matrix of the distances

        Returns
        -------
        X_new : array, shape = [n_samples, out_dim]

        Notes
        -----
        Because of scaling performed by this method, it is discouraged to use
        it together with methods that are not scale-invariant (like SVMs)
        """
        self._set_params(**params)
        X = np.atleast_2d(X)
        if not hasattr(self, 'ball_tree'):
            raise ValueError('The model is not fitted')
        ind = self.ball_tree.query(X, k=self.n_neighbors,
                                   return_distance=False)
        weights = barycenter_weights(X, self.ball_tree.data[ind], reg=reg)
        X_new = np.empty((X.shape[0], self.out_dim))
        for i in range(X.shape[0]):
            X_new[i] = np.dot(self.embedding_[ind[i]].T, weights[i])
        return X_new
