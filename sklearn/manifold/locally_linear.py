"""Locally Linear Embedding"""

# Author: Fabian Pedregosa -- <fabian.pedregosa@inria.fr>
#         Jake Vanderplas  -- <vanderplas@astro.washington.edu>
# License: BSD, (C) INRIA 2011

import numpy as np
from scipy.linalg import eigh, svd, qr, solve
from scipy.sparse import eye, csr_matrix
from ..base import BaseEstimator
from ..utils import array2d
from ..utils.arpack import eigsh
from ..neighbors import NearestNeighbors


def barycenter_weights(X, Z, reg=1e-3):
    """Compute barycenter weights of X from Y along the first axis

    We estimate the weights to assign to each point in Y[i] to recover
    the point X[i]. The barycenter weights sum to 1.

    Parameters
    ----------
    X : array-like, shape (n_samples, n_dim)

    Z : array-like, shape (n_samples, n_neighbors, n_dim)

    reg: float, optional
        amount of regularization to add for the problem to be
        well-posed in the case of n_neighbors > n_dim

    Returns
    -------
    B : array-like, shape (n_samples, n_neighbors)

    Notes
    -----
    See developers note for more information.
    """
    X = np.asarray(X)
    Z = np.asarray(Z)

    n_samples, n_neighbors = X.shape[0], Z.shape[1]
    if X.dtype.kind == 'i':
        X = X.astype(np.float)
    if Z.dtype.kind == 'i':
        Z = Z.astype(np.float)
    B = np.empty((n_samples, n_neighbors), dtype=X.dtype)
    v = np.ones(n_neighbors, dtype=X.dtype)

    # this might raise a LinalgError if G is singular and has trace
    # zero
    for i, A in enumerate(Z.transpose(0, 2, 1)):
        C = A.T - X[i]  # broadcasting
        G = np.dot(C, C.T)
        trace = np.trace(G)
        if trace > 0:
            R = reg * trace
        else:
            R = reg
        G.flat[::Z.shape[1] + 1] += R
        w = solve(G, v, sym_pos=True)
        B[i, :] = w / np.sum(w)
    return B


def barycenter_kneighbors_graph(X, n_neighbors, reg=1e-3):
    """Computes the barycenter weighted graph of k-Neighbors for points in X

    Parameters
    ----------
    X : {array-like, sparse matrix, BallTree, cKDTree, NearestNeighbors}
        Sample data, shape = (n_samples, n_features), in the form of a
        numpy array, sparse array, precomputed tree, or NearestNeighbors
        object.

    n_neighbors : int
        Number of neighbors for each sample.

    reg : float, optional
        Amount of regularization when solving the least-squares
        problem. Only relevant if mode='barycenter'. If None, use the
        default.

    Returns
    -------
    A : sparse matrix in CSR format, shape = [n_samples, n_samples]
        A[i, j] is assigned the weight of edge that connects i to j.

    See also
    --------
    sklearn.neighbors.kneighbors_graph
    sklearn.neighbors.radius_neighbors_graph
    """
    knn = NearestNeighbors(n_neighbors + 1).fit(X)
    X = knn._fit_X
    n_samples = X.shape[0]
    ind = knn.kneighbors(X, return_distance=False)[:, 1:]
    data = barycenter_weights(X, X[ind], reg=reg)
    indptr = np.arange(0, n_samples * n_neighbors + 1, n_neighbors)
    return csr_matrix((data.ravel(), ind.ravel(), indptr),
                      shape=(n_samples, n_samples))


def null_space(M, k, k_skip=1, eigen_solver='arpack', tol=1E-6, max_iter=100):
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
        dense  : use standard dense matrix operations for the eigenvalue
                    decomposition.  For this method, M must be an array
                    or matrix type.  This method should be avoided for
                    large problems.

    tol : float, optional
        Tolerance for 'arpack' method.
        Not used if eigen_solver=='dense'.

    max_iter : maximum number of iterations for 'arpack' method
        not used if eigen_solver=='dense'
    """
    if eigen_solver == 'auto':
        if M.shape[0] > 200 and k + k_skip < 10:
            eigen_solver = 'arpack'
        else:
            eigen_solver = 'dense'

    if eigen_solver == 'arpack':
        try:
            eigen_values, eigen_vectors = eigsh(M, k + k_skip, sigma=0.0,
                                                tol=tol, maxiter=max_iter)
        except RuntimeError as msg:
            raise ValueError("Error in determining null-space with ARPACK. "
                             "Error message: '%s'. "
                             "Note that method='arpack' can fail when the "
                             "weight matrix is singular or otherwise "
                             "ill-behaved.  method='dense' is recommended. "
                             "See online documentation for more information."
                             % msg)
        except:
            #let other errors pass through
            raise

        return eigen_vectors[:, k_skip:], np.sum(eigen_values[k_skip:])
    elif eigen_solver == 'dense':
        if hasattr(M, 'toarray'):
            M = M.toarray()
        eigen_values, eigen_vectors = eigh(
            M, eigvals=(k_skip, k + k_skip - 1), overwrite_a=True)
        index = np.argsort(np.abs(eigen_values))
        return eigen_vectors[:, index], np.sum(eigen_values)
    else:
        raise ValueError("Unrecognized eigen_solver '%s'" % eigen_solver)


def locally_linear_embedding(
    X, n_neighbors, out_dim, reg=1e-3, eigen_solver='auto',
    tol=1e-6, max_iter=100, method='standard',
    hessian_tol=1E-4, modified_tol=1E-12):
    """Perform a Locally Linear Embedding analysis on the data.

    Parameters
    ----------
    X : {array-like, sparse matrix, BallTree, cKDTree, NearestNeighbors}
        Sample data, shape = (n_samples, n_features), in the form of a
        numpy array, sparse array, precomputed tree, or NearestNeighbors
        object.

    n_neighbors : integer
        number of neighbors to consider for each point.

    out_dim : integer
        number of coordinates for the manifold.

    reg : float
        regularization constant, multiplies the trace of the local covariance
        matrix of the distances.

    eigen_solver : string, {'auto', 'arpack', 'dense'}
        auto : algorithm will attempt to choose the best method for input data

        arpack : use arnoldi iteration in shift-invert mode.
                    For this method, M may be a dense matrix, sparse matrix,
                    or general linear operator.

        dense  : use standard dense matrix operations for the eigenvalue
                    decomposition.  For this method, M must be an array
                    or matrix type.  This method should be avoided for
                    large problems.

    tol : float, optional
        Tolerance for 'arpack' method
        Not used if eigen_solver=='dense'.

    max_iter : integer
        maximum number of iterations for the arpack solver.

    method : {'standard', 'hessian', 'modified', 'ltsa'}
        standard : use the standard locally linear embedding algorithm.
                   see reference [Roweis2000]_
        hessian  : use the Hessian eigenmap method.  This method requires
                   n_neighbors > out_dim * (1 + (out_dim + 1) / 2.
                   see reference [Donoho2003]_
        modified : use the modified locally linear embedding algorithm.
                   see reference [Zhang2007]_
        ltsa     : use local tangent space alignment algorithm
                   see reference [Zhang2004]_

    hessian_tol : float, optional
        Tolerance for Hessian eigenmapping method.
        Only used if method == 'hessian'

    modified_tol : float, optional
        Tolerance for modified LLE method.
        Only used if method == 'modified'

    Returns
    -------
    Y : array-like, shape [n_samples, out_dim]
        Embedding vectors.

    squared_error : float
        Reconstruction error for the embedding vectors. Equivalent to
        ``norm(Y - W Y, 'fro')**2``, where W are the reconstruction weights.

    Notes
    -----
    **References**:

      .. [Roweis2000] `Roweis, S. & Saul, L. Nonlinear dimensionality reduction by
          locally linear embedding.  Science 290:2323 (2000).`
      .. [Donoho2003] `Donoho, D. & Grimes, C. Hessian eigenmaps: Locally linear embedding
          techniques for high-dimensional data. Proc Natl Acad Sci U S A.
          100:5591 (2003).`
      .. [Zhang2007] `Zhang, Z. & Wang, J. MLLE: Modified Locally Linear Embedding
          Using Multiple Weights.`
          http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.70.382
      .. [Zhang2004] `Zhang, Z. & Zha, H. Principal manifolds and nonlinear dimensionality
          reduction via tangent space alignment. Journal of Shanghai Univ.
          8:406 (2004)`
    """
    if eigen_solver not in ('auto', 'arpack', 'dense'):
        raise ValueError("unrecognized eigen_solver '%s'" % eigen_solver)

    if method not in ('standard', 'hessian', 'modified', 'ltsa'):
        raise ValueError("unrecognized method '%s'" % method)

    nbrs = NearestNeighbors(n_neighbors=n_neighbors + 1)
    nbrs.fit(X)
    X = nbrs._fit_X

    N, d_in = X.shape

    if out_dim > d_in:
        raise ValueError("output dimension must be less than or equal "
                         "to input dimension")
    if n_neighbors >= N:
        raise ValueError("n_neighbors must be less than number of points")

    if n_neighbors <= 0:
        raise ValueError("n_neighbors must be positive")

    M_sparse = (eigen_solver != 'dense')

    if method == 'standard':
        W = barycenter_kneighbors_graph(
            nbrs, n_neighbors=n_neighbors, reg=reg)

        # we'll compute M = (I-W)'(I-W)
        # depending on the solver, we'll do this differently
        if M_sparse:
            M = eye(*W.shape, format=W.format) - W
            M = (M.T * M).tocsr()
        else:
            M = (W.T * W - W.T - W).toarray()
            M.flat[::M.shape[0] + 1] += 1  # W = W - I = W - I

    elif method == 'hessian':
        dp = out_dim * (out_dim + 1) / 2

        if n_neighbors <= out_dim + dp:
            raise ValueError("for method='hessian', n_neighbors must be "
                             "greater than [out_dim * (out_dim + 3) / 2]")

        neighbors = nbrs.kneighbors(X, n_neighbors=n_neighbors + 1,
                                    return_distance=False)
        neighbors = neighbors[:, 1:]

        Yi = np.empty((n_neighbors, 1 + out_dim + dp), dtype=np.float)
        Yi[:, 0] = 1

        M = np.zeros((N, N), dtype=np.float)

        use_svd = (n_neighbors > d_in)

        for i in range(N):
            Gi = X[neighbors[i]]
            Gi -= Gi.mean(0)

            #build Hessian estimator
            if use_svd:
                U = svd(Gi, full_matrices=0)[0]
            else:
                Ci = np.dot(Gi, Gi.T)
                U = eigh(Ci)[1][:, ::-1]

            Yi[:, 1:1 + out_dim] = U[:, :out_dim]

            j = 1 + out_dim
            for k in range(out_dim):
                Yi[:, j:j + out_dim - k] = U[:, k:k + 1] * U[:, k:out_dim]
                j += out_dim - k

            Q, R = qr(Yi)

            w = Q[:, out_dim + 1:]
            S = w.sum(0)

            S[np.where(abs(S) < hessian_tol)] = 1
            w /= S

            nbrs_x, nbrs_y = np.meshgrid(neighbors[i], neighbors[i])
            M[nbrs_x, nbrs_y] += np.dot(w, w.T)

        if M_sparse:
            M = csr_matrix(M)

    elif method == 'modified':
        if n_neighbors < out_dim:
            raise ValueError("modified LLE requires n_neighbors >= out_dim")

        neighbors = nbrs.kneighbors(X, n_neighbors=n_neighbors + 1,
                                    return_distance=False)
        neighbors = neighbors[:, 1:]

        #find the eigenvectors and eigenvalues of each local covariance
        # matrix. We want V[i] to be a [n_neighbors x n_neighbors] matrix,
        # where the columns are eigenvectors
        V = np.zeros((N, n_neighbors, n_neighbors))
        nev = min(d_in, n_neighbors)
        evals = np.zeros([N, nev])

        #choose the most efficient way to find the eigenvectors
        use_svd = (n_neighbors > d_in)

        if use_svd:
            for i in range(N):
                X_nbrs = X[neighbors[i]] - X[i]
                V[i], evals[i], _ = svd(X_nbrs,
                                        full_matrices=True)
            evals **= 2
        else:
            for i in range(N):
                X_nbrs = X[neighbors[i]] - X[i]
                C_nbrs = np.dot(X_nbrs, X_nbrs.T)
                evi, vi = eigh(C_nbrs)
                evals[i] = evi[::-1]
                V[i] = vi[:, ::-1]

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
            if norm_h < modified_tol:
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
            # >> M += np.dot(W_hat,W_hat.T)
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
        neighbors = nbrs.kneighbors(X, n_neighbors=n_neighbors + 1,
                                    return_distance=False)
        neighbors = neighbors[:, 1:]

        M = np.zeros((N, N))

        use_svd = (n_neighbors > d_in)

        for i in range(N):
            Xi = X[neighbors[i]]
            Xi -= Xi.mean(0)

            # compute out_dim largest eigenvalues of Xi * Xi^T
            if use_svd:
                v = svd(Xi, full_matrices=True)[0]
            else:
                Ci = np.dot(Xi, Xi.T)
                v = eigh(Ci)[1][:, ::-1]

            Gi = np.zeros((n_neighbors, out_dim + 1))
            Gi[:, 1:] = v[:, :out_dim]
            Gi[:, 0] = 1. / np.sqrt(n_neighbors)

            GiGiT = np.dot(Gi, Gi.T)

            nbrs_x, nbrs_y = np.meshgrid(neighbors[i], neighbors[i])
            M[nbrs_x, nbrs_y] -= GiGiT
            M[neighbors[i], neighbors[i]] += 1

    return null_space(M, out_dim, k_skip=1, eigen_solver=eigen_solver,
                      tol=tol, max_iter=max_iter)


class LocallyLinearEmbedding(BaseEstimator):
    """Locally Linear Embedding

    Parameters
    ----------
    n_neighbors : integer
        number of neighbors to consider for each point.

    out_dim : integer
        number of coordinates for the manifold

    reg : float
        regularization constant, multiplies the trace of the local covariance
        matrix of the distances.

    eigen_solver : string, {'auto', 'arpack', 'dense'}
        auto : algorithm will attempt to choose the best method for input data

        arpack : use arnoldi iteration in shift-invert mode.
                    For this method, M may be a dense matrix, sparse matrix,
                    or general linear operator.

        dense  : use standard dense matrix operations for the eigenvalue
                    decomposition.  For this method, M must be an array
                    or matrix type.  This method should be avoided for
                    large problems.

    tol : float, optional
        Tolerance for 'arpack' method
        Not used if eigen_solver=='dense'.

    max_iter : integer
        maximum number of iterations for the arpack solver.
        Not used if eigen_solver=='dense'.

    method : string ['standard' | 'hessian' | 'modified']
        standard : use the standard locally linear embedding algorithm.
                   see reference [1]
        hessian  : use the Hessian eigenmap method.  This method requires
                   n_neighbors > out_dim * (1 + (out_dim + 1) / 2.
                   see reference [2]
        modified : use the modified locally linear embedding algorithm.
                   see reference [3]
        ltsa     : use local tangent space alignment algorithm
                   see reference [4]

    hessian_tol : float, optional
        Tolerance for Hessian eigenmapping method.
        Only used if method == 'hessian'

    modified_tol : float, optional
        Tolerance for modified LLE method.
        Only used if method == 'modified'

    neighbors_algorithm : string ['auto'|'brute'|'kd_tree'|'ball_tree']
        algorithm to use for nearest neighbors search,
        passed to neighbors.NearestNeighbors instance

    Attributes
    ----------
    `embedding_vectors_` : array-like, shape [out_dim, n_samples]
        Stores the embedding vectors

    `reconstruction_error_` : float
        Reconstruction error associated with `embedding_vectors_`

    `nbrs_` : NearestNeighbors object
        Stores nearest neighbors instance, including BallTree or KDtree
        if applicable.
    """

    def __init__(self, n_neighbors=5, out_dim=2, reg=1E-3,
                 eigen_solver='auto', tol=1E-6, max_iter=100,
                 method='standard', hessian_tol=1E-4, modified_tol=1E-12,
                 neighbors_algorithm='auto'):
        self.n_neighbors = n_neighbors
        self.out_dim = out_dim
        self.reg = reg
        self.eigen_solver = eigen_solver
        self.tol = tol
        self.max_iter = max_iter
        self.method = method
        self.hessian_tol = hessian_tol
        self.modified_tol = modified_tol
        self.nbrs_ = NearestNeighbors(n_neighbors,
                                      algorithm=neighbors_algorithm)

    def _fit_transform(self, X):
        self.nbrs_.fit(X)
        self.embedding_, self.reconstruction_error_ = \
            locally_linear_embedding(
                self.nbrs_, self.n_neighbors, self.out_dim,
                eigen_solver=self.eigen_solver, tol=self.tol,
                max_iter=self.max_iter, method=self.method,
                hessian_tol=self.hessian_tol, modified_tol=self.modified_tol)

    def fit(self, X, y=None):
        """Compute the embedding vectors for data X

        Parameters
        ----------
        X : array-like of shape [n_samples, n_features]
            training set.

        Returns
        -------
        self : returns an instance of self.
        """
        self._fit_transform(X)
        return self

    def fit_transform(self, X, y=None):
        """Compute the embedding vectors for data X and transform X.

        Parameters
        ----------
        X : array-like of shape [n_samples, n_features]
            training set.

        Returns
        -------
        X_new: array-like, shape (n_samples, out_dim)
        """
        self._fit_transform(X)
        return self.embedding_

    def transform(self, X):
        """
        Transform new points into embedding space.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]

        Returns
        -------
        X_new : array, shape = [n_samples, out_dim]

        Notes
        -----
        Because of scaling performed by this method, it is discouraged to use
        it together with methods that are not scale-invariant (like SVMs)
        """
        X = array2d(X)
        ind = self.nbrs_.kneighbors(X, n_neighbors=self.n_neighbors,
                                    return_distance=False)
        weights = barycenter_weights(X, self.nbrs_._fit_X[ind],
                                     reg=self.reg)
        X_new = np.empty((X.shape[0], self.out_dim))
        for i in range(X.shape[0]):
            X_new[i] = np.dot(self.embedding_[ind[i]].T, weights[i])
        return X_new
