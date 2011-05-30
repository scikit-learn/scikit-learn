"""Locally Linear Embedding"""

# Author: Fabian Pedregosa -- <fabian.pedregosa@inria.fr>
# License: BSD, (C) INRIA 2011

import numpy as np
from ..base import BaseEstimator
from ..neighbors import kneighbors_graph, BallTree, barycenter_weights


def locally_linear_embedding(
    X, n_neighbors, out_dim, reg=1e-3, eigen_solver='lobpcg', tol=1e-6,
    max_iter=100):
    """
    Perform a Locally Linear Embedding analysis on the data.

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

    eigen_solver : {'lobpcg', 'dense'}
        use the lobpcg eigensolver or a dense eigensolver based on LAPACK
        routines. The lobpcg solver is usually faster but depends on PyAMG.

    max_iter : integer
        maximum number of iterations for the lobpcg solver.

    Returns
    -------
    Y : array-like, shape [n_samples, out_dim]
        Embedding vectors.

    squared_error : float
       Reconstruction error for the embedding vectors. Equivalent to
       norm(Y - W Y, 'fro')**2, where W are the reconstruction weights.

    References
    ----------
    "An Introduction to Locally Linear Embedding", Lawrence Saul & Sam Roweis.
    """

    if eigen_solver == 'lobpcg':
        try:
            import pyamg
        except ImportError:
            import warnings
            warnings.warn('amg was not found. Using (slow) dense eigensolver')
            eigen_solver = 'dense'

    W = kneighbors_graph(
        X, n_neighbors=n_neighbors, mode='barycenter', reg=reg)

    if eigen_solver == 'dense':
        import scipy.linalg
        M = (np.dot(W.T, W) - (W.T + W)).todense()
        M.flat[::M.shape[0] + 1] += 1  # W = W - I

        eigen_values, eigen_vectors = scipy.linalg.eigh(
            M, eigvals=(1, out_dim + 1), overwrite_a=True)
        index = np.argsort(np.abs(eigen_values))
        return eigen_vectors[:, index], np.sum(eigen_values)

    elif eigen_solver == 'lobpcg':
        from scipy.sparse import linalg, eye
        # M = (I-W)' (I-W)
        A = eye(*W.shape, format=W.format) - W
        A = np.dot(A.T, A).tocsr()

        # initial approximation to the eigenvectors
        X = np.random.rand(W.shape[0], out_dim)
        try:
            ml = pyamg.smoothed_aggregation_solver(A, symmetry='symmetric')
        except TypeError:
            ml = pyamg.smoothed_aggregation_solver(A, mat_flag='symmetric')
        prec = ml.aspreconditioner()

        # compute eigenvalues and eigenvectors with LOBPCG
        eigen_values, eigen_vectors = linalg.lobpcg(
            A, X, M=prec, largest=False, tol=tol, maxiter=max_iter)

        index = np.argsort(eigen_values)
        return eigen_vectors[:, index], np.sum(eigen_values)

    else:
        raise NotImplementedError('Method %s not implemented' % eigen_solver)


class LocallyLinearEmbedding(BaseEstimator):
    """
    Locally Linear Embedding

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

    def __init__(self, n_neighbors=5, out_dim=2):
        self.n_neighbors = n_neighbors
        self.out_dim = out_dim

    def fit(self, X, Y=None, reg=1e-3, eigen_solver='lobpcg', tol=1e-6,
            max_iter=100, **params):
        """
        Compute the embedding vectors for data X.

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
        self._set_params(**params)
        self.ball_tree = BallTree(X)
        self.embedding_, self.reconstruction_error_ = \
            locally_linear_embedding(
                self.ball_tree, self.n_neighbors, self.out_dim, reg=reg,\
                eigen_solver=eigen_solver, tol=tol, max_iter=max_iter)
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
        ind = self.ball_tree.query(X, k=self.n_neighbors, return_distance=False)
        weights = barycenter_weights(X, self.ball_tree.data[ind], reg=reg)
        X_new = np.empty((X.shape[0], self.out_dim))
        for i in range(X.shape[0]):
            X_new[i] = np.dot(weights[i], self.embedding_[ind[i]])
        return X_new
