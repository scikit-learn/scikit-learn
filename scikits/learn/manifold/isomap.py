"""Isomap for manifold learning"""

# Author: Jake Vanderplas  -- <vanderplas@astro.washington.edu>
# License: BSD, (C) 2011

import numpy as np
from ..base import BaseEstimator
from ..neighbors import kneighbors_graph, NeighborsClassifier
from ..utils.graph_shortest_path import graph_shortest_path

from ..decomposition import KernelPCA


class Isomap(BaseEstimator):
    """Isomap Embedding

    Non-linear dimensionality reduction through Isometric Mapping

    Parameters
    ----------
    n_neighbors : integer
        number of neighbors to consider for each point.

    out_dim : integer
        number of coordinates for the manifold

    eigen_solver : ['auto'|'arpack'|'dense']
        'auto' : attempt to choose the most efficient solver
            for the given problem.
        'arpack' : use Arnoldi decomposition to find the eigenvalues
            and eigenvectors.  Note that arpack can handle both dense
            and sparse data efficiently
        'dense' : use a direct solver (i.e. LAPACK)
            for the eigenvalue decomposition.

    tol : float
        convergence tolerance passed to arpack or lobpcg.
        not used if eigen_solver == 'dense'

    max_iter : integer
        maximum number of iterations for the arpack solver.
        not used if eigen_solver == 'dense'

    path_method : string ['auto'|'FW'|'D']
        method to use in finding shortest path.
        'auto' : attempt to choose the best algorithm automatically
        'FW' : Floyd-Warshall algorithm
        'D' : Dijkstra algorithm with Fibonacci Heaps

    Attributes
    ----------
    `embedding_` : array-like, shape (n_samples, out_dim)
        Stores the embedding vectors

    `kernel_pca_` : `KernelPCA` object used to implement the embedding

    `training_data_` : array-like, shape (n_samples, n_features)
        Stores the training data

    `dist_matrix_` : array-like, shape (n_samples, n_samples)
        Stores the geodesic distance matrix of training data

    References
    ----------
    [1] Tenenbaum, J.B.; De Silva, V.; & Langford, J.C. A global geometric
        framework for nonlinear dimensionality reduction. Science 290 (5500)
    """

    def __init__(self, n_neighbors=5, out_dim=2,
                 eigen_solver='auto', tol=0,
                 max_iter=None, path_method='auto'):
        self.n_neighbors = n_neighbors
        self.out_dim = out_dim
        self.eigen_solver = eigen_solver
        self.tol = tol
        self.max_iter = max_iter
        self.path_method = path_method

    def _fit_transform(self, X):
        self.training_data_ = X
        self.kernel_pca_ = KernelPCA(n_components=self.out_dim,
                                     kernel="precomputed",
                                     eigen_solver=self.eigen_solver,
                                     tol=self.tol, max_iter=self.max_iter)

        # it would be best to store the BallTree of X here, but there's no
        # good way to do this without duplicating kneighbors_graph code.
        kng = kneighbors_graph(X, self.n_neighbors,
                               mode='distance')

        self.dist_matrix_ = graph_shortest_path(kng,
                                                method=self.path_method,
                                                directed=False)
        G = self.dist_matrix_ ** 2
        G *= -0.5

        self.embedding_ = self.kernel_pca_.fit_transform(G)

    def fit(self, X, Y=None, **params):
        """Compute the embedding vectors for data X

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            training set.

        Returns
        -------
        self : returns an instance of self.
        """
        self._set_params(**params)
        self._fit_transform(X)
        return self

    def fit_transform(self, X, Y=None, **params):
        """Fit the model from data in X and transform X.

        Parameters
        ----------
        X: array-like, shape (n_samples, n_features)
            Training vector, where n_samples in the number of samples
            and n_features is the number of features.

        Returns
        -------
        X_new: array-like, shape (n_samples, out_dim)
        """
        self._set_params(**params)
        self._fit_transform(X)
        return self.embedding_

    def transform(self, X, **params):
        """Transform X.

        Parameters
        ----------
        X: array-like, shape (n_samples, n_features)

        Returns
        -------
        X_new: array-like, shape (n_samples, out_dim)
        """
        neighbors = NeighborsClassifier(self.n_neighbors)
        neighbors.fit(self.training_data_, 0)
        distances, indices = neighbors.kneighbors(X)

        #Create the graph of shortest distances from X to self.training_data_
        # via the nearest neighbors of X.
        #This can be done as a single array operation, but it potentially
        # takes a lot of memory.  To avoid that, use a loop:
        G_X = np.zeros((X.shape[0], self.training_data_.shape[0]))
        for i in range(X.shape[0]):
            G_X[i] = np.min((self.dist_matrix_[indices[i]]
                             + distances[i][:, None]), 0)

        G_X **= 2
        G_X *= -0.5

        return self.kernel_pca_.transform(G_X)
