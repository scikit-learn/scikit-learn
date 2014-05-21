"""Incremental Isomap for Manifold Learning"""

# Author: Peter Fischer <peter.fischer@fau.de>
# License: BSD, (C) 2014
from __future__ import division, print_function
import collections

import numpy as np
import scipy.sparse

from ..base import BaseEstimator, TransformerMixin
from ..neighbors import kneighbors_graph
from ..utils import check_arrays
from ..decomposition.kernel_pca import KernelPCA
from ..preprocessing import KernelCenterer
from ..metrics.pairwise import euclidean_distances


from ..utils.graph import graph_shortest_path
from ..utils.graph import modified_graph_shortest_path
from ..utils.incremental_isomap_utilities import _update_insert, _update_edge
from ..utils.incremental_isomap_utilities import _construct_f, _determine_edge_changes


class IncrementalIsomap(BaseEstimator, TransformerMixin):
    """Isomap Embedding

    Non-linear dimensionality reduction through Isomap. Incremental using
    fit_transform_incremental and fit_incremental.
    Constraint: incremental fitting with 1 new sample at a time

    Parameters
    ----------
    n_neighbors : integer
        number of neighbors to consider for each point.

    n_components : integer
        number of coordinates for the manifold

    eigen_solver : string ['auto'|'arpack'|'dense']
        'auto' : Attempt to choose the most efficient solver
            for the given problem.
        'arpack' : Use Arnoldi decomposition to find the eigenvalues
            and eigenvectors.
        'dense' : Use a direct solver (i.e. LAPACK)
            for the eigenvalue decomposition.
        Only for initialization. In incremental phase, subspace iterations are
        used.

    tol : float
        Convergence tolerance passed to arpack or lobpcg.
        not used if eigen_solver == 'dense'.

    max_iter : integer
        Maximum number of iterations for the arpack solver.
        not used if eigen_solver == 'dense'.

    path_method : string ['auto'|'FW'|'D']
        Method to use in finding shortest path.
        'auto' : attempt to choose the best algorithm automatically
        'FW'   : Floyd-Warshall algorithm
        'D'    : Dijkstra algorithm with Fibonacci Heaps
        Only for initialization. In incremental phase, a modified Dijkstra
        algorithm is used.

    n_iter : integer
        Number of iterations of the incremental eigendecomposition method.
        Important for accuracy of the embedding.

    n_max_samples : integer
        Number of samples that can maximally be stored. Memory for that many
        samples is allocated

    overflow_mode : string ['extend'| 'adapt' | 'embed']
        Selects what happens after more than n_max_samples where observed
        'extend' : n_max_samples is increased
        'adapt'  : delete oldest samples
        'embed'  : just transform, do not learn anymore

    Attributes
    ----------
    `embedding_` : array-like, shape (n_samples, n_components)
        Stores the embedding vectors.

    `kernel_pca_` : object
        `KernelPCA` object used to implement the embedding.

    `train_data_` : array-like, shape (n_samples, n_features)
        Stores the training data.

    `train_data_norm2_` : array-like, shape (1, n_samples)
        Stores the norm of each training sample. Avoids recomputing it in each
        iteration

    `dist_matrix_` : array-like, shape (n_samples, n_samples)
        Stores the geodesic distance matrix of training data.

    `predecessor_matrix_` : array-like, shape (n_samples, n_samples)
        Stores the shortest paths used to compute the geodesic distances.
        Necessary for updating geodesic distances incrementally

    `kng_` : array-like, shape (n_samples, n_samples)
        Stores the k nearest neighbor graph, encoding the most similar images.

    References
    ----------

    [1] Tenenbaum, J.B.; De Silva, V.; & Langford, J.C. A global geometric
        framework for nonlinear dimensionality reduction. Science 290 (5500)
    [2] Law, M.; Jain, A. Incremental Nonlinear Dimensionality Reduction
        by Manifold Learning. IEEE PAMI 28(3), 377-391, 2006
    """

    def __init__(self, n_neighbors=5, n_components=2, eigen_solver='auto',
                 tol=0, max_iter=None, n_iter=10, path_method='auto',
                 n_max_samples=600, overflow_mode="adapt"):

        self.n_neighbors = n_neighbors
        self.n_components = n_components
        self.eigen_solver = eigen_solver
        self.tol = tol
        self.max_iter = max_iter
        self.n_iter = n_iter
        self.path_method = path_method
        self.n_max_samples = n_max_samples
        self.n = 0  # number of processed samples
        self.i = 0  # current index of sample to delete/insert
        self.overflow_mode = overflow_mode
        return

    def _fit_transform(self, X):
        X, = check_arrays(X, sparse_format='dense')
        self.n = X.shape[0]
        self.n_max_samples = max(self.n_max_samples, self.n)
        self.i = self.n - 1
        self.train_data_ = X
        self.train_data_norm2_ = np.sum(X * X, axis=1)[np.newaxis, :]
        self.kernel_pca_ = KernelPCA(n_components=self.n_components,
                                     kernel="precomputed",
                                     eigen_solver=self.eigen_solver,
                                     tol=self.tol, max_iter=self.max_iter)

        self.kng_ = kneighbors_graph(X, self.n_neighbors,
                               mode='distance').tolil()

        # Create symmetric version of self.kng_
        self.kng_symmetric_ = self.kng_.copy()
        for i in xrange(self.kng_.shape[0]):
            for j in self.kng_.getrowview(i).rows[0]:
                if self.kng_[i, j] > 0:
                    self.kng_symmetric_[j, i] = self.kng_[i, j]
                    self.kng_symmetric_[i, j] = self.kng_[i, j]

        # Compute geodesic distances
        self.dist_matrix_, self.predecessor_matrix_ = graph_shortest_path(
                                                self.kng_,
                                                method=self.path_method,
                                                directed=False)

        # Embed to low dimensions
        G = self.dist_matrix_ ** 2
        G *= -0.5
        self.embedding_ = self.kernel_pca_.fit_transform(G)
        self.embedding_ *= np.sqrt(self.n) / np.sqrt(self.kernel_pca_.lambdas_)
        self._resize_data_structures()
        return

    def _resize_data_structures(self):
        """Resize matrices to n_max_samples.

        All matrices necessary are increased to size n_max_samples.
        Old values are copied into the new matrices.
        """
        # Create matrices of self.n_max_samples size
        shape = (self.n_max_samples, self.n_max_samples)
        new = scipy.sparse.lil_matrix(shape, dtype=self.kng_.dtype)
        new.rows[:self.n] = self.kng_.rows
        new.data[:self.n] = self.kng_.data
        self.kng_ = new

        new = scipy.sparse.lil_matrix(shape, dtype=self.kng_symmetric_.dtype)
        new.rows[:self.n] = self.kng_symmetric_.rows
        new.data[:self.n] = self.kng_symmetric_.data
        self.kng_symmetric_ = new

        new = -np.ones(shape, dtype=self.predecessor_matrix_.dtype)
        new[:self.n, :self.n] = self.predecessor_matrix_
        self.predecessor_matrix_ = new

        new = np.inf * np.ones(shape, dtype=self.dist_matrix_.dtype)
        new[:self.n, :self.n] = self.dist_matrix_
        np.fill_diagonal(new, 0)
        self.dist_matrix_ = new

        new = np.zeros((self.n_max_samples, self.n_components),
                       dtype=self.embedding_.dtype)
        new[:self.n, :] = self.embedding_
        self.embedding_ = new

        new = np.zeros((self.n_max_samples, self.train_data_.shape[1]),
                       np.float)
        new[:self.n, :] = self.train_data_
        self.train_data_ = new

        new = np.zeros((1, self.n_max_samples), np.float)
        new[:, :self.n] = self.train_data_norm2_
        self.train_data_norm2_ = new
        return

    def _fit_transform_incremental(self, X):
        # Update data structures for this iteration
        if self.n < self.n_max_samples:
            self.n += 1
            self.i += 1
        elif self.overflow_mode == "embed":
            return self.transform(X)
        elif self.overflow_mode == "adapt":
            # Forgetting of old points
            # Calculate next point with wrapping
            self.i = (self.i + 1) % self.n_max_samples

            ## Build set of deleted points: all points incident on oldest point
            deleted_edges = []
            # delete edges for nearest neighbors of self.i
            for e in self.kng_.rows[self.i]:
                deleted_edges.append((self.i, e))

            # delete edges where self.i is nearest neighbor
            affected_indices = self.kng_.tocsc()[:, self.i].indices
            added_edges = []
            for e in affected_indices:
                deleted_edges.append((e, self.i))
                # TODO: avoid recomputation of distances -> very expensive
                knn_dist, knn_point = self._find_nearest_neighbors(
                                    self.train_data_[e, :].reshape(1, -1),
                                    exclude_sample=self.i)
                # One point is deleted, next larger distance might be added
                dist = knn_dist[0, self.n_neighbors]
                idx = knn_point[0, self.n_neighbors]
                if self.kng_[e, idx] < np.spacing(0):
                    added_edges.append(((e, idx), dist))

            ## Actually add and remove edges
            self._delete_edges(deleted_edges)

            # Compute add new edges to kng_ and to geodesics for this graph
            for ((j1, j2), dist) in added_edges:
                # add to neighborhood graph
                self.kng_[j1, j2] = dist
                self.kng_symmetric_[j1, j2] = dist
                self.kng_symmetric_[j2, j1] = dist
                # compute geodesic distances
                self.dist_matrix_[j1, j2] = dist
                self.dist_matrix_[j2, j1] = dist
                # compute predecessors
                self.predecessor_matrix_[j1, j2] = j1
                self.predecessor_matrix_[j2, j1] = j2

            # Shorten geodesic distances using new edges
            for ((a, b), dist) in added_edges:
                a_rowview = self.kng_symmetric_.getrowview(a)
                b_rowview = self.kng_symmetric_.getrowview(b)
                _update_edge(a,
                             np.array(a_rowview.rows[0], np.int),
                             b,
                             np.array(b_rowview.rows[0], np.int),
                             self.n + 1,
                             self.dist_matrix_,
                             self.predecessor_matrix_)

            # self.n does not increase due to deletion
        elif self.overflow_mode == "extend":
            # extend and larger than n_max_samples: double size of all matrices
            self.n_max_samples *= 2
            self._resize_data_structures()
            self.n += 1
            self.i += 1
        else:
            raise NotImplementedError("Overflow_mode " +
                                      str(self.overflow_mode) +
                                      " not implemented.")

        # Initialization
        X, = check_arrays(X, sparse_format='dense')

        # store training data
        self.train_data_[self.i, :] = X
        self.train_data_norm2_[:, self.i] = \
            np.sum(X * X, axis=1)[np.newaxis, :]

        # Find kNN of new point
        knn_dist, knn_point = self._find_nearest_neighbors(X,
                                            exclude_sample=self.i)

        # Find added and deleted edges caused by the new point
        knn_dist = np.ascontiguousarray(knn_dist, dtype=np.float64)
        knn_point = np.ascontiguousarray(knn_point, dtype=np.int)
        added_edges, deleted_edges = _determine_edge_changes(
                                                self.i,
                                                self.n,
                                                self.n_neighbors,
                                                knn_dist,
                                                knn_point,
                                                self.kng_,
                                                self.train_data_)

        ## Process deleted edges
        self._delete_edges(deleted_edges)

        ## Process added edges
        # Compute geodesic distances for new point
        self._add_point(added_edges)

        # Shorten geodesic distances using new point
        _update_insert(self.i, self.n + 1,
                       np.array(self.kng_symmetric_.getrowview(self.i).rows[0],
                                np.int),
                       self.dist_matrix_,
                       self.predecessor_matrix_)

        ## Embed new point
        # Convert to kernel
        G = self.dist_matrix_[:self.n, :self.n] ** 2
        G *= -0.5
        B = KernelCenterer().fit_transform(G)
        K = B[np.arange(self.n) != self.i, self.i]
        idx = np.arange(self.n) != self.i
        old_embedding = self.embedding_[:self.n, :][idx]
        new_embedding = np.dot(K, old_embedding / self.kernel_pca_.lambdas_)
        self.embedding_[self.i, :] = new_embedding.reshape(1, -1)

        ## Update all coordinates
        # incremental eigenvalue problem
        # converges to solution of linalg.eigh in few iterations
        for _ in xrange(self.n_iter):
            Z = np.dot(B, self.embedding_[:self.n, :] / np.sqrt(self.n))
            Q, R = np.linalg.qr(Z)
            Q *= np.sign(np.diag(R))
            Zstar = np.dot(Q.T, np.dot(B, Q))
            lambdas, alphas = np.linalg.eigh(Zstar)
            indices = lambdas.argsort()[::-1]
            lambdas = lambdas[indices]
            alphas = alphas[:, indices]
            # New eigenvalues and eigenvectors
            self.kernel_pca_.lambdas_ = lambdas
            self.kernel_pca_.alphas_ = np.dot(Q, alphas)
            # New embedding
            self.embedding_[:self.n, :] = np.sqrt(self.n) * \
                                            self.kernel_pca_.alphas_
        return self.embedding_[self.i, :].reshape(1, -1)

    def reconstruction_error(self):
        """Compute the reconstruction error for the embedding.

        Returns
        -------
        reconstruction_error : float

        Notes
        -------
        The cost function of an isomap embedding is

        ``E = frobenius_norm[K(D) - K(D_fit)] / n_samples``

        Where D is the matrix of distances for the input data X,
        D_fit is the matrix of distances for the output embedding X_fit,
        and K is the operator to transform distances to a kernel:

        ``K(D) = -0.5 * (I - 1/n_samples) * D^2 * (I - 1/n_samples)``
        """
        G = -0.5 * self.dist_matrix_[:self.n, :self.n] ** 2
        G_center = KernelCenterer().fit_transform(G)
        evals = self.embedding_[:self.n, :] / np.sqrt(self.n) * \
                    np.sqrt(self.kernel_pca_.lambdas_)
        euclidean_matrix = euclidean_distances(evals, evals, squared=True)
        K_center = KernelCenterer().fit_transform(-0.5 * euclidean_matrix)
        temp = np.sqrt(np.sum(G_center ** 2) - np.sum(K_center ** 2))
        return temp / G.shape[0]

    def fit(self, X, y=None):
        """Compute the embedding vectors for data X

        Parameters
        ----------
        X : {array-like, sparse matrix, BallTree, cKDTree, NearestNeighbors}
            Sample data, shape = (n_samples, n_features), in the form of a
            numpy array, precomputed tree, or NearestNeighbors
            object.

        Returns
        -------
        self : returns an instance of self.
        """
        self._fit_transform(X)
        return self

    def fit_transform(self, X, y=None):
        """Fit the model from data in X and transform X.

        Parameters
        ----------
        X: {array-like, sparse matrix, BallTree, cKDTree}
            Training vector, where n_samples in the number of samples
            and n_features is the number of features.

        Returns
        -------
        X_new: array-like, shape (n_samples, n_components)
        """
        self._fit_transform(X)
        return self.embedding_[:self.n, :]

    def fit_incremental(self, X, y=None):
        """Update the model from data in X

        Parameters
        ----------
        X : {array-like, sparse matrix}
            SNew data point, shape = (1, n_features)

        Returns
        -------
        self : returns an instance of self.
        """
        self._fit_transform_incremental(X)
        return self

    def fit_transform_incremental(self, X, y=None):
        """Update the model from data in X and transform X.

        Parameters
        ----------
        X: {array-like, sparse matrix}
            New data point, shape = (1, n_features)

        Returns
        -------
        X_new: array-like, shape (1, n_components)
        """
        return self._fit_transform_incremental(X)

    def transform(self, X):
        """Transform X.

        This is implemented by linking the points X into the graph of geodesic
        distances of the training data. First the `n_neighbors` nearest
        neighbors of X are found in the training data, and from these the
        shortest geodesic distances from each point in X to each point in
        the training data are computed in order to construct the kernel.
        The embedding of X is the projection of this kernel onto the
        embedding vectors of the training set.

        Parameters
        ----------
        X: array-like, shape (n_samples, n_features)

        Returns
        -------
        X_new: array-like, shape (n_samples, n_components)
        """
        # fit kernel centerer to train_data_
        G = self.dist_matrix_[:self.n, :self.n] ** 2
        G *= -0.5
        kernelCenterer = KernelCenterer().fit(G)

        # Find nearest neighbors of X in train_data_
        distances, indices = self._find_nearest_neighbors(X)
        indices = indices[:, :self.n_neighbors]
        distances = distances[:, :self.n_neighbors]

        #Create the graph of shortest distances from X to self.train_data_
        # via the nearest neighbors of X.
        #This can be done as a single array operation, but it potentially
        # takes a lot of memory.  To avoid that, use a loop:
        G_X = np.zeros((X.shape[0], self.n))
        for i in xrange(X.shape[0]):
            G_X[i] = np.min((self.dist_matrix_[indices[i], :self.n]
                             + distances[i][:, None]), 0)
        G_X **= 2
        G_X *= -0.5

        K = kernelCenterer.transform(G_X)
        temp = np.dot(K, self.kernel_pca_.alphas_ / self.kernel_pca_.lambdas_)
        return  temp * np.sqrt(self.n)

    def _find_nearest_neighbors(self, X, exclude_sample=None):
        '''Find nearest neighbors in the train_data_

        Adapted from sklearn pairwise and kneighbors.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            samples to find neighbors for
        exclude_sample : int
            which element of train_data_ is ignored for distance computation

        Returns
        -------
        distances : array-like, shape (n_samples, self.n)
            distances to all training points (sorted)
        indices : array-like, shape (n_samples, self.n)
            indices to which training sample the distance corresponds
        '''
        # use pairwise.euclidean_distances with precomputed norm (MUCH faster)
        distances = euclidean_distances(X, self.train_data_[:self.n, :],
                                        Y_norm_squared=self.train_data_norm2_[:, :self.n],
                                        squared=True)
        # adapted from KNeighborsMixin.kneighbors
        neigh_ind = distances.argsort(axis=1)
        if exclude_sample is not None:
            neigh_ind = neigh_ind[neigh_ind != exclude_sample].reshape(
                                                                    X.shape[0],
                                                                    self.n - 1)
        # return sorted distances sorted (necessary for later)
        j = np.arange(X.shape[0])[:, None]
        distances = np.sqrt(distances[j, neigh_ind].reshape(neigh_ind.shape))
        return distances, neigh_ind

    def _add_point(self, added_edges):
        '''Add edges to graphs.

        Add a point to the Isomap problem by computing the geodesic distances
        to all points using the added_edges.

        Parameters
        ----------
        added_edges : list
            list of edges to be added,
            format ((edge_start, edge_end), distance)
        '''
        # TODO: slow loop --> speed up with cython
        for ((j1, j2), dist) in added_edges:
            # add to neighborhood graph
            self.kng_[j1, j2] = dist
            self.kng_symmetric_[j1, j2] = dist
            self.kng_symmetric_[j2, j1] = dist
            # compute geodesic distances
            if j1 == self.i:  # select index of point != added point
                j = j2
            else:
                j = j1
            idx = self.dist_matrix_[self.i, :self.n] > \
                    (self.dist_matrix_[j, :self.n] + dist)
            self.dist_matrix_[self.i, :self.n][idx] = \
                                    self.dist_matrix_[j, :self.n][idx] + dist
            # compute predecessors
            self.predecessor_matrix_[idx, self.i] = j
            self.predecessor_matrix_[self.i, idx] = \
                                    self.predecessor_matrix_[j, idx]
            if idx[j]:
                self.predecessor_matrix_[self.i, j] = self.i
        self.dist_matrix_[:self.n, self.i] = \
                                self.dist_matrix_[self.i, :self.n].T
        return

    def _delete_edges(self, deleted_edges):
        '''Delete edges from graph.

        Remove edges from Isomap problem and recompute geodesic distances that
        changed.

        Parameters
        ----------
        deleted_edges : list
            list of edges to be removed, format (edge_start, edge_end)
        '''
        # Construct set of all influenced vertex pairs
        F = _construct_f(deleted_edges, self.n,
                         self.kng_symmetric_, self.predecessor_matrix_)

        # Delete edges from self.kng after _construct_f
        for e in deleted_edges:
            self.kng_[e[0], e[1]] = 0
            if self.kng_[e[1], e[0]] == 0:
                self.kng_symmetric_[e[1], e[0]] = 0
                self.kng_symmetric_[e[0], e[1]] = 0

        # Reconstruct shortest paths
        B_auxiliary = F
        l = [collections.deque([]) for _ in xrange(self.n)]
        l_index = np.zeros((self.n,), np.int)
        for i in xrange(self.n):
            # -1 to convert to 0-based indexing
            i_deg = np.count_nonzero(B_auxiliary[i, :]) - 1
            if i_deg >= 0:
                l[i_deg].append(i)
                l_index[i] = i_deg
        pos = 0
        self.kng_symmetric_ = self.kng_symmetric_.tocsr()
        for i in xrange(self.n):
            try:
                while not l[pos]:
                    pos += 1
            except IndexError:
                break
            u = l[pos].pop()
            l_index[u] = 0
            C_u = np.ascontiguousarray(np.nonzero(B_auxiliary[u, :])[0],
                                       dtype=np.int)
            self.dist_matrix_, self.predecessor_matrix_ = \
                modified_graph_shortest_path(u, self.n, C_u,
                                             self.kng_symmetric_,
                                             self.dist_matrix_,
                                             self.predecessor_matrix_,
                                             directed=False)
            # Decrease degree of elements in C_u
            for j in C_u:
                l[l_index[j]].remove(j)
                l_index[j] -= 1
                if not l_index[j] < 0:
                    l[l_index[j]].append(j)
                pos = min(pos, max(l_index[j], 0))
            B_auxiliary[u, :] = 0
            B_auxiliary[:, u] = 0
        self.kng_symmetric_ = self.kng_symmetric_.tolil()
        return
