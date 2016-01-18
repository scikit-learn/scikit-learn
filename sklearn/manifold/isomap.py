"""Isomap for manifold learning"""

# Author: Jake Vanderplas  -- <vanderplas@astro.washington.edu>
#         Lucas David      -- <ld492@drexel.edu>
# License: BSD 3 clause (C) 2011

import warnings

import numpy as np
from scipy.spatial import distance

from ..base import BaseEstimator, TransformerMixin
from ..decomposition import KernelPCA
from ..neighbors import NearestNeighbors, kneighbors_graph
from ..preprocessing import KernelCenterer
from ..utils import check_array, check_random_state
from ..utils.graph import shortest_path


class Isomap(BaseEstimator, TransformerMixin):
    """Isomap Embedding

    Non-linear dimensionality reduction through Isometric Feature Mapping.

    Read more in the :ref:`User Guide <isomap>`.

    Parameters
    ----------
    n_neighbors : integer
        number of neighbors to consider for each point.

    n_components : integer
        number of coordinates for the manifold

    eigen_solver : ['auto'|'arpack'|'dense']
        'auto' : Attempt to choose the most efficient solver
        for the given problem.

        'arpack' : Use Arnoldi decomposition to find the eigenvalues
        and eigenvectors.

        'dense' : Use a direct solver (i.e. LAPACK)
        for the eigenvalue decomposition.

    tol : float
        Convergence tolerance passed to arpack or lobpcg.
        not used if eigen_solver == 'dense'.

    max_iter : integer
        Maximum number of iterations for the arpack solver.
        not used if eigen_solver == 'dense'.

    path_method : string ['auto'|'FW'|'D'], optional
        Algorithm to use for shortest paths.  Options are:

           'auto' -- (default) select the best among 'FW', 'D', 'BF', or 'J'
                     based on the input data.

           'FW'   -- Floyd-Warshall algorithm.  Computational cost is
                     approximately ``O[N^3]``.  The input csgraph will be
                     converted to a dense representation.

           'D'    -- Dijkstra's algorithm with Fibonacci heaps.  Computational
                     cost is approximately ``O[N(N*k + N*log(N))]``, where
                     ``k`` is the average number of connected edges per node.
                     The input csgraph will be converted to a csr
                     representation.

           'BF'   -- Bellman-Ford algorithm.  This algorithm can be used when
                     weights are negative.  If a negative cycle is encountered,
                     an error will be raised.  Computational cost is
                     approximately ``O[N(N^2 k)]``, where ``k`` is the average
                     number of connected edges per node. The input csgraph will
                     be converted to a csr representation.

           'J'    -- Johnson's algorithm.  Like the Bellman-Ford algorithm,
                     Johnson's algorithm is designed for use when the weights
                     are negative.  It combines the Bellman-Ford algorithm
                     with Dijkstra's algorithm for faster computation.

    neighbors_algorithm : string ['auto'|'brute'|'kd_tree'|'ball_tree']
        Algorithm to use for nearest neighbors search,
        passed to neighbors.NearestNeighbors instance.

    landmarks : [None|'auto'|int|array-like] (default = None)
        The number or list of landmarks to use.
        If this parameter was set, L-Isomap (Landmark Isomap) will execute
        instead of the original algorithm.

        Options are:

            None       -- All samples will be used.
                          The original Isomap algorithm will be executed.

            'auto'     -- Automatically infers the number of landmarks to use.

            int        -- Use exactly landmarks randomly-selected landmarks.
                          Should be a number sufficiently greater than
                          n_components + 1, for stability [2].

            array-like -- Use the landmarks passed in the array-like structure.

    landmarks_method : ['min-max'|'random'] (default = 'min-max')
        Algorithm used to select the landmarks.

        Options are:

            'random'   -- Randomly selects landmarks.

            'min-max'  -- Uses greedy optimization to find well-distributed
                          landmarks.

    n_jobs : int, optional (default = 1)
        The number of parallel jobs to run.
        If ``-1``, then the number of jobs is set to the number of CPU cores.

    random_state : int, optional (default = None)
        The seed for the random landmark selector.
        This will be ignored if landmarks=None or the landmarks
        parameter was passed.


    Attributes
    ----------
    embedding_ : array-like, shape (n_samples, n_components)
        Stores the embedding vectors.

    kernel_pca_ : object
        `KernelPCA` object used to implement the embedding.

    training_data_ : array-like, shape (n_samples, n_features)
        Stores the training data.

    nbrs_ : sklearn.neighbors.NearestNeighbors instance
        Stores nearest neighbors instance, including BallTree or KDtree
        if applicable.

    dist_matrix_ : array-like, shape (n_samples, landmarks)
        Stores the geodesic distance matrix of training data.
        If the original Isomap algorithm was executed (landmarks=None),
        dist_matrix_ has shape (n_samples, n_samples).

    landmarks_ : array-like, shape (landmarks,)
        Stores the indices of the selected landmarks.

    References
    ----------

    .. [1] Tenenbaum, J.B.; De Silva, V.; & Langford, J.C. A global geometric
           framework for nonlinear dimensionality reduction. Science 290 (5500)
    .. [2] Silva, V. D.; & Tenenbaum, J. B. (2002). Global versus local
           methods in nonlinear dimensionality reduction. In Advances
           in neural information processing systems (pp. 705-712).
    .. [3] De Silva, V. D.; & Tenenbaum, J. B. Sparse multidimensional
           scaling using landmark points. Technical report,
           Stanford University, 2004.
    """

    def __init__(self, n_neighbors=5, n_components=2, eigen_solver='auto',
                 tol=0, max_iter=None, path_method='auto',
                 neighbors_algorithm='auto',
                 landmarks=None, landmarks_method='min-max',
                 n_jobs=1, random_state=None):
        self.n_neighbors = n_neighbors
        self.n_components = n_components
        self.eigen_solver = eigen_solver
        self.tol = tol
        self.max_iter = max_iter
        self.path_method = path_method
        self.neighbors_algorithm = neighbors_algorithm
        self.landmarks = landmarks
        self.landmarks_method = landmarks_method
        self.n_jobs = n_jobs
        self.random_state = random_state

    def _compute_landmarks(self):
        """Computes the landmarks in the training data that will be used
        when fitting the data set.

        Returns
        -------
        landmarks_: array-like, shape (landmarks,)
            The array of landmarks to use, or None, if all samples should
            be used.
        """
        if self.landmarks is None or isinstance(self.landmarks, np.ndarray):
            self.landmarks_ = self.landmarks
            return self.landmarks_

        n_samples = self.training_data_.shape[0]
        n_landmarks = self.landmarks

        if n_landmarks == 'auto':
            n_landmarks = min(self.n_components + 10, n_samples)

        if not isinstance(n_landmarks, int):
            raise ValueError("unrecognized landmarks '%s'" % n_landmarks)

        random_state = check_random_state(self.random_state)

        if self.landmarks_method == 'random':
            self.landmarks_ = np.arange(n_samples)
            random_state.shuffle(self.landmarks_)
            self.landmarks_ = self.landmarks_[:n_landmarks]

        elif self.landmarks_method == 'min-max':
            seed = random_state.randint(n_landmarks)
            landmarks = [seed]

            m = distance.cdist(self.training_data_[landmarks],
                               self.training_data_).flatten()

            for i in range(1, n_landmarks):
                landmarks.append(np.argmax(m))
                e = distance.cdist(
                    self.training_data_[landmarks[i], None],
                    self.training_data_).flatten()
                m = np.minimum(m, e)

            self.landmarks_ = np.array(landmarks)
        else:
            raise ValueError("unrecognized landmark method '%s'"
                             % self.landmarks_method)

        return self.landmarks_

    def fit(self, X, y=None):
        """Compute the `n_components`-dimensional embedding of data set X.

        Parameters
        ----------
        X : [array-like|sparse matrix|BallTree|KDTree|NearestNeighbors]
            Sample data, shape = (n_samples, n_features),
            in the form of a numpy array, precomputed tree,
            or NearestNeighbors object.

        y : array-like, optional
            Samples' labels. This information is ignored by ISOMAP.

        Returns
        -------
        self : returns an instance of self.
        """
        X = check_array(X)

        self.nbrs_ = NearestNeighbors(n_neighbors=self.n_neighbors,
                                      algorithm=self.neighbors_algorithm,
                                      n_jobs=self.n_jobs)
        self.nbrs_.fit(X)
        self.training_data_ = self.nbrs_._fit_X

        self.kernel_pca_ = KernelPCA(n_components=self.n_components,
                                     kernel='precomputed',
                                     eigen_solver=self.eigen_solver,
                                     tol=self.tol, max_iter=self.max_iter,
                                     n_jobs=self.n_jobs)

        kng = kneighbors_graph(self.nbrs_, self.n_neighbors,
                               mode='distance', n_jobs=self.n_jobs)

        landmarks = self._compute_landmarks()
        self.dist_matrix_ = shortest_path(kng,
                                          method=self.path_method,
                                          directed=False,
                                          indices=landmarks).T
        infinities = np.isinf(self.dist_matrix_)
        if infinities.any():
            warnings.warn('Neighborhood graph is disconnected, ISOMAP may '
                          'result in a bad lower-dimensional embedding. '
                          'Try to increase n_neighbors parameter.')
            self.dist_matrix_[infinities] = 0

        G = (self.dist_matrix_ if landmarks is None
             else self.dist_matrix_[landmarks])
        G = G ** 2
        G *= -.5

        self.embedding_ = self.kernel_pca_.fit_transform(G)

        if landmarks is not None:
            # Project the remaining samples.
            self.embedding_ = self.transform(X)

        return self

    def fit_transform(self, X, y=None, **fit_params):
        """Fit the model from data contained in X and transform X.

        Parameters
        ----------
        X : [array-like|sparse matrix|BallTree|KDTree|NearestNeighbors]
            Sample data, shape = (n_samples, n_features),
            in the form of a numpy array, precomputed tree,
            or NearestNeighbors object.

        y : array-like, optional
            Samples' labels. This information is ignored by ISOMAP.

        Returns
        -------
        X_new : array-like, shape (n_samples, n_components)
        """
        return self.fit(X, y).embedding_

    def transform(self, X):
        """Transform X.

        This is implemented by linking the points X into the graph of geodesic
        distances of the training data. First the `n_neighbors` nearest
        neighbors of X are found in the training data, and from these the
        shortest geodesic distances from each point in X to each landmark in
        the training data are computed in order to construct the kernel.
        The embedding of X is the projection of this kernel onto the
        embedding vectors of the training set.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)

        Returns
        -------
        X_new : array-like, shape (n_samples, n_components)
        """
        X = check_array(X)
        distances, indices = self.nbrs_.kneighbors(X, return_distance=True)

        # Create the graph of shortest distances from X to self.training_data_
        # (or to the landmarks, when executing L-Isomap) via the nearest
        # neighbors of X.
        # This can be done as a single array operation, but it potentially
        # takes a lot of memory.  To avoid that, use a loop:
        columns = (self.training_data_.shape[0] if self.landmarks_ is None
                   else self.landmarks_.shape[0])

        G_X = np.zeros((X.shape[0], columns))
        for i in range(X.shape[0]):
            G_X[i] = np.min((self.dist_matrix_[indices[i]] +
                             distances[i][:, None]), axis=0)

        G_X **= 2
        G_X *= -0.5

        return self.kernel_pca_.transform(G_X)

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
        and K is the isomap kernel:

        ``K(D) = -0.5 * (I - 1/n_samples) * D^2 * (I - 1/n_samples)``
        """
        G = (self.dist_matrix_
             if self.landmarks_ is None
             else self.dist_matrix_[self.landmarks_])
        G = -0.5 * G ** 2
        G_center = KernelCenterer().fit_transform(G)
        evals = self.kernel_pca_.lambdas_
        return np.sqrt(np.sum(G_center ** 2) - np.sum(evals ** 2)) / G.shape[0]
