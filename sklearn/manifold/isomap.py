"""Isomap for manifold learning"""

# Author: Jake Vanderplas  -- <vanderplas@astro.washington.edu>
# License: BSD 3 clause (C) 2011
import numpy as np

from ..base import BaseEstimator, TransformerMixin
from ..decomposition import KernelPCA
from ..neighbors import NearestNeighbors, kneighbors_graph
from ..preprocessing import KernelCenterer
from ..utils import check_array, check_random_state
from ..utils.graph import graph_shortest_path


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

    path_method : string ['auto'|'FW'|'D']
        Method to use in finding shortest path.

        'auto' : attempt to choose the best algorithm automatically.

        'FW' : Floyd-Warshall algorithm.

        'D' : Dijkstra's algorithm.

    neighbors_algorithm : string ['auto'|'brute'|'kd_tree'|'ball_tree']
        Algorithm to use for nearest neighbors search,
        passed to neighbors.NearestNeighbors instance.

    n_landmarks : [None|int|'auto'] (default = None)
        The number of landmarks. If this parameter was set,
        L-Isomap (Landmark Isomap) will execute instead of the
        original algorithm.

        Should be a number sufficiently greater than n_components + 1, for
        stability [2].

        None : all samples will be used. The original Isomap algorithm will be
        executed.

        'auto' : automatically computes the number of landmarks to use and
        randomly-selects them.

        int :  use exactly n_landmarks randomly-selected landmarks.

    landmarks : array-like, optional (default = None)
        The list of landmarks to use. If none was passed,
        infers it from n_landmarks.

    n_jobs : int, optional (default = 1)
        The number of parallel jobs to run.
        If ``-1``, then the number of jobs is set to the number of CPU cores.

    random_state : int, optional (default = None)
        The seed for the random landmark selector.
        This will be ignored if n_landmarks=None or the landmarks
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

    dist_matrix_ : array-like, shape (n_samples, n_landmarks)
        Stores the geodesic distance matrix of training data.
        If the original Isomap algorithm was executed (n_landmarks=None),
        dist_matrix_ has shape (n_samples, n_samples).

    landmarks_ : array-like, shape (n_landmarks,)
        Stores the indices of the selected landmarks.

    References
    ----------

    .. [1] Tenenbaum, J.B.; De Silva, V.; & Langford, J.C. A global geometric
           framework for nonlinear dimensionality reduction. Science 290 (5500)
    .. [2] Silva, V. D., & Tenenbaum, J. B. (2002). Global versus local
           methods in nonlinear dimensionality reduction. In Advances
           in neural information processing systems (pp. 705-712).
    """

    def __init__(self, n_neighbors=5, n_components=2, eigen_solver='auto',
                 tol=0, max_iter=None, path_method='auto',
                 neighbors_algorithm='auto', n_landmarks=None, landmarks=None,
                 n_jobs=1, random_state=None):
        self.n_neighbors = n_neighbors
        self.n_components = n_components
        self.eigen_solver = eigen_solver
        self.tol = tol
        self.max_iter = max_iter
        self.path_method = path_method
        self.neighbors_algorithm = neighbors_algorithm
        self.n_landmarks = n_landmarks
        self.landmarks_ = landmarks
        self.n_jobs = n_jobs
        self.random_state = random_state

    def _compute_landmarks(self, X):
        """Computes the landmarks from a data set X.

        If n_landmarks and landmarks parameters were passed,
        all samples are used as landmarks (original Isomap algorithm).
        If the landmarks were passed in the constructor or they were
        already computed, returns them.
        Otherwise, randomly selects n_landmarks (contained in the
        interval [0, |X|)).

        If n_landmarks is 'auto', selects n_components + 1 (the necessary
        landmarks required to triangulate a samples' position in the
        n_components-dimensional embedding plus a security margin.

        Parameters
        ----------
        X
            The data set containing the samples from where the landmarks
            should be inferred.

        Returns
        -------
        landmarks: array-like, shape (n_landmarks,)
        """
        if self.landmarks_ is None:
            n_samples = X.shape[0]
            n_landmarks = self.n_landmarks

            if n_landmarks is None:
                self.landmarks_ = np.arange(X.shape[0])
            else:
                if n_landmarks == 'auto':
                    # To guarantee some stability, the number of landmarks
                    # should be strictly superior to the number of dimensions
                    # of the reduced data set added by some safety margin.
                    n_landmarks = (self.n_components + 1 +
                                   min(self.n_components // 2, 100))

                    # Limit number of landmarks by the number of samples.
                    n_landmarks = min(n_landmarks, n_samples)

                self.landmarks_ = np.arange(n_samples)
                check_random_state(self.random_state).shuffle(self.landmarks_)
                self.landmarks_ = self.landmarks_[:n_landmarks]

        return self.landmarks_
>>>>>>> Merge Isomap and LandmarkIsomap classes.

    def _fit_transform(self, X):
        X = check_array(X)
        landmarks = self._compute_landmarks(X)
        self.nbrs_ = NearestNeighbors(n_neighbors=self.n_neighbors,
                                      algorithm=self.neighbors_algorithm,
                                      n_jobs=self.n_jobs)
        self.nbrs_.fit(X)
        self.training_data_ = self.nbrs_._fit_X[landmarks, :]
        self.kernel_pca_ = KernelPCA(n_components=self.n_components,
                                     kernel="precomputed",
                                     eigen_solver=self.eigen_solver,
                                     tol=self.tol, max_iter=self.max_iter,
                                     n_jobs=self.n_jobs)

        kng = kneighbors_graph(self.nbrs_, self.n_neighbors,
                               mode='distance', n_jobs=self.n_jobs)

        self.dist_matrix_ = graph_shortest_path(kng,
                                                method=self.path_method,
                                                directed=False,
                                                only_vertices=landmarks).T

        G = self.dist_matrix_[landmarks] ** 2
        G *= -.5

        # Selectively replaces embedding_ rows.
        # This preserves the order of the samples in X.
        self.embedding_ = np.zeros((X.shape[0], self.n_components))
        self.embedding_[landmarks] = self.kernel_pca_.fit_transform(G)

        if landmarks.shape[0] < X.shape[0]:
            # Embed the samples that are not landmarks.
            others = np.ones(X.shape[0], dtype=bool)
            others[landmarks] = 0
            self.embedding_[others] = self.transform(X[others])

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
        G = -0.5 * self.dist_matrix_[self.landmarks_] ** 2
        G_center = KernelCenterer().fit_transform(G)
        evals = self.kernel_pca_.lambdas_
        return np.sqrt(np.sum(G_center ** 2) - np.sum(evals ** 2)) / G.shape[0]

    def fit(self, X, y=None):
        """Compute the embedding vectors for data X

        Parameters
        ----------
        X : {array-like, sparse matrix, BallTree, KDTree, NearestNeighbors}
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
        X : {array-like, sparse matrix, BallTree, KDTree}
            Training vector, where n_samples in the number of samples
            and n_features is the number of features.

        Returns
        -------
        X_new : array-like, shape (n_samples, n_components)
        """
        self._fit_transform(X)
        return self.embedding_

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
        # via the nearest neighbors of X.
        # This can be done as a single array operation, but it potentially
        # takes a lot of memory.  To avoid that, use a loop:
        G_X = np.zeros((X.shape[0], self.training_data_.shape[0]))
        for i in range(X.shape[0]):
            G_X[i] = np.min((self.dist_matrix_[indices[i]] +
                             distances[i][:, None]), axis=0)

        G_X **= 2
        G_X *= -0.5

        return self.kernel_pca_.transform(G_X)
