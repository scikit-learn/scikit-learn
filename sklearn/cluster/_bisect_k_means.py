# Author: Michal Krawczyk <mkrwczyk.1@gmail.com>

import warnings

import numpy as np
import scipy.sparse as sp
from threadpoolctl import threadpool_limits

from ..exceptions import ConvergenceWarning
from ..exceptions import EfficiencyWarning

from ._kmeans import KMeans
from ._kmeans import _kmeans_single_elkan
from ._kmeans import _kmeans_single_lloyd

from ._k_means_common import _inertia_dense
from ._k_means_common import _inertia_sparse

from ._k_means_lloyd import lloyd_iter_chunked_dense
from ._k_means_lloyd import lloyd_iter_chunked_sparse

from ..utils.extmath import row_norms
from ..utils._openmp_helpers import _openmp_effective_n_threads

from ..utils.validation import check_array
from ..utils.validation import _check_sample_weight
from ..utils.validation import check_random_state


def _check_labels(X, sample_weight, x_squared_norms, centers, n_threads=1):
    """ Compute the labels of the given samples and centers.

    Parameters
    ----------
    X : {ndarray, sparse matrix} of shape (n_samples, n_features)
        The input samples to assign to the labels. If sparse matrix, must
        be in CSR format.

    sample_weight : ndarray of shape (n_samples,)
        The weights for each observation in X.

    x_squared_norms : ndarray of shape (n_samples,)
        Precomputed squared euclidean norm of each data point, to speed up
        computations.

    centers : ndarray of shape (n_clusters, n_features)
        The cluster centers.

    n_threads : int, default=1
        The number of OpenMP threads to use for the computation. Parallelism is
        sample-wise on the main cython loop which assigns each sample to its
        closest center.

    Returns
    -------
    labels : ndarray of shape (n_samples,)
        The resulting assignment (Labels of each point).
    """
    n_samples = X.shape[0]
    n_clusters = centers.shape[0]

    labels = np.full(n_samples, -1, dtype=np.int32)
    weight_in_clusters = np.zeros(n_clusters, dtype=centers.dtype)
    center_shift = np.zeros_like(weight_in_clusters)

    if sp.issparse(X):
        _labels = lloyd_iter_chunked_sparse
    else:
        _labels = lloyd_iter_chunked_dense

    _labels(X, sample_weight, x_squared_norms, centers, centers,
            weight_in_clusters, labels, center_shift, n_threads,
            update_centers=False)

    return labels


def _check_labels_threadpool_limit(X, sample_weight, x_squared_norms,
                                   centers, n_threads=1):
    """Same as _check_labels but in a threadpool_limits context."""
    with threadpool_limits(limits=1, user_api="blas"):
        labels = _check_labels(X, sample_weight, x_squared_norms,
                               centers, n_threads)

    return labels


class BisectKMeans(KMeans):
    """ Bisecting K-Means clustering
    K-Means variant that splits consecutively data with two centroids.
    Centroids with lower SSE (inertia) are kept as new cluster centers.
    Centroids with higher SSE are further split until the desired
    number of cluster is reached.

    That algorithm can produce partitional/hierarchical clustering and
    should be able to recognize clusters of any shape and size.

    That approach is also preferable to agglomerative clustering
    if the number of clusters is small, compared to the number of data points.

    Parameters
    ----------
    n_clusters : int, default=8
        The number of clusters to form as well as the number of
        centroids to generate.

    init : {'k-means++', 'random'} or callable, default='k-means++'
        Method for initialization:

        'k-means++' : selects initial cluster centers for k-mean
        clustering in a smart way to speed up convergence. See section
        Notes in k_init for more details.

        'random': choose `n_clusters` observations (rows) at random from data
        for the initial centroids.

        If a callable is passed, it should take arguments X, n_clusters and a
        random state and return an initialization.

    n_init : int, default=10
        Number of time the k-means algorithm will be run with different
        centroid seeds in each bisection.
        That will result producing for each bisection best output of n_init
        consecutive runs in terms of inertia.

    random_state : int, RandomState instance or None, default=None
        Determines random number generation for centroid initialization. Use
        an int to make the randomness deterministic.
        See :term:`Glossary <random_state>`.

    max_iter : int, default=30
        Maximum number of iterations of the inner k-means algorithm at each
        bisection.

    verbose : int, default=0
        Verbosity mode.

    tol : float, default=1e-4
        Relative tolerance with regards to Frobenius norm of the difference
        in the cluster centers of two consecutive iterations  to declare
        convergence. Used in inner k-means algorithm at each bisection to pick
        best possible clusters.

    copy_x : bool, default=True
        When pre-computing distances it is more numerically accurate to center
        the data first. If copy_x is True (default), then the original data is
        not modified. If False, the original data is modified, and put back
        before the function returns, but small numerical differences may be
        introduced by subtracting and then adding the data mean. Note that if
        the original data is not C-contiguous, a copy will be made even if
        copy_x is False. If the original data is sparse, but not in CSR format,
        a copy will be made even if copy_x is False.

    algorithm : {"auto", "full", "elkan"}, default="auto"
        K-means algorithm to use. The classical EM-style algorithm is "full".
        The "elkan" variation is more efficient on data with well-defined
        clusters, by using the triangle inequality. However it's more memory
        intensive due to the allocation of an extra array of shape
        (n_samples, n_clusters).
        For now "auto" (kept for backward compatibiliy) chooses "elkan" but it
        might change in the future for a better heuristic.

    bisect_strategy : {"biggest_sse", "largest_cluster"},
        default="biggest_sse"
        Defines how should bisection by performed:
        - "biggest_sse" means that Bisect K-Means will always check
        all calculated cluster for cluster with biggest SSE
        (Sum of squared errors) and bisect it. That way calculated clusters
        will be more balanced.
        - "largest_cluster" - Bisect K-Means will always split cluster with
        largest amount of points assigned to it from all clusters
        previously calculated. That should work faster than picking by SSE
        ('biggest_sse') and may produce similar results in most cases.


    Attributes
    ----------
    cluster_centers_ : ndarray of shape (n_clusters, n_features)
        Coordinates of cluster centers. If the algorithm stops before fully
        converging (see ``tol`` and ``max_iter``), these will not be
        consistent with ``labels_``.

    labels_ : ndarray of shape (n_samples,)
        Labels of each point.

    inertia_ : float
        Sum of squared distances of samples to their closest cluster center,
        weighted by the sample weights if provided.

    n_iter_ : int
        Number of iterations run.

    Notes
    -----
    Bisection cannot be performed if n_cluster < 2.

    Also it might be inefficient when n_cluster is equal to 2

    Examples
    --------
    >>> from sklearn.cluster import BisectKMeans
    >>> import numpy as np
    >>> X = np.array([[1, 2], [1, 4], [1, 0],
    ...               [10, 2], [10, 4], [10, 0],
    ...               [10, 6], [10, 8], [10, 10]])
    >>> bisect_means = BisectKMeans(n_clusters=3, random_state=0).fit(X)
    >>> bisect_means.labels_
    array([0, 0, 0, 1, 1, 1, 2, 2, 2], dtype=int32)
    >>> bisect_means.predict([[0, 0], [12, 3]])
    array([0, 1], dtype=int32)
    >>> bisect_means.cluster_centers_
    array([[ 1.,  2.],
           [10.,  2.],
           [10.,  8.]])
    """
    def __init__(self,  n_clusters=8, init='k-means++', n_init=10,
                 random_state=None, max_iter=30, verbose=0,
                 tol=1e-4, copy_x=True, algorithm='auto',
                 bisect_strategy='biggest_sse'):

        super().__init__(
            n_clusters=n_clusters, init=init, max_iter=max_iter,
            verbose=verbose, random_state=random_state, tol=tol,
            n_init=n_init, copy_x=copy_x, algorithm=algorithm)

        self.bisect_strategy = bisect_strategy

    def _compute_bisect_errors(self, X, centers, labels, sample_weight):
        """
        Calculate the squared error of each sample and group them by label.

        Parameters
        ----------
        X : {ndarray, sparse matrix} of shape (n_samples, n_features)
            The input samples.

        centers : ndarray of shape (n_clusters, n_features)
            The cluster centers.

        labels : ndarray of shape (n_samples,)
            Index of the cluster each sample belongs to.

        sample_weight : array-like of shape (n_samples,), default=None
            The weights for each observation in X.

        Returns
        -------
        errors_by_label : dict
            dictionary containing squared error of each point by label
            as ndarray.
        """
        errors_by_label = {}

        _inertia = _inertia_sparse if sp.issparse(X) else _inertia_dense

        for value in range(centers.shape[0]):
            indexes = (labels == value)

            data = X[indexes]
            weights = sample_weight[indexes]
            center = centers[value][np.newaxis, :]
            label = np.zeros(data.shape[0], dtype=np.intc)

            errors_by_label[value] = _inertia(data, weights, center,
                                              label, self._n_threads)

        return errors_by_label

    def _check_params(self, X):
        super()._check_params(X)

        # bisect_strategy
        if self.bisect_strategy not in \
                ["biggest_sse", "largest_cluster"]:
            raise ValueError(f"Bisect Strategy must be 'biggest_sse', "
                             f"or 'largest_cluster' "
                             f"got {self.bisect_strategy} instead.")

        # Regular K-Means should do less computations when there are only
        # less than 3 clusters
        if self.n_clusters < 3:
            warnings.warn("BisectKMeans might be inefficient for n_cluster "
                          "smaller than 3  "
                          "- Use Normal KMeans from sklearn.cluster instead.",
                          EfficiencyWarning)

        if X.shape[0] <= 1:
            raise ValueError("BisectKMeans needs more than one sample "
                             "to perform bisection.")

        if hasattr(self.init, '__array__'):
            raise ValueError("BisectKMeans does not support "
                             "init as array.")

    def _bisect(self, X, sample_weight, random_state):
        """ Bisection of data
        Attempts to get best bisection of data by performing regular K-Means
        for different pairs of centroids

        .. note:: Number of attempts is specified by self.n_init

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)

            Training instances to cluster.

            .. note:: The data will be converted to C ordering,
                which will cause a memory copy
                if the given data is not C-contiguous.

        sample_weight : array-like of shape (n_samples,)
            The weights for each observation in X. If None, all observations
            are assigned equal weight.

        random_state : int, RandomState instance
            Determines random number generation for centroid initialization.

        Returns
        -------
        self
            Fitted estimator.
        """
        init = self.init
        x_squared_norms = row_norms(X, squared=True)

        best_inertia = None

        for i in range(self.n_init):
            centers_init = self._init_centroids(X, x_squared_norms, init,
                                                random_state, n_centroids=2)

            labels, inertia, centers, _ = self._kmeans_single(
                X, sample_weight, centers_init, max_iter=self.max_iter,
                verbose=self.verbose, tol=self.tol,
                x_squared_norms=x_squared_norms, n_threads=self._n_threads
            )

            if best_inertia is None or inertia < best_inertia:
                best_labels = labels
                best_centers = centers
                best_inertia = inertia

        distinct_clusters = len(set(best_labels))
        if distinct_clusters != 2:
            warnings.warn(
                "Number of distinct clusters ({}) found smaller than "
                "n_clusters ({}). Possibly due to duplicate points "
                "in X.".format(distinct_clusters, 2),
                ConvergenceWarning, stacklevel=2)

        return best_centers, best_labels

    def fit(self, X, y=None, sample_weight=None):
        """Compute bisecting k-means clustering.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)

            Training instances to cluster.

            .. note:: The data will be converted to C ordering,
                which will cause a memory copy
                if the given data is not C-contiguous.

        y : Ignored
            Not used, present here for API consistency by convention.

        sample_weight : array-like of shape (n_samples,), default=None
            The weights for each observation in X. If None, all observations
            are assigned equal weight.

        Returns
        -------
        self
            Fitted estimator.
        """
        X = self._validate_data(X, accept_sparse='csr',
                                dtype=[np.float64, np.float32],
                                order='C', copy=self.copy_x,
                                accept_large_sparse=False)

        self._check_params(X)
        random_state = check_random_state(self.random_state)
        sample_weight = _check_sample_weight(sample_weight, X, dtype=X.dtype)
        self._n_threads = _openmp_effective_n_threads()

        # Validate init array
        init = self.init

        if hasattr(init, '__array__'):
            init = check_array(init, dtype=X.dtype, copy=True, order='C')
            self._validate_center_shape(X, init)

        if self._algorithm == "full":
            self._kmeans_single = _kmeans_single_lloyd
            self._check_mkl_vcomp(X, X.shape[0])
        else:
            self._kmeans_single = _kmeans_single_elkan

        _inertia = _inertia_sparse if sp.issparse(X) else _inertia_dense

        if self.verbose:
            print("Running Bisecting K-Means with parameters:")
            print(f"-> number of clusters: {self.n_clusters}")
            print(f"-> number of centroid initializations: {self.n_init}")
            print("-> relative tolerance: {:.4e}".format(self.tol))
            print(f"-> bisect strategy: {self.bisect_strategy} \n")

        # Subtract of mean of X for more accurate distance computations
        if not sp.issparse(X):
            X_mean = X.mean(axis=0)
            X -= X_mean

        # Only assign to created centroid when n_clusters == 1
        if self.n_clusters == 1:
            x_squared_norms = row_norms(X, squared=True)

            clusters = self._init_centroids(X, x_squared_norms,
                                            init, random_state,
                                            n_centroids=1)
            warnings.warn("Bisection won't be performed - "
                          "needs at least two clusters to run.")

            self.labels_ = np.zeros(X.shape[0], dtype=np.intc)

        else:
            # Run proper bisection to gather
            # self.cluster_centers_ and self.labels_
            clusters, self.labels_ = self._run_bisect_kmeans(X, sample_weight,
                                                             random_state)

        # Restore original data
        if not sp.issparse(X):
            X += X_mean
            clusters += X_mean

        self.cluster_centers_ = clusters

        self.inertia_ = _inertia(X, sample_weight, self.cluster_centers_,
                                 self.labels_, self._n_threads)

        # number of iterations will always be equal to
        # (number of clusters - 1)
        self.n_iter_ = self.n_clusters - 1

        return self

    def _run_bisect_kmeans(self, X, sample_weight, random_state):
        """ Performs Bisecting K-Means, which splits always cluster depending
        on 'bisect_strategy' attribute:

        - "biggest sse": Picks cluster with biggest SSE (Sum of Squared Errors)
         from all calculated.

        - "largest_cluster": Picks always cluster with largest number of
         points assigned from all calculated. That method will perform faster
          than picking by SSE methods, while producing similar results.

        .. note:: All of passed parameters must be pre-calculated

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            Training instances to cluster.

        sample_weight : array-like of shape (n_samples,)
            The weights for each observation in X.

        random_state : int, RandomState instance
            Determines random number generation for centroid initialization.

        Returns
        ----------
        cluster_centers_ : ndarray of shape (n_clusters, n_features)
            Coordinates of cluster centers.

        labels_ : ndarray of shape (n_samples,)
            Labels of each point.
        """
        strategy = self.bisect_strategy
        x_squared_norms = row_norms(X, squared=True)

        # Dictionary storing calculated centroids
        centers_dict = {
            0: None,
        }
        # Boolean mask for picking data to bisect
        picked_labels = np.ones(X.shape[0], dtype=bool)

        # ID of biggest center stored in centers_dict
        biggest_id = 0

        last_center_id = 0

        for n_iter in range(self.n_clusters - 1):
            # Pick data and weights to bisect
            picked_data = X[picked_labels]
            picked_weights = sample_weight[picked_labels]

            # Perform Bisection
            _centers, _ = self._bisect(picked_data, picked_weights,
                                       random_state)

            if self.verbose:
                print(f"Centroid Found: {_centers[0]}")
                print(f"Centroid Found: {_centers[1]}")

            # Save obtained centers
            centers_dict[last_center_id + 1] = _centers[0]
            centers_dict[last_center_id + 2] = _centers[1]

            # Delete cluster that was split from dict
            del centers_dict[biggest_id]

            # Convert centers to numpy array for further calculations
            centers = np.asarray(list(centers_dict.values()))

            # Recalculate labels for all centers
            # That is made to be sure that all data is split to proper clusters
            labels = _check_labels_threadpool_limit(X, sample_weight,
                                                    x_squared_norms,
                                                    centers,
                                                    self._n_threads)

            if strategy == "biggest_sse":
                # Compute SSE (Sum of Squared Errors) for each cluster
                errors = self._compute_bisect_errors(X, centers, labels,
                                                     sample_weight)

                # Pick cluster with biggest SSE
                biggest_label, _ = max(errors.items(), key=lambda x: x[1])
            else:
                # "largest_cluster"
                # Count occurrences of each label
                labels_occurrences = np.bincount(labels)
                # Pick cluster with largest number of data points assigned
                biggest_label = np.argmax(labels_occurrences)

            # Pick indexes of data for further split
            picked_labels = (labels == biggest_label)

            # Pick id of cluster for further split
            biggest_id = list(centers_dict.keys())[biggest_label]

            last_center_id += 2

        return centers, labels
