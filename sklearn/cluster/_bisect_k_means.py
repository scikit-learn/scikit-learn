# Author: Michal Krawczyk <mkrwczyk.1@gmail.com>

import warnings

import numpy as np
import scipy.sparse as sp

from ..exceptions import ConvergenceWarning
from ..exceptions import EfficiencyWarning

from ._kmeans import KMeans
from ._kmeans import _kmeans_single_elkan
from ._kmeans import _kmeans_single_lloyd
from ._kmeans import _labels_inertia_threadpool_limit
from ._kmeans import _kmeans_plusplus

from ..utils.extmath import row_norms
from ..utils._openmp_helpers import _openmp_effective_n_threads

from ..utils.validation import check_array
from ..utils.validation import _check_sample_weight
from ..utils.validation import check_random_state


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

    init : {'k-means++', 'random'}, callable or array-like of shape \
            (n_clusters, n_features), default='k-means++'
        Method for initialization:

        'k-means++' : selects initial cluster centers for k-mean
        clustering in a smart way to speed up convergence. See section
        Notes in k_init for more details.

        'random': choose `n_clusters` observations (rows) at random from data
        for the initial centroids.

        If an array is passed, it should be of shape (n_clusters, n_features)
        and gives the initial centers.

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

    max_iter : int, default=300
        Maximum number of iterations of the k-means algorithm for a
        single run.

    verbose : int, default=0
        Verbosity mode.

    tol : float, default=1e-4
        Relative tolerance with regards to Frobenius norm of the difference
        in the cluster centers of two consecutive iterations to declare
        convergence.

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

    Attributes
    ----------
    cluster_centers_ : ndarray of shape (n_clusters, n_features)
        Coordinates of cluster centers. If the algorithm stops before fully
        converging (see ``tol`` and ``max_iter``), these will not be
        consistent with ``labels_``.

    labels_ : ndarray of shape (n_samples,)
        Labels of each point

    inertia_ : float
        Sum of squared distances of samples to their closest cluster center,
        weighted by the sample weights if provided.

    n_iter_ : int
        Number of iterations run.

    Notes
    -----
    That algorithm will not work if n_cluster is smaller than 2.

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
    array([0, 0, 0, 2, 2, 2, 1, 1, 1], dtype=int32)
    >>> bisect_means.predict([[0, 0], [12, 3]])
    array([0, 2], dtype=int32)
    >>> bisect_means.cluster_centers_
    array([[ 1.,  2.],
           [10.,  8.],
           [10.,  2.]])
    """
    def __init__(self,  n_clusters=8, init='k-means++', n_init=10,
                 random_state=None, max_iter=30, verbose=0,
                 tol=1e-4, copy_x=True, algorithm='auto'):

        super().__init__(
            n_clusters=n_clusters, init=init, max_iter=max_iter,
            verbose=verbose, random_state=random_state, tol=tol,
            n_init=n_init, copy_x=copy_x, algorithm=algorithm)

    def _compute_bisect_errors(self, X, centers, labels):
        """
        Calculate the squared error of each sample and group them by label.

        .. note:: That function works only if there are two labels (0,1) and
        may be less efficient for sparse data.

        Parameters
        ----------
        X : {ndarray, sparse matrix} of shape (n_samples, n_features)
            The input samples.

        centers : ndarray of shape (n_clusters, n_features)
            The cluster centers.

        labels : ndarray of shape (n_samples,)
            Index of the cluster each sample belongs to.

        Returns
        -------
        errors_by_label : dict
            dictionary containing squared error of each point by label
            as ndarray.
        """
        errors_by_label = {}

        if sp.issparse(X):
            X = X.toarray()

        for value in range(2):
            errors_by_label[value] = ((X[np.where(labels == value)]
                                       - centers[value]) ** 2).sum(axis=1)

        return errors_by_label

    def _check_params(self, X):
        super()._check_params(X)

        if self.n_clusters < 3:
            warnings.warn("BisectKMeans might be inefficient for n_cluster "
                          "smaller than 3  "
                          "- Use Normal KMeans from sklearn.cluster instead",
                          EfficiencyWarning)

    def _init_two_centroids(self, X, x_squared_norms, init, random_state,
                            init_size=None):
        """Compute two centroids for bisecting

        Parameters
        ----------
        X : {ndarray, sparse matrix} of shape (n_samples, n_features)
            The input samples.

        x_squared_norms : ndarray of shape (n_samples,)
            Squared euclidean norm of each data point. Pass it if you have it
            at hands already to avoid it being recomputed here.

        init : {'k-means++', 'random'}, callable or ndarray of shape \
                (n_clusters, n_features)
            Method for initialization.

        random_state : RandomState instance
            Determines random number generation for centroid initialization.
            See :term:`Glossary <random_state>`.

        init_size : int, default=None
            Number of samples to randomly sample for speeding up the
            initialization (sometimes at the expense of accuracy).
        Returns
        -------
        centers : ndarray of shape (n_clusters = 2, n_features)
        """
        n_samples = X.shape[0]
        n_clusters = 2

        if init_size is not None and init_size < n_samples:
            init_indices = random_state.randint(0, n_samples, init_size)
            X = X[init_indices]
            x_squared_norms = x_squared_norms[init_indices]
            n_samples = X.shape[0]

        if isinstance(init, str) and init == 'k-means++':
            centers, _ = _kmeans_plusplus(X, n_clusters,
                                          x_squared_norms, random_state)
        elif isinstance(init, str) and init == 'random':
            seeds = random_state.permutation(n_samples)[:n_clusters]
            centers = X[seeds]
        elif hasattr(init, '__array__'):
            centers = init
        elif callable(init):
            centers = init(X, n_clusters, random_state=random_state)
            centers = check_array(centers,
                                  dtype=X.dtype, copy=False, order='C')
            self._validate_center_shape(X, centers)

        if sp.issparse(centers):
            centers = centers.toarray()

        return centers

    def _bisect(self, X, init, sample_weight=None,
                random_state=None):
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

        init : {'k-means++', 'random'}, callable or ndarray of shape \
                (n_clusters, n_features)
            Method for initialization.

        sample_weight : array-like of shape (n_samples,), default=None
            The weights for each observation in X. If None, all observations
            are assigned equal weight.

        Returns
        -------
        self
            Fitted estimator.
        """

        x_squared_norms = row_norms(X, squared=True)

        best_inertia = None

        for i in range(self.n_init):
            centers_init = self._init_two_centroids(X, x_squared_norms,
                                                    init, random_state)

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

        # Subtract of mean of X for more accurate distance computations
        if not sp.issparse(X):
            X_mean = X.mean(axis=0)
            X -= X_mean

            if hasattr(init, '__array__'):
                init -= X_mean

        if self._algorithm == "full":
            self._kmeans_single = _kmeans_single_lloyd
            self._check_mkl_vcomp(X, X.shape[0])
        else:
            self._kmeans_single = _kmeans_single_elkan

        data_left = X
        weights_left = sample_weight

        centroids = []

        if self.verbose:
            print("Running Bisecting K-Means with parameters:")
            print(f"-> number of clusters: {self.n_clusters}")
            print(f"-> number of centroid initializations: {self.n_init}")
            print("-> relative tolerance: {:.4e} \n".format(self.tol))

        init_org = init

        for n_iter in range(self.n_clusters):

            # If init array is provided -
            # Take only part of init that is dedicated for that part of bisect
            if hasattr(init, '__array__'):
                init = init_org[n_iter: n_iter + 2].copy()

            # Perform Bisection
            centers, labels = self._bisect(data_left, init, weights_left,
                                           random_state)

            # Check SSE (Sum of Squared Errors) of each of computed centroids.
            # SSE is calculated with distances between data points
            # and assigned centroids
            errors = self._compute_bisect_errors(X, centers, labels)

            if not sp.issparse(X):
                centers += X_mean

            lower_sse_index = 0 if \
                errors[0].sum(axis=0) < errors[1].sum(axis=0) else 1

            higher_sse_index = 1 if lower_sse_index == 0 else 0

            # Centroid with lower SSE is added to result centroids.
            centroids.append(centers[lower_sse_index])

            if self.verbose:
                print(f"Centroid Found: {centers[lower_sse_index]}")

            if len(centroids) + 1 == self.n_clusters:
                # Desired number of cluster centers is reached so centroid with
                # higher SSE would not be split further
                centroids.append(centers[higher_sse_index])

                if self.verbose:
                    print(f"Centroid Found: {centers[higher_sse_index]}")
                break

            # Data with label that has higher SSE will be further processed
            indexes_left = (labels == higher_sse_index)

            data_left = data_left[indexes_left]
            weights_left = weights_left[indexes_left]

        x_squared_norms = row_norms(X, squared=True)

        _centers = np.asarray(centroids)

        # Restore X
        if not sp.issparse(X):
            X += X_mean

        self.cluster_centers_, self.n_iter_ = _centers, n_iter + 1

        # Since all clusters are calculated - label each data to valid center
        self.labels_, self.inertia_ = _labels_inertia_threadpool_limit(
            X, sample_weight, x_squared_norms,
            _centers, n_threads=self._n_threads)

        return self
