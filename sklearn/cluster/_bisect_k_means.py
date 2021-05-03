import warnings

import numpy as np

from sklearn.utils.validation import check_array, check_is_fitted, _check_sample_weight, check_random_state
from ..utils.extmath import row_norms
from ._kmeans import KMeans, _kmeans_single_elkan, _kmeans_single_lloyd, _labels_inertia_threadpool_limit

import scipy.sparse as sp
from sklearn.exceptions import ConvergenceWarning, EfficiencyWarning
from ._kmeans import _kmeans_plusplus


class BisectKMeans(KMeans):
    def __init__(self,  n_clusters: int = 8, n_init: int = 10, random_state=None,
                 max_iter: int = 300, verbose=0, tol=1e4):

        super().__init__(
            n_clusters=n_clusters, max_iter=max_iter, verbose=verbose,
            random_state=random_state, tol=tol, n_init=n_init)

    def _calc_squared_errors(self, X, centers, labels):
        """
        Calculates SSE of each point and group them by label
        .. note:: This

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
        sse - dictionary containing sse of each point by label
        """
        labels_by_index = [[], []]

        for x in range(0, len(labels)):
            if labels[x] == 0:
                labels_by_index[0].append(x)
            else:
                labels_by_index[1].append(x)

        sse = {}

        for index, labels in enumerate(labels_by_index):
            sse[index] = ((X[labels] - centers[index]) ** 2).sum(axis=1)

        return sse

    def _check_params(self, X):
        super()._check_params(X)

        if self.n_clusters < 3:
            warnings.warn("BisectKMeans is inefficient for n_cluster smaller than 3"
                          " - Use Normal KMeans from sklearn.cluster instead",
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
        centers : ndarray of shape (n_clusters, n_features)
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

    def _bisect(self, X, y=None, sample_weight=None):
        X = self._validate_data(X, accept_sparse='csr',
                                dtype=[np.float64, np.float32],
                                order='C', copy=self.copy_x,
                                accept_large_sparse=False)

        # self._check_params(X)
        # self.random_state = check_random_state(self.random_state)
        # sample_weight = _check_sample_weight(sample_weight, X, dtype=X.dtype)

        #Validate init array
        init = self.init

        if hasattr(init, '__array__'):
            init = check_array(init, dtype=X.dtype, copy=True, order='C')
            self._validate_center_shape(X, init)

        if not sp.issparse(X):
            X_mean = X.mean(axis=0)
            X -= X_mean

            if hasattr(init, '__array__'):
                init -= X_mean

        x_squared_norms = row_norms(X, squared=True)

        # if self._algorithm == "full":
        #     self._kmeans_single = _kmeans_single_lloyd
        #     self._check_mkl_vcomp(X, X.shape[0])
        # else:
        #     self._kmeans_single = _kmeans_single_elkan

        best_inertia = None

        if self.verbose:
            print("Initializing Centroids")

        for i in range(self.n_init):
            centers_init = self._init_two_centroids(X, x_squared_norms,
                                                    init, self.random_state)

            labels, inertia, centers, n_iter_ = self._kmeans_single(
                X, sample_weight, centers_init, max_iter=self.max_iter,
                verbose=self.verbose, tol=self.tol, x_squared_norms=x_squared_norms,
                n_threads=self._n_threads
            )

            if best_inertia is None or inertia < best_inertia:
                best_labels = labels
                best_centers = centers
                best_inertia = inertia
                best_n_iter = n_iter_

        if not sp.issparse(X):
            if not self.copy_x:
                X += X_mean
            best_centers += X_mean

        distinct_clusters = len(set(best_labels))
        if distinct_clusters < 2:
            warnings.warn(
                "Number of distinct clusters ({}) found smaller than "
                "n_clusters ({}). Possibly due to duplicate points "
                "in X.".format(distinct_clusters, 2),
                ConvergenceWarning, stacklevel=2)

        return best_centers, best_labels, best_inertia, best_n_iter

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
        self.random_state = check_random_state(self.random_state)
        sample_weight = _check_sample_weight(sample_weight, X, dtype=X.dtype)

        if self._algorithm == "full":
            self._kmeans_single = _kmeans_single_lloyd
            self._check_mkl_vcomp(X, X.shape[0])
        else:
            self._kmeans_single = _kmeans_single_elkan

        data_left = X
        weights_left = sample_weight

        centroids = []

        for n_iter in range(self.max_iter):

            centers, labels, _,  _ = self._bisect(data_left, y, weights_left)

            sse = self._calc_squared_errors(X, centers, labels)

            lower_sse_index = 0 if sse[0].sum(axis=0) < sse[1].sum(axis=0) else 1
            higher_sse_index = 1 if lower_sse_index == 0 else 0

            centroids.append(centers[lower_sse_index])

            if len(centroids) + 1 == self.n_clusters:
                centroids.append(centers[higher_sse_index])
                break

            indexes = [x for x in range(len(labels)) if labels[x] == higher_sse_index]
            data_left = data_left[indexes]
            weights_left = weights_left[indexes]

            n_iter += 1

        x_squared_norms = row_norms(X, squared=True)

        _centers = np.asarray(centroids)

        self.cluster_centers_, self.n_iter_ = _centers, n_iter

        self.labels_, self.inertia_ = _labels_inertia_threadpool_limit(
            X, sample_weight, x_squared_norms,
            _centers, n_threads=self._n_threads)

        return self

    def predict(self, X, sample_weight=None):
        """Predict the closest cluster each sample in X belongs to.
                In the vector quantization literature, `cluster_centers_` is called
                the code book and each value returned by `predict` is the index of
                the closest code in the code book.

                Parameters
                ----------
                X : {array-like, sparse matrix} of shape (n_samples, n_features)
                    New data to predict.

                sample_weight : array-like of shape (n_samples,), default=None
                    The weights for each observation in X. If None, all observations
                    are assigned equal weight.

                Returns
                -------
                labels : ndarray of shape (n_samples,)
                    Index of the cluster each sample belongs to.
                """
        check_is_fitted(self)

        X = self._check_test_data(X)
        x_squared_norms = row_norms(X, squared=True)
        sample_weight = _check_sample_weight(sample_weight, X, dtype=X.dtype)

        labels, _ = _labels_inertia_threadpool_limit(
            X, sample_weight, x_squared_norms,
            self.cluster_centers_, n_threads=self._n_threads)

        return labels