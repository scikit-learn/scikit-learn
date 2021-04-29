import warnings

import numpy as np
from sklearn.utils.validation import check_X_y, check_array, check_is_fitted, _check_sample_weight, check_random_state
from ..utils.extmath import row_norms
from ._kmeans import KMeans, _kmeans_single_elkan, _kmeans_single_lloyd, _labels_inertia_threadpool_limit


class BisectKMeans(KMeans):
    def __init__(self,  n_clusters: int = 8, n_init: int = 10, random_state=None,
                 max_iter: int = 300, verbose=0, tol=1e4):

        super().__init__(
            n_clusters=n_clusters, max_iter=max_iter, verbose=verbose,
            random_state=random_state, tol=tol, n_init=n_init)

    def _calc_squared_errors(self, X, centers, labels):
        """
        Calculates SSE of each point and assign them by label
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
            warnings.warn("BisectKMeans is inefficient for n_cluster smaller than "
                          " - Use Normal KMeans from sklearn.cluster instead")

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
        X = self._check_test_data(X)

        self._check_params(X)
        self.random_state = check_random_state(self.random_state)
        sample_weight = _check_sample_weight(sample_weight, X, dtype=X.dtype)

        n_iter = 0

        data_left = X
        centroids = []

        while n_iter < self.max_iter:
            child = KMeans(
                n_clusters=2, n_init=self.n_init, max_iter=self.max_iter,
                tol=self.tol, verbose=self.verbose, random_state=self.random_state,
                copy_x=self.copy_x, algorithm=self._algorithm).fit(data_left)

            sse = self._calc_squared_errors(X, child.cluster_centers_, child.labels_)

            lower_sse_index = 0 if sse[0].sum(axis=0) < sse[1].sum(axis=0) else 1
            higher_sse_index = int(not lower_sse_index)

            centroids.append(child.cluster_centers_[lower_sse_index])

            if len(centroids) + 1 == self.n_clusters:
                centroids.append(child.cluster_centers_[higher_sse_index])
                break

            data_left = data_left[[x for x in range(len(child.labels_)) if child.labels_[x] == higher_sse_index]]

            n_iter += 1

        if n_iter == self.max_iter:
            raise RuntimeError("Number of iteration exceeded maximum allowed")

        x_squared_norms = row_norms(X, squared=True)

        _centers = np.asarray(centroids)

        if self._algorithm == "full":
            self._kmeans_single = _kmeans_single_lloyd
            self._check_mkl_vcomp(X, X.shape[0])
        else:
            self._kmeans_single = _kmeans_single_elkan

        self.labels_, self.inertia_, self.cluster_centers_, self.n_iter_ = self._kmeans_single(
            X, sample_weight, _centers, self.max_iter,
            self.verbose, x_squared_norms, self.tol, self._n_threads
        )

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
