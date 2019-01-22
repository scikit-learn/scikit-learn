# -*- coding: utf-8 -*-
"""
Nearest Centroid Classification
"""

# Author: Robert Layton <robertlayton@gmail.com>
#         Olivier Grisel <olivier.grisel@ensta.org>
#
# License: BSD 3 clause

import warnings
import numpy as np
from scipy import sparse as sp

from ..base import BaseEstimator, ClassifierMixin
from ..metrics.pairwise import pairwise_distances
from ..preprocessing import LabelEncoder
from ..utils.validation import check_array, check_X_y, check_is_fitted
from ..utils.sparsefuncs import csc_median_axis_0
from ..utils.multiclass import check_classification_targets, _check_partial_fit_first_call
from copy import deepcopy


class NearestCentroid(BaseEstimator, ClassifierMixin):
    """Nearest centroid classifier.

    Each class is represented by its centroid, with test samples classified to
    the class with the nearest centroid.

    Read more in the :ref:` >`.

    Parameters
    ----------
    metric : string, or callable
        The metric to use when calculating distance between instances in a
        feature array. If metric is a string or callable, it must be one of
        the options allowed by metrics.pairwise.pairwise_distances for its
        metric parameter.
        The centroids for the samples corresponding to each class is the point
        from which the sum of the distances (according to the metric) of all
        samples that belong to that particular class are minimized.
        If the "manhattan" metric is provided, this centroid is the median and
        for all other metrics, the centroid is now set to be the mean.

    shrink_threshold : float, optional (default = None)
        Threshold for shrinking centroids to remove features.

    Attributes
    ----------
    centroids_ : array-like, shape = [n_classes, n_features]
        Centroid of each class

    Examples
    --------
    >>> from sklearn.neighbors.nearest_centroid import NearestCentroid
    >>> import numpy as np
    >>> X = np.array([[-1, -1], [-2, -1], [-3, -2], [1, 1], [2, 1], [3, 2]])
    >>> y = np.array([1, 1, 1, 2, 2, 2])
    >>> clf = NearestCentroid()
    >>> clf.fit(X, y)
    NearestCentroid(metric='euclidean', shrink_threshold=None)
    >>> print(clf.predict([[-0.8, -1]]))
    [1]

    See also
    --------
    sklearn.neighbors.KNeighborsClassifier: nearest neighbors classifier

    Notes
    -----
    When used for text classification with tf-idf vectors, this classifier is
    also known as the Rocchio classifier.

    References
    ----------
    Tibshirani, R., Hastie, T., Narasimhan, B., & Chu, G. (2002). Diagnosis of
    multiple cancer types by shrunken centroids of gene expression. Proceedings
    of the National Academy of Sciences of the United States of America,
    99(10), 6567-6572. The National Academy of Sciences.

    """

    def __init__(self, metric='euclidean', shrink_threshold=None):
        self.metric = metric
        self.shrink_threshold = shrink_threshold

    def fit(self, X, y):
        """
        Fit the NearestCentroid model according to the given training data.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape = [n_samples, n_features]
            Training vector, where n_samples is the number of samples and
            n_features is the number of features.
            Note that centroid shrinking cannot be used with sparse matrices.

        y : array, shape = [n_samples]
            Target values (integers)
        """
        return self._partial_fit(X, y, np.unique(y), _refit=True)

    def partial_fit(self, X, y, classes=None):
        """Incremental fit on a batch of samples.

        This method is expected to be called several times consecutively
        on different chunks of a dataset so as to implement out-of-core
        or online learning.

        This is especially useful when the whole dataset is too big to fit in
        memory at once.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape = [n_samples, n_features]
            Training vector, where n_samples is the number of samples and
            n_features is the number of features.
            Note that centroid shrinking cannot be used with sparse matrices.

        y : array, shape = [n_samples]
            Target values (integers)

        classes : array-like, shape (n_classes,), optional (default=None)
            List of all the classes that can possibly appear in the y vector.

            Must be provided at the first call to partial_fit, can be omitted
            in subsequent calls.
        """
        if self.metric == 'manhattan':
            raise ValueError("Partial fitting with manhattan is not supported.")
        return self._partial_fit(X, y, classes, _refit=False)

    def _partial_fit(self, X, y, classes=None, _refit=False):
        """
        Actual implementation of the Nearest Centroid fitting.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape = [n_samples, n_features]
            Training vector, where n_samples is the number of samples and
            n_features is the number of features.
            Note that centroid shrinking cannot be used with sparse matrices.

        y : array, shape = [n_samples]
            Target values (integers)

        classes : array-like, shape (n_classes,), optional (default=None)
            List of all the classes that can possibly appear in the y vector.

            Must be provided at the first call to partial_fit, can be omitted
            in subsequent calls.

        _refit : bool, optional (default=False)
            If true, act as though this were the first time we called
            _partial_fit (ie, throw away any past fitting and start over).
        """
        if self.metric == 'precomputed':
            raise ValueError("Precomputed is not supported.")
        # If X is sparse and the metric is "manhattan", store it in a csc
        # format is easier to calculate the median.
        if self.metric == 'manhattan':
            X, y = check_X_y(X, y, ['csc'])
        else:
            X, y = check_X_y(X, y, ['csr', 'csc'])
        is_X_sparse = sp.issparse(X)
        if is_X_sparse and self.shrink_threshold:
            raise ValueError("threshold shrinking not supported"
                             " for sparse input")
        check_classification_targets(y)

        n_samples, n_features = X.shape

        if _refit or _check_partial_fit_first_call(self, classes):
            self.true_classes_ = classes = np.asarray(classes)
            # Mask mapping each class to its members.
            self.true_centroids_ = np.zeros((classes.size, n_features), dtype=np.float64)
            # Number of clusters in each class.
            self.nk_ = np.zeros(classes.size)

            if self.shrink_threshold:
                self.ssd_ = np.zeros((classes.size, n_features), dtype=np.float64)
                self.dataset_centroid_ = np.mean(X, axis=0)

        le = LabelEncoder()
        le.fit(self.true_classes_)
        y_ind = le.transform(y)
        n_classes = self.true_classes_.size
        if n_classes < 2:
            raise ValueError('The number of classes has to be greater than'
                             ' one; got %d class' % (n_classes))

        old_nk = deepcopy(self.nk_)
        old_centroids = deepcopy(self.true_centroids_)

        for cur_class in range(n_classes):
            center_mask = y_ind == cur_class

            # Ignore if no data for this class
            if X[center_mask].size == 0:
                continue
            if is_X_sparse:
                center_mask = np.where(center_mask)[0]

            # XXX: Update other averaging methods according to the metrics.
            if self.metric == "manhattan":
                self.nk_[cur_class] += np.sum(center_mask)
                # NumPy does not calculate median of sparse matrices.
                if not is_X_sparse:
                    self.true_centroids_[cur_class] = np.median(X[center_mask], axis=0)
                else:
                    self.true_centroids_[cur_class] = csc_median_axis_0(X[center_mask])
            else:
                if self.metric != 'euclidean':
                    warnings.warn("Averaging for metrics other than "
                                  "euclidean and manhattan not supported. "
                                  "The average is set to be the mean."
                                  )
                # Update each centroid weighted by the number of samples
                self.true_centroids_[cur_class] = X[center_mask].mean(axis=0) * np.sum(center_mask) +\
                                                  self.true_centroids_[cur_class] * self.nk_[cur_class]
                self.nk_[cur_class] += np.sum(center_mask)
                self.true_centroids_[cur_class] /= self.nk_[cur_class]

        # Filtering out centroids without any data
        self.classes_ = self.true_classes_[self.nk_ != 0]
        self.centroids_ = self.true_centroids_[self.nk_ != 0]

        if self.shrink_threshold:
            n_total = np.sum(self.nk_)
            self.dataset_centroid_ = (self.dataset_centroid_ * old_nk.sum(axis=0) + np.sum(X, axis=0)) / n_total

            # Update sum of square distances of each class
            for cur_class in range(n_classes):
                n_old = old_nk[cur_class]
                n_new = self.nk_[cur_class] - n_old
                if n_new == 0:
                    continue
                center_mask = y_ind == cur_class
                old_ssd = self.ssd_[cur_class]
                new_ssd = ((X[center_mask] - X[center_mask].mean(axis=0))**2).sum(axis=0)
                self.ssd_[cur_class] = (old_ssd + new_ssd +
                         (n_old / float(n_new * (n_new + n_old))) *
                         (n_new * old_centroids[cur_class] - n_new * X[center_mask].mean(axis=0)) ** 2)

            # m parameter for determining deviation
            m = np.sqrt((1. / self.nk_) - (1. / np.sum(self.nk_)))

            # Calculate deviation using the standard deviation of centroids.
            ssd = self.ssd_.sum(axis=0)
            s = np.sqrt(ssd / (n_total - n_classes))
            s += np.median(s)  # To deter outliers from affecting the results.
            mm = m.reshape(len(m), 1)  # Reshape to allow broadcasting.
            ms = mm * s
            deviation = ((self.true_centroids_ - self.dataset_centroid_) / ms)

            # Soft thresholding: if the deviation crosses 0 during shrinking,
            # it becomes zero.
            signs = np.sign(deviation)
            deviation = (np.abs(deviation) - self.shrink_threshold)
            np.clip(deviation, 0, None, out=deviation)
            deviation *= signs
            # Now adjust the centroids using the deviation
            msd = ms * deviation
            self.centroids_ = self.dataset_centroid_[np.newaxis, :] + msd
        return self

    def predict(self, X):
        """Perform classification on an array of test vectors X.

        The predicted class C for each sample in X is returned.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]

        Returns
        -------
        C : array, shape = [n_samples]

        Notes
        -----
        If the metric constructor parameter is "precomputed", X is assumed to
        be the distance matrix between the data to be predicted and
        ``self.centroids_``.
        """
        check_is_fitted(self, 'centroids_')

        X = check_array(X, accept_sparse='csr')
        return self.classes_[pairwise_distances(
            X, self.centroids_, metric=self.metric).argmin(axis=1)]