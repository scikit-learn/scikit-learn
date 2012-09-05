# -*- coding: utf-8 -*-
"""
Nearest Centroid Classification
"""

# Author: Robert Layton <robertlayton@gmail.com>
#         Olivier Grisel <olivier.grisel@ensta.org>
#
# License: BSD Style.

import numpy as np
from scipy import sparse as sp

from ..base import BaseEstimator, ClassifierMixin
from ..utils.validation import check_arrays, atleast2d_or_csr
from ..metrics.pairwise import pairwise_distances


class NearestCentroid(BaseEstimator, ClassifierMixin):
    """Nearest centroid classifier.

    Each class is represented by its centroid, with test samples classified to
    the class with the nearest centroid.

    Parameters
    ----------
    metric: string, or callable
        The metric to use when calculating distance between instances in a
        feature array. If metric is a string or callable, it must be one of
        the options allowed by metrics.pairwise.pairwise_distances for its
        metric parameter.
    shrink_threshold : float, optional (default = None)
        Threshold for shrinking centroids to remove features.

    Attributes
    ----------
    `centroids_` : array-like, shape = [n_classes, n_features]
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
    >>> print clf.predict([[-0.8, -1]])
    [1]

    See also
    --------
    sklearn.neighbors.KNeighborsClassifier: nearest neighbors classifier

    Notes
    -----
    When used for text classification with tf–idf vectors, this classifier is
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
            Training vector, where n_samples in the number of samples and
            n_features is the number of features.
            Note that centroid shrinking cannot be used with sparse matrices.
        y : array, shape = [n_samples]
            Target values (integers)
        """
        X, y = check_arrays(X, y, sparse_format="csr")
        if sp.issparse(X) and self.shrink_threshold:
            raise ValueError("threshold shrinking not supported"
                             " for sparse input")

        n_samples, n_features = X.shape
        classes = np.unique(y)
        self.classes_ = classes
        n_classes = classes.size
        if n_classes < 2:
            raise ValueError('y has less than 2 classes')

        # Mask mapping each class to it's members.
        self.centroids_ = np.empty((n_classes, n_features), dtype=np.float64)
        for i, cur_class in enumerate(classes):
            center_mask = y == cur_class
            if sp.issparse(X):
                center_mask = np.where(center_mask)[0]
            self.centroids_[i] = X[center_mask].mean(axis=0)

        if self.shrink_threshold:
            dataset_centroid_ = np.array(X.mean(axis=0))[0]
            # Number of clusters in each class.
            nk = np.array([np.sum(classes == cur_class)
                           for cur_class in classes])
            # m parameter for determining deviation
            m = np.sqrt((1. / nk) + (1. / n_samples))
            # Calculate deviation using the standard deviation of centroids.
            variance = np.array(np.power(X - self.centroids_[y], 2))
            variance = variance.sum(axis=0)
            s = np.sqrt(variance / (n_samples - n_classes))
            s += np.median(s)  # To deter outliers from affecting the results.
            mm = m.reshape(len(m), 1)  # Reshape to allow broadcasting.
            ms = mm * s
            deviation = ((self.centroids_ - dataset_centroid_) / ms)
            # Soft thresholding: if the deviation crosses 0 during shrinking,
            # it becomes zero.
            signs = np.sign(deviation)
            deviation = (np.abs(deviation) - self.shrink_threshold)
            deviation[deviation < 0] = 0
            deviation = np.multiply(deviation, signs)
            # Now adjust the centroids using the deviation
            msd = np.multiply(ms, deviation)
            self.centroids_ = np.array([dataset_centroid_ + msd[i]
                                        for i in xrange(n_classes)])
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
        X = atleast2d_or_csr(X)
        if not hasattr(self, "centroids_"):
            raise AttributeError("Model has not been trained yet.")
        return self.classes_[pairwise_distances(
            X, self.centroids_, metric=self.metric).argmin(axis=1)]
