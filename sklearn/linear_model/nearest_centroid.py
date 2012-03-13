"""
Nearest Centroid Classification
"""

# Author: Robert Layton <robertlayton@gmail.com>
#         Olivier Grisel <olivier.grisel@ensta.org>
#
# License: BSD Style.

import warnings

import numpy as np
import scipy.ndimage as ndimage

from ..base import BaseEstimator, ClassifierMixin
from ..utils.validation import check_arrays
from ..base import ClassifierMixin
from ..metrics.pairwise import pairwise_distances


class NearestCentroid(BaseEstimator, ClassifierMixin):
    """
    Nearest Centroid Classification

    Each class is represented by its centroid, with test samples classified to
    the class with the nearest centroid.

    Parameters
    ----------
    metric: string, or callable
        The metric to use when calculating distance between instances in a
        feature array. If metric is a string or callable, it must be one of
        the options allowed by metrics.pairwise.calculate_distance for its
        metric parameter.
        If metric is "precomputed", X is assumed to be a distance matrix and
        must be square.
    shrink_threshold : float, optional
        Threshold for shrinking centroids to remove features.

    Attributes
    ----------
    `centroids_` : array-like, shape = [n_classes, n_features]
        Centroid of each class

    Examples
    --------
    >>> from sklearn.linear_model.nearest_centroid import NearestCentroid
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
    sklearn.neighbourds.KNeighborsClassifier: Nearest Neighbours Classifier

    Reference
    ---------
    Diagnosis of multiple cancer types by shrunken centroids of gene expression
    TODO: better referencing!

    """

    def __init__(self, metric='euclidean', shrink_threshold=None):
        self.metric = metric
        self.shrink_threshold = shrink_threshold

    def fit(self, X, y):
        """
        Fit the NearestCentroid model according to the given training data.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]
            Training vector, where n_samples in the number of samples and
            n_features is the number of features.
        y : array, shape = [n_samples]
            Target values (integers)
        """
        X, y = check_arrays(X, y)
        n_samples, n_features = X.shape
        classes = np.unique(y)
        classes.sort()
        self.classes_ = classes
        n_classes = classes.size
        if n_classes < 2:
            raise ValueError('y has less than 2 classes')
        self.centroids_ = np.array([np.mean([X[i] for i in range(n_samples)
                                             if y[i] == cur_class], axis=0)
                                    for cur_class in classes])
        self.dataset_centroid_ = np.mean(X, axis=0)
        if self.shrink_threshold:
            # Number of clusters in each class.
            nk = np.array([np.sum(classes == cur_class)
                           for cur_class in classes])
            # m parameter for determining deviation
            m = np.sqrt((1. / nk) + (1. / n_samples))
            # Calculate deviation using the standard deviation of centroids.
            variance = np.sum(np.power(X - self.centroids_[y], 2), axis=0)
            s = np.sqrt(variance / (n_samples - n_classes))
            s0 = np.median(s)  # To deter outliers from affecting the results.
            s = s + s0
            assert s.shape == (n_features,), "{},{}".format(s.shape, n_features)
            ms = np.array([[m[j] * s[i] for i in range(n_features)]
                           for j in range(n_classes)])
            deviation = ((self.centroids_ - self.dataset_centroid_) / ms)
            # Soft thresholding: if the deviation crosses 0 during shrinking,
            # it becomes zero.
            signs = np.sign(deviation)
            deviation = (np.abs(deviation) - self.shrink_threshold)
            deviation[deviation < 0] = 0
            deviation *= signs
            # Now adjust the centroids using the deviation
            self.centroids_ = self.dataset_centroid_ + (ms * deviation)
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
        """
        if not hasattr(self, "centroids_"):
            raise AttributeError("Model has not been trained yet.")
        return self.classes_[pairwise_distances(
            X, self.centroids_, metric=self.metric).argmin(axis=1)]

