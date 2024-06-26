"""
Nearest Centroid Classification
"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import warnings
from numbers import Real

import numpy as np
from scipy import sparse as sp

from ..base import BaseEstimator, ClassifierMixin, _fit_context
from ..metrics.pairwise import pairwise_distances_argmin
from ..preprocessing import LabelEncoder
from ..utils._param_validation import Interval, StrOptions
from ..utils.multiclass import (
    _check_partial_fit_first_call,
    check_classification_targets,
)
from ..utils.sparsefuncs import csc_median_axis_0
from ..utils.validation import check_is_fitted, column_or_1d


class NearestCentroid(ClassifierMixin, BaseEstimator):
    """Nearest centroid classifier.

    Each class is represented by its centroid, with test samples classified to
    the class with the nearest centroid.

    Read more in the :ref:`User Guide <nearest_centroid_classifier>`.

    Parameters
    ----------
    metric : {"euclidean", "manhattan"}, default="euclidean"
        Metric to use for distance computation.

        If `metric="euclidean"`, the centroid for the samples corresponding to each
        class is the arithmetic mean, which minimizes the sum of squared L1 distances.
        If `metric="manhattan"`, the centroid is the feature-wise median, which
        minimizes the sum of L1 distances.

        .. versionchanged:: 1.5
            All metrics but `"euclidean"` and `"manhattan"` were deprecated and
            now raise an error.

        .. versionchanged:: 0.19
            `metric='precomputed'` was deprecated and now raises an error

    shrink_threshold : float, default=None
        Threshold for shrinking centroids to remove features.

    Attributes
    ----------
    centroids_ : array-like of shape (n_classes, n_features)
        Centroid of each class.

    classes_ : array of shape (n_classes,)
        The unique classes labels.

    n_features_in_ : int
        Number of features seen during :term:`fit`.

        .. versionadded:: 0.24

    feature_names_in_ : ndarray of shape (`n_features_in_`,)
        Names of features seen during :term:`fit`. Defined only when `X`
        has feature names that are all strings.

        .. versionadded:: 1.0

    See Also
    --------
    KNeighborsClassifier : Nearest neighbors classifier.

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

    Examples
    --------
    >>> from sklearn.neighbors import NearestCentroid
    >>> import numpy as np
    >>> X = np.array([[-1, -1], [-2, -1], [-3, -2], [1, 1], [2, 1], [3, 2]])
    >>> y = np.array([1, 1, 1, 2, 2, 2])
    >>> clf = NearestCentroid()
    >>> clf.fit(X, y)
    NearestCentroid()
    >>> print(clf.predict([[-0.8, -1]]))
    [1]

    For a more detailed example see:
    :ref:`sphx_glr_auto_examples_neighbors_plot_nearest_centroid.py`
    """

    _parameter_constraints: dict = {
        "metric": [StrOptions({"manhattan", "euclidean"})],
        "shrink_threshold": [Interval(Real, 0, None, closed="neither"), None],
    }

    def __init__(self, metric="euclidean", *, shrink_threshold=None):
        self.metric = metric
        self.shrink_threshold = shrink_threshold

    @_fit_context(prefer_skip_nested_validation=True)
    def fit(self, X, y):
        """
        Fit the NearestCentroid model according to the given training data.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            Training vector, where `n_samples` is the number of samples and
            `n_features` is the number of features.
            Note that centroid shrinking cannot be used with sparse matrices.
        y : array-like of shape (n_samples,)
            Target values.

        Returns
        -------
        self : object
            Fitted estimator.
        """
        X, y = self._validate_data(X, y, accept_sparse=True)
        y = column_or_1d(y, warn=True)
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
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            Training vector, where `n_samples` is the number of samples and
            `n_features` is the number of features.
            Note that centroid shrinking cannot be used with sparse matrices.

        y : array-like of shape (n_samples,)
            Target values (integers).

        classes : array-like of shape (n_classes,), optional (default=None)
            List of all the classes that can possibly appear in the `y` vector.
            Must be provided at the first call to `partial_fit`, can be omitted
            in subsequent calls.
        """
        if self.metric == "manhattan":
            raise ValueError(
                "Partial fitting with 'manhattan' metric is not supported."
            )
        return self._partial_fit(X, y, classes, _refit=False)

    def _partial_fit(self, X, y, classes=None, _refit=False):
        """
        Actual implementation of the Nearest Centroid fitting.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            Training vector, where `n_samples` is the number of samples and
            `n_features` is the number of features.
            Note that centroid shrinking cannot be used with sparse matrices.

        y : array-like of shape (n_samples,)
            Target values.

        classes : array-like of shape (n_classes,), optional (default=None)
            List of all the classes that can possibly appear in the `y` vector.
            Must be provided at the first call to `partial_fit`, can be omitted
            in subsequent calls.

        _refit : bool, optional (default=False)
            If True, act as though this were the first time we called
            `_partial_fit` (i.e. throw away any past fitting and start over).

        Returns
        -------
        self : object
            Fitted estimator.
        """
        if self.metric == "precomputed":
            raise ValueError('Metric "precomputed" is not supported.')

        # If X is sparse and the metric is "manhattan", store it in a csc
        # format is easier to calculate the median.
        if self.metric == "manhattan":
            X, y = self._validate_data(X, y, accept_sparse=["csc"])
        else:
            X, y = self._validate_data(X, y, accept_sparse=["csr", "csc"])

        X_is_sparse = sp.issparse(X)

        if X_is_sparse and self.shrink_threshold:
            raise ValueError("Threshold shrinking not supported for sparse input")
        check_classification_targets(y)

        n_samples, n_features = X.shape

        if _refit or _check_partial_fit_first_call(self, classes):
            if self.metric != "euclidean":
                warnings.warn(
                    "Averaging for metrics other than "
                    "euclidean and manhattan not supported. "
                    "The average is set to be the mean."
                )

            self.true_classes_ = np.asarray(classes)

            # Mask mapping each class to its members.
            self.true_centroids_ = np.zeros(
                (self.true_classes_.size, n_features), dtype=np.float64
            )

            # Number of clusters in each class.
            self.nk_ = np.zeros(self.true_classes_.size)

            if self.shrink_threshold:
                self.ssd_ = np.zeros(
                    (self.true_classes_.size, n_features), dtype=np.float64
                )
                self.dataset_centroid_ = np.mean(X, axis=0)

        le = LabelEncoder()
        le.fit(self.true_classes_)
        y_ind = le.transform(y)

        self.classes_ = le.classes_
        n_classes = self.classes_.size

        if n_classes < 2:
            raise ValueError(
                f"The number of classes has to be greater than one; got {n_classes}"
            )

        if self.shrink_threshold:
            old_nk = self.nk_.copy()
            old_centroids = self.true_centroids_.copy()

        for cur_class in range(n_classes):
            center_mask = y_ind == cur_class

            if X[center_mask].size == 0:
                continue

            if X_is_sparse:
                center_mask = np.where(center_mask)[0]

            if self.metric == "manhattan":
                self.nk_[cur_class] += np.sum(center_mask)

                # NumPy does not calculate median of sparse matrices.
                if not X_is_sparse:
                    self.true_centroids_[cur_class] = np.median(X[center_mask], axis=0)
                else:
                    self.true_centroids_[cur_class] = csc_median_axis_0(X[center_mask])
            else:
                # Update each centroid weighted by the number of samples
                self.true_centroids_[cur_class] = (
                    X[center_mask].mean(axis=0) * np.sum(center_mask)
                    + self.true_centroids_[cur_class] * self.nk_[cur_class]
                )
                self.nk_[cur_class] += np.sum(center_mask)
                self.true_centroids_[cur_class] /= self.nk_[cur_class]

        # Filtering out centroids without any data
        self.classes_ = self.true_classes_[self.nk_ != 0]
        self.centroids_ = self.true_centroids_[self.nk_ != 0]

        if self.shrink_threshold:
            if np.all(np.ptp(X, axis=0) == 0):
                raise ValueError("All features have zero variance. Division by zero.")

            n_total = np.sum(self.nk_)
            self.dataset_centroid_ = (
                self.dataset_centroid_ * old_nk.sum(axis=0) + np.sum(X, axis=0)
            ) / n_total

            # Update sum of square distances of each class
            for cur_class in range(n_classes):
                n_old = old_nk[cur_class]
                n_new = self.nk_[cur_class] - n_old
                if n_new == 0:
                    continue
                center_mask = y_ind == cur_class
                old_ssd = self.ssd_[cur_class]
                new_ssd = (X[center_mask] - X[center_mask].mean(axis=0)) ** 2
                new_ssd = new_ssd.sum(axis=0)
                self.ssd_[cur_class] = (
                    old_ssd
                    + new_ssd
                    + (n_old / float(n_new * (n_new + n_old)))
                    * (
                        n_new * old_centroids[cur_class]
                        - n_new * X[center_mask].mean(axis=0)
                    )
                    ** 2
                )

            # m parameter for determining deviation
            m = np.sqrt((1.0 / self.nk_) - (1.0 / np.sum(self.nk_)))

            # Calculate deviation using the standard deviation of centroids.
            ssd = self.ssd_.sum(axis=0)
            s = np.sqrt(ssd / (n_total - n_classes))
            s += np.median(s)  # To deter outliers from affecting the results.
            ms = m.reshape(-1, 1) * s  # Reshape to allow broadcasting.
            deviation = (self.true_centroids_ - self.dataset_centroid_) / ms

            # Soft thresholding: if the deviation crosses 0 during shrinking,
            # it becomes zero.
            signs = np.sign(deviation)
            deviation = np.abs(deviation) - self.shrink_threshold
            np.clip(deviation, 0, None, out=deviation)
            deviation *= signs

            # Now adjust the centroids using the deviation
            msd = ms * deviation
            self.centroids_ = self.dataset_centroid_[np.newaxis, :] + msd

        return self

    def predict(self, X):
        """Perform classification on an array of test vectors `X`.

        The predicted class `C` for each sample in `X` is returned.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            Test samples.

        Returns
        -------
        C : ndarray of shape (n_samples,)
            The predicted classes.
        """
        check_is_fitted(self)

        X = self._validate_data(X, accept_sparse="csr", reset=False)
        return self.classes_[
            pairwise_distances_argmin(X, self.centroids_, metric=self.metric)
        ]
