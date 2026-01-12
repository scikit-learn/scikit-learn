"""
Nearest Centroid Classification
"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import warnings
from numbers import Real

import numpy as np
from scipy import sparse as sp

from sklearn.base import BaseEstimator, ClassifierMixin, _fit_context
from sklearn.discriminant_analysis import DiscriminantAnalysisPredictionMixin
from sklearn.metrics.pairwise import pairwise_distances, pairwise_distances_argmin
from sklearn.preprocessing import LabelEncoder
from sklearn.utils import get_tags
from sklearn.utils._available_if import available_if
from sklearn.utils._param_validation import Interval, StrOptions
from sklearn.utils.multiclass import (
    check_classification_targets,
    _check_partial_fit_first_call,
)
from sklearn.utils.sparsefuncs import csc_median_axis_0
from sklearn.utils.validation import check_is_fitted, validate_data


class NearestCentroid(
    DiscriminantAnalysisPredictionMixin, ClassifierMixin, BaseEstimator
):
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

    priors : {"uniform", "empirical"} or array-like of shape (n_classes,), \
        default="uniform"
        The class prior probabilities. By default, the class proportions are
        inferred from the training data.

        .. versionadded:: 1.6

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

    deviations_ : ndarray of shape (n_classes, n_features)
        Deviations (or shrinkages) of the centroids of each class from the
        overall centroid. Equal to eq. (18.4) if `shrink_threshold=None`,
        else (18.5) p. 653 of [2]. Can be used to identify features used
        for classification.

        .. versionadded:: 1.6

    within_class_std_dev_ : ndarray of shape (n_features,)
        Pooled or within-class standard deviation of input data.

        .. versionadded:: 1.6

    class_prior_ : ndarray of shape (n_classes,)
        The class prior probabilities.

        .. versionadded:: 1.6

    See Also
    --------
    KNeighborsClassifier : Nearest neighbors classifier.

    Notes
    -----
    When used for text classification with tf-idf vectors, this classifier is
    also known as the Rocchio classifier.

    References
    ----------
    [1] Tibshirani, R., Hastie, T., Narasimhan, B., & Chu, G. (2002). Diagnosis of
    multiple cancer types by shrunken centroids of gene expression. Proceedings
    of the National Academy of Sciences of the United States of America,
    99(10), 6567-6572. The National Academy of Sciences.

    [2] Hastie, T., Tibshirani, R., Friedman, J. (2009). The Elements of Statistical
    Learning Data Mining, Inference, and Prediction. 2nd Edition. New York, Springer.

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
    """

    _parameter_constraints: dict = {
        "metric": [StrOptions({"manhattan", "euclidean"})],
        "shrink_threshold": [Interval(Real, 0, None, closed="neither"), None],
        "priors": ["array-like", StrOptions({"empirical", "uniform"})],
    }

    def __init__(
        self,
        metric="euclidean",
        *,
        shrink_threshold=None,
        priors="uniform",
    ):
        self.metric = metric
        self.shrink_threshold = shrink_threshold
        self.priors = priors

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
        # If X is sparse and the metric is "manhattan", store it in a csc
        # format is easier to calculate the median.
        if self.metric == "manhattan":
            X, y = validate_data(self, X, y, accept_sparse=["csc"])
        else:
            ensure_all_finite = (
                "allow-nan" if get_tags(self).input_tags.allow_nan else True
            )
            X, y = validate_data(
                self,
                X,
                y,
                ensure_all_finite=ensure_all_finite,
                accept_sparse=["csr", "csc"],
            )
        is_X_sparse = sp.issparse(X)
        check_classification_targets(y)

        n_samples, n_features = X.shape
        le = LabelEncoder()
        y_ind = le.fit_transform(y)
        self.classes_ = classes = le.classes_
        n_classes = classes.size
        if n_classes < 2:
            raise ValueError(
                "The number of classes has to be greater than one; got %d class"
                % (n_classes)
            )

        if self.priors == "empirical":  # estimate priors from sample
            _, class_counts = np.unique(y, return_inverse=True)  # non-negative ints
            self.class_prior_ = np.bincount(class_counts) / float(len(y))
        elif self.priors == "uniform":
            self.class_prior_ = np.asarray([1 / n_classes] * n_classes)
        else:
            self.class_prior_ = np.asarray(self.priors)

        if (self.class_prior_ < 0).any():
            raise ValueError("priors must be non-negative")
        if not np.isclose(self.class_prior_.sum(), 1.0):
            warnings.warn(
                "The priors do not sum to 1. Normalizing such that it sums to one.",
                UserWarning,
            )
            self.class_prior_ = self.class_prior_ / self.class_prior_.sum()

        # Mask mapping each class to its members.
        self.centroids_ = np.empty((n_classes, n_features), dtype=np.float64)

        # Number of clusters in each class.
        nk = np.zeros(n_classes)

        for cur_class in range(n_classes):
            center_mask = y_ind == cur_class
            nk[cur_class] = np.sum(center_mask)
            if is_X_sparse:
                center_mask = np.where(center_mask)[0]

            if self.metric == "manhattan":
                # NumPy does not calculate median of sparse matrices.
                if not is_X_sparse:
                    self.centroids_[cur_class] = np.median(X[center_mask], axis=0)
                else:
                    self.centroids_[cur_class] = csc_median_axis_0(X[center_mask])
            else:  # metric == "euclidean"
                self.centroids_[cur_class] = X[center_mask].mean(axis=0)

        # Compute within-class std_dev with unshrunked centroids
        variance = np.array(X - self.centroids_[y_ind], copy=False) ** 2
        self.within_class_std_dev_ = np.array(
            np.sqrt(variance.sum(axis=0) / (n_samples - n_classes)), copy=False
        )
        if any(self.within_class_std_dev_ == 0):
            warnings.warn(
                "self.within_class_std_dev_ has at least 1 zero standard deviation."
                "Inputs within the same classes for at least 1 feature are identical."
            )

        err_msg = "All features have zero variance. Division by zero."
        if is_X_sparse and np.all((X.max(axis=0) - X.min(axis=0)).toarray() == 0):
            raise ValueError(err_msg)
        elif not is_X_sparse and np.all(np.ptp(X, axis=0) == 0):
            raise ValueError(err_msg)

        dataset_centroid_ = X.mean(axis=0)
        # m parameter for determining deviation
        m = np.sqrt((1.0 / nk) - (1.0 / n_samples))
        # Calculate deviation using the standard deviation of centroids.
        # To deter outliers from affecting the results.
        s = self.within_class_std_dev_ + np.median(self.within_class_std_dev_)
        mm = m.reshape(len(m), 1)  # Reshape to allow broadcasting.
        ms = mm * s
        self.deviations_ = np.array(
            (self.centroids_ - dataset_centroid_) / ms, copy=False
        )
        # Soft thresholding: if the deviation crosses 0 during shrinking,
        # it becomes zero.
        if self.shrink_threshold:
            signs = np.sign(self.deviations_)
            self.deviations_ = np.abs(self.deviations_) - self.shrink_threshold
            np.clip(self.deviations_, 0, None, out=self.deviations_)
            self.deviations_ *= signs
            # Now adjust the centroids using the deviation
            msd = ms * self.deviations_
            self.centroids_ = np.array(dataset_centroid_ + msd, copy=False)
        return self

    @_fit_context(prefer_skip_nested_validation=True)
    def partial_fit(self, X, y, classes=None):
        """
        Incremental fit on a batch of samples.

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
        y : array-like of shape (n_samples,)
            Target values.
        classes : array-like of shape (n_classes,), default=None
            List of all the classes that can possibly appear in the y vector.

            Must be provided at the first call to partial_fit, can be omitted
            in subsequent calls.

        Returns
        -------
        self : object
            Fitted estimator.
        """
        first_call = _check_partial_fit_first_call(self, classes)
        ensure_all_finite = "allow-nan" if get_tags(self).input_tags.allow_nan else True
        X, y = validate_data(
            self,
            X,
            y,
            reset=first_call,
            ensure_all_finite=ensure_all_finite,
            accept_sparse=["csr", "csc"],
        )
        check_classification_targets(y)

        n_samples, n_features = X.shape

        if self.metric == "manhattan":
            raise ValueError("partial_fit does not support the 'manhattan' metric.")

        if first_call:
            self.n_features_in_ = n_features
            n_classes = len(self.classes_)
            self.nk_ = np.zeros(n_classes, dtype=np.int64)
            self.true_centroids_ = np.zeros((n_classes, n_features), dtype=np.float64)
        else:
            n_classes = len(self.classes_)
            if n_features != self.n_features_in_:
                raise ValueError("Number of features does not match previous calls.")

        unique_y = np.unique(y)
        for y_i in unique_y:
            if y_i not in self.classes_:
                raise ValueError(f"Unknown class {y_i} in y.")
            i = np.searchsorted(self.classes_, y_i)
            mask = y == y_i
            n_new = mask.sum()
            if n_new == 0:
                continue
            mu_new = X[mask].mean(axis=0)
            n_old = self.nk_[i]
            mu_old = self.true_centroids_[i]
            self.true_centroids_[i] = (n_old * mu_old + n_new * mu_new) / (
                n_old + n_new
            )
            self.nk_[i] += n_new

        self.centroids_ = self.true_centroids_.copy()
        return self

    def predict(self, X):
        """Perform classification on an array of test vectors `X`.

        The predicted class `C` for each sample in `X` is returned.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            Input data.

        Returns
        -------
        y_pred : ndarray of shape (n_samples,)
            The predicted classes.
        """
        check_is_fitted(self)
        if np.isclose(self.class_prior_, 1 / len(self.classes_)).all():
            # `validate_data` is called here since we are not calling `super()`
            ensure_all_finite = (
                "allow-nan" if get_tags(self).input_tags.allow_nan else True
            )
            X = validate_data(
                self,
                X,
                ensure_all_finite=ensure_all_finite,
                accept_sparse="csr",
                reset=False,
            )
            return self.classes_[
                pairwise_distances_argmin(X, self.centroids_, metric=self.metric)
            ]
        else:
            return super().predict(X)

    def _decision_function(self, X):
        # return discriminant scores, see eq. (18.2) p. 652 of the ESL.
        check_is_fitted(self, "centroids_")

        X_normalized = validate_data(
            self, X, copy=True, reset=False, accept_sparse="csr", dtype=np.float64
        )

        discriminant_score = np.empty(
            (X_normalized.shape[0], self.classes_.size), dtype=np.float64
        )

        mask = self.within_class_std_dev_ != 0
        X_normalized[:, mask] /= self.within_class_std_dev_[mask]
        centroids_normalized = self.centroids_.copy()
        centroids_normalized[:, mask] /= self.within_class_std_dev_[mask]

        for class_idx in range(self.classes_.size):
            distances = pairwise_distances(
                X_normalized, centroids_normalized[[class_idx]], metric=self.metric
            ).ravel()
            distances **= 2
            discriminant_score[:, class_idx] = np.squeeze(
                -distances + 2.0 * np.log(self.class_prior_[class_idx])
            )

        return discriminant_score

    def _check_euclidean_metric(self):
        return self.metric == "euclidean"

    decision_function = available_if(_check_euclidean_metric)(
        DiscriminantAnalysisPredictionMixin.decision_function
    )

    predict_proba = available_if(_check_euclidean_metric)(
        DiscriminantAnalysisPredictionMixin.predict_proba
    )

    predict_log_proba = available_if(_check_euclidean_metric)(
        DiscriminantAnalysisPredictionMixin.predict_log_proba
    )

    def __sklearn_tags__(self):
        tags = super().__sklearn_tags__()
        tags.input_tags.allow_nan = self.metric == "nan_euclidean"
        tags.input_tags.sparse = True
        return tags
