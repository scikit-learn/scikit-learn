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
from sklearn.metrics.pairwise import pairwise_distances
from sklearn.preprocessing import LabelEncoder
from sklearn.utils import get_tags
from sklearn.utils._available_if import available_if
from sklearn.utils._param_validation import Interval, StrOptions
from sklearn.utils.multiclass import check_classification_targets
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

    shrink_threshold : float, default=None
        Threshold for shrinking centroids to remove features.

    priors : {"uniform", "empirical"} or array-like of shape (n_classes,), \
        default="uniform"
        The class prior probabilities. By default, the class proportions are
        inferred from the training data.

    Attributes
    ----------
    centroids_ : array-like of shape (n_classes, n_features)
        Centroid of each class.

    classes_ : array of shape (n_classes,)
        The unique classes labels.

    n_features_in_ : int
        Number of features seen during :term:`fit`.

    feature_names_in_ : ndarray of shape (`n_features_in_`,)
        Names of features seen during :term:`fit`. Defined only when `X`
        has feature names that are all strings.

    deviations_ : ndarray of shape (n_classes, n_features)
        Deviations (or shrinkages) of the centroids of each class from the
        overall centroid.

    within_class_std_dev_ : ndarray of shape (n_features,)
        Pooled or within-class standard deviation of input data.

    class_prior_ : ndarray of shape (n_classes,)
        The class prior probabilities.

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
        """Fit the NearestCentroid model according to the given training data.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            Training vector, where `n_samples` is the number of samples and
            `n_features` is the number of features.
        y : array-like of shape (n_samples,)
            Target values.

        Returns
        -------
        self : object
            Fitted estimator.
        """
        if self.metric == "euclidean":
            # Clear internal state to treat this as the first partial_fit call
            for attr in ["classes_", "_class_counts", "_true_centroids", "_variance"]:
                if hasattr(self, attr):
                    delattr(self, attr)

            # Strip Mock objects from tests to avoid array_function errors
            y = np.asarray(y)
            classes = np.unique(y)
            return self.partial_fit(X, y, classes=classes)

        # ---------- ORIGINAL BATCH MANHATTAN LOGIC ----------
        X, y = validate_data(self, X, y, accept_sparse=["csc"])
        is_X_sparse = sp.issparse(X)
        check_classification_targets(y)

        n_samples, n_features = X.shape
        le = LabelEncoder()
        y_ind = le.fit_transform(y)
        self.classes_ = classes = le.classes_
        n_classes = classes.size

        if n_classes < 2:
            raise ValueError(
                "The number of classes has to be greater than one; "
                "got {n_classes} class"
            )

        if self.priors == "empirical":
            _, class_counts = np.unique(y, return_inverse=True)
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

        self.centroids_ = np.empty((n_classes, n_features), dtype=np.float64)
        nk = np.zeros(n_classes)

        for cur_class in range(n_classes):
            center_mask = y_ind == cur_class
            nk[cur_class] = np.sum(center_mask)
            if is_X_sparse:
                center_mask = np.where(center_mask)[0]

            if not is_X_sparse:
                self.centroids_[cur_class] = np.median(X[center_mask], axis=0)
            else:
                self.centroids_[cur_class] = csc_median_axis_0(X[center_mask])

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
        if (
            is_X_sparse
            and np.all((X.max(axis=0) - X.min(axis=0)).toarray() == 0)
            or not is_X_sparse
            and np.all(np.ptp(X, axis=0) == 0)
        ):
            raise ValueError(err_msg)

        dataset_centroid_ = X.mean(axis=0)
        m = np.sqrt((1.0 / nk) - (1.0 / n_samples))
        s = self.within_class_std_dev_ + np.median(self.within_class_std_dev_)
        mm = m.reshape(len(m), 1)
        ms = mm * s
        self.deviations_ = np.array(
            (self.centroids_ - dataset_centroid_) / ms, copy=False
        )

        if self.shrink_threshold:
            signs = np.sign(self.deviations_)
            self.deviations_ = np.abs(self.deviations_) - self.shrink_threshold
            np.clip(self.deviations_, 0, None, out=self.deviations_)
            self.deviations_ *= signs
            msd = ms * self.deviations_
            self.centroids_ = np.array(dataset_centroid_ + msd, copy=False)

        return self

    @_fit_context(prefer_skip_nested_validation=True)
    def partial_fit(self, X, y, classes=None):
        """Incremental fit on a batch of samples.

        This method is expected to be called several times consecutively
        on different chunks of a dataset so as to implement out-of-core
        or online learning.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            Training vectors, where `n_samples` is the number of samples and
            `n_features` is the number of features.
        y : array-like of shape (n_samples,)
            Target values.
        classes : array-like of shape (n_classes,), default=None
            List of all the classes that can possibly appear in the `y` vector.
            Must be provided at the first call to partial_fit, can be omitted
            in subsequent calls.

        Returns
        -------
        self : object
            Returns the instance itself.
        """
        if self.metric == "manhattan":
            raise ValueError(
                "The 'manhattan' metric does not support partial_fit "
                "because the median cannot be updated incrementally."
            )

        first_call = not hasattr(self, "classes_")

        ensure_all_finite = "allow-nan" if get_tags(self).input_tags.allow_nan else True
        X, y = validate_data(
            self,
            X,
            y,
            ensure_all_finite=ensure_all_finite,
            accept_sparse=["csr", "csc"],
            reset=first_call,
        )
        check_classification_targets(y)

        # Strip Mock objects to avoid array_function errors
        y_arr = np.asarray(y)

        if first_call:
            if classes is None:
                raise ValueError(
                    "classes must be passed on the first call to partial_fit."
                )
            self.classes_ = np.array(classes)
            n_classes = self.classes_.shape[0]
            n_features = X.shape[1]

            if n_classes < 2:
                raise ValueError(
                    f"The number of classes has to be greater than one; "
                    f"got {n_classes} class"
                )

            self._class_counts = np.zeros(n_classes, dtype=np.float64)
            self._true_centroids = np.zeros((n_classes, n_features), dtype=np.float64)
            self._variance = np.zeros((n_classes, n_features), dtype=np.float64)

            if self.priors == "empirical":
                self.class_prior_ = np.zeros(n_classes, dtype=np.float64)
            elif self.priors == "uniform":
                self.class_prior_ = np.full(n_classes, 1.0 / n_classes)
            else:
                self.class_prior_ = np.asarray(self.priors)
                if (self.class_prior_ < 0).any():
                    raise ValueError("priors must be non-negative")
                if not np.isclose(self.class_prior_.sum(), 1.0):
                    warnings.warn(
                        "The priors do not sum to 1. "
                        "Normalizing such that it sums to one.",
                        UserWarning,
                    )
                    self.class_prior_ = self.class_prior_ / self.class_prior_.sum()

        n_classes = self.classes_.size
        is_X_sparse = sp.issparse(X)

        # Track global feature min/max to perfectly detect zero-variance features
        if is_X_sparse:
            batch_min = X.min(axis=0).toarray().ravel()
            batch_max = X.max(axis=0).toarray().ravel()
        else:
            batch_min = np.min(X, axis=0)
            batch_max = np.max(X, axis=0)

        if first_call:
            self._feature_min = batch_min.copy()
            self._feature_max = batch_max.copy()
        else:
            self._feature_min = np.minimum(self._feature_min, batch_min)
            self._feature_max = np.maximum(self._feature_max, batch_max)

        le = LabelEncoder()
        le.fit(self.classes_)
        if np.setdiff1d(y_arr, self.classes_).size > 0:
            raise ValueError("`y` has classes not in `classes`")
        y_ind = le.transform(y_arr)

        # Online mean and variance updates per class
        for cur_class in range(n_classes):
            center_mask = y_ind == cur_class
            if is_X_sparse:
                center_mask_idx = np.where(center_mask)[0]
                X_class = X[center_mask_idx]
            else:
                X_class = X[center_mask]

            N_new = X_class.shape[0]
            if N_new == 0:
                continue

            if is_X_sparse:
                batch_mean = np.asarray(X_class.mean(axis=0)).ravel()
                batch_var = (
                    np.asarray(X_class.multiply(X_class).mean(axis=0)).ravel()
                    - batch_mean**2
                )
            else:
                batch_mean = np.mean(X_class, axis=0)
                batch_var = np.var(X_class, axis=0)

            N_old = self._class_counts[cur_class]
            N_total = N_old + N_new

            if N_old == 0:
                self._true_centroids[cur_class] = batch_mean
                self._variance[cur_class] = batch_var
            else:
                old_mean = self._true_centroids[cur_class]
                old_var = self._variance[cur_class]

                new_mean = (N_old * old_mean + N_new * batch_mean) / N_total
                new_var = (N_old * old_var + N_new * batch_var) / N_total + (
                    N_old * N_new * (old_mean - batch_mean) ** 2
                ) / (N_total**2)

                self._true_centroids[cur_class] = new_mean
                self._variance[cur_class] = new_var

            self._class_counts[cur_class] = N_total

        total_samples = self._class_counts.sum()

        if self.priors == "empirical":
            self.class_prior_ = self._class_counts / total_samples

        self.centroids_ = self._true_centroids.copy()

        # Calculate pooled standard deviation using class variances
        total_sse = np.sum(self._variance * self._class_counts[:, np.newaxis], axis=0)

        if total_samples > n_classes:
            self.within_class_std_dev_ = np.sqrt(
                total_sse / (total_samples - n_classes)
            )
        else:
            self.within_class_std_dev_ = np.zeros(self._true_centroids.shape[1])

        if any(self.within_class_std_dev_ == 0) and total_samples > n_classes:
            warnings.warn(
                "self.within_class_std_dev_ has at least 1 zero standard deviation. "
                "Inputs within the same classes for at least 1 feature are identical."
            )

        if np.all((self._feature_max - self._feature_min) == 0):
            raise ValueError("All features have zero variance. Division by zero.")

        # Recalculate global deviations and apply soft thresholding
        dataset_centroid_ = np.average(
            self._true_centroids, axis=0, weights=self._class_counts
        )
        m = np.sqrt(
            (1.0 / np.maximum(self._class_counts, 1)) - (1.0 / max(total_samples, 1))
        )
        s = self.within_class_std_dev_ + np.median(self.within_class_std_dev_)
        ms = m.reshape(-1, 1) * s

        with np.errstate(divide="ignore", invalid="ignore"):
            self.deviations_ = np.where(
                ms == 0, 0, (self.centroids_ - dataset_centroid_) / ms
            )

        if self.shrink_threshold:
            signs = np.sign(self.deviations_)
            self.deviations_ = np.abs(self.deviations_) - self.shrink_threshold
            np.clip(self.deviations_, 0, None, out=self.deviations_)
            self.deviations_ *= signs

            msd = ms * self.deviations_
            self.centroids_ = dataset_centroid_ + msd

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
        check_is_fitted(self, ["centroids_", "classes_"])

        ensure_all_finite = "allow-nan" if get_tags(self).input_tags.allow_nan else True

        X = validate_data(
            self,
            X,
            ensure_all_finite=ensure_all_finite,
            accept_sparse="csr",
            reset=False,
        )

        centroids = np.asarray(self.centroids_, dtype=np.float64, order="C")

        # ---------- EUCLIDEAN ----------
        if self.metric == "euclidean":
            if sp.issparse(X):
                X_sq = np.asarray(X.multiply(X).sum(axis=1)).ravel()
                C_sq = np.sum(centroids * centroids, axis=1)
                cross = X @ centroids.T
                distances = X_sq[:, None] + C_sq[None, :] - 2 * cross
            else:
                X = np.asarray(X, dtype=np.float64, order="C")
                X_sq = np.sum(X * X, axis=1)
                C_sq = np.sum(centroids * centroids, axis=1)
                cross = X @ centroids.T
                distances = X_sq[:, None] + C_sq[None, :] - 2 * cross

        # ---------- MANHATTAN ----------
        elif self.metric == "manhattan":
            if sp.issparse(X):
                X = X.toarray()
            else:
                X = np.asarray(X, dtype=np.float64, order="C")

            distances = np.empty((X.shape[0], centroids.shape[0]), dtype=np.float64)
            for i, c in enumerate(centroids):
                distances[:, i] = np.sum(np.abs(X - c), axis=1)

        else:
            raise ValueError(f"Unsupported metric: {self.metric}")

        # ---------- APPLY PRIORS ----------
        if hasattr(self, "class_prior_") and self.class_prior_ is not None:
            log_priors = np.log(self.class_prior_)
            scores = -distances + 2.0 * log_priors
            return self.classes_[np.argmax(scores, axis=1)]

        return self.classes_[np.argmin(distances, axis=1)]

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
