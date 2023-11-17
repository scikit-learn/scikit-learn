"""
Nearest Centroid Classification
"""

# Author: Robert Layton <robertlayton@gmail.com>
#         Olivier Grisel <olivier.grisel@ensta.org>
#         Andreas W. Kempa-Liehr <a.kempa-liehr@auckland.ac.nz>
#         Matthew Ning <mhn@bu.edu>
#
# License: BSD 3 clause

import warnings
from numbers import Real

import numpy as np
from scipy import sparse as sp

from sklearn.metrics.pairwise import _VALID_METRICS

from ..base import BaseEstimator, ClassifierMixin, _fit_context
from ..discriminant_analysis import DiscriminantAnalysisPredictionMixin
from ..metrics.pairwise import pairwise_distances, pairwise_distances_argmin
from ..preprocessing import LabelEncoder
from ..utils._param_validation import Interval, StrOptions
from ..utils.multiclass import check_classification_targets
from ..utils.sparsefuncs import csc_median_axis_0
from ..utils.validation import check_is_fitted


class NearestCentroid(
    DiscriminantAnalysisPredictionMixin, ClassifierMixin, BaseEstimator
):
    """Nearest centroid classifier.

    Each class is represented by its centroid, with test samples classified to
    the class with the nearest centroid.

    Read more in the :ref:`User Guide <nearest_centroid_classifier>`.

    Parameters
    ----------
    metric : str or callable, default="euclidean"
        Metric to use for distance computation. See the documentation of
        `scipy.spatial.distance
        <https://docs.scipy.org/doc/scipy/reference/spatial.distance.html>`_ and
        the metrics listed in
        :class:`~sklearn.metrics.pairwise.distance_metrics` for valid metric
        values. Note that "wminkowski", "seuclidean" and "mahalanobis" are not
        supported.

        The centroids for the samples corresponding to each class is
        the point from which the sum of the distances (according to the metric)
        of all samples that belong to that particular class are minimized.
        If the `"manhattan"` metric is provided, this centroid is the median
        and for all other metrics, the centroid is now set to be the mean.

        .. deprecated:: 1.3
            Support for metrics other than `euclidean` and `manhattan` and for
            callables was deprecated in version 1.3 and will be removed in
            version 1.5.

        .. versionchanged:: 0.19
            `metric='precomputed'` was deprecated and now raises an error

    shrink_threshold : float, default=None
        Threshold for shrinking centroids to remove features.

    priors : {"uniform", "empirical"} or array-like of shape (n_classes,), \
    default="empirical"
        The class prior probabilities. By default, the class proportions are
        inferred from the training data.

        .. versionadded:: 1.4

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

    deviations_ : ndarray of shape(n_classes, n_features)
        Deviation of each class using soft thresholding.

        .. versionadded:: 1.4

    within_class_std_ : ndarray of shape(n_features,)
        Within-class standard deviation with unshrunked centroids.

        .. versionadded:: 1.4

    class_priors_ : ndarray of shape(n_classes,)
        The class prior probabilities.

        .. versionadded:: 1.4

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

    _valid_metrics = set(_VALID_METRICS) - {"mahalanobis", "seuclidean", "wminkowski"}

    _parameter_constraints: dict = {
        "metric": [
            StrOptions(
                _valid_metrics, deprecated=_valid_metrics - {"manhattan", "euclidean"}
            ),
            callable,
        ],
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
        if isinstance(self.metric, str) and self.metric not in (
            "manhattan",
            "euclidean",
        ):
            warnings.warn(
                (
                    "Support for distance metrics other than euclidean and "
                    "manhattan and for callables was deprecated in version "
                    "1.3 and will be removed in version 1.5."
                ),
                FutureWarning,
            )

        # If X is sparse and the metric is "manhattan", store it in a csc
        # format is easier to calculate the median.
        if self.metric == "manhattan":
            X, y = self._validate_data(X, y, accept_sparse=["csc"])
        else:
            X, y = self._validate_data(X, y, accept_sparse=["csr", "csc"])
        is_X_sparse = sp.issparse(X)
        if is_X_sparse and self.shrink_threshold:
            raise ValueError("threshold shrinking not supported for sparse input")
        check_classification_targets(y)

        n_samples, n_features = X.shape
        le = LabelEncoder()
        y_ind = le.fit_transform(y)
        self.classes_ = classes = le.classes_
        n_classes = classes.size
        if n_classes < 2:
            raise ValueError(
                "The number of classes has to be greater than one;got %d class"
                % (n_classes)
            )

        if self.priors == "empirical":  # estimate priors from sample
            _, class_counts = np.unique(y, return_inverse=True)  # non-negative ints
            self.class_priors_ = np.bincount(class_counts) / float(len(y))
        elif self.priors == "uniform":
            self.class_priors_ = np.asarray([1 / n_classes] * n_classes)
        else:
            self.class_priors_ = np.asarray(self.priors)

        if (self.class_priors_ < 0).any():
            raise ValueError("priors must be non-negative")
        if not np.isclose(self.class_priors_.sum(), 1.0):
            warnings.warn(
                "The priors do not sum to 1. Normalizing such that it sums to one.",
                UserWarning,
            )
            self.class_priors_ = self.class_priors_ / self.class_priors_.sum()

        # Mask mapping each class to its members.
        self.centroids_ = np.empty((n_classes, n_features), dtype=np.float64)
        self.deviations_ = np.empty((n_classes, n_features), dtype=np.float64)
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
            else:
                # TODO(1.5) remove warning when metric is only manhattan or euclidean
                if self.metric != "euclidean":
                    warnings.warn(
                        "Averaging for metrics other than "
                        "euclidean and manhattan not supported. "
                        "The average is set to be the mean."
                    )
                self.centroids_[cur_class] = X[center_mask].mean(axis=0)

        # Compute within-class std_dev with unshrunked centroids
        variance = np.square(X - self.centroids_[y_ind])
        self.within_class_std_ = np.sqrt(variance.sum(axis=0) / (n_samples - n_classes))

        if not is_X_sparse:
            if np.all(np.ptp(X, axis=0) == 0):
                raise ValueError("All features have zero variance. Division by zero.")
        else:
            if np.all((X.max(axis=0) - X.min(axis=0)).toarray() == 0):
                raise ValueError("All features have zero variance. Division by zero.")
        dataset_centroid_ = np.mean(X, axis=0)
        # m parameter for determining deviation
        m = np.sqrt((1.0 / nk) - (1.0 / n_samples))
        # Calculate deviation using the standard deviation of centroids.
        # To deter outliers from affecting the results.
        s = self.within_class_std_ + np.median(self.within_class_std_)
        mm = m.reshape(len(m), 1)  # Reshape to allow broadcasting.
        ms = mm * s
        self.deviations_ = (self.centroids_ - dataset_centroid_) / ms
        # Soft thresholding: if the deviation crosses 0 during shrinking,
        # it becomes zero.
        if self.shrink_threshold:
            signs = np.sign(self.deviations_)
            self.deviations_ = np.abs(self.deviations_) - self.shrink_threshold
            np.clip(self.deviations_, 0, None, out=self.deviations_)
            self.deviations_ *= signs
            # Now adjust the centroids using the deviation
            msd = ms * self.deviations_
            self.centroids_ = dataset_centroid_[np.newaxis, :] + msd
        return self

    # TODO(1.5) remove note about precomputed metric
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

        Notes
        -----
        If the metric constructor parameter is `"precomputed"`, `X` is assumed
        to be the distance matrix between the data to be predicted and
        `self.centroids_`.
        """
        check_is_fitted(self)

        if len(np.unique(self.class_priors_, return_counts=True)[0]) == 1:
            X = self._validate_data(X, accept_sparse="csr", reset=False)
            return self.classes_[
                pairwise_distances_argmin(X, self.centroids_, metric=self.metric)
            ]
        else:
            return super().predict(X)

    def _decision_function(self, X):
        check_is_fitted(self, "centroids_")

        X = self._validate_data(X, reset=False, accept_sparse="csr")

        # check if this works with sparse matrix formats later
        discriminant_score = np.empty(
            (X.shape[0], self.classes_.size), dtype=np.float64
        )

        for cur_class in range(self.classes_.size):
            # use pairwise_distances function instead. line 2093
            Xdist_norm = pairwise_distances(
                X / self.within_class_std_,
                (self.centroids_[cur_class, :].reshape(1, -1) / self.within_class_std_),
                metric=self.metric,
            ).reshape(-1)
            Xdist_norm = Xdist_norm**2
            discriminant_score[:, cur_class] = np.squeeze(
                -Xdist_norm + 2.0 * np.log(self.class_priors_[cur_class])
            )

        return discriminant_score

    def decision_function(self, X):
        """Apply decision function to an array of samples.

        The estimation has been implemented according to
        Hastie et al. (2009), p. 652 equation (18.2).

        Note that this function is only supported for
        ``euclidean`` metric.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            Array of samples (test vectors).

        Returns
        -------
        C : ndarray of shape (n_samples,) or (n_samples, n_classes)
            Decision function values related to each class, per sample.
            In the two-class case, the shape is (n_samples,), giving the
            log likelihood ratio of the positive class.
        """
        if self.metric == "euclidean":
            return super().decision_function(X)
        else:
            raise TypeError("decision_function is only supported for Euclidean metric")

    def predict_proba(self, X):
        """Estimate class probabilities.

        Note that this function is only supported for
        ``euclidean`` metric.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            Input data.

        Returns
        -------
        y_proba : ndarray of shape (n_samples, n_classes)
            Returns the probability estimate of the sample for each class in the
            model, where classes are ordered as they are in `self.classes_`.
        """
        if self.metric == "euclidean":
            return super().predict_proba(X)
        else:
            raise TypeError("predict_proba is only supported for Euclidean metric")

    def predict_log_proba(self, X):
        """Estimate log class probabilities.

        Note that this function is only supported for
        ``euclidean`` metric.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            Input data.

        Returns
        -------
        y_log_proba : ndarray of shape (n_samples, n_classes)
            Estimated log probabilities.
        """
        if self.metric == "euclidean":
            return super().predict_log_proba(X)
        else:
            raise TypeError("predict_log_proba is only supported for Euclidean metric")
