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
from ..preprocessing import LabelEncoder
from ..utils._param_validation import Interval, StrOptions
from ..utils.multiclass import check_classification_targets
from ..utils.sparsefuncs import csc_median_axis_0
from ..utils.validation import check_is_fitted


class NearestCentroid(ClassifierMixin, BaseEstimator):
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

    priors : array-like of shape (n_classes,), default=None
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
        "priors": ["array-like", None],
    }

    def __init__(
        self,
        metric="euclidean",
        *,
        shrink_threshold=None,
        priors=None,
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

        if self.priors is None:  # estimate priors from sample
            _, cnts = np.unique(y, return_inverse=True)  # non-negative ints
            self.priors_ = np.bincount(cnts) / float(len(y))
        else:
            self.priors_ = np.asarray(self.priors)

        if (self.priors_ < 0).any():
            raise ValueError("priors must be non-negative")
        if not np.isclose(self.priors_.sum(), 1.0):
            warnings.warn("The priors do not sum to 1. Renormalizing", UserWarning)
            self.priors_ = self.priors_ / self.priors_.sum()

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
        variance = variance.sum(axis=0)
        self.s_ = np.sqrt(variance / (n_samples - n_classes))

        if self.shrink_threshold:
            if np.all(np.ptp(X, axis=0) == 0):
                raise ValueError("All features have zero variance.Division by zero.")
            dataset_centroid_ = np.mean(X, axis=0)

            # m parameter for determining deviation
            m = np.sqrt((1.0 / nk) - (1.0 / n_samples))
            # Calculate deviation using the standard deviation of centroids.
            # To deter outliers from affecting the results.
            s = self.s_ + np.median(self.s_)
            mm = m.reshape(len(m), 1)  # Reshape to allow broadcasting.
            ms = mm * s
            deviation = (self.centroids_ - dataset_centroid_) / ms
            # Soft thresholding: if the deviation crosses 0 during shrinking,
            # it becomes zero.
            signs = np.sign(deviation)
            deviation = np.abs(deviation) - self.shrink_threshold
            np.clip(deviation, 0, None, out=deviation)
            deviation *= signs
            # Now adjust the centroids using the deviation
            msd = ms * deviation
            self.centroids_ = dataset_centroid_[np.newaxis, :] + msd
        return self

    # TODO(1.5) remove note about precomputed metric
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

        Notes
        -----
        If the metric constructor parameter is `"precomputed"`, `X` is assumed
        to be the distance matrix between the data to be predicted and
        `self.centroids_`.
        """
        check_is_fitted(self)

        d = self._decision_function(X)
        y_pred = self.classes_.take(d.argmax(1))
        return y_pred

    def _decision_function(self, X):
        check_is_fitted(self, "centroids_")

        X = self._validate_data(X, reset=False, accept_sparse="csr")

        # check if this works with sparse matrix formats later
        discriminant_score = np.empty(
            (X.shape[0], self.classes_.size), dtype=np.float64
        )

        for cur_class in range(self.classes_.size):
            Xdist = X - self.centroids_[cur_class, :]
            Xdist_norm = np.square(Xdist / self.s_)
            # Hastie et al. (2009), p. 652, Eq. (18.2)
            discriminant_score[:, cur_class] = np.squeeze(
                -np.sum(Xdist_norm, axis=1) + 2.0 * np.log(self.priors_[cur_class])
            )

        return discriminant_score

    def decision_function(self, X):
        """Apply decision function to an array of samples.

        The decision function is equal (up to a constant factor) to the
        log-posterior of the model, i.e. `log p(y = k | x)`.
        However, the decision function uses the argmax definition instead
        of argmin definition as in the predict() function, which are just
        -1 factor of each other.
        The estimation has been implemented according to
        Hastie et al. (2009), p. 652 equation (18.2)

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

        dec_func = self._decision_function(X)
        # handle special case of two classes
        if len(self.classes_) == 2:
            return dec_func[:, 1] - dec_func[:, 0]
        return dec_func

    def predict_proba(self, X):
        """Estimate Class probabilities.
        The returned estimates for all classes are ordered by the
        label of classes. The estimation has been implemented according to
        Hastie et al. (2009), p. 652 equation (18.8)

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]

        Returns
        -------
        T : array-like, shape = [n_samples, n_classes]
            Returns the probability of the sample for each class in the model,
            where classes are ordered as they are in ``self.classes_``.
        """
        coef = self._decision_function(X)
        likelihood = np.exp(coef - coef.max(axis=1)[:, np.newaxis])
        return likelihood / likelihood.sum(axis=1)[:, np.newaxis]

    def predict_log_proba(self, X):
        """Estimate log probability.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            Input data.

        Returns
        -------
        C : ndarray of shape (n_samples, n_classes)
            Estimated log probabilities.
        """
        prediction = self.predict_proba(X)
        prediction[prediction == 0.0] += np.finfo(prediction.dtype).tiny
        return np.log(prediction)
