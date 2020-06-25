# -*- coding: utf-8 -*-
"""
Nearest Centroid Classification
"""

# Author: Robert Layton <robertlayton@gmail.com>
#         Olivier Grisel <olivier.grisel@ensta.org>
#         Andreas W. Kempa-Liehr <a.kempa-liehr@auckland.ac.nz>
#
# License: BSD 3 clause

import warnings
import numpy as np
from scipy import sparse as sp

from ..base import BaseEstimator, ClassifierMixin
from ..metrics.pairwise import pairwise_distances
from ..preprocessing import LabelEncoder
from ..utils.validation import check_array, check_is_fitted
from ..utils.validation import _deprecate_positional_args
from ..utils.sparsefuncs import csc_median_axis_0
from ..utils.multiclass import check_classification_targets


class NearestCentroid(ClassifierMixin, BaseEstimator):
    """Nearest centroid classifier.

    Each class is represented by its centroid, with test samples classified to
    the class with the nearest centroid.

    Read more in the :ref:`User Guide <nearest_centroid_classifier>`.

    Parameters
    ----------
    metric : str or callable
        The metric to use when calculating distance between instances in a
        feature array. If metric is a string or callable, it must be one of
        the options allowed by metrics.pairwise.pairwise_distances for its
        metric parameter.
        The centroids for the samples corresponding to each class is the point
        from which the sum of the distances (according to the metric) of all
        samples that belong to that particular class are minimized.
        If the "manhattan" metric is provided, this centroid is the median and
        for all other metrics, the centroid is now set to be the mean.

        .. versionchanged:: 0.19
            ``metric='precomputed'`` was deprecated and now raises an error

    shrink_threshold : float, default=None
        Threshold for shrinking centroids to remove features.

    Attributes
    ----------
    centroids_ : array-like of shape (n_classes, n_features)
        Centroid of each class.

    classes_ : array of shape (n_classes,)
        The unique classes labels.

    Examples
    --------
    >>> from sklearn.neighbors import NearestCentroid
    >>> import numpy as np
    >>> X = np.array([[-1, -1], [-2, -1], [-3, -2], [1, 1], [2, 1], [3, 2]])
    >>> y = np.array([1, 1, 1, 2, 2, 2])
    >>> clf = NearestCentroid()
    >>> clf.fit(X, y)
    NearestCentroid(prior=None, prob=None)
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
    Tibshirani, R., Hastie, T., Narasimhan, B., & Chu, G. (2002a). Diagnosis of
    multiple cancer types by shrunken centroids of gene expression. Proceedings
    of the National Academy of Sciences of the United States of America,
    99(10), 6567-6572. The National Academy of Sciences.

    Tibshirani, R., Hastie, T., Narasimhan, B., & Chu, G. (2002b). Class Prediction
    by Nearest Shrunken Centroids, with Applications to DNA Microaarrays. Statistical
    Science 18(1), p. 104-117.
    """

    @_deprecate_positional_args
    def __init__(self, metric='euclidean', *, shrink_threshold=None, prob=False, prior='sample'):
        self.metric = metric
        self.shrink_threshold = shrink_threshold
        self.prob_ = prob
        if not prior in ['sample', 'uninformative']:
            raise ValueError('Prior `{}` is not implemented. Choose either `sample` or `uninformative`.'.format(prior))
        self.prior_ = prior

    def fit(self, X, y):
        """
        Fit the NearestCentroid model according to the given training data.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            Training vector, where n_samples is the number of samples and
            n_features is the number of features.
            Note that centroid shrinking cannot be used with sparse matrices.
        y : array-like of shape (n_samples,)
            Target values (integers)
        """
        if self.metric == 'precomputed':
            raise ValueError("Precomputed is not supported.")
        # If X is sparse and the metric is "manhattan", store it in a csc
        # format is easier to calculate the median.
        if self.metric == 'manhattan':
            X, y = self._validate_data(X, y, accept_sparse=['csc'])
        else:
            X, y = self._validate_data(X, y, accept_sparse=['csr', 'csc'])
        is_X_sparse = sp.issparse(X)
        if is_X_sparse and self.shrink_threshold:
            raise ValueError("threshold shrinking not supported"
                             " for sparse input")
        check_classification_targets(y)

        n_samples, n_features = X.shape
        le = LabelEncoder()
        y_ind = le.fit_transform(y)
        self.classes_ = classes = le.classes_
        self.n_classes_ = n_classes = classes.size
        if n_classes < 2:
            raise ValueError('The number of classes has to be greater than'
                             ' one; got %d class' % (n_classes))

        # Mask mapping each class to its members.
        self.centroids_ = np.empty((n_classes, n_features), dtype=np.float64)
        # Number of clusters in each class.
        nk = np.zeros(n_classes)

        for cur_class in range(n_classes):
            center_mask = y_ind == cur_class
            nk[cur_class] = np.sum(center_mask)
            if is_X_sparse:
                center_mask = np.where(center_mask)[0]

            # XXX: Update other averaging methods according to the metrics.
            if self.metric == "manhattan":
                # NumPy does not calculate median of sparse matrices.
                if not is_X_sparse:
                    self.centroids_[cur_class] = np.median(X[center_mask], axis=0)
                else:
                    self.centroids_[cur_class] = csc_median_axis_0(X[center_mask])
            else:
                if self.metric != 'euclidean':
                    warnings.warn("Averaging for metrics other than "
                                  "euclidean and manhattan not supported. "
                                  "The average is set to be the mean."
                                  )
                self.centroids_[cur_class] = X[center_mask].mean(axis=0)

        if self.shrink_threshold:
            dataset_centroid_ = np.mean(X, axis=0)

            # m parameter for determining deviation
            m = np.sqrt((1. / nk) - (1. / n_samples))
            # Calculate deviation using the standard deviation of centroids.
            variance = (X - self.centroids_[y_ind]) ** 2
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
            np.clip(deviation, 0, None, out=deviation)
            deviation *= signs
            # Now adjust the centroids using the deviation
            msd = ms * deviation
            self.centroids_ = dataset_centroid_[np.newaxis, :] + msd

        if self.prob_:
            # Calculate deviation using the standard deviation of centroids.
            variance = (X - self.centroids_[y_ind]) ** 2
            variance = variance.sum(axis=0)
            self.s_ = np.sqrt(variance / (n_samples - n_classes))
            self.sample_prior_ = nk / n_samples

        return self

    def predict(self, X):
        """Perform classification on an array of test vectors X.

        The predicted class C for each sample in X is returned.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)

        Returns
        -------
        C : ndarray of shape (n_samples,)

        Notes
        -----
        If the metric constructor parameter is "precomputed", X is assumed to
        be the distance matrix between the data to be predicted and
        ``self.centroids_``.
        """
        check_is_fitted(self)

        X = check_array(X, accept_sparse='csr')
        return self.classes_[pairwise_distances(
            X, self.centroids_, metric=self.metric).argmin(axis=1)]

    def predict_proba(self, X):
        """Probability estimates.

        The returned estimates for all classes are ordered by the
        label of classes. The estimation has been implemented according to
        Tibshirani et al. (2002b), p. 108f

        Parameters
        ----------
        X : array-like, shape = (n_samples, n_features)

        Returns
        -------
        T : ndarray of shape (n_samples, n_classes)
            Returns the probability of the sample for each class in the model,
            where classes are ordered as they are in ``self.classes_``.
        """

        check_is_fitted(self, 's_')

        # Tibshirani et al. (2002b), p. 109
        if self.prior_ == 'uninformative':
            prior = np.ones(self.n_classes_) / self.n_classes_
        else:
            prior = self.sample_prior_

        X = check_array(X, accept_sparse='csr')

        discriminant_score = np.empty((X.shape[0], self.n_classes_), dtype=np.float64)
        coef = np.empty((X.shape[0], self.n_classes_), dtype=np.float64)
        for c in range(self.n_classes_):
            Xdist = X - self.centroids_[c, :]
            Xdist_norm = (Xdist / self.s_) ** 2
            # Tibshirani et al. (2002b), p. 108, Eq. (8)
            discriminant_score[:, c] = np.sum(Xdist_norm, axis=1) - 2.0 * np.log(prior[c])
            coef[:, c] = np.exp(-0.5 * discriminant_score[:, c])
        # Tibshirani et al. (2002b), p. 109, Eq. (8)
        return (coef.transpose() / np.sum(coef, axis=1)).transpose()
