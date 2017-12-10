"""
Class for outlier detection.

This class provides a framework for outlier detection. It consists in
several methods that can be added to a covariance estimator in order to
assess the outlying-ness of the observations of a data set.
Such a "outlier detector" object is proposed constructed from a robust
covariance estimator (the Minimum Covariance Determinant).

"""
# Author: Virgile Fritsch <virgile.fritsch@inria.fr>
#
# License: BSD 3 clause

import numpy as np
import scipy as sp
from . import MinCovDet
from ..utils.validation import check_is_fitted, check_array
from ..metrics import accuracy_score


class EllipticEnvelope(MinCovDet):
    """An object for detecting outliers in a Gaussian distributed dataset.

    Read more in the :ref:`User Guide <outlier_detection>`.

    Parameters
    ----------
    store_precision : boolean, optional (default=True)
        Specify if the estimated precision is stored.

    assume_centered : boolean, optional (default=False)
        If True, the support of robust location and covariance estimates
        is computed, and a covariance estimate is recomputed from it,
        without centering the data.
        Useful to work with data whose mean is significantly equal to
        zero but is not exactly zero.
        If False, the robust location and covariance are directly computed
        with the FastMCD algorithm without additional treatment.

    support_fraction : float in (0., 1.), optional (default=None)
        The proportion of points to be included in the support of the raw
        MCD estimate. If None, the minimum value of support_fraction will
        be used within the algorithm: `[n_sample + n_features + 1] / 2`.

    contamination : float in (0., 0.5), optional (default=0.1)
        The amount of contamination of the data set, i.e. the proportion
        of outliers in the data set.

    random_state : int, RandomState instance or None, optional (default=None)
        The seed of the pseudo random number generator to use when shuffling
        the data.  If int, random_state is the seed used by the random number
        generator; If RandomState instance, random_state is the random number
        generator; If None, the random number generator is the RandomState
        instance used by `np.random`.

    Attributes
    ----------
    location_ : array-like, shape (n_features,)
        Estimated robust location

    covariance_ : array-like, shape (n_features, n_features)
        Estimated robust covariance matrix

    precision_ : array-like, shape (n_features, n_features)
        Estimated pseudo inverse matrix.
        (stored only if store_precision is True)

    support_ : array-like, shape (n_samples,)
        A mask of the observations that have been used to compute the
        robust estimates of location and shape.

    See Also
    --------
    EmpiricalCovariance, MinCovDet

    Notes
    -----
    Outlier detection from covariance estimation may break or not
    perform well in high-dimensional settings. In particular, one will
    always take care to work with ``n_samples > n_features ** 2``.

    References
    ----------
    ..  [1] Rousseeuw, P.J., Van Driessen, K. "A fast algorithm for the minimum
        covariance determinant estimator" Technometrics 41(3), 212 (1999)

    """
    def __init__(self, store_precision=True, assume_centered=False,
                 support_fraction=None, contamination=0.1,
                 random_state=None):
        super(EllipticEnvelope, self).__init__(
            store_precision=store_precision,
            assume_centered=assume_centered,
            support_fraction=support_fraction,
            random_state=random_state)
        self.contamination = contamination

    def fit(self, X, y=None):
        """Fit the EllipticEnvelope model with X.

        Parameters
        ----------
        X : numpy array or sparse matrix of shape [n_samples, n_features]
            Training data
        y : (ignored)
        """
        super(EllipticEnvelope, self).fit(X)
        self.threshold_ = sp.stats.scoreatpercentile(
            self.dist_, 100. * (1. - self.contamination))
        return self

    def decision_function(self, X, raw_values=False):
        """Compute the decision function of the given observations.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)

        raw_values : bool
            Whether or not to consider raw Mahalanobis distances as the
            decision function. Must be False (default) for compatibility
            with the others outlier detection tools.

        Returns
        -------
        decision : array-like, shape (n_samples, )
            Decision function of the samples.
            It is equal to the Mahalanobis distances if `raw_values`
            is True. By default (``raw_values=False``), it is equal
            to the cubic root of the shifted Mahalanobis distances.
            In that case, the threshold for being an outlier is 0, which
            ensures a compatibility with other outlier detection tools
            such as the One-Class SVM.

        """
        check_is_fitted(self, 'threshold_')
        X = check_array(X)
        mahal_dist = self.mahalanobis(X)
        if raw_values:
            decision = mahal_dist
        else:
            transformed_mahal_dist = mahal_dist ** 0.33
            decision = self.threshold_ ** 0.33 - transformed_mahal_dist

        return decision

    def predict(self, X):
        """Outlyingness of observations in X according to the fitted model.

        Parameters
        ----------
        X : array-like, shape = (n_samples, n_features)

        Returns
        -------
        is_outliers : array, shape = (n_samples, ), dtype = bool
            For each observation, tells whether or not it should be considered
            as an outlier according to the fitted model.

        threshold : float,
            The values of the less outlying point's decision function.

        """
        check_is_fitted(self, 'threshold_')
        X = check_array(X)
        is_inlier = -np.ones(X.shape[0], dtype=int)
        if self.contamination is not None:
            values = self.decision_function(X, raw_values=True)
            is_inlier[values <= self.threshold_] = 1
        else:
            raise NotImplementedError("You must provide a contamination rate.")

        return is_inlier

    def score(self, X, y, sample_weight=None):
        """Returns the mean accuracy on the given test data and labels.

        In multi-label classification, this is the subset accuracy
        which is a harsh metric since you require for each sample that
        each label set be correctly predicted.

        Parameters
        ----------
        X : array-like, shape = (n_samples, n_features)
            Test samples.

        y : array-like, shape = (n_samples,) or (n_samples, n_outputs)
            True labels for X.

        sample_weight : array-like, shape = (n_samples,), optional
            Sample weights.

        Returns
        -------
        score : float
            Mean accuracy of self.predict(X) wrt. y.

        """
        return accuracy_score(y, self.predict(X), sample_weight=sample_weight)
