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
# License: BSD Style.

import numpy as np
import  scipy as sp
from . import MinCovDet
from ..utils import deprecated
from ..base import ClassifierMixin


class OutlierDetectionMixin(object):
    """Set of methods for outliers detection with covariance estimators.

    Parameters
    ----------
    contamination: float, 0. < contamination < 0.5
        The amount of contamination of the data set, i.e. the proportion
        of outliers in the data set.

    Notes
    -----
    Outlier detection from covariance estimation may break or not
    perform well in high-dimensional settings. In particular, one will
    always take care to work with ``n_samples > n_features ** 2``.

    """
    def __init__(self, contamination=0.1):
        self.contamination = contamination
        self.threshold = None

    def decision_function(self, X, raw_values=False):
        """Compute the decision function of the given observations.

        Parameters
        ----------
        X: array-like, shape (n_samples, n_features)

        raw_values: bool
            Whether or not to consider raw Mahalanobis distances as the
            decision function. Must be False (default) for compatibility
            with the others outlier detection tools.

        Returns
        -------
        decision: array-like, shape (n_samples, )
            The values of the decision function for each observations.
            It is equal to the Mahalanobis distances if `raw_values`
            is True. By default (``raw_values=True``), it is equal
            to the cubic root of the shifted Mahalanobis distances.
            In that case, the threshold for being an outlier is 0, which
            ensures a compatibility with other outlier detection tools
            such as the One-Class SVM.

        """
        mahal_dist = self.mahalanobis(X)
        if raw_values:
            decision = mahal_dist
        else:
            if self.threshold is None:
                raise Exception("Please fit data before predicting")
            transformed_mahal_dist = mahal_dist ** 0.33
            decision = self.threshold ** 0.33 - transformed_mahal_dist

        return decision

    def predict(self, X):
        """Outlyingness of observations in X according to the fitted model.

        Parameters
        ----------
        X: array-like, shape = (n_samples, n_features)

        Returns
        -------
        is_outliers: array, shape = (n_samples, ), dtype = bool
            For each observations, tells whether or not it should be considered
            as an outlier according to the fitted model.
        threshold: float,
            The values of the less outlying point's decision function.

        """
        if self.threshold is None:
            raise Exception("Please fit data before predicting")
        is_inlier = -np.ones(X.shape[0], dtype=int)
        if self.contamination is not None:
            values = self.decision_function(X, raw_values=True)
            is_inlier[values <= self.threshold] = 1
        else:
            raise NotImplemented("You must provide a contamination rate.")

        return is_inlier


class EllipticEnvelope(ClassifierMixin, OutlierDetectionMixin, MinCovDet):
    """An object for detecting outliers in a Gaussian distributed dataset.

    Attributes
    ----------
    `contamination`: float, 0. < contamination < 0.5
      The amount of contamination of the data set, i.e. the proportion of \
      outliers in the data set.

    `location_`: array-like, shape (n_features,)
        Estimated robust location

    `covariance_`: array-like, shape (n_features, n_features)
        Estimated robust covariance matrix

    `precision_`: array-like, shape (n_features, n_features)
        Estimated pseudo inverse matrix.
        (stored only if store_precision is True)

    `support_`: array-like, shape (n_samples,)
        A mask of the observations that have been used to compute the
        robust estimates of location and shape.

    Parameters
    ----------
    store_precision: bool
        Specify if the estimated precision is stored
    assume_centered: Boolean
        If True, the support of robust location and covariance estimates
        is computed, and a covariance estimate is recomputed from it,
        without centering the data.
        Useful to work with data whose mean is significantly equal to
        zero but is not exactly zero.
        If False, the robust location and covariance are directly computed
        with the FastMCD algorithm without additional treatment.
    support_fraction: float, 0 < support_fraction < 1
        The proportion of points to be included in the support of the raw
        MCD estimate. Default is ``None``, which implies that the minimum
        value of support_fraction will be used within the algorithm:
        [n_sample + n_features + 1] / 2
    contamination: float, 0. < contamination < 0.5
        The amount of contamination of the data set, i.e. the proportion
        of outliers in the data set.

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
                 support_fraction=None, contamination=0.1):
        MinCovDet.__init__(self, store_precision=store_precision,
                           assume_centered=assume_centered,
                           support_fraction=support_fraction)
        OutlierDetectionMixin.__init__(self, contamination=contamination)

    def fit(self, X):
        """
        """
        MinCovDet.fit(self, X)
        self.threshold = sp.stats.scoreatpercentile(
            self.dist_, 100. * (1. - self.contamination))

        return self


# Deprecated classes
@deprecated("Use EllipticEnvelope instead. To be removed in 0.13.")
class EllipticEnvelop(EllipticEnvelope):
    pass
