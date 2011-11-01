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
from sklearn.covariance import MinCovDet


class OutlierDetection():
    """Set of methods for outliers detection with covariance estimators.

    """
    def __init__(self, contamination=0.1):
        """

        Parameters
        ----------
        contamination: float, 0. < contamination < 0.5
          The amount of contamination of the data set, i.e. the proportion
          of outliers in the data set.

        """
        self.contamination = contamination

    def decision_function(self, X):
        """Compute the decision function of the given observations.

        Parameters
        ----------
        X: array-like, shape = (n_samples, n_features)

        Returns
        -------
        mahal_dist: array-like, shape = (n_samples, )
          The values of the decision function for each observations.

        """
        if self.covariance_ is None:
            raise Exception("No decision function defined. " \
                                "Please fit some data first")
        X_centered = X - self.location_
        mahal_dist = self.mahalanobis(X_centered)

        return mahal_dist

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
        is_outlier = np.zeros(X.shape[0], dtype=bool)
        if self.contamination is not None:
            X_centered = X - self.location_
            values = self.decision_function(X_centered)
            threshold = sp.stats.scoreatpercentile(
                values, 100. * (1. - self.contamination))
            is_outlier[values > threshold] = True
        else:
            raise NotImplemented("You must provide a contamination rate.")

        return is_outlier, threshold


class EllipticData(OutlierDetection, MinCovDet):
    """An object for detecting outliers in a Gaussian distributed dataset.

    Attributes
    ----------
    `contamination`: float, 0. < contamination < 0.5
      The amount of contamination of the data set, i.e. the proportion
      of outliers in the data set.

    `location_`: array-like, shape (n_features,)
        Estimated robust location

    `covariance_`: array-like, shape (n_features, n_features)
        Estimated robust covariance matrix

    `precision_`: array-like, shape (n_features, n_features)
        Estimated pseudo inverse matrix.
        (stored only if store_precision is True)

    `support_`: array-like, shape (n_samples,)
        A mask of the observations that have been used to compute
        the robust estimates of location and shape.

    """
    def __init__(self, store_precision=True, assume_centered=False,
                 h=None, correction="empirical", reweighting=None,
                 contamination=0.1):
        """

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
        h: float, 0 < h < 1
          The proportion of points to be included in the support of the raw
          MCD estimate. Default is None, which implies that the minimum
          value of h will be used within the algorithm:
          [n_sample + n_features + 1] / 2
        correction: str
          Improve the covariance estimator consistency at Gaussian models
            - "empirical" (default): correction using the empirical correction
              factor suggested by Rousseeuw and Van Driessen in [1]
            - "theoretical": correction using the theoretical correction factor
              derived in [2]
            - else: no correction
        reweighting: str
          Computation of a reweighted estimator:
            - "rousseeuw" (default): Reweight observations using Rousseeuw's
              method (equivalent to deleting outlying observations from the
              data set before computing location and covariance estimates)
            - else: no re-weighting
        contamination: float, 0. < contamination < 0.5
          The amount of contamination of the data set, i.e. the proportion
          of outliers in the data set.

        """
        MinCovDet.__init__(self, store_precision=store_precision,
                           assume_centered=assume_centered, h=h,
                           correction=correction, reweighting=reweighting)
        OutlierDetection.__init__(self, contamination=contamination)
