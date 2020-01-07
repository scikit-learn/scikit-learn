# Author: Virgile Fritsch <virgile.fritsch@inria.fr>
#
# License: BSD 3 clause

import numpy as np
from . import MinCovDet
from ..utils.validation import check_is_fitted, check_array
from ..metrics import accuracy_score
from ..base import OutlierMixin


class EllipticEnvelope(OutlierMixin, MinCovDet):
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

    offset_ : float
        Offset used to define the decision function from the raw scores.
        We have the relation: ``decision_function = score_samples - offset_``.
        The offset depends on the contamination parameter and is defined in
        such a way we obtain the expected number of outliers (samples with
        decision function < 0) in training.

    Examples
    --------
    >>> import numpy as np
    >>> from sklearn.covariance import EllipticEnvelope
    >>> true_cov = np.array([[.8, .3],
    ...                      [.3, .4]])
    >>> X = np.random.RandomState(0).multivariate_normal(mean=[0, 0],
    ...                                                  cov=true_cov,
    ...                                                  size=500)
    >>> cov = EllipticEnvelope(random_state=0).fit(X)
    >>> # predict returns 1 for an inlier and -1 for an outlier
    >>> cov.predict([[0, 0],
    ...              [3, 3]])
    array([ 1, -1])
    >>> cov.covariance_
    array([[0.7411..., 0.2535...],
           [0.2535..., 0.3053...]])
    >>> cov.location_
    array([0.0813... , 0.0427...])

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
    .. [1] Rousseeuw, P.J., Van Driessen, K. "A fast algorithm for the
       minimum covariance determinant estimator" Technometrics 41(3), 212
       (1999)

    """
    def __init__(self, store_precision=True, assume_centered=False,
                 support_fraction=None, contamination=0.1,
                 random_state=None):
        super().__init__(
            store_precision=store_precision,
            assume_centered=assume_centered,
            support_fraction=support_fraction,
            random_state=random_state)
        self.contamination = contamination

    def fit(self, X, y=None):
        """Fit the EllipticEnvelope model.

        Parameters
        ----------
        X : numpy array or sparse matrix, shape (n_samples, n_features).
            Training data

        y : Ignored
            not used, present for API consistency by convention.

        """
        super().fit(X)
        self.offset_ = np.percentile(-self.dist_, 100. * self.contamination)
        return self

    def decision_function(self, X):
        """Compute the decision function of the given observations.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)

        Returns
        -------

        decision : array-like, shape (n_samples, )
            Decision function of the samples.
            It is equal to the shifted Mahalanobis distances.
            The threshold for being an outlier is 0, which ensures a
            compatibility with other outlier detection algorithms.

        """
        check_is_fitted(self)
        negative_mahal_dist = self.score_samples(X)
        return negative_mahal_dist - self.offset_

    def score_samples(self, X):
        """Compute the negative Mahalanobis distances.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)

        Returns
        -------
        negative_mahal_distances : array-like, shape (n_samples, )
            Opposite of the Mahalanobis distances.
        """
        check_is_fitted(self)
        return -self.mahalanobis(X)

    def predict(self, X):
        """
        Predict the labels (1 inlier, -1 outlier) of X according to the
        fitted model.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)

        Returns
        -------
        is_inlier : array, shape (n_samples,)
            Returns -1 for anomalies/outliers and +1 for inliers.
        """
        X = check_array(X)
        is_inlier = np.full(X.shape[0], -1, dtype=int)
        values = self.decision_function(X)
        is_inlier[values >= 0] = 1

        return is_inlier

    def score(self, X, y, sample_weight=None):
        """Returns the mean accuracy on the given test data and labels.

        In multi-label classification, this is the subset accuracy
        which is a harsh metric since you require for each sample that
        each label set be correctly predicted.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            Test samples.

        y : array-like, shape (n_samples,) or (n_samples, n_outputs)
            True labels for X.

        sample_weight : array-like, shape (n_samples,), optional
            Sample weights.

        Returns
        -------
        score : float
            Mean accuracy of self.predict(X) wrt. y.

        """
        return accuracy_score(y, self.predict(X), sample_weight=sample_weight)
