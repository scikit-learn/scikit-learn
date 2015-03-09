# coding: utf-8

# Author: Johannes Schönberger
#
# License: BSD 3 clause

import numpy as np

from ..base import BaseEstimator, MetaEstimatorMixin, RegressorMixin, clone
from ..utils import check_random_state, check_array, check_consistent_length
from ..utils.random import sample_without_replacement
from ..utils.validation import check_is_fitted
from .base import LinearRegression


_EPSILON = np.spacing(1)


def _dynamic_max_trials(n_inliers, n_samples, min_samples, probability):
    """Determine number trials such that at least one outlier-free subset is
    sampled for the given inlier/outlier ratio.

    Parameters
    ----------
    n_inliers : int
        Number of inliers in the data.

    n_samples : int
        Total number of samples in the data.

    min_samples : int
        Minimum number of samples chosen randomly from original data.

    probability : float
        Probability (confidence) that one outlier-free sample is generated.

    Returns
    -------
    trials : int
        Number of trials.

    """
    inlier_ratio = n_inliers / float(n_samples)
    nom = max(_EPSILON, 1 - probability)
    denom = max(_EPSILON, 1 - inlier_ratio ** min_samples)
    if nom == 1:
        return 0
    if denom == 1:
        return float('inf')
    return abs(float(np.ceil(np.log(nom) / np.log(denom))))


class RANSACRegressor(BaseEstimator, MetaEstimatorMixin, RegressorMixin):
    """RANSAC (RANdom SAmple Consensus) algorithm.

    RANSAC is an iterative algorithm for the robust estimation of parameters
    from a subset of inliers from the complete data set. More information can
    be found in the general documentation of linear models.

    A detailed description of the algorithm can be found in the documentation
    of the ``linear_model`` sub-package.

    Parameters
    ----------
    base_estimator : object, optional
        Base estimator object which implements the following methods:

         * `fit(X, y)`: Fit model to given training data and target values.
         * `score(X, y)`: Returns the mean accuracy on the given test data,
           which is used for the stop criterion defined by `stop_score`.
           Additionally, the score is used to decide which of two equally
           large consensus sets is chosen as the better one.

        If `base_estimator` is None, then
        ``base_estimator=sklearn.linear_model.LinearRegression()`` is used for
        target values of dtype float.

        Note that the current implementation only supports regression
        estimators.

    min_samples : int (>= 1) or float ([0, 1]), optional
        Minimum number of samples chosen randomly from original data. Treated
        as an absolute number of samples for `min_samples >= 1`, treated as a
        relative number `ceil(min_samples * X.shape[0]`) for
        `min_samples < 1`. This is typically chosen as the minimal number of
        samples necessary to estimate the given `base_estimator`. By default a
        ``sklearn.linear_model.LinearRegression()`` estimator is assumed and
        `min_samples` is chosen as ``X.shape[1] + 1``.

    residual_threshold : float, optional
        Maximum residual for a data sample to be classified as an inlier.
        By default the threshold is chosen as the MAD (median absolute
        deviation) of the target values `y`.

    is_data_valid : callable, optional
        This function is called with the randomly selected data before the
        model is fitted to it: `is_data_valid(X, y)`. If its return value is
        False the current randomly chosen sub-sample is skipped.

    is_model_valid : callable, optional
        This function is called with the estimated model and the randomly
        selected data: `is_model_valid(model, X, y)`. If its return value is
        False the current randomly chosen sub-sample is skipped.
        Rejecting samples with this function is computationally costlier than
        with `is_data_valid`. `is_model_valid` should therefore only be used if
        the estimated model is needed for making the rejection decision.

    max_trials : int, optional
        Maximum number of iterations for random sample selection.

    stop_n_inliers : int, optional
        Stop iteration if at least this number of inliers are found.

    stop_score : float, optional
        Stop iteration if score is greater equal than this threshold.

    stop_probability : float in range [0, 1], optional
        RANSAC iteration stops if at least one outlier-free set of the training
        data is sampled in RANSAC. This requires to generate at least N
        samples (iterations)::

            N >= log(1 - probability) / log(1 - e**m)

        where the probability (confidence) is typically set to high value such
        as 0.99 (the default) and e is the current fraction of inliers w.r.t.
        the total number of samples.

    residual_metric : callable, optional
        Metric to reduce the dimensionality of the residuals to 1 for
        multi-dimensional target values ``y.shape[1] > 1``. By default the sum
        of absolute differences is used::

            lambda dy: np.sum(np.abs(dy), axis=1)

    random_state : integer or numpy.RandomState, optional
        The generator used to initialize the centers. If an integer is
        given, it fixes the seed. Defaults to the global numpy random
        number generator.

    Attributes
    ----------
    estimator_ : object
        Best fitted model (copy of the `base_estimator` object).

    n_trials_ : int
        Number of random selection trials until one of the stop criteria is
        met. It is always ``<= max_trials``.

    inlier_mask_ : bool array of shape [n_samples]
        Boolean mask of inliers classified as ``True``.

    References
    ----------
    .. [1] http://en.wikipedia.org/wiki/RANSAC
    .. [2] http://www.cs.columbia.edu/~belhumeur/courses/compPhoto/ransac.pdf
    .. [3] http://www.bmva.org/bmvc/2009/Papers/Paper355/Paper355.pdf
    """

    def __init__(self, base_estimator=None, min_samples=None,
                 residual_threshold=None, is_data_valid=None,
                 is_model_valid=None, max_trials=100,
                 stop_n_inliers=np.inf, stop_score=np.inf,
                 stop_probability=0.99, residual_metric=None,
                 random_state=None):

        self.base_estimator = base_estimator
        self.min_samples = min_samples
        self.residual_threshold = residual_threshold
        self.is_data_valid = is_data_valid
        self.is_model_valid = is_model_valid
        self.max_trials = max_trials
        self.stop_n_inliers = stop_n_inliers
        self.stop_score = stop_score
        self.stop_probability = stop_probability
        self.residual_metric = residual_metric
        self.random_state = random_state

    def fit(self, X, y):
        """Fit estimator using RANSAC algorithm.

        Parameters
        ----------
        X : array-like or sparse matrix, shape [n_samples, n_features]
            Training data.

        y : array-like, shape = [n_samples] or [n_samples, n_targets]
            Target values.

        Raises
        ------
        ValueError
            If no valid consensus set could be found. This occurs if
            `is_data_valid` and `is_model_valid` return False for all
            `max_trials` randomly chosen sub-samples.

        """
        X = check_array(X, accept_sparse='csr')
        y = check_array(y, ensure_2d=False)
        if y.ndim == 1:
            y = y.reshape(-1, 1)
        check_consistent_length(X, y)

        if self.base_estimator is not None:
            base_estimator = clone(self.base_estimator)
        else:
            base_estimator = LinearRegression()

        if self.min_samples is None:
            # assume linear model by default
            min_samples = X.shape[1] + 1
        elif 0 < self.min_samples < 1:
            min_samples = np.ceil(self.min_samples * X.shape[0])
        elif self.min_samples >= 1:
            if self.min_samples % 1 != 0:
                raise ValueError("Absolute number of samples must be an "
                                 "integer value.")
            min_samples = self.min_samples
        else:
            raise ValueError("Value for `min_samples` must be scalar and "
                             "positive.")
        if min_samples > X.shape[0]:
            raise ValueError("`min_samples` may not be larger than number "
                             "of samples ``X.shape[0]``.")

        if self.stop_probability < 0 or self.stop_probability > 1:
            raise ValueError("`stop_probability` must be in range [0, 1].")

        if self.residual_threshold is None:
            # MAD (median absolute deviation)
            residual_threshold = np.median(np.abs(y - np.median(y)))
        else:
            residual_threshold = self.residual_threshold

        if self.residual_metric is None:
            residual_metric = lambda dy: np.sum(np.abs(dy), axis=1)
        else:
            residual_metric = self.residual_metric

        random_state = check_random_state(self.random_state)

        try:  # Not all estimator accept a random_state
            base_estimator.set_params(random_state=random_state)
        except ValueError:
            pass

        n_inliers_best = 0
        score_best = np.inf
        inlier_mask_best = None
        X_inlier_best = None
        y_inlier_best = None

        # number of data samples
        n_samples = X.shape[0]
        sample_idxs = np.arange(n_samples)

        n_samples, _ = X.shape

        for self.n_trials_ in range(1, self.max_trials + 1):

            # choose random sample set
            subset_idxs = sample_without_replacement(n_samples, min_samples,
                                                     random_state=random_state)
            X_subset = X[subset_idxs]
            y_subset = y[subset_idxs]

            # check if random sample set is valid
            if (self.is_data_valid is not None
                    and not self.is_data_valid(X_subset, y_subset)):
                continue

            # fit model for current random sample set
            base_estimator.fit(X_subset, y_subset)

            # check if estimated model is valid
            if (self.is_model_valid is not None and not
                    self.is_model_valid(base_estimator, X_subset, y_subset)):
                continue

            # residuals of all data for current random sample model
            y_pred = base_estimator.predict(X)
            if y_pred.ndim == 1:
                y_pred = y_pred[:, None]

            residuals_subset = residual_metric(y_pred - y)

            # classify data into inliers and outliers
            inlier_mask_subset = residuals_subset < residual_threshold
            n_inliers_subset = np.sum(inlier_mask_subset)

            # less inliers -> skip current random sample
            if n_inliers_subset < n_inliers_best:
                continue

            # extract inlier data set
            inlier_idxs_subset = sample_idxs[inlier_mask_subset]
            X_inlier_subset = X[inlier_idxs_subset]
            y_inlier_subset = y[inlier_idxs_subset]

            # score of inlier data set
            score_subset = base_estimator.score(X_inlier_subset,
                                                y_inlier_subset)

            # same number of inliers but worse score -> skip current random
            # sample
            if (n_inliers_subset == n_inliers_best
                    and score_subset < score_best):
                continue

            # save current random sample as best sample
            n_inliers_best = n_inliers_subset
            score_best = score_subset
            inlier_mask_best = inlier_mask_subset
            X_inlier_best = X_inlier_subset
            y_inlier_best = y_inlier_subset

            # break if sufficient number of inliers or score is reached
            if (n_inliers_best >= self.stop_n_inliers
                    or score_best >= self.stop_score
                    or self.n_trials_
                       >= _dynamic_max_trials(n_inliers_best, n_samples,
                                              min_samples,
                                              self.stop_probability)):
                break

        # if none of the iterations met the required criteria
        if inlier_mask_best is None:
            raise ValueError(
                "RANSAC could not find valid consensus set, because"
                " either the `residual_threshold` rejected all the samples or"
                " `is_data_valid` and `is_model_valid` returned False for all"
                " `max_trials` randomly ""chosen sub-samples. Consider "
                "relaxing the ""constraints.")

        # estimate final model using all inliers
        base_estimator.fit(X_inlier_best, y_inlier_best)

        self.estimator_ = base_estimator
        self.inlier_mask_ = inlier_mask_best
        return self

    def predict(self, X):
        """Predict using the estimated model.

        This is a wrapper for `estimator_.predict(X)`.

        Parameters
        ----------
        X : numpy array of shape [n_samples, n_features]

        Returns
        -------
        y : array, shape = [n_samples] or [n_samples, n_targets]
            Returns predicted values.
        """
        check_is_fitted(self, 'estimator_')

        return self.estimator_.predict(X)

    def score(self, X, y):
        """Returns the score of the prediction.

        This is a wrapper for `estimator_.score(X, y)`.

        Parameters
        ----------
        X : numpy array or sparse matrix of shape [n_samples, n_features]
            Training data.

        y : array, shape = [n_samples] or [n_samples, n_targets]
            Target values.

        Returns
        -------
        z : float
            Score of the prediction.
        """
        check_is_fitted(self, 'estimator_')

        return self.estimator_.score(X, y)
