# coding: utf-8

# Author: Johannes Schönberger
#
# License: BSD 3 clause

import numpy as np

from ..base import BaseEstimator, clone
from ..utils import check_random_state, atleast2d_or_csr
from .base import LinearRegression


class RANSAC(BaseEstimator):
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
         * `score(X)`: Returns the mean accuracy on the given test data, which
           is used for the stop criterion defined by `stop_score`.
           Additionally, the score is used to decide which of two equally
           large consensus sets is chosen as the better one.

        If no base estimator is specified, by default
        ``sklearn.linear_model.LinearRegression`` is used for target values of
        dtype float.

        Note that the current implementation only supports regression
        estimators.

    min_n_samples : int (>= 1) or float ([0, 1]), optional
        Minimum number of samples chosen randomly from original data. Treated
        as an absolute number of samples for `min_n_samples >= 1`, treated as a
        relative number `ceil(min_n_samples * X.shape[0]`) for
        `min_n_samples < 1`. By default a
        ``sklearn.linear_model.LinearRegression`` estimator is assumed and
        `min_n_samples` is chosen as ``X.shape[1] + 1``.

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

    max_trials : int, optional
        Maximum number of iterations for random sample selection.

    stop_n_inliers : int, optional
        Stop iteration if at least this number of inliers are found.

    stop_score : float, optional
        Stop iteration if score is greater equal than this threshold.

    residual_metric : callable, optional
        Metric to reduce the dimensionality of the residuals to 1 for
        multi-dimensional target values ``y.shape[1] > 1``. By default the sum
        of absolute differences is used::

            lambda dy: np.sum(np.abs(dx), axis=1)

    random_state : integer or numpy.RandomState, optional
        The generator used to initialize the centers. If an integer is
        given, it fixes the seed. Defaults to the global numpy random
        number generator.

    Attributes
    ----------
    estimator_ : object
        Best fitted model (copy of the `base_estimator` object).

    n_trials_ : int
        Number of random selection trials.

    inlier_mask_ : bool array of shape [n_samples]
        Boolean mask of inliers classified as ``True``.

    References
    ----------
    .. [1] http://en.wikipedia.org/wiki/RANSAC
    .. [2] http://www.cs.columbia.edu/~belhumeur/courses/compPhoto/ransac.pdf
    .. [3] http://www.bmva.org/bmvc/2009/Papers/Paper355/Paper355.pdf
    """

    def __init__(self, base_estimator=None, min_n_samples=None,
                 residual_threshold=None, is_data_valid=None,
                 is_model_valid=None, max_trials=100,
                 stop_n_inliers=np.inf, stop_score=np.inf,
                 residual_metric=None, random_state=None):

        self.base_estimator = base_estimator
        self.min_n_samples = min_n_samples
        self.residual_threshold = residual_threshold
        self.is_data_valid = is_data_valid
        self.is_model_valid = is_model_valid
        self.max_trials = max_trials
        self.stop_n_inliers = stop_n_inliers
        self.stop_score = stop_score
        self.residual_metric = residual_metric
        self.random_state = random_state

    def fit(self, X, y):
        """Fit estimator using RANSAC algorithm.

        Parameters
        ----------
        X : numpy array or sparse matrix of shape [n_samples, n_features]
            Training data.

        y : numpy array of shape [n_samples, n_targets]
            Target values.

        Raises
        ------
        ValueError
            If no valid consensus set could be found. This occurs if
            `is_data_valid` and `is_model_valid` return False for all
            `max_trials` randomly chosen sub-samples.

        """
        if self.base_estimator is not None:
            base_estimator = clone(self.base_estimator)
        elif y.dtype.kind == 'f':
            base_estimator = LinearRegression()
        else:
            raise ValueError("`base_estimator` not specified.")

        if self.min_n_samples is None:
            # assume linear model by default
            min_n_samples = X.shape[1] + 1
        elif 0 < self.min_n_samples < 1:
            min_n_samples = np.ceil(self.min_n_samples * X.shape[0])
        elif self.min_n_samples >= 1:
            min_n_samples = self.min_n_samples
        else:
            raise ValueError("Value for `min_n_samples` must be scalar and "
                             "positive.")

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

        best_n_inliers = 0
        best_score = np.inf
        best_inlier_mask = None
        best_inlier_X = None
        best_inlier_y = None

        # number of data samples
        n_samples = X.shape[0]
        sample_idxs = np.arange(n_samples)

        X = atleast2d_or_csr(X)
        y = np.asarray(y)

        if y.ndim == 1:
            y = y[:, None]

        for n_trials in range(self.max_trials):

            # choose random sample set
            random_idxs = random_state.randint(0, n_samples, min_n_samples)
            rsample_X = X[random_idxs]
            rsample_y = y[random_idxs]

            # check if random sample set is valid
            if (self.is_data_valid is not None
                    and not self.is_data_valid(rsample_X, rsample_y)):
                continue

            # fit model for current random sample set
            base_estimator.fit(rsample_X, rsample_y)

            # check if estimated model is valid
            if (self.is_model_valid is not None and not
                    self.is_model_valid(base_estimator, rsample_X, rsample_y)):
                continue

            # residuals of all data for current random sample model
            rsample_residuals = residual_metric(base_estimator.predict(X) - y)

            # classify data into inliers and outliers
            rsample_inlier_mask = rsample_residuals < residual_threshold
            rsample_n_inliers = np.sum(rsample_inlier_mask)

            # less inliers -> skip current random sample
            if rsample_n_inliers < best_n_inliers:
                continue

            # extract inlier data set
            rsample_inlier_idxs = sample_idxs[rsample_inlier_mask]
            rsample_inlier_X = X[rsample_inlier_idxs]
            rsample_inlier_y = y[rsample_inlier_idxs]

            # score of inlier data set
            rsample_score = base_estimator.score(rsample_inlier_X,
                                                 rsample_inlier_y)

            # same number of inliers but worse score -> skip current random
            # sample
            if (rsample_n_inliers == best_n_inliers
                    and rsample_score < best_score):
                continue

            # save current random sample as best sample
            best_n_inliers = rsample_n_inliers
            best_score = rsample_score
            best_inlier_mask = rsample_inlier_mask
            best_inlier_X = rsample_inlier_X
            best_inlier_y = rsample_inlier_y

            # break if sufficient number of inliers or score is reached
            if (best_n_inliers >= self.stop_n_inliers
                    or best_score >= self.stop_score):
                break

        # if none of the iterations met the required criteria
        if best_inlier_mask is None:
            raise ValueError("RANSAC could not find valid consensus set, "
                             "because `is_data_valid` and `is_model_valid` "
                             "returned False for all `max_trials` randomly "
                             "chosen sub-samples. Consider relaxing the "
                             "constraints.")

        # estimate final model using all inliers
        base_estimator.fit(best_inlier_X, best_inlier_y)

        self.estimator_ = base_estimator
        self.n_trials_ = n_trials + 1
        self.inlier_mask_ = best_inlier_mask

    def predict(self, X):
        """Predict using the estimated model.

        This is a wrapper for `estimator_.predict(X)`.

        Parameters
        ----------
        X : numpy array of shape [n_samples, n_features]

        Returns
        -------
        C : array, shape = [n_samples]
            Returns predicted values.
        """
        return self.estimator_.predict(X)

    def score(self, X, y):
        """Returns the score of the prediction.

        This is a wrapper for `estimator_.score(X, y)`.

        Parameters
        ----------
        X : numpy array or sparse matrix of shape [n_samples, n_features]
            Training data.

        y : numpy array of shape [n_samples, n_targets]
            Target values.

        Returns
        -------
        z : float
            Score of the prediction.
        """
        return self.estimator_.score(X, y)
