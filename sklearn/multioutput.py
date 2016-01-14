"""
This module implements multioutput regression.

The estimator provided in this module is a meta-estimator: it requires
a base estimator to be provided in their constructor. The meta-estimator
extends single output regressors to support multioutput regression.
"""

# Author: Tim Head <betatim@gmail.com>
#
# License: BSD 3 clause

import numpy as np
import scipy.sparse as sp

from .base import BaseEstimator, clone, MetaEstimatorMixin, RegressorMixin
from .utils import check_array, check_random_state, check_X_y
from .utils.fixes import parallel_helper
from .utils.validation import check_is_fitted, has_fit_parameter
from .externals.joblib import Parallel
from .externals.joblib import delayed

__all__ = ["MultiOutputRegressor"]


def _fit_regression(estimator, X, y, sample_weight):
    estimator = clone(estimator)
    if sample_weight is not None:
        estimator.fit(X, y, sample_weight=sample_weight)
    else:
        estimator.fit(X, y)
    return estimator


class MultiOutputRegressor(BaseEstimator, RegressorMixin, MetaEstimatorMixin):
    """Multi target regression

    This strategy consits of fitting one regressor per target. This is a
    simple strategy for extending regressors that do not natively support
    multi-target regression.

    Parameters
    ----------
    estimator : estimator object
        An estimator object implementing `fit` and `predict`.

    n_jobs : int, optional, default: 1
        The number of jobs to run in parallel for `fit`. If -1,
        then the number of jobs is set to the number of cores.
        When individual estimators are fast to train or predict
        using `n_jobs>1` can result in slower performance due
        to the overhead of spawning processes.
    """
    def __init__(self, estimator, n_jobs=1):
        self.estimator = estimator
        self.n_jobs = n_jobs

    def fit(self, X, y, sample_weight=None):
        """Fit underlying regressors

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]
            Data.

        y : array-like, shape = [n_samples, n_targets]
            The multiple target values.

        sample_weight : array-like, shape = [n_samples] or None
            Sample weights. If None, then samples are equally weighted.
            Only supported if the underlying regressor supports sample
            weights.

        Returns
        -------
        self
        """
        X, y = check_X_y(X, y,
                         multi_output=True,
                         accept_sparse=True)

        if y.ndim == 1:
            raise ValueError("y must have at least two dimensions for "
                             "multi target regression but has only one.")

        if (sample_weight is not None and
            not has_fit_parameter(self.estimator, 'sample_weight')):
            raise ValueError("Underlying regressor does not support"
                             " sample weights.")

        self.estimators_ = Parallel(n_jobs=self.n_jobs)(delayed(_fit_regression)(
            self.estimator, X, y[:, i], sample_weight) for i in range(y.shape[1]))
        return self

    def predict(self, X):
        """Predict regression values for X.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]
            Data.

        Returns
        -------
        y : array of shape [n_samples, n_targets]
            Predicted values.
        """
        check_is_fitted(self, 'estimators_')

        X = check_array(X, accept_sparse=True)

        pred = Parallel(n_jobs=self.n_jobs)(delayed(parallel_helper)(e, 'predict', X)
                                            for e in self.estimators_)

        return np.asarray(pred).T

    def score(self, X, y, sample_weight=None):
        """Returns the coefficient of determination R^2 of the prediction.

        The coefficient R^2 is defined as (1 - u/v), where u is the regression
        sum of squares ((y_true - y_pred) ** 2).sum() and v is the residual
        sum of squares ((y_true - y_true.mean()) ** 2).sum().
        Best possible score is 1.0 and it can be negative (because the
        model can be arbitrarily worse). A constant model that always
        predicts the expected value of y, disregarding the input features,
        would get a R^2 score of 0.0.

        Note
        ----
        R^2 is calculated by weighting all the targets equally using
        `multioutput='uniform_average'`.

        Parameters
        ----------
        X : array-like, shape = (n_samples, n_features)
            Test samples.

        y : array-like, shape = (n_samples) or (n_samples, n_outputs)
            True values for X.

        sample_weight : array-like, shape = [n_samples], optional
            Sample weights.

        Returns
        -------
        score : float
            R^2 of self.predict(X) wrt. y.
        """
        # XXX remove in 0.19 when r2_score default for multioutput changes
        from .metrics import r2_score
        return r2_score(y, self.predict(X), sample_weight=sample_weight,
                        multioutput='uniform_average')
