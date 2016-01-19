"""
Multioutput regression strategies
=================================

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
from .utils.validation import check_is_fitted
from .externals.joblib import Parallel
from .externals.joblib import delayed

__all__ = ["MultiOutputRegressor"]


def _fit_regression(estimator, X, y):
    estimator = clone(estimator)
    estimator.fit(X, y)
    return estimator


def _parallel_helper(obj, methodname, *args, **kwargs):
    """Private helper to workaround Python 2 pickle limitations"""
    return getattr(obj, methodname)(*args, **kwargs)


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
    """
    def __init__(self, estimator, n_jobs=1):
        self.estimator = estimator
        self.n_jobs = n_jobs

    def fit(self, X, y):
        """Fit underlying regressors

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]
            Data.

        y : array-like, shape = [n_samples, n_targets]
            Multiple targets.

        Returns
        -------
        self
        """
        X, y = check_X_y(X, y,
                         multi_output=True,
                         accept_sparse=['csr', 'csc', 'coo', 'dok', 'lil'])

        if y.ndim == 1:
            raise ValueError("y must have at least two dimensions for "
                             "multi target regression but has only one.")

        self.estimators_ = Parallel(n_jobs=self.n_jobs)(delayed(_fit_regression)(
            self.estimator, X, y[:, i]) for i in range(y.shape[1]))
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

        X = check_array(X, accept_sparse=['csr', 'csc', 'coo', 'dok', 'lil'])

        pred = Parallel(n_jobs=self.n_jobs)(delayed(_parallel_helper)(e, 'predict', X)
                                            for e in self.estimators_)

        return np.asarray(pred).T
