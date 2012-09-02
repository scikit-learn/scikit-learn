"""
Generalized Linear models.
"""

# Author: Alexandre Gramfort <alexandre.gramfort@inria.fr>
#         Fabian Pedregosa <fabian.pedregosa@inria.fr>
#         Olivier Grisel <olivier.grisel@ensta.org>
#         Vincent Michel <vincent.michel@inria.fr>
#         Peter Prettenhofer <peter.prettenhofer@gmail.com>
#         Mathieu Blondel <mathieu@mblondel.org>
#
# License: BSD Style.

from abc import ABCMeta, abstractmethod

import numpy as np
import scipy.sparse as sp
from scipy import linalg

from ..externals.joblib import Parallel, delayed
from ..base import BaseEstimator
from ..base import RegressorMixin
from ..utils.extmath import safe_sparse_dot
from ..utils import as_float_array, safe_asarray
from ..utils.fixes import lsqr
from ..utils.sparsefuncs import csc_mean_variance_axis0, \
                                inplace_csc_column_scale
from cd_fast import sparse_std


###
### TODO: intercept for all models
### We should define a common function to center data instead of
### repeating the same code inside each fit method.

### TODO: bayesian_ridge_regression and bayesian_regression_ard
### should be squashed into its respective objects.


def sparse_center_data(X, y, fit_intercept, normalize=False):
    """
    Compute informations needed to center data to have mean zero along
    axis 0. Be aware that X will not be centered since it would break
    the sparsity, but will be normalized if asked so.
    """
    X_data = np.array(X.data, np.float64)
    if fit_intercept:
        # copy if 'normalize' is True or X is not a csc matrix
        X = sp.csc_matrix(X, copy=normalize)
        X_mean, X_std = csc_mean_variance_axis0(X)
        if normalize:
            X_std = sparse_std(
                X.shape[0], X.shape[1],
                X_data, X.indices, X.indptr, X_mean)
            X_std[X_std == 0] = 1
            inplace_csc_column_scale(X, X_std)
        else:
            X_std = np.ones(X.shape[1])
        y_mean = y.mean(axis=0)
        y = y - y_mean
    else:
        X_mean = np.zeros(X.shape[1])
        X_std = np.ones(X.shape[1])
        y_mean = 0. if y.ndim == 1 else np.zeros(y.shape[1], dtype=X.dtype)

    X_data = np.array(X.data, np.float64)
    return X_data, y, X_mean, y_mean, X_std


def center_data(X, y, fit_intercept, normalize=False, copy=True):
    """
    Centers data to have mean zero along axis 0. This is here because
    nearly all linear models will want their data to be centered.
    """
    X = as_float_array(X, copy)

    if fit_intercept:
        if sp.issparse(X):
            X_mean = np.zeros(X.shape[1])
            X_std = np.ones(X.shape[1])
        else:
            X_mean = X.mean(axis=0)
            X -= X_mean
            if normalize:
                X_std = np.sqrt(np.sum(X ** 2, axis=0))
                X_std[X_std == 0] = 1
                X /= X_std
            else:
                X_std = np.ones(X.shape[1])
        y_mean = y.mean(axis=0)
        y = y - y_mean
    else:
        X_mean = np.zeros(X.shape[1])
        X_std = np.ones(X.shape[1])
        y_mean = 0. if y.ndim == 1 else np.zeros(y.shape[1], dtype=X.dtype)
    return X, y, X_mean, y_mean, X_std


class LinearModel(BaseEstimator):
    """Base class for Linear Models"""
    __metaclass__ = ABCMeta

    @abstractmethod
    def fit(self, X, y):
        """Fit model."""

    def decision_function(self, X):
        """Decision function of the linear model

        Parameters
        ----------
        X : numpy array of shape [n_samples, n_features]

        Returns
        -------
        C : array, shape = [n_samples]
            Returns predicted values.
        """
        X = safe_asarray(X)
        return safe_sparse_dot(X, self.coef_.T) + self.intercept_

    def predict(self, X):
        """Predict using the linear model

        Parameters
        ----------
        X : numpy array of shape [n_samples, n_features]

        Returns
        -------
        C : array, shape = [n_samples]
            Returns predicted values.
        """
        return self.decision_function(X)

    _center_data = staticmethod(center_data)

    def _set_intercept(self, X_mean, y_mean, X_std):
        """Set the intercept_
        """
        if self.fit_intercept:
            self.coef_ = self.coef_ / X_std
            self.intercept_ = y_mean - np.dot(X_mean, self.coef_.T)
        else:
            self.intercept_ = 0.


class LinearRegression(LinearModel, RegressorMixin):
    """
    Ordinary least squares Linear Regression.

    Attributes
    ----------
    `coef_` : array
        Estimated coefficients for the linear regression problem.

    `intercept_` : array
        Independent term in the linear model.

    Parameters
    ----------
    fit_intercept : boolean, optional
        wether to calculate the intercept for this model. If set
        to false, no intercept will be used in calculations
        (e.g. data is expected to be already centered).
    normalize : boolean, optional
        If True, the regressors X are normalized

    Notes
    -----
    From the implementation point of view, this is just plain Ordinary
    Least Squares (numpy.linalg.lstsq) wrapped as a predictor object.

    """

    def __init__(self, fit_intercept=True, normalize=False, copy_X=True):
        self.fit_intercept = fit_intercept
        self.normalize = normalize
        self.copy_X = copy_X

    def fit(self, X, y, n_jobs=1):
        """
        Fit linear model.

        Parameters
        ----------
        X : numpy array or sparse matrix of shape [n_samples,n_features]
            Training data
        y : numpy array of shape [n_samples, n_responses]
            Target values
        n_jobs : The number of jobs to use for the computation.
            If -1 all CPUs are used. This will only provide speedup for
            n_response > 1 and sufficient large problems

        Returns
        -------
        self : returns an instance of self.
        """
        X = safe_asarray(X)
        y = np.asarray(y)

        X, y, X_mean, y_mean, X_std = self._center_data(X, y,
                self.fit_intercept, self.normalize, self.copy_X)

        if sp.issparse(X):
            if y.ndim < 2:
                out = lsqr(X, y)
                self.coef_ = out[0]
                self.residues_ = out[3]
            else:
                # sparse_lstsq cannot handle y with shape (M, K)
                outs = Parallel(n_jobs=n_jobs)(delayed(lsqr)
                        (X, y[:, j].ravel()) for j in range(y.shape[1]))
                self.coef_ = np.vstack(out[0] for out in outs)
                self.residues_ = np.vstack(out[3] for out in outs)
        else:
            self.coef_, self.residues_, self.rank_, self.singular_ = \
                    linalg.lstsq(X, y)
            self.coef_ = self.coef_.T

        if y.ndim == 1:
            self.coef_ = np.ravel(self.coef_)
        self._set_intercept(X_mean, y_mean, X_std)
        return self
