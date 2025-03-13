"""Base class for Multiple Kernel Learning."""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from abc import ABCMeta, abstractmethod
from numbers import Integral, Real

import numpy as np

from ..base import (
    BaseEstimator,
    MetaEstimatorMixin,
    TransformerMixin,
    _fit_context,
)
from ..utils._param_validation import Interval, StrOptions
from ..utils.validation import check_is_fitted
from ._algo import (
    _average_mkl,
    _simple_mkl,
    _sum_mkl,
)

ALGORITHMS = {
    "average": _average_mkl,
    "simple": _simple_mkl,
    "sum": _sum_mkl,
}


class BaseMKL(BaseEstimator, MetaEstimatorMixin, TransformerMixin, metaclass=ABCMeta):
    """Multiple Kernel Learning base class."""

    # TODO: DOC: kernels_params list of ({"single", "all"}, dict)

    _parameter_constraints: dict = {
        "kernels": [list, StrOptions({"precomputed"})],
        "kernels_params": [list, None],
        "algo": [StrOptions({algo for algo in ALGORITHMS})],
        "epsilon": [Interval(Real, 0.0, None, closed="neither"), None],
        "tol": [Interval(Real, 0.0, None, closed="neither")],
        "verbose": ["verbose"],
        "max_iter": [Interval(Integral, -1, None, closed="left")],
        "random_state": ["random_state"],
    }

    @abstractmethod
    def __init__(
        self,
        kernels,
        kernels_params,
        algo,
        epsilon,
        tol,
        verbose,
        max_iter,
        random_state,
    ):
        self.kernels = kernels
        self.kernels_params = kernels_params
        self.algo = algo
        self.epsilon = epsilon
        self.tol = tol
        self.verbose = verbose
        self.max_iter = max_iter
        self.random_state = random_state

    def __call__(self, X, Xfit=None):
        """
        Computes the combined kernel matrix of X.

        This method allows the instance to be used as a kernel argument in an estimator.
        The instance must be fitted using `fit` before calling this method.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            The input data to transform.

        Xfit : None, default=None
            This parameter is not used. It is only present to allow the instance to be
            passed as a kernel argument to an estimator.

        Returns
        -------
        K : array-like of shape (n_samples, n_fitted_samples)
            The computed kernel matrix.

        Raises
        ------
        NotFittedError
            If the instance has not been fitted yet.
        """
        check_is_fitted(
            self,
            msg=(
                "This %(name)s instance needs to be fitted before calling it. "
                "Call 'fit' with appropriate arguments before calling this instance."
            ),
        )
        return self.transform(X)

    @_fit_context(prefer_skip_nested_validation=True)
    def fit(self, X, y=None):
        # TODO: DOC: X : list of kernels matrices (n, n) or array-like of shape (n, m)
        # TODO: Use utils check functions
        X_fit = np.asarray(X, dtype=np.float64)
        y_fit = y

        if self.epsilon is None:
            if self._svm.__class__.__name__ == "SVC" and np.unique(y).shape[0] > 2:
                self.epsilon = 1e-1
            else:
                self.epsilon = 1e-2

        # TODO: Manage kernels and kernels_params
        self.weights_, self._svm = ALGORITHMS[self.algo].learn(
            self._svm,
            X_fit,
            y_fit,
            self.epsilon,
            self.tol,
            self.max_iter,
            self.verbose,
        )

        return self
