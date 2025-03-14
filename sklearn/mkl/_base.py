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
from ..metrics.pairwise import (
    laplacian_kernel,
    linear_kernel,
    polynomial_kernel,
    rbf_kernel,
    sigmoid_kernel,
)
from ..utils._param_validation import Interval, StrOptions
from ..utils.validation import check_is_fitted, check_X_y
from ._algo import (
    _average_mkl,
    _simple_mkl,
    _sum_mkl,
)
from ._utils import kernel_generator, number_of_kernels

ALGORITHMS = {
    "average": _average_mkl,
    "simple": _simple_mkl,
    "sum": _sum_mkl,
}


KERNELS = {
    "laplace": laplacian_kernel,
    "laplacian": laplacian_kernel,
    "linear": linear_kernel,
    "poly": polynomial_kernel,
    "polynomial": polynomial_kernel,
    "rbf": rbf_kernel,
    "gaussian": rbf_kernel,
    "sigmoid": sigmoid_kernel,
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
        precompute_kernels,
        algo,
        epsilon,
        tol,
        verbose,
        max_iter,
        random_state,
    ):
        self.kernels = self._check_and_prepare_kernels(kernels)
        self.kernels_params = kernels_params
        self.precompute_kernels = precompute_kernels
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
                "Please call 'fit' with appropriate arguments before using this "
                "instance as a callable."
            ),
        )
        return self.transform(X)

    @_fit_context(prefer_skip_nested_validation=True)
    def fit(self, X, y=None):
        # TODO: DOC: X : list of kernels matrices (n, n) or array-like of shape (n, m)
        self.n_kernels_ = number_of_kernels(X, self.kernels, self.kernels_params)
        X_fit, y_fit = self._check_and_prepare_X_y(X, y)

        if self.epsilon is None:
            if self._svm.__class__.__name__ == "SVC" and np.unique(y).shape[0] > 2:
                self.epsilon = 1e-1
            else:
                self.epsilon = 1e-2

        # TODO: Manage kernels and kernels_params
        self.weights_, self._svm = ALGORITHMS[self.algo].learn(
            X=X_fit,
            y=y_fit,
            svm=self._svm,
            kernels=self.kernels,
            kernels_params=self.kernels_params,
            precompute_kernels=self.precompute_kernels,
            n_kernels=self.n_kernels_,
            epsilon=self.epsilon,
            tol=self.tol,
            verbose=self.verbose,
            max_iter=self.max_iter,
        )

        return self

    @staticmethod
    def _check_and_prepare_kernels(kernels):  # TODO: Take kernels_params as input
        if kernels == "precomputed":
            return kernels

        if isinstance(kernels, list):
            kernel_list = []
            for kernel in kernels:
                if callable(kernel):
                    kernel_list.append(kernel)
                elif isinstance(kernel, str):
                    if kernel in KERNELS:
                        kernel_list.append(KERNELS[kernel])
                    else:
                        raise ValueError(
                            f"Invalid kernel '{kernel}'. "
                            f"Valid kernels are {set(KERNELS)}."
                        )
                else:
                    raise ValueError(
                        "Invalid kernel. Kernels must be callables or strings."
                    )
            return kernel_list

        raise ValueError(
            "Invalid kernels. Kernels must be 'precomputed' "
            "or a list of callables and/or strings."
        )

    def _check_and_prepare_X_y(self, X, y):
        X = np.asarray(X, dtype=np.float64)

        if isinstance(self.kernels, str) and self.kernels == "precomputed":
            if len(X.shape) != 3 or X.shape[1] != X.shape[2]:
                raise ValueError(
                    "X must be a 3D array of shape (n_kernels, n_samples, n_samples) "
                    "when using precomputed kernels."
                )
            for i in range(X.shape[0]):
                X[i, :, :], y = check_X_y(X[i, :, :], y)
        else:
            X, y = check_X_y(X, y)
            if self.precompute_kernels:
                new_X = np.empty(
                    (self.n_kernels_, X.shape[0], X.shape[0]), dtype=np.float64
                )
                for i, kernel in enumerate(
                    kernel_generator(X, self.kernels, self.kernels_params)
                ):
                    new_X[i, :, :] = kernel
                return new_X, y

        return X, y
