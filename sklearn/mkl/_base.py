"""Base class for Multiple Kernel Learning."""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import warnings
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
        "kernels": [list, None],
        "kernels_scope": [list, None],
        "kernels_params": [list, None],
        "precompute_kernels": [bool, None],
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
        kernels_scope,
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
        self.kernels_scope = self._check_and_prepare_kernels_scope(kernels_scope)
        self.kernels_params = self._check_and_prepare_kernels_params(kernels_params)
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
        self.n_kernels_ = number_of_kernels(
            X=X,
            kernels=self.kernels,
            kernels_scope=self.kernels_scope,
            kernels_params=self.kernels_params,
            precomputed_kernels=False,  # Kernels are not precomputed yet
        )
        self.n_samples_ = len(y)
        X, y = self._check_and_prepare_X_y(X, y)

        if self.epsilon is None:
            if self._svm.__class__.__name__ == "SVC" and np.unique(y).shape[0] > 2:
                self.epsilon = 1e-1
            else:
                self.epsilon = 1e-2

        # TODO: Manage kernels and kernels_params
        self.weights_, self._svm = ALGORITHMS[self.algo].learn(
            X=X,
            y=y,
            svm=self._svm,
            kernels=self.kernels,
            kernels_scope=self.kernels_scope,
            kernels_params=self.kernels_params,
            precomputed_kernels=self._precomputed_kernels,
            n_kernels=self.n_kernels_,
            n_samples=self.n_samples_,
            epsilon=self.epsilon,
            tol=self.tol,
            verbose=self.verbose,
            max_iter=self.max_iter,
        )

        return self

    def _check_and_prepare_kernels(self, kernels):
        if kernels is None:
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
            "Invalid kernels. Kernels must be a list of "
            "callables and/or strings, or None."
        )

    def _check_and_prepare_kernels_scope(self, kernels_scope):
        if self.kernels is None:
            if kernels_scope is not None:
                warnings.warn(
                    "Attribute 'kernels_scope' is not None while 'kernels' is None. "
                    "'kernels_scope' will be ignored."
                )
            return kernels_scope

        if kernels_scope is None:
            return ["all" for _ in range(len(self.kernels))]

        if len(kernels_scope) == len(self.kernels) and all(
            scope in {"single", "all"} for scope in kernels_scope
        ):
            return kernels_scope

        raise ValueError(
            "Invalid 'kernels_scope'. Attribute 'kernels_scope' must be a list of "
            "{'all', 'single'} with the same size as 'kernels', or None."
        )

    def _check_and_prepare_kernels_params(self, kernels_params):
        if self.kernels is None:
            if kernels_params is not None:
                warnings.warn(
                    "Attribute 'kernels_params' is not None while 'kernels' is None. "
                    "'kernels_params' will be ignored."
                )
            return kernels_params

        if kernels_params is None:
            return [{} for _ in range(len(self.kernels))]

        if len(kernels_params) == len(self.kernels) and all(
            isinstance(params, dict) for params in kernels_params
        ):
            for params in kernels_params:
                prev_size = -1
                for key, value in params.items():
                    if not isinstance(key, str):
                        raise ValueError(
                            "Invalid 'kernels_params'. Keys must be strings "
                            "identifying the parameter name."
                        )
                    if not isinstance(value, list):
                        raise ValueError(
                            "Invalid 'kernels_params'. Values must be lists of the "
                            "values the parameter (identified by the key) will take."
                        )
                    if prev_size == -1:
                        prev_size = len(value)
                    elif prev_size != len(value):
                        raise ValueError(
                            "Invalid 'kernels_params'. All lists of "
                            "parameter values must have the same size."
                        )
            return kernels_params

        raise ValueError(
            "Invalid 'kernels_params'. Attribute 'kernels_params' must be a list of "
            "dictionaries with the same size as 'kernels', or None."
        )

    def _check_and_prepare_X_y(self, X, y):
        X = np.asarray(X, dtype=np.float64)
        self._precomputed_kernels = self.precompute_kernels

        if self.kernels is None:
            if len(X.shape) != 3 or X.shape[1] != X.shape[2]:
                raise ValueError(
                    "X must be a 3D array of shape (n_kernels, n_samples, n_samples) "
                    "when using precomputed kernels."
                )
            for i in range(X.shape[0]):
                X[i, :, :], y = check_X_y(X[i, :, :], y)
        else:
            X, y = check_X_y(X, y)
            if self.precompute_kernels or self.precompute_kernels is None:
                try:
                    new_X = np.empty(
                        (self.n_kernels_, X.shape[0], X.shape[0]),
                        dtype=np.float64,
                    )
                    kernel = np.empty(
                        (X.shape[0], X.shape[0]),
                        dtype=np.float64,
                    )  # Used to check if enough memory is available for computing
                    self._precomputed_kernels = True

                    if self.verbose:
                        print(f"[{self.__class__.__name__}] Precomputing kernels...")

                    for i, kernel in enumerate(
                        kernel_generator(
                            X=X,
                            kernels=self.kernels,
                            kernels_scope=self.kernels_scope,
                            kernels_params=self.kernels_params,
                            precomputed_kernels=False,  # We are precomputing kernels
                        )
                    ):
                        new_X[i, :, :] = kernel

                    return new_X, y
                except MemoryError:
                    if self.precompute_kernels:
                        raise MemoryError(
                            "Memory error occurred while precomputing kernels. "
                            "Try setting 'precompute_kernels=False'."
                        )
                    if self.verbose:
                        print(
                            f"[{self.__class__.__name__}] Not enough memory to "
                            "precompute kernels. Kernels will be computed on-the-fly. "
                            "This may significantly slow down the computation."
                        )
                    self._precomputed_kernels = False

        return X, y
