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
from ..metrics.pairwise import PAIRWISE_KERNEL_FUNCTIONS
from ..utils._param_validation import Interval, StrOptions
from ..utils.validation import (
    _check_n_features,
    check_consistent_length,
    check_is_fitted,
    validate_data,
)
from ._algo import (
    _average_mkl,
    _simple_mkl,
    _sum_mkl,
)
from ._algo._utils import combine_kernels_nonsym  # type: ignore
from ._utils import kernel_generator, number_of_kernels

MKL_ALGORITHMS = {
    "average": _average_mkl,
    "simple": _simple_mkl,
    "sum": _sum_mkl,
}


class BaseMKL(TransformerMixin, MetaEstimatorMixin, BaseEstimator, metaclass=ABCMeta):
    """Base class for Multiple Kernel Learning (MKL) estimators.

    This class provides a framework for learning a combination of multiple kernels
    for classification and regression tasks. It implements different MKL algorithms
    to optimize kernel-based learning models.
    """

    _parameter_constraints: dict = {
        "kernels": ["array-like", None],
        "kernels_scopes": ["array-like", None],
        "kernels_param_grids": ["array-like", None],
        "precompute_kernels": ["boolean", None],
        "algo": [StrOptions({algo for algo in MKL_ALGORITHMS})],
        "tol": [Interval(Real, 0.0, None, closed="neither"), None],
        "numeric_tol": [Interval(Real, 0.0, None, closed="neither")],
        "verbose": ["verbose"],
        "max_iter": [Interval(Integral, -1, None, closed="left")],
        "random_state": ["random_state"],
    }

    @abstractmethod
    def __init__(
        self,
        algo,
        kernels,
        kernels_scopes,
        kernels_param_grids,
        precompute_kernels,
        tol,
        numeric_tol,
        verbose,
        max_iter,
        random_state,
    ):
        self.algo = algo
        self.kernels = kernels
        self.kernels_scopes = kernels_scopes
        self.kernels_param_grids = kernels_param_grids
        self.precompute_kernels = precompute_kernels
        self.tol = tol
        self.numeric_tol = numeric_tol
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

        Notes
        -----
        Most estimators check the consistency of the number of samples between `X`
        and `y`. This method is primarily designed for dynamically computed kernels
        and may not work correctly when using precomputed kernels (i.e., when `kernels`
        is set to None). In such cases, shape mismatches can occur, making it
        impractical to use as a kernel argument.
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
        X, y = self._prepare_X_y_for_learning(X, y)
        self._set_svm()  # Defined in subclasses
        self._learn(X, y)
        self._post_learning_processing()
        return self

    def transform(self, X):
        check_is_fitted(self)
        X, _ = self._validate_data(X, fit=False)

        if self._kernels is None:
            return combine_kernels_nonsym(
                weights=self.weights_,
                kernels=X,
                n_samples=X[0].shape[0],
                m_samples=X[0].shape[1],
            ).base
        else:
            return combine_kernels_nonsym(
                weights=self.weights_,
                kernels=kernel_generator(
                    X=X,
                    Y=self.__Xfit,
                    kernels=self._kernels,
                    kernels_scopes=self._kernels_scopes,
                    kernels_param_grids=self._kernels_param_grids,
                    precomputed_kernels=False,  # We need to compute the kernels
                ),
                n_samples=(X.shape[0] if hasattr(X, "shape") else len(X)),
                m_samples=(
                    self.__Xfit.shape[0]
                    if hasattr(self.__Xfit, "shape")
                    else len(self.__Xfit)
                ),
            ).base

    def predict(self, X):
        kernel = self.transform(X)
        return self._svm.predict(kernel)

    def _learn(self, X, y=None):
        if self.algo in ["average", "sum"]:
            self.weights_, self._svm, self.n_iter_ = MKL_ALGORITHMS[self.algo].learn(
                X=X,
                y=y,
                svm=self._svm,
                kernels=self._kernels,
                kernels_scopes=self._kernels_scopes,
                kernels_param_grids=self._kernels_param_grids,
                precomputed_kernels=self._precomputed_kernels,
                n_kernels=self.n_kernels_,
                n_samples=self.n_samples_in_,
                verbose=self.verbose,
            )
        elif self.algo == "simple":
            if self.tol is None:
                if self._svm.__class__.__name__ == "SVC" and np.unique(y).shape[0] > 2:
                    self._tol = 1e-1
                else:
                    self._tol = 1e-2
            else:
                self._tol = self.tol

            self.weights_, self._svm, self.n_iter_ = MKL_ALGORITHMS[self.algo].learn(
                X=X,
                y=y,
                svm=self._svm,
                kernels=self._kernels,
                kernels_scopes=self._kernels_scopes,
                kernels_param_grids=self._kernels_param_grids,
                precomputed_kernels=self._precomputed_kernels,
                n_kernels=self.n_kernels_,
                n_samples=self.n_samples_in_,
                tol=self._tol,
                numeric_tol=self.numeric_tol,
                verbose=self.verbose,
                max_iter=self.max_iter,
            )

    def _prepare_X_y_for_learning(self, X, y=None):
        self._kernels, self._kernels_scopes, self._kernels_param_grids = (
            self._validate_kernels(self.kernels),
            self._validate_kernels_scopes(self.kernels_scopes),
            self._validate_kernels_param_grids(self.kernels_param_grids),
        )
        X, y = self._validate_data(X, y)

        self.n_kernels_ = number_of_kernels(
            X=X,
            kernels=self._kernels,
            kernels_scopes=self._kernels_scopes,
            kernels_param_grids=self._kernels_param_grids,
            precomputed_kernels=False,  # Kernels are not precomputed yet
        )
        self.n_samples_in_ = (
            X[0].shape[0]
            if self._kernels is None
            else X.shape[0] if hasattr(X, "shape") else len(X)
        )

        if self._kernels is not None:
            # Reference to X to compute the kernel in predict/transform
            self.__Xfit = X

        X = self._precompute_kernels_if_needed(X)

        return X, y

    def _post_learning_processing(self):
        if hasattr(self._svm, "alpha_init_"):
            del self._svm.alpha_init_
        if hasattr(self._svm, "alpha_raw_"):
            del self._svm.alpha_raw_

    def _validate_kernels(self, kernels):
        if kernels is None:
            return kernels

        if isinstance(kernels, list):
            kernel_list = []
            for kernel in kernels:
                if callable(kernel):
                    kernel_list.append(kernel)
                elif isinstance(kernel, str):
                    if kernel in PAIRWISE_KERNEL_FUNCTIONS:
                        kernel_list.append(PAIRWISE_KERNEL_FUNCTIONS[kernel])
                    else:
                        raise ValueError(
                            f"Invalid kernel '{kernel}'. "
                            f"Valid kernels are {set(PAIRWISE_KERNEL_FUNCTIONS)}."
                        )
                else:
                    raise ValueError(
                        "Invalid kernel. Kernels must be callables or strings, "
                        f"not `{type(kernel)}`."
                    )
            return kernel_list

        raise ValueError(
            "Invalid kernels. Kernels must be a list of "
            f"callables and/or strings, or None, not `{type(kernels)}`."
        )

    def _validate_kernels_scopes(self, kernels_scopes):
        if self.kernels is None:
            if kernels_scopes is not None:
                warnings.warn(
                    "Attribute 'kernels_scopes' is not None while 'kernels' is None. "
                    "'kernels_scopes' will be ignored."
                )
            return kernels_scopes

        if kernels_scopes is None:
            return ["all" for _ in range(len(self.kernels))]

        if len(kernels_scopes) == len(self.kernels) and all(
            scope in {"single", "all"} for scope in kernels_scopes
        ):
            return kernels_scopes

        raise ValueError(
            "Invalid 'kernels_scopes'. Attribute 'kernels_scopes' must be a list of "
            "{'all', 'single'} with the same size as 'kernels', or None."
        )

    def _validate_kernels_param_grids(self, kernels_param_grids):
        if self.kernels is None:
            if kernels_param_grids is not None:
                warnings.warn(
                    "Attribute 'kernels_param_grids' is not None while 'kernels' is "
                    "None. 'kernels_param_grids' will be ignored."
                )
            return kernels_param_grids

        if kernels_param_grids is None:
            return [{} for _ in range(len(self.kernels))]

        if len(kernels_param_grids) == len(self.kernels) and all(
            isinstance(params, dict) for params in kernels_param_grids
        ):
            for params in kernels_param_grids:
                prev_size = -1
                for key, value in params.items():
                    if not isinstance(key, str):
                        raise ValueError(
                            "Invalid 'kernels_param_grids'. Keys must be strings "
                            "identifying the parameter name."
                        )
                    if not isinstance(value, (list, tuple, np.ndarray)):
                        raise ValueError(
                            "Invalid 'kernels_param_grids'. "
                            "Values must be lists of the values the "
                            "parameter (identified by the key) will take."
                        )
                    if prev_size == -1:
                        prev_size = len(value)
                    elif prev_size != len(value):
                        raise ValueError(
                            "Invalid 'kernels_param_grids'. All lists of "
                            "parameter values must have the same size."
                        )
            return kernels_param_grids

        raise ValueError(
            "Invalid 'kernels_param_grids'. Attribute 'kernels_param_grids' must be a "
            "list of dictionaries with the same size as 'kernels', or None."
        )

    def _validate_data(self, X, y=None, fit=True):
        if not fit:
            check_is_fitted(self)

        if self._kernels is None:
            if fit and (
                len(np.shape(X[0])) != 2 or np.shape(X[0])[0] != np.shape(X[0])[1]
            ):
                raise ValueError(
                    "X must be a 3D array of shape (n_kernels, n_samples, "
                    "n_samples) when using precomputed kernels."
                )
            elif not fit and (
                len(np.shape(X[0])) != 2 or np.shape(X[0])[1] != self.n_samples_in_
            ):
                raise ValueError(
                    "X must be a 3D array of shape (n_kernels, n_samples, "
                    "n_fit_samples) when using precomputed kernels."
                )
            for i in range(np.shape(X)[0]):
                check_consistent_length(X[i], y)
        else:
            try:
                if y is None:
                    X = validate_data(
                        self,
                        X,
                        reset=fit,
                        dtype=np.float64,
                        order="C",
                        accept_sparse="csr",
                        accept_large_sparse=False,
                        ensure_min_samples=(2 if fit else 1),
                    )
                else:
                    X, y = validate_data(
                        self,
                        X,
                        y,
                        reset=fit,
                        dtype=np.float64,
                        order="C",
                        accept_sparse="csr",
                        accept_large_sparse=False,
                        ensure_min_samples=(2 if fit else 1),
                    )
            except ValueError:
                check_consistent_length(X, y)
                _check_n_features(self, X, reset=fit)

        return X, y

    def _precompute_kernels_if_needed(self, X):
        self._precomputed_kernels = self.precompute_kernels
        if self._kernels is not None and (
            self.precompute_kernels or self.precompute_kernels is None
        ):
            try:
                shape = X.shape[0] if hasattr(X, "shape") else len(X)
                new_X = np.empty(
                    (self.n_kernels_, shape, shape),
                    dtype=np.float64,
                )
                kernel = np.empty(
                    (shape, shape),
                    dtype=np.float64,
                )  # Used to check if enough memory is available for computing
                self._precomputed_kernels = True

                if self.verbose:
                    print(f"[{self.__class__.__name__}] Precomputing kernels...")

                for i, kernel in enumerate(
                    kernel_generator(
                        X=X,
                        kernels=self._kernels,
                        kernels_scopes=self._kernels_scopes,
                        kernels_param_grids=self._kernels_param_grids,
                        precomputed_kernels=False,  # We are precomputing kernels
                    )
                ):
                    new_X[i, :, :] = kernel

                return new_X
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
        return X

    def __sklearn_tags__(self):
        tags = super().__sklearn_tags__()
        tags.input_tags.sparse = True
        return tags
