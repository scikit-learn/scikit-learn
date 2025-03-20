# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import warnings

from ..base import (
    ClassifierMixin,
    OutlierMixin,
    RegressorMixin,
)
from ._base import BaseMKL
from ._svm import SVC, SVR, OneClassSVM


class MKLC(ClassifierMixin, BaseMKL):
    """Multiple Kernel Learning for Classification."""

    _parameter_constraints: dict = {
        **BaseMKL._parameter_constraints,
        "svm_params": [dict, None],
        # "algo": [StrOptions({"algo1", "algo2", ...})] # Algorithms supported
    }

    def __init__(
        self,
        *,
        kernels=None,  # None or list of functions or list of strings
        kernels_scopes=None,  # None or list of {"single", "all"}
        kernels_param_grids=None,  # None or list of (str, dict)
        algo="simple",
        svm_params=None,
        precompute_kernels=None,  # If none, it tries to compute the kernels
        tol=None,  # DOC: auto depending on algo
        numeric_tol=1e-8,
        verbose=False,
        max_iter=-1,
        random_state=None,
    ):
        super().__init__(
            kernels=kernels,
            kernels_scopes=kernels_scopes,
            kernels_param_grids=kernels_param_grids,
            precompute_kernels=precompute_kernels,
            algo=algo,
            tol=tol,
            numeric_tol=numeric_tol,
            verbose=verbose,
            max_iter=max_iter,
            random_state=random_state,
        )
        _warn_svm_params(svm_params)
        self.svm_params = svm_params

    def fit(self, X, y=None):
        self._svm = SVC(
            kernel="precomputed",
            random_state=self.random_state,
            **(self.svm_params if self.svm_params is not None else {}),
        )

        super().fit(X, y)

        self.classes_ = self._svm.classes_

        return self

    def decision_function(self, X):
        kernel = self.transform(X)
        return self._svm.decision_function(kernel)


class MKLR(RegressorMixin, BaseMKL):
    """Multiple Kernel Learning for Regression."""

    _parameter_constraints: dict = {
        **BaseMKL._parameter_constraints,
        "svm_params": [dict, None],
        # "algo": [StrOptions({"algo1", "algo2", ...})] # Algorithms supported
    }

    def __init__(
        self,
        *,
        kernels=None,
        kernels_scopes=None,
        kernels_param_grids=None,
        algo="simple",
        svm_params=None,
        precompute_kernels=None,
        tol=None,
        numeric_tol=1e-8,
        verbose=False,
        max_iter=-1,
        random_state=None,
    ):
        super().__init__(
            kernels=kernels,
            kernels_scopes=kernels_scopes,
            kernels_param_grids=kernels_param_grids,
            precompute_kernels=precompute_kernels,
            algo=algo,
            tol=tol,
            numeric_tol=numeric_tol,
            verbose=verbose,
            max_iter=max_iter,
            random_state=random_state,
        )
        _warn_svm_params(svm_params)
        self.svm_params = svm_params

    def fit(self, X, y=None):
        self._svm = SVR(
            kernel="precomputed",
            **(self.svm_params if self.svm_params is not None else {}),
        )

        super().fit(X, y)

        return self


class OneClassMKL(OutlierMixin, BaseMKL):
    """Multiple Kernel Learning for Outlier Detection."""

    _parameter_constraints: dict = {
        **BaseMKL._parameter_constraints,
        "svm_params": [dict, None],
        # "algo": [StrOptions({"algo1", "algo2", ...})] # Algorithms supported
    }

    def __init__(
        self,
        *,
        kernels=None,
        kernels_scopes=None,
        kernels_param_grids=None,
        algo="simple",
        svm_params=None,
        precompute_kernels=None,
        tol=None,
        numeric_tol=1e-8,
        verbose=False,
        max_iter=-1,
        random_state=None,
    ):
        super().__init__(
            kernels=kernels,
            kernels_scopes=kernels_scopes,
            kernels_param_grids=kernels_param_grids,
            precompute_kernels=precompute_kernels,
            algo=algo,
            tol=tol,
            numeric_tol=numeric_tol,
            verbose=verbose,
            max_iter=max_iter,
            random_state=random_state,
        )
        _warn_svm_params(svm_params)
        self.svm_params = svm_params

    def fit(self, X, y=None):
        self._svm = OneClassSVM(
            kernel="precomputed",
            **(self.svm_params if self.svm_params is not None else {}),
        )

        super().fit(X, y)

        self.offset_ = self._svm.offset_

        return self

    def decision_function(self, X):
        kernel = self.transform(X)
        return self._svm.decision_function(kernel)

    def score_samples(self, X):
        kernel = self.transform(X)
        return self._svm.score_samples(kernel)


def _warn_svm_params(svm_params):
    if isinstance(svm_params, dict):
        if "kernel" in svm_params:
            warnings.warn("Internal SVM kernel should not be set with MKL.")
        if "random_state" in svm_params:
            warnings.warn("Internal SVM random state can be set with MKL random state.")
