"""Multiple Kernel Learning (MKL) classes."""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

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
        C=1.0,
        algo="simple",
        svm_params=None,
        precompute_kernels=None,  # If none, it tries to compute the kernels
        tol=None,  # DOC: auto depending on algo
        numeric_tol=1e-8,
        verbose=False,
        max_iter=200,
        random_state=None,
    ):
        super().__init__(
            algo=algo,
            kernels=kernels,
            kernels_scopes=kernels_scopes,
            kernels_param_grids=kernels_param_grids,
            precompute_kernels=precompute_kernels,
            tol=tol,
            numeric_tol=numeric_tol,
            verbose=verbose,
            max_iter=max_iter,
            random_state=random_state,
        )
        self.C = C
        self._warn_svm_params(svm_params, ["C", "random_state"])
        self.svm_params = svm_params

    def decision_function(self, X):
        kernel = self.transform(X)
        return self._svm.decision_function(kernel)

    def _set_svm(self):
        self._svm = SVC(
            kernel="precomputed",
            C=self.C,
            random_state=self.random_state,
            **({} if self.svm_params is None else self.svm_params),
        )

    def _post_learning_processing(self):
        super()._post_learning_processing()
        self.classes_ = self._svm.classes_


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
        C=1.0,
        epsilon=0.1,
        algo="simple",
        svm_params=None,
        precompute_kernels=None,
        tol=None,
        numeric_tol=1e-8,
        verbose=False,
        max_iter=200,
        random_state=None,
    ):
        super().__init__(
            algo=algo,
            kernels=kernels,
            kernels_scopes=kernels_scopes,
            kernels_param_grids=kernels_param_grids,
            precompute_kernels=precompute_kernels,
            tol=tol,
            numeric_tol=numeric_tol,
            verbose=verbose,
            max_iter=max_iter,
            random_state=random_state,
        )
        self.C = C
        self.epsilon = epsilon
        self._warn_svm_params(svm_params, ["C", "epsilon"])
        self.svm_params = svm_params

    def _set_svm(self):
        self._svm = SVR(
            kernel="precomputed",
            C=self.C,
            epsilon=self.epsilon,
            **({} if self.svm_params is None else self.svm_params),
        )


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
        nu=0.5,
        algo="simple",
        svm_params=None,
        precompute_kernels=None,
        tol=None,
        numeric_tol=1e-8,
        verbose=False,
        max_iter=200,
        random_state=None,
    ):
        super().__init__(
            algo=algo,
            kernels=kernels,
            kernels_scopes=kernels_scopes,
            kernels_param_grids=kernels_param_grids,
            precompute_kernels=precompute_kernels,
            tol=tol,
            numeric_tol=numeric_tol,
            verbose=verbose,
            max_iter=max_iter,
            random_state=random_state,
        )
        self.nu = nu
        self._warn_svm_params(svm_params, ["nu"])
        self.svm_params = svm_params

    def decision_function(self, X):
        kernel = self.transform(X)
        return self._svm.decision_function(kernel)

    def score_samples(self, X):
        kernel = self.transform(X)
        return self._svm.score_samples(kernel)

    def _set_svm(self):
        svm_params = {} if self.svm_params is None else self.svm_params
        self._svm = OneClassSVM(
            kernel="precomputed",
            nu=self.nu,
            **svm_params,
        )

    def _post_learning_processing(self):
        super()._post_learning_processing()
        self.offset_ = self._svm.offset_
