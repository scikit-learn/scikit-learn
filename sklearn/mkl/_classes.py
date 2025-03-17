# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from ..base import (
    ClassifierMixin,
    OutlierMixin,
    RegressorMixin,
)
from ._base import BaseMKL
from ._svm import SVC, SVR, OneClassSVM


class MKLC(BaseMKL, ClassifierMixin):
    """Multiple Kernel Learning for Classification."""

    def __init__(
        self,
        *,
        kernels=None,  # None or list of functions or list of strings
        kernels_scope=None,  # None or list of {"single", "all"}
        kernels_params=None,  # None or list of (str, dict)
        precompute_kernels=None,  # If none, it tries to compute the kernels
        algo="simple",
        epsilon=None,  # DOC: auto depending on algo
        tol=1e-8,
        verbose=False,
        max_iter=-1,
        random_state=None,
        **svm_params,
    ):
        super().__init__(
            kernels=kernels,
            kernels_scope=kernels_scope,
            kernels_params=kernels_params,
            precompute_kernels=precompute_kernels,
            algo=algo,
            epsilon=epsilon,
            tol=tol,
            verbose=verbose,
            max_iter=max_iter,
            random_state=random_state,
        )
        self._svm = SVC(
            kernel="precomputed",
            random_state=random_state,
            **_check_and_prepare_svm_params(svm_params),
        )


class MKLR(BaseMKL, RegressorMixin):
    """Multiple Kernel Learning for Regression."""

    def __init__(
        self,
        *,
        kernels=None,
        kernels_scope=None,
        kernels_params=None,
        precompute_kernels=None,
        algo="simple",
        epsilon=None,
        tol=1e-8,
        verbose=False,
        max_iter=-1,
        random_state=None,
        **svm_params,
    ):
        super().__init__(
            kernels=kernels,
            kernels_scope=kernels_scope,
            kernels_params=kernels_params,
            precompute_kernels=precompute_kernels,
            algo=algo,
            epsilon=epsilon,
            tol=tol,
            verbose=verbose,
            max_iter=max_iter,
            random_state=random_state,
        )
        self._svm = SVR(
            kernel="precomputed",
            **_check_and_prepare_svm_params(svm_params),
        )


class OneClassMKL(BaseMKL, OutlierMixin):
    """Multiple Kernel Learning for Outlier Detection."""

    def __init__(
        self,
        *,
        kernels=None,
        kernels_scope=None,
        kernels_params=None,
        precompute_kernels=None,
        algo="simple",
        epsilon=None,
        tol=1e-8,
        verbose=False,
        max_iter=-1,
        random_state=None,
        **svm_params,
    ):
        super().__init__(
            kernels=kernels,
            kernels_scope=kernels_scope,
            kernels_params=kernels_params,
            precompute_kernels=precompute_kernels,
            algo=algo,
            epsilon=epsilon,
            tol=tol,
            verbose=verbose,
            max_iter=max_iter,
            random_state=random_state,
        )
        self._svm = OneClassSVM(
            kernel="precomputed",
            **_check_and_prepare_svm_params(svm_params),
        )


def _check_and_prepare_svm_params(svm_params):
    assert all(
        k.startswith("svm__") for k in svm_params
    ), "Internal SVM parameters must start with 'svm__'."
    assert "svm__kernel" not in svm_params, "Internal SVM kernel cannot be set."
    assert (
        "svm__random_state" not in svm_params
    ), "Internal SVM random state can be set with MKL random state."

    return {k.removeprefix("svm__"): v for k, v in svm_params.items()}
