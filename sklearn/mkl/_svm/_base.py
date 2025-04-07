# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import numpy as np

from ...svm._base import BaseLibSVM
from . import _libsvm as libsvm  # type: ignore


class BaseLibSVMforMKL(BaseLibSVM):
    """
    Base class for Multiple Kernel Learning (MKL) with SVM. Inherits from
    `BaseLibSVM` and adds functionality for working with support vector machines in a
    multiple kernel learning context.

    This class is responsible for fitting the SVM model using MKL, with various
    attributes to handle the optimization process.

    Added Attributes
    ----------------
    objective_val_ : float
        The value of the objective function from the maximization problem during
        optimization. This value reflects the current objective at the end of the
        fitting process.

    alpha_raw_ : ndarray of shape (any, any)
        The raw dual coefficients returned by the LibSVM solver.

    alpha_raw_lengths_ : ndarray (any,)
        The lengths of the arrays of `alpha_raw_`.

    alpha_init_ : ndarray (any, any)
        Initializes the dual coefficients before the fitting process. If not provided,
        the coefficients are initialized automatically in LivSVM.
    """

    def _dense_fit(self, X, y, sample_weight, solver_type, kernel, random_seed):
        if callable(self.kernel):
            # you must store a reference to X to compute the kernel in predict
            # TODO: add keyword copy to copy on demand
            self._BaseLibSVM__Xfit = X
            X = self._compute_kernel(X)

            if X.shape[0] != X.shape[1]:
                raise ValueError("X.shape[0] should be equal to X.shape[1]")

        libsvm.set_verbosity_wrap(self.verbose)
        # we don't pass **self.get_params() to allow subclasses to
        # add other parameters to __init__
        (
            self.support_,
            self.support_vectors_,
            self._n_support,
            self.dual_coef_,
            self.intercept_,
            self._probA,
            self._probB,
            self.alpha_raw_,
            self.alpha_raw_lengths_,
            self.objective_val_,
            self.fit_status_,
            self._num_iter,
        ) = libsvm.fit(
            X,
            y,
            svm_type=solver_type,
            sample_weight=sample_weight,
            class_weight=getattr(self, "class_weight_", np.empty(0)),
            alpha_init=getattr(self, "alpha_init_", np.empty((0, 0))),
            kernel=kernel,
            C=self.C,
            nu=self.nu,
            probability=self.probability,
            degree=self.degree,
            shrinking=self.shrinking,
            tol=self.tol,
            cache_size=self.cache_size,
            coef0=self.coef0,
            gamma=self._gamma,
            epsilon=self.epsilon,
            max_iter=self.max_iter,
            random_seed=random_seed,
        )
        self._warn_from_fit_status()
