"""Module :mod:`sklearn.kernel_ridge` implements kernel ridge regression."""

# Authors: Mathieu Blondel <mathieu@mblondel.org>
#          Jan Hendrik Metzen <jhm@informatik.uni-bremen.de>
#          Carlos Perales <sir.perales@gmail.com>
# License: BSD 3 clause

import numpy as np

from ..base import BaseEstimator, RegressorMixin
from ..metrics.pairwise import pairwise_kernels
from ..linear_model.ridge import _solve_cholesky_kernel
from ..utils import check_array, check_X_y
from ..utils.validation import check_is_fitted
from ..utils.multiclass import unique_labels


class KernelRidge(BaseEstimator, RegressorMixin):
    """Kernel ridge regression and Kernel ridge classification
    (also known as Kernel Extreme Learning Machine).
    """
    # By default, it works as a regressor
    regression = True

    def __init__(self, alpha=1, kernel="rbf", gamma=None, degree=3, coef0=1,
                 kernel_params=None):
        self.alpha = alpha
        self.kernel = kernel
        self.gamma = gamma
        self.degree = degree
        self.coef0 = coef0
        self.kernel_params = kernel_params

    def _get_kernel(self, X, Y=None):
        if callable(self.kernel):
            params = self.kernel_params or {}
        else:
            params = {"gamma": self.gamma,
                      "degree": self.degree,
                      "coef0": self.coef0}
        return pairwise_kernels(X, Y, metric=self.kernel,
                                filter_params=True, **params)

    @property
    def _pairwise(self):
        return self.kernel == "precomputed"

    def fit(self, X, y=None, sample_weight=None):
        """Fit Kernel Ridge regression model

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape = [n_samples, n_features]
            Training data

        y : array-like, shape = [n_samples] or [n_samples, n_targets]
            Target values

        sample_weight : float or array-like of shape [n_samples]
            Individual weights for each sample, ignored if None is passed.

        Returns
        -------
        self : returns an instance of self.
        """
        # Convert data
        if sample_weight is not None and not isinstance(sample_weight, float):
            sample_weight = check_array(sample_weight, ensure_2d=False)

        K = self._get_kernel(X)
        alpha = np.atleast_1d(self.alpha)

        copy = self.kernel == "precomputed"

        if self.regression is True:
            X, y = check_X_y(X, y, accept_sparse=("csr", "csc"), multi_output=True,
                             y_numeric=True)
            ravel = False
            if len(y.shape) == 1:
                y = y.reshape(-1, 1)
                ravel = True

            self.dual_coef_ = _solve_cholesky_kernel(K, y, alpha,
                                                     sample_weight,
                                                     copy)
            if ravel:
                self.dual_coef_ = self.dual_coef_.ravel()

        else:
            # Check alpha
            if self.alpha <= 0:
                raise ValueError('C must be positive')
            X, y = check_X_y(X, y)
            self.classes_ = unique_labels(y)

            # y must be encoded by zero arrays with 1 in a position
            # assigned to a label
            n = X.shape[0]
            if self.classes_.shape[0] == 1:
                self.n_classes_ = int(self.classes_[0] + 1)
                self.classes_ = np.arange(self.n_classes_)
            else:
                self.n_classes_ = len(self.classes_)

            T = np.zeros((n, self.n_classes_), dtype=np.float64)
            # It is essential in order to adapt it to string labels
            self.class_corr_ = {}
            for i in range(n):
                row = [y[i] == self.classes_]
                T[i] = np.array(row, dtype=np.float64)
                pos = np.argmax(row)
                self.class_corr_[str(pos)] = self.classes_[pos]
                # T is already encoded

            # # ELM
            # alpha = np.eye(n) * 2 * self.alpha + K
            # self.dual_coef_ = np.linalg.solve(alpha, T)

            # Ridge
            self.dual_coef_ = _solve_cholesky_kernel(K, T, alpha,
                                                     sample_weight,
                                                     copy)

            # This value is similar to number of neurons in hidden
            # layer in neural version of Extreme Learning Machine
            self.h_ = self.dual_coef_.shape[0]

        self.X_fit_ = X

        return self

    def predict(self, X):
        """Predict using the kernel ridge model

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape = [n_samples, n_features]
            Samples.

        Returns
        -------
        C : array, shape = [n_samples] or [n_samples, n_targets]
            Returns predicted values.
        """
        check_is_fitted(self, ["X_fit_", "dual_coef_"])
        # Input validation
        K = self._get_kernel(X, self.X_fit_)
        indicator = np.dot(K, self.dual_coef_)

        if self.regression is True:
            y = indicator
        else:
            try:
                X = check_array(X)
            except TypeError:
                raise ValueError('Predict with sparse input '
                                 'when trained with dense')
            # Decoding
            y = []
            for y_i in indicator:
                pos = np.argmax(y_i)
                y.append(self.class_corr_[str(pos)])
            y = np.array(y)

        return y
