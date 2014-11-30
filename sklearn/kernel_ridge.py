import numpy as np
from scipy import linalg

from sklearn.base import BaseEstimator, RegressorMixin
from sklearn.metrics.pairwise import pairwise_kernels


def _solve_dense_cholesky_kernel(K, y, alpha, sample_weight=1.0, copy=True):
    if copy:
        K = K.copy()

    # dual_coef = inv(X X^t + alpha*Id) y
    n_samples = K.shape[0]
    n_targets = y.shape[1]

    one_alpha = np.array_equal(alpha, len(alpha) * [alpha[0]])
    has_sw = isinstance(sample_weight, np.ndarray) or sample_weight != 1.0

    if has_sw:
        sw = np.sqrt(sample_weight)
        y = y * sw[:, np.newaxis]
        K *= np.outer(sw, sw)

    if one_alpha:
        # Only one penalty, we can solve multi-target problems in one time.
        K.flat[::n_samples + 1] += alpha[0]

        try:
            dual_coef = linalg.solve(K, y, sym_pos=True, overwrite_a=True)
        except:
            print "Warning: Singular Matrix"
            Kinv = linalg.pinv2(K)
            dual_coef = np.dot(Kinv, y)

        # K is expensive to compute and store in memory so change it back in
        # case it was user-given.
        K.flat[::n_samples + 1] -= alpha[0]

        if has_sw:
            dual_coef *= sw[:, np.newaxis]

        return dual_coef
    else:
        # One penalty per target. We need to solve each target separately.
        dual_coefs = np.empty([n_targets, n_samples])

        for dual_coef, target, current_alpha in zip(dual_coefs, y.T, alpha):
            K.flat[::n_samples + 1] += current_alpha

            dual_coef[:] = linalg.solve(K, target, sym_pos=True,
                                        overwrite_a=False).ravel()

            K.flat[::n_samples + 1] -= current_alpha

        if has_sw:
            dual_coefs *= sw[np.newaxis, :]

        return dual_coefs.T


class KernelRidge(BaseEstimator, RegressorMixin):

    def __init__(self, alpha=1, kernel="linear", gamma=None, degree=3, coef0=1,
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

    def fit(self, X, y=None, sample_weight=1.0):
        n_samples = X.shape[0]
        K = self._get_kernel(X)
        alpha = np.array([self.alpha])

        ravel = False
        if len(y.shape) == 1:
            y = y.reshape(-1, 1)
            ravel = True


        copy = self.kernel == "precomputed"
        self.dual_coef_ = _solve_dense_cholesky_kernel(K, y, alpha,
                                                       sample_weight,
                                                       copy)
        if ravel:
            self.dual_coef_ = self.dual_coef_.ravel()

        self.X_fit_ = X

        return self

    def predict(self, X):
        K = self._get_kernel(X, self.X_fit_)
        return np.dot(K, self.dual_coef_)
