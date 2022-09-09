import numpy as np
from sklearn.gaussian_process.kernels import Hyperparameter, Kernel


class CustomMinT(Kernel):
    """
    A custom kernel that has a diag method that returns the first column of the
    input matrix X. This is a helper for the test to check that the input
    matrix X is not mutated.
    """

    def __init__(self, sigma_0=1.0, sigma_0_bounds=(0.01, 10)):
        self.sigma_0 = sigma_0
        self.sigma_0_bounds = sigma_0_bounds

    @property
    def hyperparameter_sigma_0(self):
        return Hyperparameter("sigma_0", "numeric", self.sigma_0_bounds)

    def __call__(self, X, Y=None, eval_gradient=False):
        if Y is not None and eval_gradient:
            raise ValueError("Gradient can only be evaluated when Y is None.")

        X = np.atleast_2d(X)
        ones_x = np.ones_like(X)

        if Y is None:
            kc, kr = X * ones_x.T, ones_x * X.T
        else:
            ones_y = np.ones_like(Y)
            kc, kr = X * ones_y.T, ones_x * Y.T

        kcr = np.concatenate((kc[..., None], kr[..., None]), axis=-1)
        k = np.min(kcr, axis=-1)

        if eval_gradient:
            if not self.hyperparameter_sigma_0.fixed:
                k_gradient = np.empty((k.shape[0], k.shape[1], 1))
                k_gradient[..., 0] = self.sigma_0
                return k, k_gradient
            else:
                return k, np.empty((X.shape[0], X.shape[0], 0))
        else:
            return k

    def diag(self, X):
        return X[:, 0]

    def is_stationary(self):
        return False
