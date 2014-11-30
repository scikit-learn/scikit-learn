"""The :mod:`sklearn.kernel_ridge` module implements kernel ridge regression."""

# Authors: Mathieu Blondel <mathieu@mblondel.org>
#          Jan Hendrik Metzen <jhm@informatik.uni-bremen.de>
# License: BSD 3 clause

import warnings

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
            # XXX: Check if overwrite_a=True is ok when copy is False
            dual_coef = linalg.solve(K, y, sym_pos=True, overwrite_a=True)
        except np.linalg.LinAlgError:
            warnings.warn("Singular matrix in solving dual problem. Using "
                          "pseudo inverse instead.")
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
    """Kernelized ridge regression.

    Kernelized ridge regression (KRR) combines ridge regression (linear least
    squares plus l2-norm  regularization) with the kernel trick. It thus learns
    a linear function in the Reproducing kernel Hilbert space induced by the
    respective kernel. This corresponds to a non-linear in the original space.

    The model learned by KRR is identical to support vector regression (SVR).
    However, different loss functions are used (ridge versus  epsilon-
    insensitive loss). In contrast to SVR, fitting a KRR can be done in closed-
    form and is typically faster for medium-sized datasets. On the other hand,
    the learned model is non-sparse and thus slower than SVR at prediction-time.

    Using KRR with a linear kernel can be advisable over non-kernelized Ridge
    in situations when the number of feature dimensions D is considerably larger
    than the number of training datapoints N since the computational cost of 
    solving the dual (kernelized) problem is O(N^3) while the standard Ridge
    requires O(D^3).

    This estimator has built-in support for multi-variate regression
    (i.e., when y is a 2d-array of shape [n_samples, n_targets]).

    Parameters
    ----------
    alpha : {float, array-like}, shape = [n_targets]
        Small positive values of alpha improve the conditioning of the problem
        and reduce the variance of the estimates.  Alpha corresponds to
        ``(2*C)^-1`` in other linear models such as LogisticRegression or
        LinearSVC. If an array is passed, penalties are assumed to be specific
        to the targets. Hence they must correspond in number.

    kernel : string or callable, default="linear"
        Kernel mapping used internally. A callable should accept two arguments
        and the keyword arguments passed to this object as kernel_params, and
        should return a floating point number.

    gamma : float, default=None
        Gamma parameter for the RBF, polynomial, exponential chi2 and
        sigmoid kernels. Interpretation of the default value is left to
        the kernel; see the documentation for sklearn.metrics.pairwise.
        Ignored by other kernels.

    degree : float, default=3
        Degree of the polynomial kernel. Ignored by other kernels.

    coef0 : float, default=1
        Zero coefficient for polynomial and sigmoid kernels.
        Ignored by other kernels.

    kernel_params : mapping of string to any, optional
        Additional parameters (keyword arguments) for kernel function passed
        as callable object.

    Attributes
    ----------
    dual_coef_ : array, shape = [n_features] or [n_targets, n_features]
        Weight vector(s) in kernel space

    X_fit_ : {array-like, sparse matrix}, shape = [n_samples, n_features]
        Training data, which is also required for prediction 

    References
    ----------
    * Kevin P. Murphy
      "Machine Learning: A Probabilistic Perspective", The MIT Press
      chapter 14.4.3, pp. 492-493

    See also
    --------
    Ridge
        Linear, non-kernelized least squares with l2 regularization.
    SVR
        Support Vector Regression implemented using libsvm.

    Examples
    --------
    >>> from sklearn.kernel_ridge import KernelRidge
    >>> import numpy as np
    >>> n_samples, n_features = 10, 5
    >>> np.random.seed(0)
    >>> y = np.random.randn(n_samples)
    >>> X = np.random.randn(n_samples, n_features)
    >>> clf = KernelRidge(alpha=1.0)
    >>> clf.fit(X, y) # doctest: +NORMALIZE_WHITESPACE
    KernelRidge(alpha=1.0, coef0=1, degree=3, gamma=None, kernel='linear',
                kernel_params=None)
    """
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
        """Fit Kernel Ridge regression model

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape = [n_samples, n_features]
            Training data

        y : array-like, shape = [n_samples] or [n_samples, n_targets]
            Target values

        sample_weight : float or numpy array of shape [n_samples]
            Individual weights for each sample

        Returns
        -------
        self : returns an instance of self.
        """
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
        """Predict using the the kernel ridge model

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape = [n_samples, n_features]
            Samples.

        Returns
        -------
        C : array, shape = [n_samples] or [n_samples, n_targets]
            Returns predicted values.
        """
        K = self._get_kernel(X, self.X_fit_)
        return np.dot(K, self.dual_coef_)
