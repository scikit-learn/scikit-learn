"""Module :mod:`sklearn.kernel_ridge` implements kernel ridge regression."""

# Authors: Mathieu Blondel <mathieu@mblondel.org>
#          Jan Hendrik Metzen <jhm@informatik.uni-bremen.de>
# License: BSD 3 clause

import numpy as np

from .base import BaseEstimator, RegressorMixin, ClassifierMixin
from .metrics.pairwise import pairwise_kernels
from .linear_model.ridge import _solve_cholesky_kernel
from .utils import check_array, check_X_y
from .utils.validation import check_is_fitted
from .utils.multiclass import unique_labels


class KernelRidge(BaseEstimator, RegressorMixin, ClassifierMixin):
    """Kernel ridge regression and Kernel ridge classification
    (also known as Kernel Extreme Learning Machine).

    Kernel ridge regression (KRR) combines ridge regression (linear least
    squares with l2-norm regularization) with the kernel trick. It thus
    learns a linear function in the space induced by the respective kernel and
    the data. For non-linear kernels, this corresponds to a non-linear
    function in the original space.

    The form of the model learned by KRR is identical to support vector
    regression (SVR). However, different loss functions are used: KRR uses
    squared error loss while support vector regression uses epsilon-insensitive
    loss, both combined with l2 regularization. In contrast to SVR, fitting a
    KRR model can be done in closed-form and is typically faster for
    medium-sized datasets. On the other  hand, the learned model is non-sparse
    and thus slower than SVR, which learns a sparse model for epsilon > 0, at
    prediction-time.

    This estimator has built-in support for multi-variate regression
    (i.e., when y is a 2d-array of shape [n_samples, n_targets]).

    Read more in the :ref:`User Guide <kernel_ridge>`.

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
        Gamma parameter for the RBF, laplacian, polynomial, exponential chi2
        and sigmoid kernels. Interpretation of the default value is left to
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

    regression : bool, default=True
         Behaviour parameter, True if this object works as a
         regressor, False if as a classifier.
         When it is working as a classifier, behaviour is identical
         to Kernel Extreme Learning Machine.

    Attributes
    ----------
    F : array, shape = [n_samples] or [n_samples, n_targets]
        Representation of weight vector(s) in kernel space

    X_fit_ : {array-like, sparse matrix}, shape = [n_samples, n_features]
        Training data, which is also required for prediction

    References
    ----------
    * Kevin P. Murphy
      "Machine Learning: A Probabilistic Perspective", The MIT Press
      chapter 14.4.3, pp. 492-493

    See also
    --------
    sklearn.linear_model.Ridge:
        Linear ridge regression.
    sklearn.svm.SVR:
        Support Vector Regression implemented using libsvm.

    Examples
    --------
    >>> from sklearn.kernel_ridge import KernelRidge
    >>> import numpy as np
    >>> n_samples, n_features = 10, 5
    >>> rng = np.random.RandomState(0)
    >>> y = rng.randn(n_samples)
    >>> X = rng.randn(n_samples, n_features)
    >>> clf = KernelRidge(alpha=1.0)
    >>> clf.fit(X, y) # doctest: +NORMALIZE_WHITESPACE
    KernelRidge(alpha=1.0, coef0=1, degree=3, gamma=None, kernel='linear',
                kernel_params=None)
    """
    def __init__(self, alpha=1, kernel="linear", gamma=None, degree=3, coef0=1,
                 kernel_params=None, regression=True):
        self.alpha = alpha
        self.kernel = kernel
        self.gamma = gamma
        self.degree = degree
        self.coef0 = coef0
        self.kernel_params = kernel_params
        self.regression = regression

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
        try:
            X = check_array(X)
        except TypeError:
            raise ValueError('Predict with sparse input '
                             'when trained with dense')
        K = self._get_kernel(X, self.X_fit_)
        indicator = np.dot(K, self.dual_coef_)

        if self.regression is True:
            y = indicator
        else:
            # Decoding
            y = []
            for y_i in indicator:
                pos = np.argmax(y_i)
                y.append(self.class_corr_[str(pos)])
            y = np.array(y)

        return y
