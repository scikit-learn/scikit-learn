"""Module :mod:`sklearn.kernel_ridge` implements kernel ridge regression."""

# Authors: Mathieu Blondel <mathieu@mblondel.org>
#          Jan Hendrik Metzen <jhm@informatik.uni-bremen.de>
#          Carlos Perales <sir.perales@gmail.com>
# License: BSD 3 clause

import numpy as np

from ..metrics.pairwise import pairwise_kernels
from ..linear_model.ridge import _solve_cholesky_kernel
from ..utils import check_array, check_X_y
from ..utils.validation import check_is_fitted
from ..preprocessing import LabelBinarizer


class _BaseKernelRidge(object):
    """Kernel ridge regression and Kernel Ridge Classifier.

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

    Attributes
    ----------
    dual_coef_ : array, shape = [n_samples] or [n_samples, n_targets]
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
    sklearn.kernel_ridge.KernelRidgeClassifier:
        Kernel Ridge implemented for multiclass clasifications.

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

        ravel = False
        if self._estimator_type == "regressor":
            X, y_ = check_X_y(X, y, accept_sparse=("csr", "csc"),
                              multi_output=True, y_numeric=True)
            if len(y_.shape) == 1:
                y_ = y_.reshape(-1, 1)
                ravel = True
        elif self._estimator_type == "classifier":  # Multiclass classification
            X, y = check_X_y(X, y)
            # Check alpha
            if self.alpha <= 0:
                raise ValueError('alpha must be positive')

            # y is going to be encoded by zero arrays
            # with 1 in a position assigned to a label
            self.label_encoder_ = LabelBinarizer()
            y_ = self.label_encoder_.fit_transform(y)
            self.classes_ = self.label_encoder_.classes_
        else:
            raise ValueError('BaseKernelRidge cannot be used as'
                             ' a estimator')

        self.dual_coef_ = _solve_cholesky_kernel(K, y_, alpha,
                                                 sample_weight,
                                                 copy)

        if ravel:
            self.dual_coef_ = self.dual_coef_.ravel()

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
        y = self.decision_function(X=X)

        if self._estimator_type == "classifier":
            # Decoding
            y = self.label_encoder_.inverse_transform(y)

        return y

    def decision_function(self, X):
        """Predict confidence scores for samples.

        The confidence score for a sample is the signed distance of that
        sample to the hyperplane.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape = (n_samples, n_features)
            Samples.

        Returns
        -------
        array, shape=(n_samples,) if n_classes == 2 else (n_samples, n_classes)
            Confidence scores per (sample, class) combination. In the binary
            case, confidence score for self.classes_[1] where >0 means this
            class would be predicted.
        """
        check_is_fitted(self, ["X_fit_", "dual_coef_"])
        # Input validation
        K = self._get_kernel(X, self.X_fit_)
        y = np.dot(K, self.dual_coef_)
        return y
