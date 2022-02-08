"""Module :mod:`sklearn.kernel_ridge` implements kernel ridge regression."""

# Authors: Mathieu Blondel <mathieu@mblondel.org>
#          Jan Hendrik Metzen <jhm@informatik.uni-bremen.de>
# License: BSD 3 clause

import numpy as np
from ..kernel_ridge.krr import KernelRidge
from ..linear_model._ridge import RidgeClassifier
from ..metrics.pairwise import pairwise_kernels
from ..linear_model._ridge import _solve_cholesky_kernel
from ..utils.validation import check_is_fitted, _check_sample_weight
from ..utils.deprecation import deprecated
from ..preprocessing import LabelBinarizer
from ..utils import column_or_1d
from ..utils import compute_sample_weight


class KernelRidgeClassifier(KernelRidge, RidgeClassifier):
    """Kernel ridge classification.

    Kernel ridge clasification (KRR) combines ridge regression (linear least
    squares with l2-norm regularization) with the kernel trick. It thus
    learns a linear function in the space induced by the respective kernel and
    the data. For non-linear kernels, this corresponds to a non-linear
    function in the original space.

    The form of the model learned by KRR is identical to support vector
    regression (SVR). However, different loss functions are used: KRR uses
    squared error loss while support vector regression uses epsilon-insensitive
    loss, both combined with l2 regularization. In contrast to SVR, fitting a
    KRR model can be done in closed-form and is typically faster for
    medium-sized datasets. On the other hand, the learned model is non-sparse
    and thus slower than SVR, which learns a sparse model for epsilon > 0, at
    prediction-time.

    This estimator has built-in support for multi-variate regression
    (i.e., when y is a 2d-array of shape [n_samples, n_targets]).

    Read more in the :ref:`User Guide <kernel_ridge>`.

    Parameters
    ----------
    alpha : float or array-like of shape (n_targets,), default=1.0
        Regularization strength; must be a positive float. Regularization
        improves the conditioning of the problem and reduces the variance of
        the estimates. Larger values specify stronger regularization.
        Alpha corresponds to ``1 / (2C)`` in other linear models such as
        :class:`~sklearn.linear_model.LogisticRegression` or
        :class:`~sklearn.svm.LinearSVC`. If an array is passed, penalties are
        assumed to be specific to the targets. Hence they must correspond in
        number. See :ref:`ridge_regression` for formula.

    kernel : str or callable, default="linear"
        Kernel mapping used internally. This parameter is directly passed to
        :class:`~sklearn.metrics.pairwise.pairwise_kernel`.
        If `kernel` is a string, it must be one of the metrics
        in `pairwise.PAIRWISE_KERNEL_FUNCTIONS` or "precomputed".
        If `kernel` is "precomputed", X is assumed to be a kernel matrix.
        Alternatively, if `kernel` is a callable function, it is called on
        each pair of instances (rows) and the resulting value recorded. The
        callable should take two rows from X as input and return the
        corresponding kernel value as a single number. This means that
        callables from :mod:`sklearn.metrics.pairwise` are not allowed, as
        they operate on matrices, not single samples. Use the string
        identifying the kernel instead.

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

    kernel_params : mapping of str to any, default=None
        Additional parameters (keyword arguments) for kernel function passed
        as callable object.

    Attributes
    ----------
    dual_coef_ : ndarray of shape (n_samples,) or (n_samples, n_targets)
        Representation of weight vector(s) in kernel space

    X_fit_ : {ndarray, sparse matrix} of shape (n_samples, n_features)
        Training data, which is also required for prediction. If
        kernel == "precomputed" this is instead the precomputed
        training matrix, of shape (n_samples, n_samples).

    n_features_in_ : int
        Number of features seen during :term:`fit`.

        .. versionadded:: 0.24

    feature_names_in_ : ndarray of shape (`n_features_in_`,)
        Names of features seen during :term:`fit`. Defined only when `X`
        has feature names that are all strings.

        .. versionadded:: 1.0

    See Also
    --------
    sklearn.gaussian_process.GaussianProcessRegressor : Gaussian
        Process regressor providing automatic kernel hyperparameters
        tuning and predictions uncertainty.
    sklearn.linear_model.Ridge : Linear ridge regression.
    sklearn.linear_model.RidgeCV : Ridge regression with built-in
        cross-validation.
    sklearn.svm.SVR : Support Vector Regression accepting a large variety
        of kernels.

    References
    ----------
    * Kevin P. Murphy
      "Machine Learning: A Probabilistic Perspective", The MIT Press
      chapter 14.4.3, pp. 492-493

    Examples
    --------
    >>> from sklearn.kernel_ridge import KernelRidgeClassifier
    >>> from sklearn.datasets import load_iris
    >>> iris = load_iris()
    >>> X, y = iris.data, iris.target
    >>> n_samples, n_features = 10, 5
    >>> krc = KernelRidgeClassifier(alpha=1.0)
    >>> krc.fit(X, y)
    KernelRidge(alpha=1.0)
    """
    def _prepare_data(self, X, y, sample_weight):
        """Validate `X` and `y` and binarize `y`.

        Parameters
        ----------
        X : {ndarray, sparse matrix} of shape (n_samples, n_features)
            Training data.

        y : ndarray of shape (n_samples,)
            Target values.

        sample_weight : float or ndarray of shape (n_samples,), default=None
            Individual weights for each sample. If given a float, every sample
            will have the same weight.

        Returns
        -------
        X : {ndarray, sparse matrix} of shape (n_samples, n_features)
            Validated training data.

        y : ndarray of shape (n_samples,)
            Validated target values.

        sample_weight : ndarray of shape (n_samples,)
            Validated sample weights.

        Y : ndarray of shape (n_samples, n_classes)
            The binarized version of `y`.
        """
        X, y = self._validate_data(
            X, y, accept_sparse=("csr", "csc"), multi_output=True, y_numeric=True
        )
        if sample_weight is not None and not isinstance(sample_weight, float):
            sample_weight = _check_sample_weight(sample_weight, X)
    
        self._label_binarizer = LabelBinarizer(pos_label=1, neg_label=-1)
        Y = self._label_binarizer.fit_transform(y)
        if not self._label_binarizer.y_type_.startswith("multilabel"):
            y = column_or_1d(y, warn=True)

        sample_weight = _check_sample_weight(sample_weight, X, dtype=X.dtype)
        return X, y, sample_weight, Y

    def fit(self, X, y, sample_weight=None):
        """Fit Kernel Ridge regression model.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            Training data. If kernel == "precomputed" this is instead
            a precomputed kernel matrix, of shape (n_samples, n_samples).

        y : array-like of shape (n_samples,) or (n_samples, n_targets)
            Target values.

        sample_weight : float or array-like of shape (n_samples,), default=None
            Individual weights for each sample, ignored if None is passed.

        Returns
        -------
        self : object
            Returns the instance itself.
        """
        # Convert data
        X, y, sample_weight, Y = self._prepare_data(X, y, sample_weight)

        if sample_weight is not None and not isinstance(sample_weight, float):
            sample_weight = _check_sample_weight(sample_weight, X)

        K = self._get_kernel(X)
        alpha = np.atleast_1d(self.alpha)

        ravel = False
        if len(y.shape) == 1:
            y = y.reshape(-1, 1)
            ravel = True

        copy = self.kernel == "precomputed"
        self.dual_coef_ = _solve_cholesky_kernel(K, y, alpha, sample_weight, copy)
        if ravel:
            self.dual_coef_ = self.dual_coef_.ravel()

        self.X_fit_ = X

        return self

    def predict(self, X):
        """Predict using the kernel ridge model.

        

        Returns
        -------
        C : ndarray of shape (n_samples,) or (n_samples, n_targets)
            Returns predicted values.
        """
        check_is_fitted(self)
        X = self._validate_data(X, accept_sparse=("csr", "csc"), reset=False)
        K = self._get_kernel(X, self.X_fit_)
        return np.dot(K, self.dual_coef_)

    def predict(self, X):
        """Predict class labels for samples in `X`.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            Samples. If kernel == "precomputed" this is instead a
            precomputed kernel matrix, shape = [n_samples,
            n_samples_fitted], where n_samples_fitted is the number of
            samples used in the fitting for this estimator.

        Returns
        -------
        y_pred : ndarray of shape (n_samples,) or (n_samples, n_outputs)
            Vector or matrix containing the predictions. In binary and
            multiclass problems, this is a vector containing `n_samples`. In
            a multilabel problem, it returns a matrix of shape
            `(n_samples, n_outputs)`.
        """
        check_is_fitted(self, attributes=["_label_binarizer"])
        X = self._validate_data(X, accept_sparse=("csr", "csc"), reset=False)
        K = self._get_kernel(X, self.X_fit_)
        if self._label_binarizer.y_type_.startswith("multilabel"):
            # Threshold such that the negative label is -1 and positive label
            # is 1 to use the inverse transform of the label binarizer fitted
            # during fit.
            scores = 2 * (self.decision_function(X) > 0) - 1
            return self._label_binarizer.inverse_transform(scores)
        return super(RidgeClassifier).predict(X)
