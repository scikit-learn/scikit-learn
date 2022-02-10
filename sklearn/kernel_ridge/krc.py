"""Module :mod:`sklearn.kernel_ridge.krc` implements kernel ridge classification."""

# Author: Carlos Perales <cperalesg@outlook.es>
# License: BSD 3 clause

import numpy as np
from ..kernel_ridge.krr import KernelRidge
from ..linear_model._ridge import RidgeClassifier
from ..linear_model._ridge import _solve_cholesky_kernel
from ..utils.validation import check_is_fitted, _check_sample_weight


class KernelRidgeClassifier(KernelRidge, RidgeClassifier):
    """Kernel ridge classification.

    Kernel ridge clasification (KRR) combines ridge regression (linear least
    squares with l2-norm regularization) with the kernel trick and the label
    binarization. It thus learns a linear function in the space induced
    y the respective kernel and the data. For non-linear kernels,
    this corresponds to a non-linear function in the original space.

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
    coef_ : ndarray of shape (1, n_features) or (n_classes, n_features)
        Coefficient of the features in the decision function.

        ``coef_`` is of shape (1, n_features) when the given problem is binary.

    classes_ : ndarray of shape (n_classes,)
        The classes labels.

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
    >>> krc = KernelRidgeClassifier()
    >>> krc.fit(X, y)
    KernelRidgeClassifier()
    """

    _estimator_type = "classifier"
    normalize = "deprecated"
    fit_intercept = False
    class_weight = None

    def __init__(
        self,
        alpha=1.0,
        *,
        kernel="rbf",
        gamma=None,
        degree=3,
        coef0=1,
        kernel_params=None,
    ):
        self.alpha = alpha
        self.kernel = kernel
        self.gamma = gamma
        self.degree = degree
        self.coef0 = coef0
        self.kernel_params = kernel_params

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
        X, y, sample_weight, Y = super(RidgeClassifier, self)._prepare_data(
            X, y, sample_weight, "auto"
        )

        if sample_weight is not None and not isinstance(sample_weight, float):
            sample_weight = _check_sample_weight(sample_weight, X)

        K = self._get_kernel(X)
        alpha = np.atleast_1d(self.alpha)

        copy = self.kernel == "precomputed"
        self.coef_ = _solve_cholesky_kernel(K, Y, alpha, sample_weight, copy)
        self.X_fit_ = X.copy()
        return self

    def decision_function(self, X):
        """
        Predict confidence scores for samples.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            The data matrix for which we want to get the confidence scores.

        Returns
        -------
        scores : ndarray of shape (n_samples,) or (n_samples, n_classes)
            Confidence scores per `(n_samples, n_classes)` combination. In the
            binary case, confidence score for `self.classes_[1]` where >0 means
            this class would be predicted.
        """
        check_is_fitted(self)
        X = self._validate_data(X, accept_sparse=("csr", "csc"), reset=False)
        K = self._get_kernel(X, self.X_fit_)
        scores = np.dot(K, self.coef_)
        return scores.ravel() if scores.shape[1] == 1 else scores

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
        return super(RidgeClassifier, self).predict(X)
