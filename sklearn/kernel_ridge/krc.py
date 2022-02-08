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
    >>> krc = KernelRidgeClassifier(kernel='rbf', alpha=1.0, gamma=0.01)
    >>> krc.fit(X, y)
    KernelRidgeClassifier(alpha=0.1)
    """

    class_weight = None

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
        # Binarize and process as a classifier, using RidgeClassifier method
        X, y, sample_weight, Y = super(RidgeClassifier, self)._prepare_data(
            X, y, sample_weight, "auto"
        )
        # Fit with kernel from KernelRidge
        super().fit(X, Y, sample_weight)
        self.coef_ = self.dual_coef_
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
