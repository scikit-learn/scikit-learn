"""
Generalized Linear Models.
"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import warnings
from abc import ABCMeta, abstractmethod
from numbers import Integral, Real

import numpy as np
import scipy.sparse as sp
from scipy import linalg, optimize, sparse
from scipy.sparse.linalg import lsqr
from scipy.special import expit

from sklearn.base import (
    BaseEstimator,
    ClassifierMixin,
    MultiOutputMixin,
    RegressorMixin,
    _fit_context,
)
from sklearn.utils import check_array, check_random_state
from sklearn.utils._array_api import (
    _asarray_with_order,
    _average,
    get_namespace,
    get_namespace_and_device,
    indexing_dtype,
    supported_float_dtypes,
)
from sklearn.utils._param_validation import Interval
from sklearn.utils._seq_dataset import (
    ArrayDataset32,
    ArrayDataset64,
    CSRDataset32,
    CSRDataset64,
)
from sklearn.utils.extmath import safe_sparse_dot
from sklearn.utils.parallel import Parallel, delayed
from sklearn.utils.sparsefuncs import mean_variance_axis
from sklearn.utils.validation import (
    _check_sample_weight,
    check_is_fitted,
    validate_data,
)

# TODO: bayesian_ridge_regression and bayesian_regression_ard
# should be squashed into its respective objects.

SPARSE_INTERCEPT_DECAY = 0.01
# For sparse data intercept updates are scaled by this decay factor to avoid
# intercept oscillation.


def make_dataset(X, y, sample_weight, random_state=None):
    """Create ``Dataset`` abstraction for sparse and dense inputs.

    This also returns the ``intercept_decay`` which is different
    for sparse datasets.

    Parameters
    ----------
    X : array-like, shape (n_samples, n_features)
        Training data

    y : array-like, shape (n_samples, )
        Target values.

    sample_weight : numpy array of shape (n_samples,)
        The weight of each sample

    random_state : int, RandomState instance or None (default)
        Determines random number generation for dataset random sampling. It is not
        used for dataset shuffling.
        Pass an int for reproducible output across multiple function calls.
        See :term:`Glossary <random_state>`.

    Returns
    -------
    dataset
        The ``Dataset`` abstraction
    intercept_decay
        The intercept decay
    """

    rng = check_random_state(random_state)
    # seed should never be 0 in SequentialDataset64
    seed = rng.randint(1, np.iinfo(np.int32).max)

    if X.dtype == np.float32:
        CSRData = CSRDataset32
        ArrayData = ArrayDataset32
    else:
        CSRData = CSRDataset64
        ArrayData = ArrayDataset64

    if sp.issparse(X):
        dataset = CSRData(X.data, X.indptr, X.indices, y, sample_weight, seed=seed)
        intercept_decay = SPARSE_INTERCEPT_DECAY
    else:
        X = np.ascontiguousarray(X)
        dataset = ArrayData(X, y, sample_weight, seed=seed)
        intercept_decay = 1.0

    return dataset, intercept_decay


def _preprocess_data(
    X,
    y,
    *,
    fit_intercept,
    copy=True,
    copy_y=True,
    sample_weight=None,
    check_input=True,
    rescale_with_sw=True,
):
    """Common data preprocessing for fitting linear models.

    This helper is in charge of the following steps:

    - `sample_weight` is assumed to be `None` or a validated array with same dtype as
      `X`.
    - If `check_input=True`, perform standard input validation of `X`, `y`.
    - Perform copies if requested to avoid side-effects in case of inplace
      modifications of the input.

    Then, if `fit_intercept=True` this preprocessing centers both `X` and `y` as
    follows:
        - if `X` is dense, center the data and
        store the mean vector in `X_offset`.
        - if `X` is sparse, store the mean in `X_offset`
        without centering `X`. The centering is expected to be handled by the
        linear solver where appropriate.
        - in either case, always center `y` and store the mean in `y_offset`.
        - both `X_offset` and `y_offset` are always weighted by `sample_weight`
          if not set to `None`.

    If `fit_intercept=False`, no centering is performed and `X_offset`, `y_offset`
    are set to zero.

    If `rescale_with_sw` is True, then X and y are rescaled with the square root of
    sample weights.

    Returns
    -------
    X_out : {ndarray, sparse matrix} of shape (n_samples, n_features)
        If copy=True a copy of the input X is triggered, otherwise operations are
        inplace.
        If input X is dense, then X_out is centered.
    y_out : {ndarray, sparse matrix} of shape (n_samples,) or (n_samples, n_targets)
        Centered version of y. Possibly performed inplace on input y depending
        on the copy_y parameter.
    X_offset : ndarray of shape (n_features,)
        The mean per column of input X.
    y_offset : float or ndarray of shape (n_features,)
    X_scale : ndarray of shape (n_features,)
        Always an array of ones. TODO: refactor the code base to make it
        possible to remove this unused variable.
    sample_weight_sqrt : ndarray of shape (n_samples, ) or None
        `np.sqrt(sample_weight)`
    """
    xp, _, device_ = get_namespace_and_device(X, y, sample_weight)
    n_samples, n_features = X.shape
    X_is_sparse = sp.issparse(X)

    if check_input:
        X = check_array(
            X, copy=copy, accept_sparse=["csr", "csc"], dtype=supported_float_dtypes(xp)
        )
        y = check_array(y, dtype=X.dtype, copy=copy_y, ensure_2d=False)
    else:
        y = xp.astype(y, X.dtype, copy=copy_y)
        if copy:
            if X_is_sparse:
                X = X.copy()
            else:
                X = _asarray_with_order(X, order="K", copy=True, xp=xp)

    dtype_ = X.dtype

    if fit_intercept:
        if X_is_sparse:
            X_offset, X_var = mean_variance_axis(X, axis=0, weights=sample_weight)
        else:
            X_offset = _average(X, axis=0, weights=sample_weight, xp=xp)

            X_offset = xp.astype(X_offset, X.dtype, copy=False)
            X -= X_offset

        y_offset = _average(y, axis=0, weights=sample_weight, xp=xp)
        y -= y_offset
    else:
        X_offset = xp.zeros(n_features, dtype=X.dtype, device=device_)
        if y.ndim == 1:
            y_offset = xp.asarray(0.0, dtype=dtype_, device=device_)
        else:
            y_offset = xp.zeros(y.shape[1], dtype=dtype_, device=device_)

    # X_scale is no longer needed. It is a historic artifact from the
    # time where linear model exposed the normalize parameter.
    X_scale = xp.ones(n_features, dtype=X.dtype, device=device_)

    if sample_weight is not None and rescale_with_sw:
        # Sample weight can be implemented via a simple rescaling.
        # For sparse X and y, it triggers copies anyway.
        # For dense X and y that already have been copied, we safely do inplace
        # rescaling.
        X, y, sample_weight_sqrt = _rescale_data(X, y, sample_weight, inplace=copy)
    else:
        sample_weight_sqrt = None
    return X, y, X_offset, y_offset, X_scale, sample_weight_sqrt


def _rescale_data(X, y, sample_weight, inplace=False):
    """Rescale data sample-wise by square root of sample_weight.

    For many linear models, this enables easy support for sample_weight because

        (y - X w)' S (y - X w)

    with S = diag(sample_weight) becomes

        ||y_rescaled - X_rescaled w||_2^2

    when setting

        y_rescaled = sqrt(S) y
        X_rescaled = sqrt(S) X

    The parameter `inplace` only takes effect for dense X and dense y.

    Returns
    -------
    X_rescaled : {array-like, sparse matrix}

    y_rescaled : {array-like, sparse matrix}

    sample_weight_sqrt : array-like of shape (n_samples,)
    """
    # Assume that _validate_data and _check_sample_weight have been called by
    # the caller.
    xp, _ = get_namespace(X, y, sample_weight)
    n_samples = X.shape[0]
    sample_weight_sqrt = xp.sqrt(sample_weight)

    if sp.issparse(X) or sp.issparse(y):
        sw_matrix = sparse.dia_matrix(
            (sample_weight_sqrt, 0), shape=(n_samples, n_samples)
        )

    if sp.issparse(X):
        X = safe_sparse_dot(sw_matrix, X)
    else:
        if inplace:
            X *= sample_weight_sqrt[:, None]
        else:
            X = X * sample_weight_sqrt[:, None]

    if sp.issparse(y):
        y = safe_sparse_dot(sw_matrix, y)
    else:
        if inplace:
            if y.ndim == 1:
                y *= sample_weight_sqrt
            else:
                y *= sample_weight_sqrt[:, None]
        else:
            if y.ndim == 1:
                y = y * sample_weight_sqrt
            else:
                y = y * sample_weight_sqrt[:, None]
    return X, y, sample_weight_sqrt


class LinearModel(BaseEstimator, metaclass=ABCMeta):
    """Base class for Linear Models"""

    @abstractmethod
    def fit(self, X, y):
        """Fit model."""

    def _decision_function(self, X):
        check_is_fitted(self)

        X = validate_data(self, X, accept_sparse=["csr", "csc", "coo"], reset=False)
        coef_ = self.coef_
        if coef_.ndim == 1:
            return X @ coef_ + self.intercept_
        else:
            return X @ coef_.T + self.intercept_

    def predict(self, X):
        """
        Predict using the linear model.

        Parameters
        ----------
        X : array-like or sparse matrix, shape (n_samples, n_features)
            Samples.

        Returns
        -------
        C : array, shape (n_samples,)
            Returns predicted values.
        """
        return self._decision_function(X)

    def _set_intercept(self, X_offset, y_offset, X_scale=None):
        """Set the intercept_"""
        xp, _ = get_namespace(X_offset, y_offset, X_scale)

        if self.fit_intercept:
            # We always want coef_.dtype=X.dtype. For instance, X.dtype can differ from
            # coef_.dtype if warm_start=True.
            self.coef_ = xp.astype(self.coef_, X_offset.dtype, copy=False)
            if X_scale is not None:
                self.coef_ = xp.divide(self.coef_, X_scale)

            if self.coef_.ndim == 1:
                self.intercept_ = y_offset - X_offset @ self.coef_
            else:
                self.intercept_ = y_offset - X_offset @ self.coef_.T

        else:
            self.intercept_ = 0.0


# XXX Should this derive from LinearModel? It should be a mixin, not an ABC.
# Maybe the n_features checking can be moved to LinearModel.
class LinearClassifierMixin(ClassifierMixin):
    """Mixin for linear classifiers.

    Handles prediction for sparse and dense X.
    """

    def decision_function(self, X):
        """
        Predict confidence scores for samples.

        The confidence score for a sample is proportional to the signed
        distance of that sample to the hyperplane.

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
        xp, _ = get_namespace(X)

        X = validate_data(self, X, accept_sparse="csr", reset=False)
        coef_T = self.coef_.T if self.coef_.ndim == 2 else self.coef_
        scores = safe_sparse_dot(X, coef_T, dense_output=True) + self.intercept_
        return (
            xp.reshape(scores, (-1,))
            if (scores.ndim > 1 and scores.shape[1] == 1)
            else scores
        )

    def predict(self, X):
        """
        Predict class labels for samples in X.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            The data matrix for which we want to get the predictions.

        Returns
        -------
        y_pred : ndarray of shape (n_samples,)
            Vector containing the class labels for each sample.
        """
        xp, _ = get_namespace(X)
        scores = self.decision_function(X)
        if len(scores.shape) == 1:
            indices = xp.astype(scores > 0, indexing_dtype(xp))
        else:
            indices = xp.argmax(scores, axis=1)

        return xp.take(self.classes_, indices, axis=0)

    def _predict_proba_lr(self, X):
        """Probability estimation for OvR logistic regression.

        Positive class probabilities are computed as
        1. / (1. + np.exp(-self.decision_function(X)));
        multiclass is handled by normalizing that over all classes.
        """
        prob = self.decision_function(X)
        expit(prob, out=prob)
        if prob.ndim == 1:
            return np.vstack([1 - prob, prob]).T
        else:
            # OvR normalization, like LibLinear's predict_probability
            prob /= prob.sum(axis=1).reshape((prob.shape[0], -1))
            return prob


class SparseCoefMixin:
    """Mixin for converting coef_ to and from CSR format.

    L1-regularizing estimators should inherit this.
    """

    def densify(self):
        """
        Convert coefficient matrix to dense array format.

        Converts the ``coef_`` member (back) to a numpy.ndarray. This is the
        default format of ``coef_`` and is required for fitting, so calling
        this method is only required on models that have previously been
        sparsified; otherwise, it is a no-op.

        Returns
        -------
        self
            Fitted estimator.
        """
        msg = "Estimator, %(name)s, must be fitted before densifying."
        check_is_fitted(self, msg=msg)
        if sp.issparse(self.coef_):
            self.coef_ = self.coef_.toarray()
        return self

    def sparsify(self):
        """
        Convert coefficient matrix to sparse format.

        Converts the ``coef_`` member to a scipy.sparse matrix, which for
        L1-regularized models can be much more memory- and storage-efficient
        than the usual numpy.ndarray representation.

        The ``intercept_`` member is not converted.

        Returns
        -------
        self
            Fitted estimator.

        Notes
        -----
        For non-sparse models, i.e. when there are not many zeros in ``coef_``,
        this may actually *increase* memory usage, so use this method with
        care. A rule of thumb is that the number of zero elements, which can
        be computed with ``(coef_ == 0).sum()``, must be more than 50% for this
        to provide significant benefits.

        After calling this method, further fitting with the partial_fit
        method (if any) will not work until you call densify.
        """
        msg = "Estimator, %(name)s, must be fitted before sparsifying."
        check_is_fitted(self, msg=msg)
        self.coef_ = sp.csr_matrix(self.coef_)
        return self


class LinearRegression(MultiOutputMixin, RegressorMixin, LinearModel):
    """
    Ordinary least squares Linear Regression.

    LinearRegression fits a linear model with coefficients w = (w1, ..., wp)
    to minimize the residual sum of squares between the observed targets in
    the dataset, and the targets predicted by the linear approximation.

    Parameters
    ----------
    fit_intercept : bool, default=True
        Whether to calculate the intercept for this model. If set
        to False, no intercept will be used in calculations
        (i.e. data is expected to be centered).

    copy_X : bool, default=True
        If True, X will be copied; else, it may be overwritten.

    tol : float, default=1e-6
        The precision of the solution (`coef_`) is determined by `tol` which
        specifies a different convergence criterion for the `lsqr` solver.
        `tol` is set as `atol` and `btol` of :func:`scipy.sparse.linalg.lsqr` when
        fitting on sparse training data. This parameter has no effect when fitting
        on dense data.

        .. versionadded:: 1.7

    n_jobs : int, default=None
        The number of jobs to use for the computation. This will only provide
        speedup in case of sufficiently large problems, that is if firstly
        `n_targets > 1` and secondly `X` is sparse or if `positive` is set
        to `True`. ``None`` means 1 unless in a
        :obj:`joblib.parallel_backend` context. ``-1`` means using all
        processors. See :term:`Glossary <n_jobs>` for more details.

    positive : bool, default=False
        When set to ``True``, forces the coefficients to be positive. This
        option is only supported for dense arrays.

        For a comparison between a linear regression model with positive constraints
        on the regression coefficients and a linear regression without such constraints,
        see :ref:`sphx_glr_auto_examples_linear_model_plot_nnls.py`.

        .. versionadded:: 0.24

    Attributes
    ----------
    coef_ : array of shape (n_features, ) or (n_targets, n_features)
        Estimated coefficients for the linear regression problem.
        If multiple targets are passed during the fit (y 2D), this
        is a 2D array of shape (n_targets, n_features), while if only
        one target is passed, this is a 1D array of length n_features.

    rank_ : int
        Rank of matrix `X`. Only available when `X` is dense.

    singular_ : array of shape (min(X, y),)
        Singular values of `X`. Only available when `X` is dense.

    intercept_ : float or array of shape (n_targets,)
        Independent term in the linear model. Set to 0.0 if
        `fit_intercept = False`.

    n_features_in_ : int
        Number of features seen during :term:`fit`.

        .. versionadded:: 0.24

    feature_names_in_ : ndarray of shape (`n_features_in_`,)
        Names of features seen during :term:`fit`. Defined only when `X`
        has feature names that are all strings.

        .. versionadded:: 1.0

    See Also
    --------
    Ridge : Ridge regression addresses some of the
        problems of Ordinary Least Squares by imposing a penalty on the
        size of the coefficients with l2 regularization.
    Lasso : The Lasso is a linear model that estimates
        sparse coefficients with l1 regularization.
    ElasticNet : Elastic-Net is a linear regression
        model trained with both l1 and l2 -norm regularization of the
        coefficients.

    Notes
    -----
    From the implementation point of view, this is just plain Ordinary
    Least Squares (:func:`scipy.linalg.lstsq`) or Non Negative Least Squares
    (:func:`scipy.optimize.nnls`) wrapped as a predictor object.

    Examples
    --------
    >>> import numpy as np
    >>> from sklearn.linear_model import LinearRegression
    >>> X = np.array([[1, 1], [1, 2], [2, 2], [2, 3]])
    >>> # y = 1 * x_0 + 2 * x_1 + 3
    >>> y = np.dot(X, np.array([1, 2])) + 3
    >>> reg = LinearRegression().fit(X, y)
    >>> reg.score(X, y)
    1.0
    >>> reg.coef_
    array([1., 2.])
    >>> reg.intercept_
    np.float64(3.0)
    >>> reg.predict(np.array([[3, 5]]))
    array([16.])
    """

    _parameter_constraints: dict = {
        "fit_intercept": ["boolean"],
        "copy_X": ["boolean"],
        "n_jobs": [None, Integral],
        "positive": ["boolean"],
        "tol": [Interval(Real, 0, None, closed="left")],
    }

    def __init__(
        self,
        *,
        fit_intercept=True,
        copy_X=True,
        tol=1e-6,
        n_jobs=None,
        positive=False,
    ):
        self.fit_intercept = fit_intercept
        self.copy_X = copy_X
        self.tol = tol
        self.n_jobs = n_jobs
        self.positive = positive

    @_fit_context(prefer_skip_nested_validation=True)
    def fit(self, X, y, sample_weight=None):
        """
        Fit linear model.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            Training data.

        y : array-like of shape (n_samples,) or (n_samples, n_targets)
            Target values. Will be cast to X's dtype if necessary.

        sample_weight : array-like of shape (n_samples,), default=None
            Individual weights for each sample.

            .. versionadded:: 0.17
               parameter *sample_weight* support to LinearRegression.

        Returns
        -------
        self : object
            Fitted Estimator.
        """
        n_jobs_ = self.n_jobs

        accept_sparse = False if self.positive else ["csr", "csc", "coo"]

        X, y = validate_data(
            self,
            X,
            y,
            accept_sparse=accept_sparse,
            y_numeric=True,
            multi_output=True,
            force_writeable=True,
        )

        has_sw = sample_weight is not None
        if has_sw:
            sample_weight = _check_sample_weight(
                sample_weight, X, dtype=X.dtype, ensure_non_negative=True
            )

        # Note that neither _rescale_data nor the rest of the fit method of
        # LinearRegression can benefit from in-place operations when X is a
        # sparse matrix. Therefore, let's not copy X when it is sparse.
        copy_X_in_preprocess_data = self.copy_X and not sp.issparse(X)

        X, y, X_offset, y_offset, _, sample_weight_sqrt = _preprocess_data(
            X,
            y,
            fit_intercept=self.fit_intercept,
            copy=copy_X_in_preprocess_data,
            sample_weight=sample_weight,
        )

        if self.positive:
            if y.ndim < 2:
                self.coef_ = optimize.nnls(X, y)[0]
            else:
                # scipy.optimize.nnls cannot handle y with shape (M, K)
                outs = Parallel(n_jobs=n_jobs_)(
                    delayed(optimize.nnls)(X, y[:, j]) for j in range(y.shape[1])
                )
                self.coef_ = np.vstack([out[0] for out in outs])
        elif sp.issparse(X):
            if has_sw:

                def matvec(b):
                    return X.dot(b) - sample_weight_sqrt * b.dot(X_offset)

                def rmatvec(b):
                    return X.T.dot(b) - X_offset * b.dot(sample_weight_sqrt)

            else:

                def matvec(b):
                    return X.dot(b) - b.dot(X_offset)

                def rmatvec(b):
                    return X.T.dot(b) - X_offset * b.sum()

            X_centered = sparse.linalg.LinearOperator(
                shape=X.shape, matvec=matvec, rmatvec=rmatvec
            )

            if y.ndim < 2:
                self.coef_ = lsqr(X_centered, y, atol=self.tol, btol=self.tol)[0]
            else:
                # sparse_lstsq cannot handle y with shape (M, K)
                outs = Parallel(n_jobs=n_jobs_)(
                    delayed(lsqr)(
                        X_centered, y[:, j].ravel(), atol=self.tol, btol=self.tol
                    )
                    for j in range(y.shape[1])
                )
                self.coef_ = np.vstack([out[0] for out in outs])
        else:
            # cut-off ratio for small singular values
            cond = max(X.shape) * np.finfo(X.dtype).eps
            self.coef_, _, self.rank_, self.singular_ = linalg.lstsq(X, y, cond=cond)
            self.coef_ = self.coef_.T

        if y.ndim == 1:
            self.coef_ = np.ravel(self.coef_)
        self._set_intercept(X_offset, y_offset)
        return self

    def __sklearn_tags__(self):
        tags = super().__sklearn_tags__()
        tags.input_tags.sparse = not self.positive
        return tags


def _check_precomputed_gram_matrix(
    X, precompute, X_offset, X_scale, rtol=None, atol=1e-5
):
    """Computes a single element of the gram matrix and compares it to
    the corresponding element of the user supplied gram matrix.

    If the values do not match a ValueError will be thrown.

    Parameters
    ----------
    X : ndarray of shape (n_samples, n_features)
        Data array.

    precompute : array-like of shape (n_features, n_features)
        User-supplied gram matrix.

    X_offset : ndarray of shape (n_features,)
        Array of feature means used to center design matrix.

    X_scale : ndarray of shape (n_features,)
        Array of feature scale factors used to normalize design matrix.

    rtol : float, default=None
        Relative tolerance; see numpy.allclose
        If None, it is set to 1e-4 for arrays of dtype numpy.float32 and 1e-7
        otherwise.

    atol : float, default=1e-5
        absolute tolerance; see :func`numpy.allclose`. Note that the default
        here is more tolerant than the default for
        :func:`numpy.testing.assert_allclose`, where `atol=0`.

    Raises
    ------
    ValueError
        Raised when the provided Gram matrix is not consistent.
    """

    n_features = X.shape[1]
    f1 = n_features // 2
    f2 = min(f1 + 1, n_features - 1)

    v1 = (X[:, f1] - X_offset[f1]) * X_scale[f1]
    v2 = (X[:, f2] - X_offset[f2]) * X_scale[f2]

    expected = np.dot(v1, v2)
    actual = precompute[f1, f2]

    dtypes = [precompute.dtype, expected.dtype]
    if rtol is None:
        rtols = [1e-4 if dtype == np.float32 else 1e-7 for dtype in dtypes]
        rtol = max(rtols)

    if not np.isclose(expected, actual, rtol=rtol, atol=atol):
        raise ValueError(
            "Gram matrix passed in via 'precompute' parameter "
            "did not pass validation when a single element was "
            "checked - please check that it was computed "
            f"properly. For element ({f1},{f2}) we computed "
            f"{expected} but the user-supplied value was "
            f"{actual}."
        )


def _pre_fit(
    X,
    y,
    Xy,
    precompute,
    fit_intercept,
    copy,
    check_gram=True,
    sample_weight=None,
):
    """Function used at beginning of fit in linear models with L1 or L0 penalty.

    This function applies _preprocess_data and additionally computes the gram matrix
    `precompute` as needed as well as `Xy`.

    It is assumed that X, y and sample_weight are already validated.

    Returns
    -------
    X
    y
    X_offset
    y_offset
    X_scale
    precompute
    Xy
    """
    n_samples, n_features = X.shape

    if sparse.issparse(X):
        # copy is not needed here as X is not modified inplace when X is sparse
        copy = False
        precompute = False
        # Rescale X and y only in dense case. Sparse cd solver directly deals with
        # sample_weight.
        rescale_with_sw = False
    else:
        # copy was done in fit if necessary
        rescale_with_sw = True

    X, y, X_offset, y_offset, X_scale, _ = _preprocess_data(
        X,
        y,
        fit_intercept=fit_intercept,
        copy=copy,
        sample_weight=sample_weight,
        check_input=False,
        rescale_with_sw=rescale_with_sw,
    )

    if hasattr(precompute, "__array__"):
        if fit_intercept and not np.allclose(X_offset, np.zeros(n_features)):
            warnings.warn(
                (
                    "Gram matrix was provided but X was centered to fit "
                    "intercept: recomputing Gram matrix."
                ),
                UserWarning,
            )
            # TODO: instead of warning and recomputing, we could just center
            # the user provided Gram matrix a-posteriori (after making a copy
            # when `copy=True`).
            # recompute Gram
            precompute = "auto"
            Xy = None
        elif check_gram:
            # If we're going to use the user's precomputed gram matrix, we
            # do a quick check to make sure its not totally bogus.
            _check_precomputed_gram_matrix(X, precompute, X_offset, X_scale)

    # precompute if n_samples > n_features
    if isinstance(precompute, str) and precompute == "auto":
        precompute = n_samples > n_features

    if precompute is True:
        # make sure that the 'precompute' array is contiguous.
        precompute = np.empty(shape=(n_features, n_features), dtype=X.dtype, order="C")
        np.dot(X.T, X, out=precompute)

    if not hasattr(precompute, "__array__"):
        Xy = None  # cannot use Xy if precompute is not Gram

    if hasattr(precompute, "__array__") and Xy is None:
        common_dtype = np.result_type(X.dtype, y.dtype)
        if y.ndim == 1:
            # Xy is 1d, make sure it is contiguous.
            Xy = np.empty(shape=n_features, dtype=common_dtype, order="C")
            np.dot(X.T, y, out=Xy)
        else:
            # Make sure that Xy is always F contiguous even if X or y are not
            # contiguous: the goal is to make it fast to extract the data for a
            # specific target.
            n_targets = y.shape[1]
            Xy = np.empty(shape=(n_features, n_targets), dtype=common_dtype, order="F")
            np.dot(y.T, X, out=Xy.T)

    return X, y, X_offset, y_offset, X_scale, precompute, Xy
