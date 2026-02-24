"""
Ridge regression
"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import numbers
import warnings
from abc import ABCMeta, abstractmethod
from functools import partial
from numbers import Integral, Real

import numpy as np
from scipy import linalg, optimize, sparse
from scipy.sparse import linalg as sp_linalg

from sklearn.base import (
    BaseEstimator,
    MultiOutputMixin,
    RegressorMixin,
    _fit_context,
    is_classifier,
)
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model._base import (
    LinearClassifierMixin,
    LinearModel,
    _preprocess_data,
    _rescale_data,
)
from sklearn.linear_model._sag import sag_solver
from sklearn.metrics import check_scoring, get_scorer, get_scorer_names
from sklearn.model_selection import GridSearchCV
from sklearn.preprocessing import LabelBinarizer
from sklearn.utils import (
    Bunch,
    check_array,
    check_consistent_length,
    check_scalar,
    column_or_1d,
    compute_sample_weight,
)
from sklearn.utils._array_api import (
    _is_numpy_namespace,
    _max_precision_float_dtype,
    _ravel,
    device,
    get_namespace,
    get_namespace_and_device,
    move_to,
)
from sklearn.utils._param_validation import Interval, StrOptions, validate_params
from sklearn.utils.extmath import row_norms, safe_sparse_dot
from sklearn.utils.fixes import _sparse_linalg_cg
from sklearn.utils.metadata_routing import (
    MetadataRouter,
    MethodMapping,
    _raise_for_params,
    _routing_enabled,
    process_routing,
)
from sklearn.utils.validation import (
    _check_sample_weight,
    check_is_fitted,
    validate_data,
)


def _get_rescaled_operator(X, X_offset, sample_weight_sqrt):
    """Create LinearOperator for matrix products with implicit centering.

    Matrix product `LinearOperator @ coef` returns `(X - X_offset) @ coef`.
    """

    def matvec(b):
        return X.dot(b) - sample_weight_sqrt * b.dot(X_offset)

    def rmatvec(b):
        return X.T.dot(b) - X_offset * b.dot(sample_weight_sqrt)

    X1 = sparse.linalg.LinearOperator(shape=X.shape, matvec=matvec, rmatvec=rmatvec)
    return X1


def _solve_sparse_cg(
    X,
    y,
    alpha,
    max_iter=None,
    tol=1e-4,
    verbose=0,
    X_offset=None,
    X_scale=None,
    sample_weight_sqrt=None,
):
    if sample_weight_sqrt is None:
        sample_weight_sqrt = np.ones(X.shape[0], dtype=X.dtype)

    n_samples, n_features = X.shape

    if X_offset is None or X_scale is None:
        X1 = sp_linalg.aslinearoperator(X)
    else:
        X_offset_scale = X_offset / X_scale
        X1 = _get_rescaled_operator(X, X_offset_scale, sample_weight_sqrt)

    coefs = np.empty((y.shape[1], n_features), dtype=X.dtype)

    if n_features > n_samples:

        def create_mv(curr_alpha):
            def _mv(x):
                return X1.matvec(X1.rmatvec(x)) + curr_alpha * x

            return _mv

    else:

        def create_mv(curr_alpha):
            def _mv(x):
                return X1.rmatvec(X1.matvec(x)) + curr_alpha * x

            return _mv

    for i in range(y.shape[1]):
        y_column = y[:, i]

        mv = create_mv(alpha[i])
        if n_features > n_samples:
            # kernel ridge
            # w = X.T * inv(X X^t + alpha*Id) y
            C = sp_linalg.LinearOperator(
                (n_samples, n_samples), matvec=mv, dtype=X.dtype
            )
            coef, info = _sparse_linalg_cg(C, y_column, rtol=tol)
            coefs[i] = X1.rmatvec(coef)
        else:
            # linear ridge
            # w = inv(X^t X + alpha*Id) * X.T y
            y_column = X1.rmatvec(y_column)
            C = sp_linalg.LinearOperator(
                (n_features, n_features), matvec=mv, dtype=X.dtype
            )
            coefs[i], info = _sparse_linalg_cg(C, y_column, maxiter=max_iter, rtol=tol)

        if info < 0:
            raise ValueError("Failed with error code %d" % info)

        if max_iter is None and info > 0 and verbose:
            warnings.warn(
                "sparse_cg did not converge after %d iterations." % info,
                ConvergenceWarning,
            )

    return coefs


def _solve_lsqr(
    X,
    y,
    *,
    alpha,
    fit_intercept=True,
    max_iter=None,
    tol=1e-4,
    X_offset=None,
    X_scale=None,
    sample_weight_sqrt=None,
):
    """Solve Ridge regression via LSQR.

    We expect that y is always mean centered.
    If X is dense, we expect it to be mean centered such that we can solve
        ||y - Xw||_2^2 + alpha * ||w||_2^2

    If X is sparse, we expect X_offset to be given such that we can solve
        ||y - (X - X_offset)w||_2^2 + alpha * ||w||_2^2

    With sample weights S=diag(sample_weight), this becomes
        ||sqrt(S) (y - (X - X_offset) w)||_2^2 + alpha * ||w||_2^2
    and we expect y and X to already be rescaled, i.e. sqrt(S) @ y, sqrt(S) @ X. In
    this case, X_offset is the sample_weight weighted mean of X before scaling by
    sqrt(S). The objective then reads
       ||y - (X - sqrt(S) X_offset) w)||_2^2 + alpha * ||w||_2^2
    """
    if sample_weight_sqrt is None:
        sample_weight_sqrt = np.ones(X.shape[0], dtype=X.dtype)

    if sparse.issparse(X) and fit_intercept:
        X_offset_scale = X_offset / X_scale
        X1 = _get_rescaled_operator(X, X_offset_scale, sample_weight_sqrt)
    else:
        # No need to touch anything
        X1 = X

    n_samples, n_features = X.shape
    coefs = np.empty((y.shape[1], n_features), dtype=X.dtype)
    n_iter = np.empty(y.shape[1], dtype=np.int32)

    # According to the lsqr documentation, alpha = damp^2.
    sqrt_alpha = np.sqrt(alpha)

    for i in range(y.shape[1]):
        y_column = y[:, i]
        info = sp_linalg.lsqr(
            X1, y_column, damp=sqrt_alpha[i], atol=tol, btol=tol, iter_lim=max_iter
        )
        coefs[i] = info[0]
        n_iter[i] = info[2]

    return coefs, n_iter


def _solve_cholesky(X, y, alpha):
    # w = inv(X^t X + alpha*Id) * X.T y
    n_features = X.shape[1]
    n_targets = y.shape[1]

    A = safe_sparse_dot(X.T, X, dense_output=True)
    Xy = safe_sparse_dot(X.T, y, dense_output=True)

    one_alpha = np.array_equal(alpha, len(alpha) * [alpha[0]])

    if one_alpha:
        A.flat[:: n_features + 1] += alpha[0]
        return linalg.solve(A, Xy, assume_a="pos", overwrite_a=True).T
    else:
        coefs = np.empty([n_targets, n_features], dtype=X.dtype)
        for coef, target, current_alpha in zip(coefs, Xy.T, alpha):
            A.flat[:: n_features + 1] += current_alpha
            coef[:] = linalg.solve(A, target, assume_a="pos", overwrite_a=False).ravel()
            A.flat[:: n_features + 1] -= current_alpha
        return coefs


def _solve_cholesky_kernel(K, y, alpha, sample_weight=None, copy=False):
    # dual_coef = inv(X X^t + alpha*Id) y
    n_samples = K.shape[0]
    n_targets = y.shape[1]

    if copy:
        K = K.copy()

    alpha = np.atleast_1d(alpha)
    one_alpha = (alpha == alpha[0]).all()
    has_sw = isinstance(sample_weight, np.ndarray) or sample_weight not in [1.0, None]

    if has_sw:
        # Unlike other solvers, we need to support sample_weight directly
        # because K might be a pre-computed kernel.
        sw = np.sqrt(np.atleast_1d(sample_weight))
        y = y * sw[:, np.newaxis]
        K *= np.outer(sw, sw)

    if one_alpha:
        # Only one penalty, we can solve multi-target problems in one time.
        K.flat[:: n_samples + 1] += alpha[0]

        try:
            # Note: we must use overwrite_a=False in order to be able to
            #       use the fall-back solution below in case a LinAlgError
            #       is raised
            dual_coef = linalg.solve(K, y, assume_a="pos", overwrite_a=False)
        except np.linalg.LinAlgError:
            warnings.warn(
                "Singular matrix in solving dual problem. Using "
                "least-squares solution instead."
            )
            dual_coef = linalg.lstsq(K, y)[0]

        # K is expensive to compute and store in memory so change it back in
        # case it was user-given.
        K.flat[:: n_samples + 1] -= alpha[0]

        if has_sw:
            dual_coef *= sw[:, np.newaxis]

        return dual_coef
    else:
        # One penalty per target. We need to solve each target separately.
        dual_coefs = np.empty([n_targets, n_samples], K.dtype)

        for dual_coef, target, current_alpha in zip(dual_coefs, y.T, alpha):
            K.flat[:: n_samples + 1] += current_alpha

            dual_coef[:] = linalg.solve(
                K, target, assume_a="pos", overwrite_a=False
            ).ravel()

            K.flat[:: n_samples + 1] -= current_alpha

        if has_sw:
            dual_coefs *= sw[np.newaxis, :]

        return dual_coefs.T


def _solve_svd(X, y, alpha, xp=None):
    xp, _ = get_namespace(X, xp=xp)
    U, s, Vt = xp.linalg.svd(X, full_matrices=False)
    idx = s > 1e-15  # same default value as scipy.linalg.pinv
    s_nnz = s[idx][:, None]
    UTy = U.T @ y
    d = xp.zeros((s.shape[0], alpha.shape[0]), dtype=X.dtype, device=device(X))
    d[idx] = s_nnz / (s_nnz**2 + alpha)
    d_UT_y = d * UTy
    return (Vt.T @ d_UT_y).T


def _solve_lbfgs(
    X,
    y,
    alpha,
    positive=True,
    max_iter=None,
    tol=1e-4,
    X_offset=None,
    X_scale=None,
    sample_weight_sqrt=None,
):
    """Solve ridge regression with LBFGS.

    The main purpose is fitting with forcing coefficients to be positive.
    For unconstrained ridge regression, there are faster dedicated solver methods.
    Note that with positive bounds on the coefficients, LBFGS seems faster
    than scipy.optimize.lsq_linear.
    """
    n_samples, n_features = X.shape

    options = {}
    if max_iter is not None:
        options["maxiter"] = max_iter
    config = {
        "method": "L-BFGS-B",
        "tol": tol,
        "jac": True,
        "options": options,
    }

    if isinstance(positive, bool):
        if positive:
            config["bounds"] = [(0, np.inf)] * n_features
    elif isinstance(positive, list):
        if len(positive) != n_features:
            raise ValueError(
                "Length of 'positive' list must be equal to the number of features"
            )
        elif not all(isinstance(p, bool) for p in positive):
            raise ValueError(
                "'positive' must be either a boolean or a list of booleans"
            )
        config["bounds"] = [(0, np.inf) if p else (None, None) for p in positive]
    else:
        raise ValueError("'positive' must be either a boolean or a list of booleans")

    if X_offset is not None and X_scale is not None:
        X_offset_scale = X_offset / X_scale
    else:
        X_offset_scale = None

    if sample_weight_sqrt is None:
        sample_weight_sqrt = np.ones(X.shape[0], dtype=X.dtype)

    coefs = np.empty((y.shape[1], n_features), dtype=X.dtype)

    for i in range(y.shape[1]):
        x0 = np.zeros((n_features,))
        y_column = y[:, i]

        def func(w):
            residual = X.dot(w) - y_column
            if X_offset_scale is not None:
                residual -= sample_weight_sqrt * w.dot(X_offset_scale)
            f = 0.5 * residual.dot(residual) + 0.5 * alpha[i] * w.dot(w)
            grad = X.T @ residual + alpha[i] * w
            if X_offset_scale is not None:
                grad -= X_offset_scale * residual.dot(sample_weight_sqrt)

            return f, grad

        result = optimize.minimize(func, x0, **config)
        if not result["success"]:
            warnings.warn(
                (
                    "The lbfgs solver did not converge. Try increasing max_iter "
                    f"or tol. Currently: max_iter={max_iter} and tol={tol}"
                ),
                ConvergenceWarning,
            )
        coefs[i] = result["x"]

    return coefs


def _get_valid_accept_sparse(is_X_sparse, solver):
    if is_X_sparse and solver in ["auto", "sag", "saga"]:
        return "csr"
    else:
        return ["csr", "csc", "coo"]


@validate_params(
    {
        "X": ["array-like", "sparse matrix", sp_linalg.LinearOperator],
        "y": ["array-like"],
        "alpha": [Interval(Real, 0, None, closed="left"), "array-like"],
        "sample_weight": [
            Interval(Real, None, None, closed="neither"),
            "array-like",
            None,
        ],
        "solver": [
            StrOptions(
                {"auto", "svd", "cholesky", "lsqr", "sparse_cg", "sag", "saga", "lbfgs"}
            )
        ],
        "max_iter": [Interval(Integral, 0, None, closed="left"), None],
        "tol": [Interval(Real, 0, None, closed="left")],
        "verbose": ["verbose"],
        "positive": ["boolean", "array-like"],
        "random_state": ["random_state"],
        "return_n_iter": ["boolean"],
        "return_intercept": ["boolean"],
        "check_input": ["boolean"],
    },
    prefer_skip_nested_validation=True,
)
def ridge_regression(
    X,
    y,
    alpha,
    *,
    sample_weight=None,
    solver="auto",
    max_iter=None,
    tol=1e-4,
    verbose=0,
    positive=False,
    random_state=None,
    return_n_iter=False,
    return_intercept=False,
    check_input=True,
):
    """Solve the ridge equation by the method of normal equations.

    Read more in the :ref:`User Guide <ridge_regression>`.

    Parameters
    ----------
    X : {array-like, sparse matrix, LinearOperator} of shape \
        (n_samples, n_features)
        Training data.

    y : array-like of shape (n_samples,) or (n_samples, n_targets)
        Target values.

    alpha : float or array-like of shape (n_targets,)
        Constant that multiplies the L2 term, controlling regularization
        strength. `alpha` must be a non-negative float i.e. in `[0, inf)`.

        When `alpha = 0`, the objective is equivalent to ordinary least
        squares, solved by the :class:`LinearRegression` object. For numerical
        reasons, using `alpha = 0` with the `Ridge` object is not advised.
        Instead, you should use the :class:`LinearRegression` object.

        If an array is passed, penalties are assumed to be specific to the
        targets. Hence they must correspond in number.

        For an illustration of the effect of alpha on the model coefficients, see
        :ref:`sphx_glr_auto_examples_linear_model_plot_ridge_coeffs.py`.

    sample_weight : float or array-like of shape (n_samples,), default=None
        Individual weights for each sample. If given a float, every sample
        will have the same weight. If sample_weight is not None and
        solver='auto', the solver will be set to 'cholesky'.

        .. versionadded:: 0.17

    solver : {'auto', 'svd', 'cholesky', 'lsqr', 'sparse_cg', \
            'sag', 'saga', 'lbfgs'}, default='auto'
        Solver to use in the computational routines:

        - 'auto' chooses the solver automatically based on the type of data.

        - 'svd' uses a Singular Value Decomposition of X to compute the Ridge
          coefficients. It is the most stable solver, in particular more stable
          for singular matrices than 'cholesky' at the cost of being slower.

        - 'cholesky' uses the standard scipy.linalg.solve function to
          obtain a closed-form solution via a Cholesky decomposition of
          dot(X.T, X)

        - 'sparse_cg' uses the conjugate gradient solver as found in
          scipy.sparse.linalg.cg. As an iterative algorithm, this solver is
          more appropriate than 'cholesky' for large-scale data
          (possibility to set `tol` and `max_iter`).

        - 'lsqr' uses the dedicated regularized least-squares routine
          scipy.sparse.linalg.lsqr. It is the fastest and uses an iterative
          procedure.

        - 'sag' uses a Stochastic Average Gradient descent, and 'saga' uses
          its improved, unbiased version named SAGA. Both methods also use an
          iterative procedure, and are often faster than other solvers when
          both n_samples and n_features are large. Note that 'sag' and
          'saga' fast convergence is only guaranteed on features with
          approximately the same scale. You can preprocess the data with a
          scaler from sklearn.preprocessing.

        - 'lbfgs' uses L-BFGS-B algorithm implemented in
          `scipy.optimize.minimize`. It can be used only when `positive`
          is True.

        All solvers except 'svd' support both dense and sparse data. However, only
        'lsqr', 'sag', 'sparse_cg', and 'lbfgs' support sparse input when
        `fit_intercept` is True.

        .. versionadded:: 0.17
           Stochastic Average Gradient descent solver.
        .. versionadded:: 0.19
           SAGA solver.

    max_iter : int, default=None
        Maximum number of iterations for conjugate gradient solver.
        For the 'sparse_cg' and 'lsqr' solvers, the default value is determined
        by scipy.sparse.linalg. For 'sag' and saga solver, the default value is
        1000. For 'lbfgs' solver, the default value is 15000.

    tol : float, default=1e-4
        Precision of the solution. Note that `tol` has no effect for solvers 'svd' and
        'cholesky'.

        .. versionchanged:: 1.2
           Default value changed from 1e-3 to 1e-4 for consistency with other linear
           models.

    verbose : int, default=0
        Verbosity level. Setting verbose > 0 will display additional
        information depending on the solver used.

    positive : bool or list of bool, default=False
        When set to ``True``, forces the coefficients to be positive.
        Positivity constraint on the coefficients. 
        If True, forces all coefficients to be positive.
        If False, no constraint is enforced.
        If list of length `n_features`, the positivity constraint is specified 
        for each feature, forcing the corresponding coefficient to be positive.
        Only 'lbfgs' solver is supported in this case.

    random_state : int, RandomState instance, default=None
        Used when ``solver`` == 'sag' or 'saga' to shuffle the data.
        See :term:`Glossary <random_state>` for details.

    return_n_iter : bool, default=False
        If True, the method also returns `n_iter`, the actual number of
        iteration performed by the solver.

        .. versionadded:: 0.17

    return_intercept : bool, default=False
        If True and if X is sparse, the method also returns the intercept,
        and the solver is automatically changed to 'sag'. This is only a
        temporary fix for fitting the intercept with sparse data. For dense
        data, use sklearn.linear_model._preprocess_data before your regression.

        .. versionadded:: 0.17

    check_input : bool, default=True
        If False, the input arrays X and y will not be checked.

        .. versionadded:: 0.21

    Returns
    -------
    coef : ndarray of shape (n_features,) or (n_targets, n_features)
        Weight vector(s).

    n_iter : int, optional
        The actual number of iteration performed by the solver.
        Only returned if `return_n_iter` is True.

    intercept : float or ndarray of shape (n_targets,)
        The intercept of the model. Only returned if `return_intercept`
        is True and if X is a scipy sparse array.

    Notes
    -----
    This function won't compute the intercept.

    Regularization improves the conditioning of the problem and
    reduces the variance of the estimates. Larger values specify stronger
    regularization. Alpha corresponds to ``1 / (2C)`` in other linear
    models such as :class:`~sklearn.linear_model.LogisticRegression` or
    :class:`~sklearn.svm.LinearSVC`. If an array is passed, penalties are
    assumed to be specific to the targets. Hence they must correspond in
    number.

    Examples
    --------
    >>> import numpy as np
    >>> from sklearn.datasets import make_regression
    >>> from sklearn.linear_model import ridge_regression
    >>> rng = np.random.RandomState(0)
    >>> X = rng.randn(100, 4)
    >>> y = 2.0 * X[:, 0] - 1.0 * X[:, 1] + 0.1 * rng.standard_normal(100)
    >>> coef, intercept = ridge_regression(X, y, alpha=1.0, return_intercept=True,
    ...                                    random_state=0)
    >>> coef
    array([ 1.97, -1., -2.69e-3, -9.27e-4 ])
    >>> intercept
    np.float64(-.0012)
    """
    return _ridge_regression(
        X,
        y,
        alpha,
        sample_weight=sample_weight,
        solver=solver,
        max_iter=max_iter,
        tol=tol,
        verbose=verbose,
        positive=positive,
        random_state=random_state,
        return_n_iter=return_n_iter,
        return_intercept=return_intercept,
        X_scale=None,
        X_offset=None,
        check_input=check_input,
    )


def _ridge_regression(
    X,
    y,
    alpha,
    sample_weight=None,
    solver="auto",
    max_iter=None,
    tol=1e-4,
    verbose=0,
    positive=False,
    random_state=None,
    return_n_iter=False,
    return_intercept=False,
    return_solver=False,
    X_scale=None,
    X_offset=None,
    check_input=True,
    fit_intercept=False,
):
    xp, is_array_api_compliant, device_ = get_namespace_and_device(
        X, y, sample_weight, X_scale, X_offset
    )
    is_numpy_namespace = _is_numpy_namespace(xp)
    X_is_sparse = sparse.issparse(X)

    has_sw = sample_weight is not None

    if isinstance(positive, bool):
        _positive = positive
    elif isinstance(positive, list):
        if len(positive) != X.shape[1]:
            raise ValueError(
                "Length of 'positive' list must be equal to the number of features"
            )
        if not all(isinstance(p, bool) for p in positive):
            raise ValueError(
                "'positive' must be either a boolean or a list of booleans"
            )
        _positive = True

    solver = resolve_solver(solver, _positive, return_intercept, X_is_sparse, xp)

    if is_numpy_namespace and not X_is_sparse:
        X = np.asarray(X)

    if not is_numpy_namespace and solver != "svd":
        raise ValueError(
            f"Array API dispatch to namespace {xp.__name__} only supports "
            f"solver 'svd'. Got '{solver}'."
        )

    if _positive and solver != "lbfgs":
        raise ValueError(
            "When positive=True, only 'lbfgs' solver can be used. "
            f"Please change solver {solver} to 'lbfgs' "
            "or set positive=False."
        )

    if solver == "lbfgs" and not _positive:
        raise ValueError(
            "'lbfgs' solver can be used only when positive=True. "
            "Please use another solver."
        )

    if return_intercept and solver != "sag":
        raise ValueError(
            "In Ridge, only 'sag' solver can directly fit the "
            "intercept. Please change solver to 'sag' or set "
            "return_intercept=False."
        )

    if check_input:
        _dtype = [xp.float64, xp.float32]
        _accept_sparse = _get_valid_accept_sparse(X_is_sparse, solver)
        X = check_array(X, accept_sparse=_accept_sparse, dtype=_dtype, order="C")
        y = check_array(y, dtype=X.dtype, ensure_2d=False, order=None)
    check_consistent_length(X, y)

    n_samples, n_features = X.shape

    if y.ndim > 2:
        raise ValueError("Target y has the wrong shape %s" % str(y.shape))

    if y.ndim == 1:
        y = xp.reshape(y, (-1, 1))

    n_samples_, n_targets = y.shape

    if n_samples != n_samples_:
        raise ValueError(
            "Number of samples in X and y does not correspond: %d != %d"
            % (n_samples, n_samples_)
        )

    if has_sw:
        sample_weight = _check_sample_weight(sample_weight, X, dtype=X.dtype)

        if solver not in ["sag", "saga"]:
            # SAG supports sample_weight directly. For other solvers,
            # we implement sample_weight via a simple rescaling.
            X, y, sample_weight_sqrt = _rescale_data(X, y, sample_weight)

    # Some callers of this method might pass alpha as single
    # element array which already has been validated.
    if alpha is not None and not isinstance(alpha, type(xp.asarray([0.0]))):
        alpha = check_scalar(
            alpha,
            "alpha",
            target_type=numbers.Real,
            min_val=0.0,
            include_boundaries="left",
        )

    # There should be either 1 or n_targets penalties
    alpha = _ravel(xp.asarray(alpha, device=device_, dtype=X.dtype), xp=xp)
    if alpha.shape[0] not in [1, n_targets]:
        raise ValueError(
            "Number of targets and number of penalties do not correspond: %d != %d"
            % (alpha.shape[0], n_targets)
        )

    if alpha.shape[0] == 1 and n_targets > 1:
        alpha = xp.full(
            shape=(n_targets,), fill_value=alpha[0], dtype=alpha.dtype, device=device_
        )

    n_iter = None
    if solver == "sparse_cg":
        coef = _solve_sparse_cg(
            X,
            y,
            alpha,
            max_iter=max_iter,
            tol=tol,
            verbose=verbose,
            X_offset=X_offset,
            X_scale=X_scale,
            sample_weight_sqrt=sample_weight_sqrt if has_sw else None,
        )

    elif solver == "lsqr":
        coef, n_iter = _solve_lsqr(
            X,
            y,
            alpha=alpha,
            fit_intercept=fit_intercept,
            max_iter=max_iter,
            tol=tol,
            X_offset=X_offset,
            X_scale=X_scale,
            sample_weight_sqrt=sample_weight_sqrt if has_sw else None,
        )

    elif solver == "cholesky":
        if n_features > n_samples:
            K = safe_sparse_dot(X, X.T, dense_output=True)
            try:
                dual_coef = _solve_cholesky_kernel(K, y, alpha)

                coef = safe_sparse_dot(X.T, dual_coef, dense_output=True).T
            except linalg.LinAlgError:
                # use SVD solver if matrix is singular
                solver = "svd"
        else:
            try:
                coef = _solve_cholesky(X, y, alpha)
            except linalg.LinAlgError:
                # use SVD solver if matrix is singular
                solver = "svd"

    elif solver in ["sag", "saga"]:
        # precompute max_squared_sum for all targets
        max_squared_sum = row_norms(X, squared=True).max()

        coef = np.empty((y.shape[1], n_features), dtype=X.dtype)
        n_iter = np.empty(y.shape[1], dtype=np.int32)
        intercept = np.zeros((y.shape[1],), dtype=X.dtype)
        for i, (alpha_i, target) in enumerate(zip(alpha, y.T)):
            init = {
                "coef": np.zeros((n_features + int(return_intercept), 1), dtype=X.dtype)
            }
            coef_, n_iter_, _ = sag_solver(
                X,
                target.ravel(),
                sample_weight,
                "squared",
                alpha_i,
                0,
                max_iter,
                tol,
                verbose,
                random_state,
                False,
                max_squared_sum,
                init,
                is_saga=solver == "saga",
            )
            if return_intercept:
                coef[i] = coef_[:-1]
                intercept[i] = coef_[-1]
            else:
                coef[i] = coef_
            n_iter[i] = n_iter_

        if intercept.shape[0] == 1:
            intercept = intercept[0]

    elif solver == "lbfgs":
        coef = _solve_lbfgs(
            X,
            y,
            alpha,
            positive=positive,
            tol=tol,
            max_iter=max_iter,
            X_offset=X_offset,
            X_scale=X_scale,
            sample_weight_sqrt=sample_weight_sqrt if has_sw else None,
        )

    if solver == "svd":
        if X_is_sparse:
            raise TypeError("SVD solver does not support sparse inputs currently")
        coef = _solve_svd(X, y, alpha, xp)

    if n_targets == 1:
        coef = _ravel(coef)

    coef = xp.asarray(coef)

    if return_n_iter and return_intercept:
        res = coef, n_iter, intercept
    elif return_intercept:
        res = coef, intercept
    elif return_n_iter:
        res = coef, n_iter
    else:
        res = coef

    return (*res, solver) if return_solver else res


def resolve_solver(solver, positive, return_intercept, is_sparse, xp):
    if solver != "auto":
        return solver

    is_numpy_namespace = _is_numpy_namespace(xp)

    auto_solver_np = resolve_solver_for_numpy(positive, return_intercept, is_sparse)
    if is_numpy_namespace:
        return auto_solver_np

    if positive:
        raise ValueError(
            "The solvers that support positive fitting do not support "
            f"Array API dispatch to namespace {xp.__name__}. Please "
            "either disable Array API dispatch, or use a numpy-like "
            "namespace, or set `positive=False`."
        )

    # At the moment, Array API dispatch only supports the "svd" solver.
    solver = "svd"
    if solver != auto_solver_np:
        warnings.warn(
            f"Using Array API dispatch to namespace {xp.__name__} with "
            f"`solver='auto'` will result in using the solver '{solver}'. "
            "The results may differ from those when using a Numpy array, "
            f"because in that case the preferred solver would be {auto_solver_np}. "
            f"Set `solver='{solver}'` to suppress this warning."
        )

    return solver


def resolve_solver_for_numpy(positive, return_intercept, is_sparse):
    if positive:
        return "lbfgs"

    if return_intercept:
        # sag supports fitting intercept directly
        return "sag"

    if not is_sparse:
        return "cholesky"

    return "sparse_cg"


class _BaseRidge(LinearModel, metaclass=ABCMeta):
    _parameter_constraints: dict = {
        "alpha": [Interval(Real, 0, None, closed="left"), np.ndarray],
        "fit_intercept": ["boolean"],
        "copy_X": ["boolean"],
        "max_iter": [Interval(Integral, 1, None, closed="left"), None],
        "tol": [Interval(Real, 0, None, closed="left")],
        "solver": [
            StrOptions(
                {"auto", "svd", "cholesky", "lsqr", "sparse_cg", "sag", "saga", "lbfgs"}
            )
        ],
        "positive": ["boolean", "array-like"],
        "random_state": ["random_state"],
    }

    @abstractmethod
    def __init__(
        self,
        alpha=1.0,
        *,
        fit_intercept=True,
        copy_X=True,
        max_iter=None,
        tol=1e-4,
        solver="auto",
        positive=False,
        random_state=None,
    ):
        self.alpha = alpha
        self.fit_intercept = fit_intercept
        self.copy_X = copy_X
        self.max_iter = max_iter
        self.tol = tol
        self.solver = solver
        self.positive = positive
        self.random_state = random_state

    def fit(self, X, y, sample_weight=None):
        xp, is_array_api_compliant = get_namespace(X, y, sample_weight)

        if self.solver == "lbfgs" and not self.positive:
            raise ValueError(
                "'lbfgs' solver can be used only when positive=True. "
                "Please use another solver."
            )

        if self.positive:
            if self.solver not in ["auto", "lbfgs"]:
                raise ValueError(
                    f"solver='{self.solver}' does not support positive fitting. Please"
                    " set the solver to 'auto' or 'lbfgs', or set `positive=False`"
                )
            else:
                solver = self.solver
        elif sparse.issparse(X) and self.fit_intercept:
            if self.solver not in ["auto", "lbfgs", "lsqr", "sag", "sparse_cg"]:
                raise ValueError(
                    "solver='{}' does not support fitting the intercept "
                    "on sparse data. Please set the solver to 'auto' or "
                    "'lsqr', 'sparse_cg', 'sag', 'lbfgs' "
                    "or set `fit_intercept=False`".format(self.solver)
                )
            if self.solver in ["lsqr", "lbfgs"]:
                solver = self.solver
            elif self.solver == "sag" and self.max_iter is None and self.tol > 1e-4:
                warnings.warn(
                    '"sag" solver requires many iterations to fit '
                    "an intercept with sparse inputs. Either set the "
                    'solver to "auto" or "sparse_cg", or set a low '
                    '"tol" and a high "max_iter" (especially if inputs are '
                    "not standardized)."
                )
                solver = "sag"
            else:
                solver = "sparse_cg"
        else:
            solver = self.solver

        if sample_weight is not None:
            sample_weight = _check_sample_weight(sample_weight, X, dtype=X.dtype)

        # when X is sparse we only remove offset from y
        X, y, X_offset, y_offset, X_scale, _ = _preprocess_data(
            X,
            y,
            fit_intercept=self.fit_intercept,
            copy=self.copy_X,
            sample_weight=sample_weight,
            rescale_with_sw=False,
        )

        if solver == "sag" and sparse.issparse(X) and self.fit_intercept:
            self.coef_, self.n_iter_, self.intercept_, self.solver_ = _ridge_regression(
                X,
                y,
                alpha=self.alpha,
                sample_weight=sample_weight,
                max_iter=self.max_iter,
                tol=self.tol,
                solver="sag",
                positive=self.positive,
                random_state=self.random_state,
                return_n_iter=True,
                return_intercept=True,
                return_solver=True,
                check_input=False,
            )
            # add the offset which was subtracted by _preprocess_data
            self.intercept_ += y_offset

        else:
            if sparse.issparse(X) and self.fit_intercept:
                # required to fit intercept with sparse_cg and lbfgs solver
                params = {"X_offset": X_offset, "X_scale": X_scale}
            else:
                # for dense matrices or when intercept is set to 0
                params = {}

            self.coef_, self.n_iter_, self.solver_ = _ridge_regression(
                X,
                y,
                alpha=self.alpha,
                sample_weight=sample_weight,
                max_iter=self.max_iter,
                tol=self.tol,
                solver=solver,
                positive=self.positive,
                random_state=self.random_state,
                return_n_iter=True,
                return_intercept=False,
                return_solver=True,
                check_input=False,
                fit_intercept=self.fit_intercept,
                **params,
            )
            self._set_intercept(X_offset, y_offset, X_scale)

        return self


class Ridge(MultiOutputMixin, RegressorMixin, _BaseRidge):
    """Linear least squares with l2 regularization.

    Minimizes the objective function::

    ||y - Xw||^2_2 + alpha * ||w||^2_2

    This model solves a regression model where the loss function is
    the linear least squares function and regularization is given by
    the l2-norm. Also known as Ridge Regression or Tikhonov regularization.
    This estimator has built-in support for multi-variate regression
    (i.e., when y is a 2d-array of shape (n_samples, n_targets)).

    Read more in the :ref:`User Guide <ridge_regression>`.

    Parameters
    ----------
    alpha : {float, ndarray of shape (n_targets,)}, default=1.0
        Constant that multiplies the L2 term, controlling regularization
        strength. `alpha` must be a non-negative float i.e. in `[0, inf)`.

        When `alpha = 0`, the objective is equivalent to ordinary least
        squares, solved by the :class:`LinearRegression` object. For numerical
        reasons, using `alpha = 0` with the `Ridge` object is not advised.
        Instead, you should use the :class:`LinearRegression` object.

        If an array is passed, penalties are assumed to be specific to the
        targets. Hence they must correspond in number.

        See :ref:`sphx_glr_auto_examples_linear_model_plot_ridge_coeffs.py`
        for an illustration of the effect of alpha on the model coefficients.

    fit_intercept : bool, default=True
        Whether to fit the intercept for this model. If set
        to false, no intercept will be used in calculations
        (i.e. ``X`` and ``y`` are expected to be centered).

    copy_X : bool, default=True
        If True, X will be copied; else, it may be overwritten.

    max_iter : int, default=None
        Maximum number of iterations for conjugate gradient solver.
        For 'sparse_cg' and 'lsqr' solvers, the default value is determined
        by scipy.sparse.linalg. For 'sag' solver, the default value is 1000.
        For 'lbfgs' solver, the default value is 15000.

    tol : float, default=1e-4
        The precision of the solution (`coef_`) is determined by `tol` which
        specifies a different convergence criterion for each solver:

        - 'svd': `tol` has no impact.

        - 'cholesky': `tol` has no impact.

        - 'sparse_cg': norm of residuals smaller than `tol`.

        - 'lsqr': `tol` is set as atol and btol of scipy.sparse.linalg.lsqr,
          which control the norm of the residual vector in terms of the norms of
          matrix and coefficients.

        - 'sag' and 'saga': relative change of coef smaller than `tol`.

        - 'lbfgs': maximum of the absolute (projected) gradient=max|residuals|
          smaller than `tol`.

        .. versionchanged:: 1.2
           Default value changed from 1e-3 to 1e-4 for consistency with other linear
           models.

    solver : {'auto', 'svd', 'cholesky', 'lsqr', 'sparse_cg', \
            'sag', 'saga', 'lbfgs'}, default='auto'
        Solver to use in the computational routines:

        - 'auto' chooses the solver automatically based on the type of data.

        - 'svd' uses a Singular Value Decomposition of X to compute the Ridge
          coefficients. It is the most stable solver, in particular more stable
          for singular matrices than 'cholesky' at the cost of being slower.

        - 'cholesky' uses the standard :func:`scipy.linalg.solve` function to
          obtain a closed-form solution.

        - 'sparse_cg' uses the conjugate gradient solver as found in
          :func:`scipy.sparse.linalg.cg`. As an iterative algorithm, this solver is
          more appropriate than 'cholesky' for large-scale data
          (possibility to set `tol` and `max_iter`).

        - 'lsqr' uses the dedicated regularized least-squares routine
          :func:`scipy.sparse.linalg.lsqr`. It is the fastest and uses an iterative
          procedure.

        - 'sag' uses a Stochastic Average Gradient descent, and 'saga' uses
          its improved, unbiased version named SAGA. Both methods also use an
          iterative procedure, and are often faster than other solvers when
          both n_samples and n_features are large. Note that 'sag' and
          'saga' fast convergence is only guaranteed on features with
          approximately the same scale. You can preprocess the data with a
          scaler from :mod:`sklearn.preprocessing`.

        - 'lbfgs' uses L-BFGS-B algorithm implemented in
          :func:`scipy.optimize.minimize`. It can be used only when `positive`
          is True.

        All solvers except 'svd' support both dense and sparse data. However, only
        'lsqr', 'sag', 'sparse_cg', and 'lbfgs' support sparse input when
        `fit_intercept` is True.

        .. versionadded:: 0.17
           Stochastic Average Gradient descent solver.
        .. versionadded:: 0.19
           SAGA solver.

    positive : bool or list of bool, default=False
        When set to ``True``, forces the coefficients to be positive.
        When set to ``[True, True, False]``, forces the first two
        coefficients to be positive and the third one to be negative.
        When set to ``False``, no specific constraint is enforced.
        When positive is list, the list must have the same length as the number
        of expected coefficients. (X.shape[1])
        Only 'lbfgs' solver is supported in this case.

    random_state : int, RandomState instance, default=None
        Used when ``solver`` == 'sag' or 'saga' to shuffle the data.
        See :term:`Glossary <random_state>` for details.

        .. versionadded:: 0.17
           `random_state` to support Stochastic Average Gradient.

    Attributes
    ----------
    coef_ : ndarray of shape (n_features,) or (n_targets, n_features)
        Weight vector(s).

    intercept_ : float or ndarray of shape (n_targets,)
        Independent term in decision function. Set to 0.0 if
        ``fit_intercept = False``.

    n_iter_ : None or ndarray of shape (n_targets,)
        Actual number of iterations for each target. Available only for
        'sag' and 'lsqr' solvers. Other solvers will return None.

        .. versionadded:: 0.17

    n_features_in_ : int
        Number of features seen during :term:`fit`.

        .. versionadded:: 0.24

    feature_names_in_ : ndarray of shape (`n_features_in_`,)
        Names of features seen during :term:`fit`. Defined only when `X`
        has feature names that are all strings.

        .. versionadded:: 1.0

    solver_ : str
        The solver that was used at fit time by the computational
        routines.

        .. versionadded:: 1.5

    See Also
    --------
    RidgeClassifier : Ridge classifier.
    RidgeCV : Ridge regression with built-in cross validation.
    :class:`~sklearn.kernel_ridge.KernelRidge` : Kernel ridge regression
        combines ridge regression with the kernel trick.

    Notes
    -----
    Regularization improves the conditioning of the problem and
    reduces the variance of the estimates. Larger values specify stronger
    regularization. Alpha corresponds to ``1 / (2C)`` in other linear
    models such as :class:`~sklearn.linear_model.LogisticRegression` or
    :class:`~sklearn.svm.LinearSVC`.

    Examples
    --------
    >>> from sklearn.linear_model import Ridge
    >>> import numpy as np
    >>> n_samples, n_features = 10, 5
    >>> rng = np.random.RandomState(0)
    >>> y = rng.randn(n_samples)
    >>> X = rng.randn(n_samples, n_features)
    >>> clf = Ridge(alpha=1.0)
    >>> clf.fit(X, y)
    Ridge()
    """

    def __init__(
        self,
        alpha=1.0,
        *,
        fit_intercept=True,
        copy_X=True,
        max_iter=None,
        tol=1e-4,
        solver="auto",
        positive=False,
        random_state=None,
    ):
        super().__init__(
            alpha=alpha,
            fit_intercept=fit_intercept,
            copy_X=copy_X,
            max_iter=max_iter,
            tol=tol,
            solver=solver,
            positive=positive,
            random_state=random_state,
        )

    @_fit_context(prefer_skip_nested_validation=True)
    def fit(self, X, y, sample_weight=None):
        """Fit Ridge regression model.

        Parameters
        ----------
        X : {ndarray, sparse matrix} of shape (n_samples, n_features)
            Training data.

        y : ndarray of shape (n_samples,) or (n_samples, n_targets)
            Target values.

        sample_weight : float or ndarray of shape (n_samples,), default=None
            Individual weights for each sample. If given a float, every sample
            will have the same weight.

        Returns
        -------
        self : object
            Fitted estimator.
        """
        _accept_sparse = _get_valid_accept_sparse(sparse.issparse(X), self.solver)
        xp, _ = get_namespace(X, y, sample_weight)
        X, y = validate_data(
            self,
            X,
            y,
            accept_sparse=_accept_sparse,
            dtype=[xp.float64, xp.float32],
            force_writeable=True,
            multi_output=True,
            y_numeric=True,
        )
        return super().fit(X, y, sample_weight=sample_weight)

    def __sklearn_tags__(self):
        tags = super().__sklearn_tags__()
        tags.array_api_support = True
        tags.input_tags.sparse = (self.solver != "svd") and (
            self.solver != "cholesky" or not self.fit_intercept
        )
        return tags


class _RidgeClassifierMixin(LinearClassifierMixin):
    def _prepare_data(self, X, y, sample_weight, solver):
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

        solver : str
            The solver used in `Ridge` to know which sparse format to support.

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
        accept_sparse = _get_valid_accept_sparse(sparse.issparse(X), solver)
        xp, _, device_ = get_namespace_and_device(X)
        sample_weight = move_to(sample_weight, xp=xp, device=device_)
        X, y = validate_data(
            self,
            X,
            y,
            accept_sparse=accept_sparse,
            multi_output=True,
            y_numeric=False,
            force_writeable=True,
        )

        self._label_binarizer = LabelBinarizer(pos_label=1, neg_label=-1)
        xp_y, y_is_array_api = get_namespace(y)
        Y = self._label_binarizer.fit_transform(y)
        Y = move_to(Y, xp=xp, device=device_)
        if y_is_array_api and xp_y.isdtype(y.dtype, "numeric"):
            self.classes_ = move_to(
                self._label_binarizer.classes_, xp=xp, device=device_
            )
        else:
            self.classes_ = self._label_binarizer.classes_
        if not self._label_binarizer.y_type_.startswith("multilabel"):
            y = column_or_1d(y, warn=True)

        sample_weight = _check_sample_weight(sample_weight, X, dtype=X.dtype)
        if self.class_weight:
            reweighting = compute_sample_weight(self.class_weight, y)
            reweighting = move_to(reweighting, xp=xp, device=device_)
            sample_weight = sample_weight * reweighting
        return X, y, sample_weight, Y

    def predict(self, X):
        """Predict class labels for samples in `X`.

        Parameters
        ----------
        X : {array-like, spare matrix} of shape (n_samples, n_features)
            The data matrix for which we want to predict the targets.

        Returns
        -------
        y_pred : ndarray of shape (n_samples,) or (n_samples, n_outputs)
            Vector or matrix containing the predictions. In binary and
            multiclass problems, this is a vector containing `n_samples`. In
            a multilabel problem, it returns a matrix of shape
            `(n_samples, n_outputs)`.
        """
        check_is_fitted(self, attributes=["_label_binarizer"])
        if self._label_binarizer.y_type_.startswith("multilabel"):
            # Threshold such that the negative label is -1 and positive label
            # is 1 to use the inverse transform of the label binarizer fitted
            # during fit.
            decision = self.decision_function(X)
            xp, _ = get_namespace(decision)
            scores = 2.0 * xp.astype(decision > 0, decision.dtype) - 1.0
            return self._label_binarizer.inverse_transform(scores)
        return super().predict(X)

    def __sklearn_tags__(self):
        tags = super().__sklearn_tags__()
        tags.classifier_tags.multi_label = True
        return tags

    def _get_scorer_instance(self):
        """Return a scorer which corresponds to what's defined in ClassiferMixin
        parent class. This is used for routing `sample_weight`.
        """
        return get_scorer("accuracy")


class RidgeClassifier(_RidgeClassifierMixin, _BaseRidge):
    """Classifier using Ridge regression.

    This classifier first converts the target values into ``{-1, 1}`` and
    then treats the problem as a regression task (multi-output regression in
    the multiclass case).

    Read more in the :ref:`User Guide <ridge_regression>`.

    Parameters
    ----------
    alpha : float, default=1.0
        Regularization strength; must be a positive float. Regularization
        improves the conditioning of the problem and reduces the variance of
        the estimates. Larger values specify stronger regularization.
        Alpha corresponds to ``1 / (2C)`` in other linear models such as
        :class:`~sklearn.linear_model.LogisticRegression` or
        :class:`~sklearn.svm.LinearSVC`.

        For an illustration of the effect of alpha on the model coefficients, see
        :ref:`sphx_glr_auto_examples_linear_model_plot_ridge_coeffs.py`.

    fit_intercept : bool, default=True
        Whether to calculate the intercept for this model. If set to false, no
        intercept will be used in calculations (e.g. data is expected to be
        already centered).

    copy_X : bool, default=True
        If True, X will be copied; else, it may be overwritten.

    max_iter : int, default=None
        Maximum number of iterations for conjugate gradient solver.
        The default value is determined by scipy.sparse.linalg.

    tol : float, default=1e-4
        The precision of the solution (`coef_`) is determined by `tol` which
        specifies a different convergence criterion for each solver:

        - 'svd': `tol` has no impact.

        - 'cholesky': `tol` has no impact.

        - 'sparse_cg': norm of residuals smaller than `tol`.

        - 'lsqr': `tol` is set as atol and btol of scipy.sparse.linalg.lsqr,
          which control the norm of the residual vector in terms of the norms of
          matrix and coefficients.

        - 'sag' and 'saga': relative change of coef smaller than `tol`.

        - 'lbfgs': maximum of the absolute (projected) gradient=max|residuals|
          smaller than `tol`.

        .. versionchanged:: 1.2
           Default value changed from 1e-3 to 1e-4 for consistency with other linear
           models.

    class_weight : dict or 'balanced', default=None
        Weights associated with classes in the form ``{class_label: weight}``.
        If not given, all classes are supposed to have weight one.

        The "balanced" mode uses the values of y to automatically adjust
        weights inversely proportional to class frequencies in the input data
        as ``n_samples / (n_classes * np.bincount(y))``.

    solver : {'auto', 'svd', 'cholesky', 'lsqr', 'sparse_cg', \
            'sag', 'saga', 'lbfgs'}, default='auto'
        Solver to use in the computational routines:

        - 'auto' chooses the solver automatically based on the type of data.

        - 'svd' uses a Singular Value Decomposition of X to compute the Ridge
          coefficients. It is the most stable solver, in particular more stable
          for singular matrices than 'cholesky' at the cost of being slower.

        - 'cholesky' uses the standard scipy.linalg.solve function to
          obtain a closed-form solution.

        - 'sparse_cg' uses the conjugate gradient solver as found in
          scipy.sparse.linalg.cg. As an iterative algorithm, this solver is
          more appropriate than 'cholesky' for large-scale data
          (possibility to set `tol` and `max_iter`).

        - 'lsqr' uses the dedicated regularized least-squares routine
          scipy.sparse.linalg.lsqr. It is the fastest and uses an iterative
          procedure.

        - 'sag' uses a Stochastic Average Gradient descent, and 'saga' uses
          its unbiased and more flexible version named SAGA. Both methods
          use an iterative procedure, and are often faster than other solvers
          when both n_samples and n_features are large. Note that 'sag' and
          'saga' fast convergence is only guaranteed on features with
          approximately the same scale. You can preprocess the data with a
          scaler from sklearn.preprocessing.

          .. versionadded:: 0.17
             Stochastic Average Gradient descent solver.
          .. versionadded:: 0.19
             SAGA solver.

        - 'lbfgs' uses L-BFGS-B algorithm implemented in
          `scipy.optimize.minimize`. It can be used only when `positive`
          is True.

    positive : bool, default=False
        When set to ``True``, forces the coefficients to be positive.
        Only 'lbfgs' solver is supported in this case.

    random_state : int, RandomState instance, default=None
        Used when ``solver`` == 'sag' or 'saga' to shuffle the data.
        See :term:`Glossary <random_state>` for details.

    Attributes
    ----------
    coef_ : ndarray of shape (1, n_features) or (n_classes, n_features)
        Coefficient of the features in the decision function.

        ``coef_`` is of shape (1, n_features) when the given problem is binary.

    intercept_ : float or ndarray of shape (n_targets,)
        Independent term in decision function. Set to 0.0 if
        ``fit_intercept = False``.

    n_iter_ : None or ndarray of shape (n_targets,)
        Actual number of iterations for each target. Available only for
        sag and lsqr solvers. Other solvers will return None.

    classes_ : ndarray of shape (n_classes,)
        The classes labels.

    n_features_in_ : int
        Number of features seen during :term:`fit`.

        .. versionadded:: 0.24

    feature_names_in_ : ndarray of shape (`n_features_in_`,)
        Names of features seen during :term:`fit`. Defined only when `X`
        has feature names that are all strings.

        .. versionadded:: 1.0

    solver_ : str
        The solver that was used at fit time by the computational
        routines.

        .. versionadded:: 1.5

    See Also
    --------
    Ridge : Ridge regression.
    RidgeClassifierCV :  Ridge classifier with built-in cross validation.

    Notes
    -----
    For multi-class classification, n_class classifiers are trained in
    a one-versus-all approach. Concretely, this is implemented by taking
    advantage of the multi-variate response support in Ridge.

    Examples
    --------
    >>> from sklearn.datasets import load_breast_cancer
    >>> from sklearn.linear_model import RidgeClassifier
    >>> X, y = load_breast_cancer(return_X_y=True)
    >>> clf = RidgeClassifier().fit(X, y)
    >>> clf.score(X, y)
    0.9595...
    """

    _parameter_constraints: dict = {
        **_BaseRidge._parameter_constraints,
        "class_weight": [dict, StrOptions({"balanced"}), None],
    }

    def __init__(
        self,
        alpha=1.0,
        *,
        fit_intercept=True,
        copy_X=True,
        max_iter=None,
        tol=1e-4,
        class_weight=None,
        solver="auto",
        positive=False,
        random_state=None,
    ):
        super().__init__(
            alpha=alpha,
            fit_intercept=fit_intercept,
            copy_X=copy_X,
            max_iter=max_iter,
            tol=tol,
            solver=solver,
            positive=positive,
            random_state=random_state,
        )
        self.class_weight = class_weight

    @_fit_context(prefer_skip_nested_validation=True)
    def fit(self, X, y, sample_weight=None):
        """Fit Ridge classifier model.

        Parameters
        ----------
        X : {ndarray, sparse matrix} of shape (n_samples, n_features)
            Training data.

        y : ndarray of shape (n_samples,)
            Target values.

        sample_weight : float or ndarray of shape (n_samples,), default=None
            Individual weights for each sample. If given a float, every sample
            will have the same weight.

            .. versionadded:: 0.17
               *sample_weight* support to RidgeClassifier.

        Returns
        -------
        self : object
            Instance of the estimator.
        """
        X, y, sample_weight, Y = self._prepare_data(X, y, sample_weight, self.solver)

        super().fit(X, Y, sample_weight=sample_weight)
        return self

    def __sklearn_tags__(self):
        tags = super().__sklearn_tags__()
        tags.input_tags.sparse = (self.solver != "svd") and (
            self.solver != "cholesky" or not self.fit_intercept
        )
        return tags


def _check_gcv_mode(X, gcv_mode):
    # svd only implemented for dense X
    if gcv_mode == "svd":
        if sparse.issparse(X):
            # TODO(1.11) raise ValueError
            msg = (
                "The 'svd' mode is not supported for sparse X, we fallback to "
                "`gcv_mode='eigen'`. Passing `gcv_mode='svd'` on sparse X will raise "
                "an error in 1.11, use the default or pass `gcv_mode='eigen'` to "
                "suppress this warning."
            )
            warnings.warn(msg, FutureWarning)
        else:
            return "svd"

    # All other cases ("auto", "eigen")
    # fallbacks to gram (n <= p) or cov (p < n)
    n, p = X.shape
    return "gram" if n <= p else "cov"


def _find_smallest_angle(query, vectors):
    """Find the column of vectors that is most aligned with the query.

    Both query and the columns of vectors must have their l2 norm equal to 1.

    Parameters
    ----------
    query : ndarray of shape (n_samples,)
        Normalized query vector.

    vectors : ndarray of shape (n_samples, n_features)
        Vectors to which we compare query, as columns. Must be normalized.
    """
    xp, _ = get_namespace(query)
    abs_cosine = xp.abs(query @ vectors)
    index = xp.argmax(abs_cosine)
    return index


class _X_CenterStackOp(sparse.linalg.LinearOperator):
    """Behaves as centered and scaled X with an added intercept column.

    This operator behaves as
    np.hstack([X - sqrt_sw[:, None] * X_mean, sqrt_sw[:, None]])
    """

    def __init__(self, X, X_mean, sqrt_sw):
        n_samples, n_features = X.shape
        super().__init__(X.dtype, (n_samples, n_features + 1))
        self.X = X
        self.X_mean = X_mean
        self.sqrt_sw = sqrt_sw

    def _matvec(self, v):
        v = v.ravel()
        return (
            safe_sparse_dot(self.X, v[:-1], dense_output=True)
            - self.sqrt_sw * self.X_mean.dot(v[:-1])
            + v[-1] * self.sqrt_sw
        )

    def _matmat(self, v):
        return (
            safe_sparse_dot(self.X, v[:-1], dense_output=True)
            - self.sqrt_sw[:, None] * self.X_mean.dot(v[:-1])
            + v[-1] * self.sqrt_sw[:, None]
        )

    def _transpose(self):
        return _XT_CenterStackOp(self.X, self.X_mean, self.sqrt_sw)


class _XT_CenterStackOp(sparse.linalg.LinearOperator):
    """Behaves as transposed centered and scaled X with an intercept column.

    This operator behaves as
    np.hstack([X - sqrt_sw[:, None] * X_mean, sqrt_sw[:, None]]).T
    """

    def __init__(self, X, X_mean, sqrt_sw):
        n_samples, n_features = X.shape
        super().__init__(X.dtype, (n_features + 1, n_samples))
        self.X = X
        self.X_mean = X_mean
        self.sqrt_sw = sqrt_sw

    def _matvec(self, v):
        v = v.ravel()
        n_features = self.shape[0]
        res = np.empty(n_features, dtype=self.X.dtype)
        res[:-1] = safe_sparse_dot(self.X.T, v, dense_output=True) - (
            self.X_mean * self.sqrt_sw.dot(v)
        )
        res[-1] = np.dot(v, self.sqrt_sw)
        return res

    def _matmat(self, v):
        n_features = self.shape[0]
        res = np.empty((n_features, v.shape[1]), dtype=self.X.dtype)
        res[:-1] = safe_sparse_dot(self.X.T, v, dense_output=True) - self.X_mean[
            :, None
        ] * self.sqrt_sw.dot(v)
        res[-1] = np.dot(self.sqrt_sw, v)
        return res


class _IdentityRegressor(RegressorMixin, BaseEstimator):
    """Fake regressor which will directly output the prediction."""

    def decision_function(self, y_predict):
        return y_predict

    def predict(self, y_predict):
        return y_predict


class _IdentityClassifier(LinearClassifierMixin, BaseEstimator):
    """Fake classifier which will directly output the prediction.

    We inherit from LinearClassifierMixin to get the proper shape for the
    output `y`.
    """

    def __init__(self, classes):
        self.classes_ = classes

    def decision_function(self, y_predict):
        return y_predict


class _RidgeGCV(LinearModel):
    """Ridge regression with built-in Leave-one-out Cross-Validation.

    This class is not intended to be used directly. Use RidgeCV instead.

    `_RidgeGCV` uses a Generalized Cross-Validation for model selection. It's an
    efficient approximation of leave-one-out cross-validation (LOO-CV), where instead of
    computing multiple models by excluding one data point at a time, it uses an
    algebraic shortcut to approximate the LOO-CV error, making it faster and
    computationally more efficient.

    Using a naive grid-search approach with a leave-one-out cross-validation in contrast
    requires to fit `n_samples` models to compute the prediction error for each sample
    and then to repeat this process for each alpha in the grid.

    Here, the prediction error for each sample is computed by solving a **single**
    linear system (in other words a single model) via a matrix factorization (i.e.
    eigendecomposition or SVD) solving the problem stated in the Notes section. Finally,
    we need to repeat this process for each alpha in the grid. The detailed complexity
    is further discussed in Sect. 4 in [1].

    This algebraic approach is only applicable for regularized least squares
    problems. It could potentially be extended to kernel ridge regression.

    See the Notes section and references for more details regarding the formulation
    and the linear system that is solved.

    Notes
    -----

    1. Unweighted and no intercept

    We start by the simplest case `fit_intercept=False` and `sample_weight=None`.
    The other cases (see below) reduce to this one after proper scaling/centering
    of the design matrix X.

    The design matrix X has shape (n, p) = (n_samples, n_features).

    Let G = (K + alpha*Id_n) where K = X X' is the Gram matrix and Id_n is the
    identity matrix of size n.

    Let H = (C + alpha*Id_p) where C = X' X is the covariance matrix and Id_p
    is the identity matrix of size p.

    The solution of the regularized least square (fitted `coef_`) is given by:
    w = H^-1 X' y = X' c where c = G^-1 y.

    Let loov (resp looe) be the leave-one-out values (resp errors), that is the
    vector of predictions (resp errors) for each single observation when the model
    was fitted with all examples but this example. As shown in [1]:
    looe = y - loov = c / d where d = diag(G^-1).

    The best score (negative mean squared error or user-provided scoring) is
    stored in the `best_score_` attribute, and the selected hyperparameter in
    `alpha_`.

    2. Leveraging a precomputed matrix decomposition

    The leave-one-out errors and coefficients can be efficiently computed for any
    alpha from the SVD of X, or the eigendecomposition of K = X X' or C = X' X.

    Reduced SVD X = U S V' when n < p (wide X)
    Let D = 1 / (S^2 + alpha)
    c = U D U' y
    d = diag(U D U')
    w = V S / (S^2 + alpha) U' y

    Eigendecomposition K = U L U'
    Let D = 1 / (L + alpha)
    c = U D U' y
    d = diag(U D U')
    w = X' c.

    Reduced SVD X = U S V' when p < n (long X)
    Let M = alpha / (S^2 + alpha) - 1
    alpha c = y + U M U' y
    alpha d = 1 + diag(U M U')
    w = V S / (S^2 + alpha) U' y

    Eigendecomposition C = V L V'
    H^-1 = V 1 / (L + alpha) V'
    alpha c = y - X H^-1 X' y
    alpha d = 1 - diag(X H^-1 X')
    w = H^-1 X' y

    3. Fitting with intercept or sample weights

    Fitting with intercept and/or sample weights reduces to the unweigthed no
    intercept case after centering and/or rescaling of X and y, as done in
    `_preprocess_data`:
    X <- sqrt(s) (X - X_mean)
    y <- sqrt(s) (y - y_mean)

    The returned looe are also rescaled by sample weights:
    looe <- sqrt(s) looe

    If we fit an intercept, there is the following correction term:
    d <- d - sqrt(s) * G^-1 sqrt(s) / sum(s)

    References
    ----------
    .. [1] R. Rifkin and R. Lippert (2007). "Notes on Regularized Least Squares."
           https://dspace.mit.edu/bitstream/handle/1721.1/37318/MIT-CSAIL-TR-2007-025.pdf
    .. [2] R. Rifkin (2007). "Regularized Least Squares."
           https://www.mit.edu/~9.520/spring07/Classes/rlsslides.pdf
    """

    def __init__(
        self,
        alphas=(0.1, 1.0, 10.0),
        *,
        fit_intercept=True,
        scoring=None,
        copy_X=True,
        gcv_mode=None,
        store_cv_results=False,
        is_clf=False,
        alpha_per_target=False,
    ):
        self.alphas = alphas
        self.fit_intercept = fit_intercept
        self.scoring = scoring
        self.copy_X = copy_X
        self.gcv_mode = gcv_mode
        self.store_cv_results = store_cv_results
        self.is_clf = is_clf
        self.alpha_per_target = alpha_per_target

    @staticmethod
    def _decomp_diag(v_prime, Q):
        # compute diagonal of the matrix: dot(Q, dot(diag(v_prime), Q.T))
        xp, _ = get_namespace(v_prime, Q)
        return xp.sum(v_prime * Q**2, axis=1)

    @staticmethod
    def _diag_dot(D, B):
        xp, _ = get_namespace(D, B)
        # compute dot(diag(D), B)
        if len(B.shape) > 1:
            # handle case where B is > 1-d
            D = D[(slice(None),) + (None,) * (len(B.shape) - 1)]
        return D * B

    def _compute_gram(self, X, X_mean, sqrt_sw):
        """Computes the Gram matrix X X' with possible centering.

        Parameters
        ----------
        X : {ndarray, sparse matrix, sparse array} of shape (n_samples, n_features)
            The preprocessed design matrix.

        X_mean : ndarray of shape (n_feature,)
            The weighted mean of X for each feature.

        sqrt_sw : ndarray of shape (n_samples,)
            Square roots of sample weights.

        Returns
        -------
        gram : ndarray of shape (n_samples, n_samples)
            The Gram matrix.

        Notes
        -----
        When self.fit_intercept is False no centering is done.

        When X is dense the centering has been done in preprocessing
        so the mean is 0 and we just compute X X'.

        When X is sparse it has not been centered in preprocessing, but
        it has been scaled by sqrt_sw. The centered X is never actually
        computed because centering would break the sparsity of X.
        """
        center = self.fit_intercept and sparse.issparse(X)
        if not center:
            # in this case centering has been done in preprocessing
            # or we are not fitting an intercept.
            return safe_sparse_dot(X, X.T, dense_output=True)
        # X is sparse and fit_intercept is True
        # centered matrix = X - sqrt_sw X_mean'
        X_Xm = safe_sparse_dot(X, X_mean, dense_output=True)
        return (
            safe_sparse_dot(X, X.T, dense_output=True)
            - X_Xm[:, None] * sqrt_sw[None, :]
            - sqrt_sw[:, None] * X_Xm[None, :]
            + (X_mean @ X_mean) * sqrt_sw[:, None] * sqrt_sw[None, :]
        )

    def _compute_covariance(self, X, X_mean, sqrt_sw):
        """Computes covariance matrix X' X with possible centering.

        Parameters
        ----------
        X : {ndarray, sparse matrix, sparse array} of shape (n_samples, n_features)
            The preprocessed design matrix.

        X_mean : ndarray of shape (n_feature,)
            The weighted mean of X for each feature.

        sqrt_sw : ndarray of shape (n_samples,)
            Square roots of sample weights.

        Returns
        -------
        covariance : ndarray of shape (n_features, n_features)
            The covariance matrix.

        Notes
        -----
        When self.fit_intercept is False no centering is done.

        When X is dense the centering has been done in preprocessing
        so the mean is 0 and we just compute X' X.

        When X is sparse it has not been centered in preprocessing, but
        it has been scaled by sqrt_sw. The centered X is never actually
        computed because centering would break the sparsity of X.
        """
        center = self.fit_intercept and sparse.issparse(X)
        if not center:
            # in this case centering has been done in preprocessing
            # or we are not fitting an intercept.
            return safe_sparse_dot(X.T, X, dense_output=True)
        # X is sparse and fit_intercept is True
        # centered matrix = X - sqrt_sw X_mean'
        sw_sum = sqrt_sw @ sqrt_sw
        return (
            safe_sparse_dot(X.T, X, dense_output=True)
            - sw_sum * X_mean[:, None] * X_mean[None, :]
        )

    def _sparse_multidot_diag(self, X, A, X_mean, sqrt_sw):
        """Compute the diagonal of X A X' with possible centering.

        Parameters
        ----------
        X : {ndarray, sparse matrix, sparse array} of shape (n_samples, n_features)
            The preprocessed design matrix.

        A : ndarray of shape (n_features, n_features)
            The inner matrix.

        X_mean : ndarray of shape (n_feature,)
            The weighted mean of X for each feature.

        sqrt_sw : ndarray of shape (n_samples,)
            Square roots of sample weights.

        Returns
        -------
        diag : np.ndarray, shape (n_samples,)
            The computed diagonal.

        Notes
        -----
        When self.fit_intercept is False no centering is done.

        When X is dense the centering has been done in preprocessing
        so the mean is 0 and we just compute diag(X A X').

        When X is sparse it has not been centered in preprocessing, but
        it has been scaled by sqrt_sw. The centered X is never actually
        computed because centering would break the sparsity of X.
        """
        xp, _ = get_namespace(X)
        XA = X @ A
        if sparse.isspmatrix(X):
            # sparse matrix use multiply for element wise multiplication
            XAX = np.ravel(X.multiply(XA).sum(axis=1))
        else:
            XAX = xp.sum(XA * X, axis=1)
        center = self.fit_intercept and sparse.issparse(X)
        if not center:
            # in this case centering has been done in preprocessing
            # or we are not fitting an intercept.
            return XAX
        # X is sparse and fit_intercept is True
        # centered matrix = X - sqrt_sw X_mean'
        XA_Xm = XA @ X_mean
        A_Xm = A @ X_mean
        sw = sqrt_sw * sqrt_sw
        return XAX - 2 * sqrt_sw * XA_Xm + sw * (X_mean @ A_Xm)

    def _eigen_decompose_gram(self, X, X_mean, y, sqrt_sw):
        """Eigendecomposition of Gram matrix X X'"""
        xp, is_array_api = get_namespace(X)
        K = self._compute_gram(X, X_mean, sqrt_sw)
        eigvals, Q = xp.linalg.eigh(K)
        QT_y = Q.T @ y
        QT_sqrt_sw = Q.T @ sqrt_sw
        XT = X.T
        return eigvals, Q, QT_y, QT_sqrt_sw, XT, X_mean

    def _solve_eigen_gram(
        self, alpha, y, sqrt_sw, eigvals, Q, QT_y, QT_sqrt_sw, XT, X_mean
    ):
        """Compute looe and coef when we have a decomposition of X X'"""
        D = 1.0 / (eigvals + alpha)
        c = Q @ self._diag_dot(D, QT_y)
        d = self._decomp_diag(D, Q)
        if self.fit_intercept:
            sw_sum = sqrt_sw @ sqrt_sw
            Ginv_sqrt_sw = Q @ self._diag_dot(D, QT_sqrt_sw)
            d -= Ginv_sqrt_sw * sqrt_sw / sw_sum
        if y.ndim == 2:
            d = d[:, None]
        XT_c = XT @ c
        if self.fit_intercept and sparse.issparse(XT):
            # centered matrix = X - sqrt_sw X_mean'
            if y.ndim == 2:
                XT_c -= X_mean[:, None] * (sqrt_sw @ c)
            else:
                XT_c -= X_mean * (sqrt_sw @ c)
        looe = c / d
        coef = XT_c
        return looe, coef

    def _eigen_decompose_covariance(self, X, X_mean, y, sqrt_sw):
        """Eigendecomposition of covariance matrix X' X"""
        xp, is_array_api = get_namespace(X)
        cov = self._compute_covariance(X, X_mean, sqrt_sw)
        eigvals, V = xp.linalg.eigh(cov)
        XT_y = safe_sparse_dot(X.T, y, dense_output=True)
        XT_sqrt_sw = safe_sparse_dot(X.T, sqrt_sw, dense_output=True)
        if self.fit_intercept and sparse.issparse(X):
            # centered matrix = X - sqrt_sw X_mean'
            if y.ndim == 2:
                XT_y -= X_mean[:, None] * (sqrt_sw @ y)
            else:
                XT_y -= X_mean * (sqrt_sw @ y)
            XT_sqrt_sw -= X_mean * (sqrt_sw @ sqrt_sw)
        return eigvals, V, X, X_mean, XT_y, XT_sqrt_sw

    def _solve_eigen_covariance(
        self, alpha, y, sqrt_sw, eigvals, V, X, X_mean, XT_y, XT_sqrt_sw
    ):
        """Compute looe and coef when we have a decomposition of X' X"""
        D = 1 / (eigvals + alpha)
        Hinv = (V * D) @ V.T
        Hinv_XT_y = Hinv @ XT_y
        Hinv_XT_sqrt_sw = Hinv @ XT_sqrt_sw
        X_Hinv_XT_y = safe_sparse_dot(X, Hinv_XT_y, dense_output=True)
        X_Hinv_XT_sqrt_sw = safe_sparse_dot(X, Hinv_XT_sqrt_sw, dense_output=True)
        if self.fit_intercept and sparse.issparse(X):
            # centered = X - sqrt_sw X_mean'
            if y.ndim == 2:
                X_Hinv_XT_y -= sqrt_sw[:, None] * (X_mean @ Hinv_XT_y)
            else:
                X_Hinv_XT_y -= sqrt_sw * (X_mean @ Hinv_XT_y)
            X_Hinv_XT_sqrt_sw -= sqrt_sw * (X_mean @ Hinv_XT_sqrt_sw)
        alpha_c = y - X_Hinv_XT_y
        alpha_d = 1 - self._sparse_multidot_diag(X, Hinv, X_mean, sqrt_sw)
        if self.fit_intercept:
            sw_sum = sqrt_sw @ sqrt_sw
            alpha_Ginv_sqrt_sw = sqrt_sw - X_Hinv_XT_sqrt_sw
            alpha_d -= alpha_Ginv_sqrt_sw * sqrt_sw / sw_sum
        if y.ndim == 2:
            alpha_d = alpha_d[:, None]
        looe = alpha_c / alpha_d
        coef = Hinv_XT_y
        return looe, coef

    def _svd_decompose_design_matrix(self, X, X_mean, y, sqrt_sw):
        """Reduced SVD decomposition of X"""
        xp, _ = get_namespace(X)
        # reduced svd
        U, singvals, VT = xp.linalg.svd(X, full_matrices=False)
        UT_y = U.T @ y
        UT_sqrt_sw = U.T @ sqrt_sw
        V = VT.T
        return singvals, U, V, UT_y, UT_sqrt_sw

    def _solve_svd_design_matrix_long(
        self, alpha, y, sqrt_sw, singvals, U, V, UT_y, UT_sqrt_sw
    ):
        """Compute looe and coef when we have an SVD decomposition of X.

        Long X case (n_features < n_samples).
        """
        M = alpha / (singvals**2 + alpha) - 1
        alpha_c = U @ self._diag_dot(M, UT_y) + y
        alpha_d = self._decomp_diag(M, U) + 1
        if self.fit_intercept:
            sw_sum = sqrt_sw @ sqrt_sw
            alpha_Ginv_sqrt_sw = U @ self._diag_dot(M, UT_sqrt_sw) + sqrt_sw
            alpha_d -= alpha_Ginv_sqrt_sw * sqrt_sw / sw_sum
        if y.ndim == 2:
            # handle case where y is 2-d
            alpha_d = alpha_d[:, None]
        looe = alpha_c / alpha_d
        coef = V @ self._diag_dot(singvals / (singvals**2 + alpha), UT_y)
        return looe, coef

    def _solve_svd_design_matrix_wide(
        self, alpha, y, sqrt_sw, singvals, U, V, UT_y, UT_sqrt_sw
    ):
        """Compute looe and coef when we have an SVD decomposition of X.

        Wide X case (n_samples < n_features).
        """
        alpha_D = alpha / (singvals**2 + alpha)
        alpha_c = U @ self._diag_dot(alpha_D, UT_y)
        alpha_d = self._decomp_diag(alpha_D, U)
        if self.fit_intercept:
            sw_sum = sqrt_sw @ sqrt_sw
            alpha_Ginv_sqrt_sw = U @ self._diag_dot(alpha_D, UT_sqrt_sw)
            alpha_d -= alpha_Ginv_sqrt_sw * sqrt_sw / sw_sum
        if y.ndim == 2:
            # handle case where y is 2-d
            alpha_d = alpha_d[:, None]
        looe = alpha_c / alpha_d
        coef = V @ self._diag_dot(singvals / (singvals**2 + alpha), UT_y)
        return looe, coef

    def _solve_svd_design_matrix(
        self, alpha, y, sqrt_sw, singvals, U, V, UT_y, UT_sqrt_sw
    ):
        """Compute looe and coef when we have an SVD decomposition of X."""
        n_samples = U.shape[0]
        n_features = V.shape[0]
        if n_samples <= n_features:
            return self._solve_svd_design_matrix_wide(
                alpha, y, sqrt_sw, singvals, U, V, UT_y, UT_sqrt_sw
            )
        else:
            return self._solve_svd_design_matrix_long(
                alpha, y, sqrt_sw, singvals, U, V, UT_y, UT_sqrt_sw
            )

    def fit(self, X, y, sample_weight=None, score_params=None):
        """Fit Ridge regression model with gcv.

        Parameters
        ----------
        X : {ndarray, sparse matrix} of shape (n_samples, n_features)
            Training data. Will be cast to float64 if necessary.

        y : ndarray of shape (n_samples,) or (n_samples, n_targets)
            Target values. Will be cast to float64 if necessary.

        sample_weight : float or ndarray of shape (n_samples,), default=None
            Individual weights for each sample. If given a float, every sample
            will have the same weight. Note that the scale of `sample_weight`
            has an impact on the loss; i.e. multiplying all weights by `k`
            is equivalent to setting `alpha / k`.

        score_params : dict, default=None
            Parameters to be passed to the underlying scorer.

            .. versionadded:: 1.5
                See :ref:`Metadata Routing User Guide <metadata_routing>` for
                more details.

        Returns
        -------
        self : object
        """
        xp, is_array_api, device_ = get_namespace_and_device(X)
        y, sample_weight = move_to(y, sample_weight, xp=xp, device=device_)
        if (is_array_api and xp.isdtype(X.dtype, "real floating")) or getattr(
            getattr(X, "dtype", None), "kind", None
        ) == "f":
            original_floating_dtype = X.dtype
        else:
            # for X that does not have a simple dtype (e.g. pandas dataframe)
            # the attributes will be stored in the dtype chosen by
            # `validate_data``, i.e. np.float64
            original_floating_dtype = None
        # Using float32 can be numerically unstable for this estimator. So if
        # the array API namespace and device allow, convert the input values
        # to float64 whenever possible before converting the results back to
        # float32.
        dtype = _max_precision_float_dtype(xp, device=device_)
        X, y = validate_data(
            self,
            X,
            y,
            accept_sparse=["csr", "csc", "coo"],
            dtype=dtype,
            multi_output=True,
            y_numeric=True,
        )

        # alpha_per_target cannot be used in classifier mode. All subclasses
        # of _RidgeGCV that are classifiers keep alpha_per_target at its
        # default value: False, so the condition below should never happen.
        assert not (self.is_clf and self.alpha_per_target)

        if sample_weight is not None:
            sample_weight = _check_sample_weight(sample_weight, X, dtype=X.dtype)

        self.alphas = np.asarray(self.alphas)

        unscaled_y = y
        X, y, X_offset, y_offset, X_scale, sqrt_sw = _preprocess_data(
            X,
            y,
            fit_intercept=self.fit_intercept,
            copy=self.copy_X,
            sample_weight=sample_weight,
            rescale_with_sw=True,
        )

        gcv_mode = _check_gcv_mode(X, self.gcv_mode)

        n_samples, n_features = X.shape
        if gcv_mode == "gram":
            decompose = self._eigen_decompose_gram
            solve = self._solve_eigen_gram
        elif gcv_mode == "cov":
            decompose = self._eigen_decompose_covariance
            solve = self._solve_eigen_covariance
        elif gcv_mode == "svd":
            decompose = self._svd_decompose_design_matrix
            solve = self._solve_svd_design_matrix
        else:
            raise ValueError(f"Unknown {gcv_mode=}")

        if sqrt_sw is None:
            sqrt_sw = xp.ones(n_samples, dtype=X.dtype, device=device_)

        decomposition = decompose(X, X_offset, y, sqrt_sw)

        n_y = 1 if y.ndim == 1 else y.shape[1]
        if (
            isinstance(self.alphas, numbers.Number)
            or getattr(self.alphas, "ndim", None) == 0
        ):
            alphas = [float(self.alphas)]
        else:
            alphas = list(map(float, self.alphas))
        n_alphas = len(alphas)

        if self.store_cv_results:
            self.cv_results_ = xp.empty(
                (n_samples * n_y, n_alphas), dtype=X.dtype, device=device_
            )

        best_coef, best_score, best_alpha = None, None, None

        for i, alpha in enumerate(alphas):
            looe, coef = solve(float(alpha), y, sqrt_sw, *decomposition)
            if self.scoring is None:
                squared_errors = looe**2
                alpha_score = self._score_without_scorer(squared_errors=squared_errors)
                if self.store_cv_results:
                    self.cv_results_[:, i] = _ravel(squared_errors)
            else:
                predictions = y - looe
                # Rescale predictions back to original scale
                if sample_weight is not None:  # avoid the unnecessary division by ones
                    if predictions.ndim > 1:
                        predictions /= sqrt_sw[:, None]
                    else:
                        predictions /= sqrt_sw
                predictions += y_offset

                if self.store_cv_results:
                    self.cv_results_[:, i] = _ravel(predictions)

                score_params = score_params or {}
                alpha_score = self._score(
                    predictions=predictions,
                    y=unscaled_y,
                    n_y=n_y,
                    scorer=self.scoring,
                    score_params=score_params,
                )

            # Keep track of the best model
            if best_score is None:
                # initialize
                if self.alpha_per_target and n_y > 1:
                    best_coef = coef
                    best_score = xp.reshape(alpha_score, shape=(-1,))
                    best_alpha = xp.full(n_y, alpha, device=device_)
                else:
                    best_coef = coef
                    best_score = alpha_score
                    best_alpha = alpha
            else:
                # update
                if self.alpha_per_target and n_y > 1:
                    to_update = alpha_score > best_score
                    best_coef[:, to_update] = coef[:, to_update]
                    best_score[to_update] = alpha_score[to_update]
                    best_alpha[to_update] = alpha
                elif alpha_score > best_score:
                    best_coef, best_score, best_alpha = coef, alpha_score, alpha

        self.alpha_ = best_alpha
        self.best_score_ = best_score
        self.coef_ = best_coef
        if y.ndim == 2:
            self.coef_ = self.coef_.T
        if y.ndim == 1 or y.shape[1] == 1:
            self.coef_ = _ravel(self.coef_)

        self._set_intercept(X_offset, y_offset, X_scale)

        if self.store_cv_results:
            if y.ndim == 1:
                cv_results_shape = n_samples, n_alphas
            else:
                cv_results_shape = n_samples, n_y, n_alphas
            self.cv_results_ = xp.reshape(self.cv_results_, shape=cv_results_shape)

        if original_floating_dtype:
            if type(self.intercept_) is not float:
                self.intercept_ = xp.astype(
                    self.intercept_, original_floating_dtype, copy=False
                )
            self.coef_ = xp.astype(self.coef_, original_floating_dtype, copy=False)
            if self.store_cv_results:
                self.cv_results_ = xp.astype(
                    self.cv_results_, original_floating_dtype, copy=False
                )

        return self

    def _score_without_scorer(self, squared_errors):
        """Performs scoring using squared errors when the scorer is None."""
        xp, _ = get_namespace(squared_errors)
        if self.alpha_per_target:
            _score = xp.mean(-squared_errors, axis=0)
        else:
            _score = xp.mean(-squared_errors)

        return _score

    def _score(self, *, predictions, y, n_y, scorer, score_params):
        """Performs scoring with the specified scorer using the
        predictions and the true y values.
        """
        xp, _, device_ = get_namespace_and_device(y)
        if self.is_clf:
            identity_estimator = _IdentityClassifier(
                classes=xp.arange(n_y, device=device_)
            )
            _score = scorer(
                identity_estimator,
                predictions,
                xp.argmax(y, axis=1),
                **score_params,
            )
        else:
            identity_estimator = _IdentityRegressor()
            if self.alpha_per_target:
                _score = xp.asarray(
                    [
                        scorer(
                            identity_estimator,
                            predictions[:, j],
                            y[:, j],
                            **score_params,
                        )
                        for j in range(n_y)
                    ],
                    device=device_,
                )
            else:
                _score = scorer(
                    identity_estimator,
                    predictions,
                    y,
                    **score_params,
                )

        return _score

    def __sklearn_tags__(self):
        tags = super().__sklearn_tags__()
        # Required since this is neither a RegressorMixin nor a ClassifierMixin
        tags.target_tags.required = True
        return tags


class _BaseRidgeCV(LinearModel):
    _parameter_constraints: dict = {
        "alphas": ["array-like", Interval(Real, 0, None, closed="neither")],
        "fit_intercept": ["boolean"],
        "scoring": [StrOptions(set(get_scorer_names())), callable, None],
        "cv": ["cv_object"],
        "gcv_mode": [StrOptions({"auto", "svd", "eigen"}), None],
        "store_cv_results": ["boolean"],
        "alpha_per_target": ["boolean"],
    }

    def __init__(
        self,
        alphas=(0.1, 1.0, 10.0),
        *,
        fit_intercept=True,
        scoring=None,
        cv=None,
        gcv_mode=None,
        store_cv_results=False,
        alpha_per_target=False,
    ):
        self.alphas = alphas
        self.fit_intercept = fit_intercept
        self.scoring = scoring
        self.cv = cv
        self.gcv_mode = gcv_mode
        self.store_cv_results = store_cv_results
        self.alpha_per_target = alpha_per_target

    def fit(self, X, y, sample_weight=None, **params):
        """Fit Ridge regression model with cv.

        Parameters
        ----------
        X : ndarray of shape (n_samples, n_features)
            Training data. If using GCV, will be cast to float64
            if necessary.

        y : ndarray of shape (n_samples,) or (n_samples, n_targets)
            Target values. Will be cast to X's dtype if necessary.

        sample_weight : float or ndarray of shape (n_samples,), default=None
            Individual weights for each sample. If given a float, every sample
            will have the same weight.

        **params : dict, default=None
            Extra parameters for the underlying scorer.

            .. versionadded:: 1.5
                Only available if `enable_metadata_routing=True`,
                which can be set by using
                ``sklearn.set_config(enable_metadata_routing=True)``.
                See :ref:`Metadata Routing User Guide <metadata_routing>` for
                more details.

        Returns
        -------
        self : object
            Fitted estimator.

        Notes
        -----
        When sample_weight is provided, the selected hyperparameter may depend
        on whether we use leave-one-out cross-validation (cv=None)
        or another form of cross-validation, because only leave-one-out
        cross-validation takes the sample weights into account when computing
        the validation score.
        """
        _raise_for_params(params, self, "fit")
        cv = self.cv
        scorer = self._get_scorer()

        # `_RidgeGCV` does not work for alpha = 0
        if cv is None:
            check_scalar_alpha = partial(
                check_scalar,
                target_type=numbers.Real,
                min_val=0.0,
                include_boundaries="neither",
            )
        else:
            check_scalar_alpha = partial(
                check_scalar,
                target_type=numbers.Real,
                min_val=0.0,
                include_boundaries="left",
            )

        if isinstance(self.alphas, (np.ndarray, list, tuple)):
            n_alphas = 1 if np.ndim(self.alphas) == 0 else len(self.alphas)
            if n_alphas != 1:
                for index, alpha in enumerate(self.alphas):
                    alpha = check_scalar_alpha(alpha, f"alphas[{index}]")
            else:
                self.alphas[0] = check_scalar_alpha(self.alphas[0], "alphas")
        alphas = np.asarray(self.alphas)

        if sample_weight is not None:
            params["sample_weight"] = sample_weight

        if cv is None:
            if _routing_enabled():
                routed_params = process_routing(
                    self,
                    "fit",
                    **params,
                )
            else:
                routed_params = Bunch(scorer=Bunch(score={}))
                if sample_weight is not None:
                    routed_params.scorer.score["sample_weight"] = sample_weight

            # reset `scorer` variable to original user-intend if no scoring is passed
            if self.scoring is None:
                scorer = None

            estimator = _RidgeGCV(
                alphas,
                fit_intercept=self.fit_intercept,
                scoring=scorer,
                gcv_mode=self.gcv_mode,
                store_cv_results=self.store_cv_results,
                is_clf=is_classifier(self),
                alpha_per_target=self.alpha_per_target,
            )
            estimator.fit(
                X,
                y,
                sample_weight=sample_weight,
                score_params=routed_params.scorer.score,
            )
            self.alpha_ = estimator.alpha_
            self.best_score_ = estimator.best_score_
            if self.store_cv_results:
                self.cv_results_ = estimator.cv_results_
        else:
            if self.store_cv_results:
                raise ValueError("cv!=None and store_cv_results=True are incompatible")
            if self.alpha_per_target:
                raise ValueError("cv!=None and alpha_per_target=True are incompatible")

            parameters = {"alpha": alphas}
            solver = "sparse_cg" if sparse.issparse(X) else "auto"
            model = RidgeClassifier if is_classifier(self) else Ridge
            estimator = model(
                fit_intercept=self.fit_intercept,
                solver=solver,
            )
            if _routing_enabled():
                estimator.set_fit_request(sample_weight=True)

            grid_search = GridSearchCV(
                estimator,
                parameters,
                cv=cv,
                scoring=scorer,
            )

            grid_search.fit(X, y, **params)
            estimator = grid_search.best_estimator_
            self.alpha_ = grid_search.best_estimator_.alpha
            self.best_score_ = grid_search.best_score_

        self.coef_ = estimator.coef_
        self.intercept_ = estimator.intercept_
        self.n_features_in_ = estimator.n_features_in_
        if hasattr(estimator, "feature_names_in_"):
            self.feature_names_in_ = estimator.feature_names_in_

        return self

    def get_metadata_routing(self):
        """Get metadata routing of this object.

        Please check :ref:`User Guide <metadata_routing>` on how the routing
        mechanism works.

        .. versionadded:: 1.5

        Returns
        -------
        routing : MetadataRouter
            A :class:`~sklearn.utils.metadata_routing.MetadataRouter` encapsulating
            routing information.
        """
        router = (
            MetadataRouter(owner=self)
            .add_self_request(self)
            .add(
                scorer=self._get_scorer(),
                method_mapping=MethodMapping().add(caller="fit", callee="score"),
            )
            .add(
                splitter=self.cv,
                method_mapping=MethodMapping().add(caller="fit", callee="split"),
            )
        )
        return router

    def _get_scorer(self):
        """Make sure the scorer is weighted if necessary.

        This uses `self._get_scorer_instance()` implemented in child objects to get the
        raw scorer instance of the estimator, which will be ignored if `self.scoring` is
        not None.
        """
        if _routing_enabled() and self.scoring is None:
            # This estimator passes an array of 1s as sample_weight even if
            # sample_weight is not provided by the user. Therefore we need to
            # always request it. But we don't set it if it's passed explicitly
            # by the user.
            return self._get_scorer_instance().set_score_request(sample_weight=True)

        return check_scoring(estimator=self, scoring=self.scoring, allow_none=True)

    def __sklearn_tags__(self):
        tags = super().__sklearn_tags__()
        tags.array_api_support = True
        tags.input_tags.sparse = True
        return tags


class RidgeCV(MultiOutputMixin, RegressorMixin, _BaseRidgeCV):
    """Ridge regression with built-in cross-validation.

    See glossary entry for :term:`cross-validation estimator`.

    By default, it performs efficient Leave-One-Out Cross-Validation.

    Read more in the :ref:`User Guide <ridge_regression>`.

    Parameters
    ----------
    alphas : array-like of shape (n_alphas,), default=(0.1, 1.0, 10.0)
        Array of alpha values to try.
        Regularization strength; must be a positive float. Regularization
        improves the conditioning of the problem and reduces the variance of
        the estimates. Larger values specify stronger regularization.
        Alpha corresponds to ``1 / (2C)`` in other linear models such as
        :class:`~sklearn.linear_model.LogisticRegression` or
        :class:`~sklearn.svm.LinearSVC`.
        If using Leave-One-Out cross-validation, alphas must be strictly positive.

        For an example on how regularization strength affects the model coefficients,
        see :ref:`sphx_glr_auto_examples_linear_model_plot_ridge_coeffs.py`.

    fit_intercept : bool, default=True
        Whether to calculate the intercept for this model. If set
        to false, no intercept will be used in calculations
        (i.e. data is expected to be centered).

    scoring : str, callable, default=None
        The scoring method to use for cross-validation. Options:

        - str: see :ref:`scoring_string_names` for options.
        - callable: a scorer callable object (e.g., function) with signature
          ``scorer(estimator, X, y)``. See :ref:`scoring_callable` for details.
        - `None`: negative :ref:`mean squared error <mean_squared_error>` if cv is
          None (i.e. when using leave-one-out cross-validation), or
          :ref:`coefficient of determination <r2_score>` (:math:`R^2`) otherwise.

    cv : int, cross-validation generator or an iterable, default=None
        Determines the cross-validation splitting strategy.
        Possible inputs for cv are:

        - None, to use the efficient Leave-One-Out cross-validation
        - integer, to specify the number of folds,
        - :term:`CV splitter`,
        - an iterable yielding (train, test) splits as arrays of indices.

        For integer/None inputs, if ``y`` is binary or multiclass,
        :class:`~sklearn.model_selection.StratifiedKFold` is used, else,
        :class:`~sklearn.model_selection.KFold` is used.

        Refer :ref:`User Guide <cross_validation>` for the various
        cross-validation strategies that can be used here.

    gcv_mode : {'auto', 'svd', 'eigen'}, default='auto'
        Flag indicating which strategy to use when performing
        Leave-One-Out Cross-Validation. Options are::

            'auto' : same as 'eigen'
            'svd' : use singular value decomposition of X when X is dense,
                fallback to 'eigen' when X is sparse
            'eigen' : use eigendecomposition of X X' when n_samples <= n_features
                or X' X when n_features < n_samples

        The 'auto' mode is the default and is intended to pick the cheaper
        option depending on the shape and sparsity of the training data.

    store_cv_results : bool, default=False
        Flag indicating if the cross-validation values corresponding to
        each alpha should be stored in the ``cv_results_`` attribute (see
        below). This flag is only compatible with ``cv=None`` (i.e. using
        Leave-One-Out Cross-Validation).

        .. versionchanged:: 1.5
            Parameter name changed from `store_cv_values` to `store_cv_results`.

    alpha_per_target : bool, default=False
        Flag indicating whether to optimize the alpha value (picked from the
        `alphas` parameter list) for each target separately (for multi-output
        settings: multiple prediction targets). When set to `True`, after
        fitting, the `alpha_` attribute will contain a value for each target.
        When set to `False`, a single alpha is used for all targets.

        .. versionadded:: 0.24

    Attributes
    ----------
    cv_results_ : ndarray of shape (n_samples, n_alphas) or \
            shape (n_samples, n_targets, n_alphas), optional
        Cross-validation values for each alpha (only available if
        ``store_cv_results=True`` and ``cv=None``). After ``fit()`` has been
        called, this attribute will contain the mean squared errors if
        `scoring is None` otherwise it will contain standardized per point
        prediction values.

        .. versionchanged:: 1.5
            `cv_values_` changed to `cv_results_`.

    coef_ : ndarray of shape (n_features) or (n_targets, n_features)
        Weight vector(s).

    intercept_ : float or ndarray of shape (n_targets,)
        Independent term in decision function. Set to 0.0 if
        ``fit_intercept = False``.

    alpha_ : float or ndarray of shape (n_targets,)
        Estimated regularization parameter, or, if ``alpha_per_target=True``,
        the estimated regularization parameter for each target.

    best_score_ : float or ndarray of shape (n_targets,)
        Score of base estimator with best alpha, or, if
        ``alpha_per_target=True``, a score for each target.

        .. versionadded:: 0.23

    n_features_in_ : int
        Number of features seen during :term:`fit`.

        .. versionadded:: 0.24

    feature_names_in_ : ndarray of shape (`n_features_in_`,)
        Names of features seen during :term:`fit`. Defined only when `X`
        has feature names that are all strings.

        .. versionadded:: 1.0

    See Also
    --------
    Ridge : Ridge regression.
    RidgeClassifier : Classifier based on ridge regression on {-1, 1} labels.
    RidgeClassifierCV : Ridge classifier with built-in cross validation.

    Examples
    --------
    >>> from sklearn.datasets import load_diabetes
    >>> from sklearn.linear_model import RidgeCV
    >>> X, y = load_diabetes(return_X_y=True)
    >>> clf = RidgeCV(alphas=[1e-3, 1e-2, 1e-1, 1]).fit(X, y)
    >>> clf.score(X, y)
    0.5166...
    """

    @_fit_context(prefer_skip_nested_validation=True)
    def fit(self, X, y, sample_weight=None, **params):
        """Fit Ridge regression model with cv.

        Parameters
        ----------
        X : ndarray of shape (n_samples, n_features)
            Training data. If using GCV, will be cast to float64
            if necessary.

        y : ndarray of shape (n_samples,) or (n_samples, n_targets)
            Target values. Will be cast to X's dtype if necessary.

        sample_weight : float or ndarray of shape (n_samples,), default=None
            Individual weights for each sample. If given a float, every sample
            will have the same weight.

        **params : dict, default=None
            Parameters to be passed to the underlying scorer.

            .. versionadded:: 1.5
                Only available if `enable_metadata_routing=True`,
                which can be set by using
                ``sklearn.set_config(enable_metadata_routing=True)``.
                See :ref:`Metadata Routing User Guide <metadata_routing>` for
                more details.

        Returns
        -------
        self : object
            Fitted estimator.

        Notes
        -----
        When sample_weight is provided, the selected hyperparameter may depend
        on whether we use leave-one-out cross-validation (cv=None)
        or another form of cross-validation, because only leave-one-out
        cross-validation takes the sample weights into account when computing
        the validation score.
        """
        super().fit(X, y, sample_weight=sample_weight, **params)
        return self

    def _get_scorer_instance(self):
        """Return a scorer which corresponds to what's defined in RegressorMixin
        parent class. This is used for routing `sample_weight`.
        """
        return get_scorer("r2")


class RidgeClassifierCV(_RidgeClassifierMixin, _BaseRidgeCV):
    """Ridge classifier with built-in cross-validation.

    See glossary entry for :term:`cross-validation estimator`.

    By default, it performs Leave-One-Out Cross-Validation. Currently,
    only the n_features > n_samples case is handled efficiently.

    Read more in the :ref:`User Guide <ridge_regression>`.

    Parameters
    ----------
    alphas : array-like of shape (n_alphas,), default=(0.1, 1.0, 10.0)
        Array of alpha values to try.
        Regularization strength; must be a positive float. Regularization
        improves the conditioning of the problem and reduces the variance of
        the estimates. Larger values specify stronger regularization.
        Alpha corresponds to ``1 / (2C)`` in other linear models such as
        :class:`~sklearn.linear_model.LogisticRegression` or
        :class:`~sklearn.svm.LinearSVC`.
        If using Leave-One-Out cross-validation, alphas must be strictly positive.

        For an example on how regularization strength affects the model coefficients,
        see :ref:`sphx_glr_auto_examples_linear_model_plot_ridge_coeffs.py`.

    fit_intercept : bool, default=True
        Whether to calculate the intercept for this model. If set
        to false, no intercept will be used in calculations
        (i.e. data is expected to be centered).

    scoring : str, callable, default=None
        The scoring method to use for cross-validation. Options:

        - str: see :ref:`scoring_string_names` for options.
        - callable: a scorer callable object (e.g., function) with signature
          ``scorer(estimator, X, y)``. See :ref:`scoring_callable` for details.
        - `None`: negative :ref:`mean squared error <mean_squared_error>` if cv is
          None (i.e. when using leave-one-out cross-validation), or
          :ref:`accuracy <accuracy_score>` otherwise.

    cv : int, cross-validation generator or an iterable, default=None
        Determines the cross-validation splitting strategy.
        Possible inputs for cv are:

        - None, to use the efficient Leave-One-Out cross-validation
        - integer, to specify the number of folds,
        - :term:`CV splitter`,
        - an iterable yielding (train, test) splits as arrays of indices.

        Refer :ref:`User Guide <cross_validation>` for the various
        cross-validation strategies that can be used here.

    class_weight : dict or 'balanced', default=None
        Weights associated with classes in the form ``{class_label: weight}``.
        If not given, all classes are supposed to have weight one.

        The "balanced" mode uses the values of y to automatically adjust
        weights inversely proportional to class frequencies in the input data
        as ``n_samples / (n_classes * np.bincount(y))``.

    store_cv_results : bool, default=False
        Flag indicating if the cross-validation results corresponding to
        each alpha should be stored in the ``cv_results_`` attribute (see
        below). This flag is only compatible with ``cv=None`` (i.e. using
        Leave-One-Out Cross-Validation).

        .. versionchanged:: 1.5
            Parameter name changed from `store_cv_values` to `store_cv_results`.

    Attributes
    ----------
    cv_results_ : ndarray of shape (n_samples, n_targets, n_alphas), optional
        Cross-validation results for each alpha (only if ``store_cv_results=True`` and
        ``cv=None``). After ``fit()`` has been called, this attribute will
        contain the mean squared errors if `scoring is None` otherwise it
        will contain standardized per point prediction values.

        .. versionchanged:: 1.5
            `cv_values_` changed to `cv_results_`.

    coef_ : ndarray of shape (1, n_features) or (n_targets, n_features)
        Coefficient of the features in the decision function.

        ``coef_`` is of shape (1, n_features) when the given problem is binary.

    intercept_ : float or ndarray of shape (n_targets,)
        Independent term in decision function. Set to 0.0 if
        ``fit_intercept = False``.

    alpha_ : float
        Estimated regularization parameter.

    best_score_ : float
        Score of base estimator with best alpha.

        .. versionadded:: 0.23

    classes_ : ndarray of shape (n_classes,)
        The classes labels.

    n_features_in_ : int
        Number of features seen during :term:`fit`.

        .. versionadded:: 0.24

    feature_names_in_ : ndarray of shape (`n_features_in_`,)
        Names of features seen during :term:`fit`. Defined only when `X`
        has feature names that are all strings.

        .. versionadded:: 1.0

    See Also
    --------
    Ridge : Ridge regression.
    RidgeClassifier : Ridge classifier.
    RidgeCV : Ridge regression with built-in cross validation.

    Notes
    -----
    For multi-class classification, n_class classifiers are trained in
    a one-versus-all approach. Concretely, this is implemented by taking
    advantage of the multi-variate response support in Ridge.

    Examples
    --------
    >>> from sklearn.datasets import load_breast_cancer
    >>> from sklearn.linear_model import RidgeClassifierCV
    >>> X, y = load_breast_cancer(return_X_y=True)
    >>> clf = RidgeClassifierCV(alphas=[1e-3, 1e-2, 1e-1, 1]).fit(X, y)
    >>> clf.score(X, y)
    0.9630...
    """

    _parameter_constraints: dict = {
        **_BaseRidgeCV._parameter_constraints,
        "class_weight": [dict, StrOptions({"balanced"}), None],
    }
    for param in ("gcv_mode", "alpha_per_target"):
        _parameter_constraints.pop(param)

    def __init__(
        self,
        alphas=(0.1, 1.0, 10.0),
        *,
        fit_intercept=True,
        scoring=None,
        cv=None,
        class_weight=None,
        store_cv_results=False,
    ):
        super().__init__(
            alphas=alphas,
            fit_intercept=fit_intercept,
            scoring=scoring,
            cv=cv,
            store_cv_results=store_cv_results,
        )
        self.class_weight = class_weight

    @_fit_context(prefer_skip_nested_validation=True)
    def fit(self, X, y, sample_weight=None, **params):
        """Fit Ridge classifier with cv.

        Parameters
        ----------
        X : ndarray of shape (n_samples, n_features)
            Training vectors, where `n_samples` is the number of samples
            and `n_features` is the number of features. When using GCV,
            will be cast to float64 if necessary.

        y : ndarray of shape (n_samples,)
            Target values. Will be cast to X's dtype if necessary.

        sample_weight : float or ndarray of shape (n_samples,), default=None
            Individual weights for each sample. If given a float, every sample
            will have the same weight.

        **params : dict, default=None
            Parameters to be passed to the underlying scorer.

            .. versionadded:: 1.5
                Only available if `enable_metadata_routing=True`,
                which can be set by using
                ``sklearn.set_config(enable_metadata_routing=True)``.
                See :ref:`Metadata Routing User Guide <metadata_routing>` for
                more details.

        Returns
        -------
        self : object
            Fitted estimator.
        """
        # `RidgeClassifier` does not accept "sag" or "saga" solver and thus support
        # csr, csc, and coo sparse matrices. By using solver="eigen" we force to accept
        # all sparse format.
        X, y, sample_weight, Y = self._prepare_data(X, y, sample_weight, solver="eigen")

        # If cv is None, gcv mode will be used and we used the binarized Y
        # since y will not be binarized in _RidgeGCV estimator.
        # If cv is not None, a GridSearchCV with some RidgeClassifier
        # estimators are used where y will be binarized. Thus, we pass y
        # instead of the binarized Y.
        target = Y if self.cv is None else y
        super().fit(X, target, sample_weight=sample_weight, **params)
        return self
