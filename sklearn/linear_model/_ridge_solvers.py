import warnings

import numpy as np
from scipy import linalg
from scipy import sparse
from scipy import optimize

from ..utils.extmath import safe_sparse_dot
from ..exceptions import ConvergenceWarning
from scipy.sparse import linalg as sp_linalg


def _cholesky_helper(X, y, alpha, n_features, n_samples):
    if n_features > n_samples:
        K = safe_sparse_dot(X, X.T, dense_output=True)
        try:
            dual_coef = _solve_cholesky_kernel(K, y, alpha)

            coef = safe_sparse_dot(X.T, dual_coef, dense_output=True).T
        except linalg.LinAlgError:
            # use SVD solver if matrix is singular
            coef = _solve_svd(X, y, alpha)
    else:
        try:
            coef = _solve_cholesky(X, y, alpha)
        except linalg.LinAlgError:
            # use SVD solver if matrix is singular
            coef = _solve_svd(X, y, alpha)
    return coef


def _solve_sparse_cg(
    X,
    y,
    alpha,
    max_iter=None,
    tol=1e-3,
    verbose=0,
    X_offset=None,
    X_scale=None,
    sample_weight_sqrt=None,
):
    if sample_weight_sqrt is None:
        sample_weight_sqrt = np.ones(X.shape[0], dtype=X.dtype)

    def _get_rescaled_operator(X):

        X_offset_scale = X_offset / X_scale

        def matvec(b):
            return X.dot(b) - sample_weight_sqrt * b.dot(X_offset_scale)

        def rmatvec(b):
            return X.T.dot(b) - X_offset_scale * b.dot(sample_weight_sqrt)

        X1 = sparse.linalg.LinearOperator(shape=X.shape, matvec=matvec, rmatvec=rmatvec)
        return X1

    n_samples, n_features = X.shape

    if X_offset is None or X_scale is None:
        X1 = sp_linalg.aslinearoperator(X)
    else:
        X1 = _get_rescaled_operator(X)

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
            # FIXME atol
            try:
                coef, info = sp_linalg.cg(C, y_column, tol=tol, atol="legacy")
            except TypeError:
                # old scipy
                coef, info = sp_linalg.cg(C, y_column, tol=tol)
            coefs[i] = X1.rmatvec(coef)
        else:
            # linear ridge
            # w = inv(X^t X + alpha*Id) * X.T y
            y_column = X1.rmatvec(y_column)
            C = sp_linalg.LinearOperator(
                (n_features, n_features), matvec=mv, dtype=X.dtype
            )
            # FIXME atol
            try:
                coefs[i], info = sp_linalg.cg(
                    C, y_column, maxiter=max_iter, tol=tol, atol="legacy"
                )
            except TypeError:
                # old scipy
                coefs[i], info = sp_linalg.cg(C, y_column, maxiter=max_iter, tol=tol)

        if info < 0:
            raise ValueError("Failed with error code %d" % info)

        if max_iter is None and info > 0 and verbose:
            warnings.warn(
                "sparse_cg did not converge after %d iterations." % info,
                ConvergenceWarning,
            )

    return coefs


def _solve_lsqr(X, y, alpha, max_iter=None, tol=1e-3):
    n_samples, n_features = X.shape
    coefs = np.empty((y.shape[1], n_features), dtype=X.dtype)
    n_iter = np.empty(y.shape[1], dtype=np.int32)

    # According to the lsqr documentation, alpha = damp^2.
    sqrt_alpha = np.sqrt(alpha)

    for i in range(y.shape[1]):
        y_column = y[:, i]
        info = sp_linalg.lsqr(
            X, y_column, damp=sqrt_alpha[i], atol=tol, btol=tol, iter_lim=max_iter
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
        return linalg.solve(A, Xy, sym_pos=True, overwrite_a=True).T
    else:
        coefs = np.empty([n_targets, n_features], dtype=X.dtype)
        for coef, target, current_alpha in zip(coefs, Xy.T, alpha):
            A.flat[:: n_features + 1] += current_alpha
            coef[:] = linalg.solve(A, target, sym_pos=True, overwrite_a=False).ravel()
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
            dual_coef = linalg.solve(K, y, sym_pos=True, overwrite_a=False)
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
                K, target, sym_pos=True, overwrite_a=False
            ).ravel()

            K.flat[:: n_samples + 1] -= current_alpha

        if has_sw:
            dual_coefs *= sw[np.newaxis, :]

        return dual_coefs.T


def _solve_svd(X, y, alpha):
    U, s, Vt = linalg.svd(X, full_matrices=False)
    idx = s > 1e-15  # same default value as scipy.linalg.pinv
    s_nnz = s[idx][:, np.newaxis]
    UTy = np.dot(U.T, y)
    d = np.zeros((s.size, alpha.size), dtype=X.dtype)
    d[idx] = s_nnz / (s_nnz**2 + alpha)
    d_UT_y = d * UTy
    return np.dot(Vt.T, d_UT_y).T


def _solve_lbfgs(
    X,
    y,
    alpha,
    positive=True,
    max_iter=None,
    tol=1e-3,
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
    if positive:
        config["bounds"] = [(0, np.inf)] * n_features

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
                "The lbfgs solver did not converge. Try increasing max_iter "
                f"or tol. Currently: max_iter={max_iter} and tol={tol}",
                ConvergenceWarning,
            )
        coefs[i] = result["x"]

    return coefs
