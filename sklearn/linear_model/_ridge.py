"""
Ridge regression
"""

# Author: Mathieu Blondel <mathieu@mblondel.org>
#         Reuben Fletcher-Costin <reuben.fletchercostin@gmail.com>
#         Fabian Pedregosa <fabian@fseoane.net>
#         Michael Eickenberg <michael.eickenberg@nsup.org>
# License: BSD 3 clause


from abc import ABCMeta, abstractmethod
import warnings

import numpy as np
from scipy import linalg
from scipy import sparse
from scipy.sparse import linalg as sp_linalg

from ._base import LinearClassifierMixin, LinearModel, _rescale_data
from ._sag import sag_solver
from ..base import RegressorMixin, MultiOutputMixin
from ..utils.extmath import safe_sparse_dot
from ..utils.extmath import row_norms
from ..utils import check_X_y
from ..utils import check_array
from ..utils import check_consistent_length
from ..utils import compute_sample_weight
from ..utils import column_or_1d
from ..utils.validation import _check_sample_weight
from ..preprocessing import LabelBinarizer
from ..model_selection import GridSearchCV
from ..metrics import check_scoring
from ..exceptions import ConvergenceWarning
from ..utils.sparsefuncs import mean_variance_axis


def _solve_sparse_cg(X, y, alpha, max_iter=None, tol=1e-3, verbose=0,
                     X_offset=None, X_scale=None):

    def _get_rescaled_operator(X):

        X_offset_scale = X_offset / X_scale

        def matvec(b):
            return X.dot(b) - b.dot(X_offset_scale)

        def rmatvec(b):
            return X.T.dot(b) - X_offset_scale * np.sum(b)

        X1 = sparse.linalg.LinearOperator(shape=X.shape,
                                          matvec=matvec,
                                          rmatvec=rmatvec)
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
                (n_samples, n_samples), matvec=mv, dtype=X.dtype)
            # FIXME atol
            try:
                coef, info = sp_linalg.cg(C, y_column, tol=tol, atol='legacy')
            except TypeError:
                # old scipy
                coef, info = sp_linalg.cg(C, y_column, tol=tol)
            coefs[i] = X1.rmatvec(coef)
        else:
            # linear ridge
            # w = inv(X^t X + alpha*Id) * X.T y
            y_column = X1.rmatvec(y_column)
            C = sp_linalg.LinearOperator(
                (n_features, n_features), matvec=mv, dtype=X.dtype)
            # FIXME atol
            try:
                coefs[i], info = sp_linalg.cg(C, y_column, maxiter=max_iter,
                                              tol=tol, atol='legacy')
            except TypeError:
                # old scipy
                coefs[i], info = sp_linalg.cg(C, y_column, maxiter=max_iter,
                                              tol=tol)

        if info < 0:
            raise ValueError("Failed with error code %d" % info)

        if max_iter is None and info > 0 and verbose:
            warnings.warn("sparse_cg did not converge after %d iterations." %
                          info, ConvergenceWarning)

    return coefs


def _solve_lsqr(X, y, alpha, max_iter=None, tol=1e-3):
    n_samples, n_features = X.shape
    coefs = np.empty((y.shape[1], n_features), dtype=X.dtype)
    n_iter = np.empty(y.shape[1], dtype=np.int32)

    # According to the lsqr documentation, alpha = damp^2.
    sqrt_alpha = np.sqrt(alpha)

    for i in range(y.shape[1]):
        y_column = y[:, i]
        info = sp_linalg.lsqr(X, y_column, damp=sqrt_alpha[i],
                              atol=tol, btol=tol, iter_lim=max_iter)
        coefs[i] = info[0]
        n_iter[i] = info[2]

    return coefs, n_iter


def _solve_cholesky(X, y, alpha):
    # w = inv(X^t X + alpha*Id) * X.T y
    n_samples, n_features = X.shape
    n_targets = y.shape[1]

    A = safe_sparse_dot(X.T, X, dense_output=True)
    Xy = safe_sparse_dot(X.T, y, dense_output=True)

    one_alpha = np.array_equal(alpha, len(alpha) * [alpha[0]])

    if one_alpha:
        A.flat[::n_features + 1] += alpha[0]
        return linalg.solve(A, Xy, sym_pos=True,
                            overwrite_a=True).T
    else:
        coefs = np.empty([n_targets, n_features], dtype=X.dtype)
        for coef, target, current_alpha in zip(coefs, Xy.T, alpha):
            A.flat[::n_features + 1] += current_alpha
            coef[:] = linalg.solve(A, target, sym_pos=True,
                                   overwrite_a=False).ravel()
            A.flat[::n_features + 1] -= current_alpha
        return coefs


def _solve_cholesky_kernel(K, y, alpha, sample_weight=None, copy=False):
    # dual_coef = inv(X X^t + alpha*Id) y
    n_samples = K.shape[0]
    n_targets = y.shape[1]

    if copy:
        K = K.copy()

    alpha = np.atleast_1d(alpha)
    one_alpha = (alpha == alpha[0]).all()
    has_sw = isinstance(sample_weight, np.ndarray) \
        or sample_weight not in [1.0, None]

    if has_sw:
        # Unlike other solvers, we need to support sample_weight directly
        # because K might be a pre-computed kernel.
        sw = np.sqrt(np.atleast_1d(sample_weight))
        y = y * sw[:, np.newaxis]
        K *= np.outer(sw, sw)

    if one_alpha:
        # Only one penalty, we can solve multi-target problems in one time.
        K.flat[::n_samples + 1] += alpha[0]

        try:
            # Note: we must use overwrite_a=False in order to be able to
            #       use the fall-back solution below in case a LinAlgError
            #       is raised
            dual_coef = linalg.solve(K, y, sym_pos=True,
                                     overwrite_a=False)
        except np.linalg.LinAlgError:
            warnings.warn("Singular matrix in solving dual problem. Using "
                          "least-squares solution instead.")
            dual_coef = linalg.lstsq(K, y)[0]

        # K is expensive to compute and store in memory so change it back in
        # case it was user-given.
        K.flat[::n_samples + 1] -= alpha[0]

        if has_sw:
            dual_coef *= sw[:, np.newaxis]

        return dual_coef
    else:
        # One penalty per target. We need to solve each target separately.
        dual_coefs = np.empty([n_targets, n_samples], K.dtype)

        for dual_coef, target, current_alpha in zip(dual_coefs, y.T, alpha):
            K.flat[::n_samples + 1] += current_alpha

            dual_coef[:] = linalg.solve(K, target, sym_pos=True,
                                        overwrite_a=False).ravel()

            K.flat[::n_samples + 1] -= current_alpha

        if has_sw:
            dual_coefs *= sw[np.newaxis, :]

        return dual_coefs.T


def _solve_svd(X, y, alpha):
    U, s, Vt = linalg.svd(X, full_matrices=False)
    idx = s > 1e-15  # same default value as scipy.linalg.pinv
    s_nnz = s[idx][:, np.newaxis]
    UTy = np.dot(U.T, y)
    d = np.zeros((s.size, alpha.size), dtype=X.dtype)
    d[idx] = s_nnz / (s_nnz ** 2 + alpha)
    d_UT_y = d * UTy
    return np.dot(Vt.T, d_UT_y).T


def _get_valid_accept_sparse(is_X_sparse, solver):
    if is_X_sparse and solver in ['auto', 'sag', 'saga']:
        return 'csr'
    else:
        return ['csr', 'csc', 'coo']


def ridge_regression(X, y, alpha, sample_weight=None, solver='auto',
                     max_iter=None, tol=1e-3, verbose=0, random_state=None,
                     return_n_iter=False, return_intercept=False,
                     check_input=True):
    """Solve the ridge equation by the method of normal equations.

    Read more in the :ref:`User Guide <ridge_regression>`.

    Parameters
    ----------
    X : {array-like, sparse matrix, LinearOperator} of shape \
        (n_samples, n_features)
        Training data

    y : array-like of shape (n_samples,) or (n_samples, n_targets)
        Target values

    alpha : float or array-like of shape (n_targets,)
        Regularization strength; must be a positive float. Regularization
        improves the conditioning of the problem and reduces the variance of
        the estimates. Larger values specify stronger regularization.
        Alpha corresponds to ``C^-1`` in other linear models such as
        LogisticRegression or LinearSVC. If an array is passed, penalties are
        assumed to be specific to the targets. Hence they must correspond in
        number.

    sample_weight : float or array-like of shape (n_samples,), default=None
        Individual weights for each sample. If given a float, every sample
        will have the same weight. If sample_weight is not None and
        solver='auto', the solver will be set to 'cholesky'.

        .. versionadded:: 0.17

    solver : {'auto', 'svd', 'cholesky', 'lsqr', 'sparse_cg', 'sag', 'saga'}
        Solver to use in the computational routines:

        - 'auto' chooses the solver automatically based on the type of data.

        - 'svd' uses a Singular Value Decomposition of X to compute the Ridge
          coefficients. More stable for singular matrices than
          'cholesky'.

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


        All last five solvers support both dense and sparse data. However, only
        'sag' and 'sparse_cg' supports sparse input when`fit_intercept` is
        True.

        .. versionadded:: 0.17
           Stochastic Average Gradient descent solver.
        .. versionadded:: 0.19
           SAGA solver.

    max_iter : int, optional
        Maximum number of iterations for conjugate gradient solver.
        For the 'sparse_cg' and 'lsqr' solvers, the default value is determined
        by scipy.sparse.linalg. For 'sag' and saga solver, the default value is
        1000.

    tol : float, default=1e-3
        Precision of the solution.

    verbose : int, default=0
        Verbosity level. Setting verbose > 0 will display additional
        information depending on the solver used.

    random_state : int, RandomState instance or None, optional, default None
        The seed of the pseudo random number generator to use when shuffling
        the data.  If int, random_state is the seed used by the random number
        generator; If RandomState instance, random_state is the random number
        generator; If None, the random number generator is the RandomState
        instance used by `np.random`. Used when ``solver`` == 'sag'.

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
    coef : array of shape (n_features,) or (n_targets, n_features)
        Weight vector(s).

    n_iter : int, optional
        The actual number of iteration performed by the solver.
        Only returned if `return_n_iter` is True.

    intercept : float or array of shape (n_targets,)
        The intercept of the model. Only returned if `return_intercept`
        is True and if X is a scipy sparse array.

    Notes
    -----
    This function won't compute the intercept.
    """
    return _ridge_regression(X, y, alpha,
                             sample_weight=sample_weight,
                             solver=solver,
                             max_iter=max_iter,
                             tol=tol,
                             verbose=verbose,
                             random_state=random_state,
                             return_n_iter=return_n_iter,
                             return_intercept=return_intercept,
                             X_scale=None,
                             X_offset=None,
                             check_input=check_input)


def _ridge_regression(X, y, alpha, sample_weight=None, solver='auto',
                      max_iter=None, tol=1e-3, verbose=0, random_state=None,
                      return_n_iter=False, return_intercept=False,
                      X_scale=None, X_offset=None, check_input=True):

    has_sw = sample_weight is not None

    if solver == 'auto':
        if return_intercept:
            # only sag supports fitting intercept directly
            solver = "sag"
        elif not sparse.issparse(X):
            solver = "cholesky"
        else:
            solver = "sparse_cg"

    if solver not in ('sparse_cg', 'cholesky', 'svd', 'lsqr', 'sag', 'saga'):
        raise ValueError("Known solvers are 'sparse_cg', 'cholesky', 'svd'"
                         " 'lsqr', 'sag' or 'saga'. Got %s." % solver)

    if return_intercept and solver != 'sag':
        raise ValueError("In Ridge, only 'sag' solver can directly fit the "
                         "intercept. Please change solver to 'sag' or set "
                         "return_intercept=False.")

    if check_input:
        _dtype = [np.float64, np.float32]
        _accept_sparse = _get_valid_accept_sparse(sparse.issparse(X), solver)
        X = check_array(X, accept_sparse=_accept_sparse, dtype=_dtype,
                        order="C")
        y = check_array(y, dtype=X.dtype, ensure_2d=False, order=None)
    check_consistent_length(X, y)

    n_samples, n_features = X.shape

    if y.ndim > 2:
        raise ValueError("Target y has the wrong shape %s" % str(y.shape))

    ravel = False
    if y.ndim == 1:
        y = y.reshape(-1, 1)
        ravel = True

    n_samples_, n_targets = y.shape

    if n_samples != n_samples_:
        raise ValueError("Number of samples in X and y does not correspond:"
                         " %d != %d" % (n_samples, n_samples_))

    if has_sw:
        sample_weight = _check_sample_weight(sample_weight, X, dtype=X.dtype)

        if solver not in ['sag', 'saga']:
            # SAG supports sample_weight directly. For other solvers,
            # we implement sample_weight via a simple rescaling.
            X, y = _rescale_data(X, y, sample_weight)

    # There should be either 1 or n_targets penalties
    alpha = np.asarray(alpha, dtype=X.dtype).ravel()
    if alpha.size not in [1, n_targets]:
        raise ValueError("Number of targets and number of penalties "
                         "do not correspond: %d != %d"
                         % (alpha.size, n_targets))

    if alpha.size == 1 and n_targets > 1:
        alpha = np.repeat(alpha, n_targets)

    n_iter = None
    if solver == 'sparse_cg':
        coef = _solve_sparse_cg(X, y, alpha,
                                max_iter=max_iter,
                                tol=tol,
                                verbose=verbose,
                                X_offset=X_offset,
                                X_scale=X_scale)

    elif solver == 'lsqr':
        coef, n_iter = _solve_lsqr(X, y, alpha, max_iter, tol)

    elif solver == 'cholesky':
        if n_features > n_samples:
            K = safe_sparse_dot(X, X.T, dense_output=True)
            try:
                dual_coef = _solve_cholesky_kernel(K, y, alpha)

                coef = safe_sparse_dot(X.T, dual_coef, dense_output=True).T
            except linalg.LinAlgError:
                # use SVD solver if matrix is singular
                solver = 'svd'
        else:
            try:
                coef = _solve_cholesky(X, y, alpha)
            except linalg.LinAlgError:
                # use SVD solver if matrix is singular
                solver = 'svd'

    elif solver in ['sag', 'saga']:
        # precompute max_squared_sum for all targets
        max_squared_sum = row_norms(X, squared=True).max()

        coef = np.empty((y.shape[1], n_features), dtype=X.dtype)
        n_iter = np.empty(y.shape[1], dtype=np.int32)
        intercept = np.zeros((y.shape[1], ), dtype=X.dtype)
        for i, (alpha_i, target) in enumerate(zip(alpha, y.T)):
            init = {'coef': np.zeros((n_features + int(return_intercept), 1),
                                     dtype=X.dtype)}
            coef_, n_iter_, _ = sag_solver(
                X, target.ravel(), sample_weight, 'squared', alpha_i, 0,
                max_iter, tol, verbose, random_state, False, max_squared_sum,
                init,
                is_saga=solver == 'saga')
            if return_intercept:
                coef[i] = coef_[:-1]
                intercept[i] = coef_[-1]
            else:
                coef[i] = coef_
            n_iter[i] = n_iter_

        if intercept.shape[0] == 1:
            intercept = intercept[0]
        coef = np.asarray(coef)

    if solver == 'svd':
        if sparse.issparse(X):
            raise TypeError('SVD solver does not support sparse'
                            ' inputs currently')
        coef = _solve_svd(X, y, alpha)

    if ravel:
        # When y was passed as a 1d-array, we flatten the coefficients.
        coef = coef.ravel()

    if return_n_iter and return_intercept:
        return coef, n_iter, intercept
    elif return_intercept:
        return coef, intercept
    elif return_n_iter:
        return coef, n_iter
    else:
        return coef


class _BaseRidge(LinearModel, metaclass=ABCMeta):
    @abstractmethod
    def __init__(self, alpha=1.0, fit_intercept=True, normalize=False,
                 copy_X=True, max_iter=None, tol=1e-3, solver="auto",
                 random_state=None):
        self.alpha = alpha
        self.fit_intercept = fit_intercept
        self.normalize = normalize
        self.copy_X = copy_X
        self.max_iter = max_iter
        self.tol = tol
        self.solver = solver
        self.random_state = random_state

    def fit(self, X, y, sample_weight=None):

        # all other solvers work at both float precision levels
        _dtype = [np.float64, np.float32]
        _accept_sparse = _get_valid_accept_sparse(sparse.issparse(X),
                                                  self.solver)
        X, y = check_X_y(X, y,
                         accept_sparse=_accept_sparse,
                         dtype=_dtype,
                         multi_output=True, y_numeric=True)
        if sparse.issparse(X) and self.fit_intercept:
            if self.solver not in ['auto', 'sparse_cg', 'sag']:
                raise ValueError(
                    "solver='{}' does not support fitting the intercept "
                    "on sparse data. Please set the solver to 'auto' or "
                    "'sparse_cg', 'sag', or set `fit_intercept=False`"
                    .format(self.solver))
            if (self.solver == 'sag' and self.max_iter is None and
                    self.tol > 1e-4):
                warnings.warn(
                    '"sag" solver requires many iterations to fit '
                    'an intercept with sparse inputs. Either set the '
                    'solver to "auto" or "sparse_cg", or set a low '
                    '"tol" and a high "max_iter" (especially if inputs are '
                    'not standardized).')
                solver = 'sag'
            else:
                solver = 'sparse_cg'
        else:
            solver = self.solver

        if sample_weight is not None:
            sample_weight = _check_sample_weight(sample_weight, X,
                                                 dtype=X.dtype)

        # when X is sparse we only remove offset from y
        X, y, X_offset, y_offset, X_scale = self._preprocess_data(
            X, y, self.fit_intercept, self.normalize, self.copy_X,
            sample_weight=sample_weight, return_mean=True)

        if solver == 'sag' and sparse.issparse(X) and self.fit_intercept:
            self.coef_, self.n_iter_, self.intercept_ = _ridge_regression(
                X, y, alpha=self.alpha, sample_weight=sample_weight,
                max_iter=self.max_iter, tol=self.tol, solver=self.solver,
                random_state=self.random_state, return_n_iter=True,
                return_intercept=True, check_input=False)
            # add the offset which was subtracted by _preprocess_data
            self.intercept_ += y_offset

        else:
            if sparse.issparse(X) and self.fit_intercept:
                # required to fit intercept with sparse_cg solver
                params = {'X_offset': X_offset, 'X_scale': X_scale}
            else:
                # for dense matrices or when intercept is set to 0
                params = {}

            self.coef_, self.n_iter_ = _ridge_regression(
                X, y, alpha=self.alpha, sample_weight=sample_weight,
                max_iter=self.max_iter, tol=self.tol, solver=solver,
                random_state=self.random_state, return_n_iter=True,
                return_intercept=False, check_input=False, **params)
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
    alpha : {float, array-like of shape (n_targets,)}, default=1.0
        Regularization strength; must be a positive float. Regularization
        improves the conditioning of the problem and reduces the variance of
        the estimates. Larger values specify stronger regularization.
        Alpha corresponds to ``C^-1`` in other linear models such as
        LogisticRegression or LinearSVC. If an array is passed, penalties are
        assumed to be specific to the targets. Hence they must correspond in
        number.

    fit_intercept : bool, default=True
        Whether to calculate the intercept for this model. If set
        to false, no intercept will be used in calculations
        (i.e. data is expected to be centered).

    normalize : bool, default=False
        This parameter is ignored when ``fit_intercept`` is set to False.
        If True, the regressors X will be normalized before regression by
        subtracting the mean and dividing by the l2-norm.
        If you wish to standardize, please use
        :class:`sklearn.preprocessing.StandardScaler` before calling ``fit``
        on an estimator with ``normalize=False``.

    copy_X : bool, default=True
        If True, X will be copied; else, it may be overwritten.

    max_iter : int, optional
        Maximum number of iterations for conjugate gradient solver.
        For 'sparse_cg' and 'lsqr' solvers, the default value is determined
        by scipy.sparse.linalg. For 'sag' solver, the default value is 1000.

    tol : float, default=1e-3
        Precision of the solution.

    solver : {'auto', 'svd', 'cholesky', 'lsqr', 'sparse_cg', 'sag', 'saga'}
        Solver to use in the computational routines:

        - 'auto' chooses the solver automatically based on the type of data.

        - 'svd' uses a Singular Value Decomposition of X to compute the Ridge
          coefficients. More stable for singular matrices than
          'cholesky'.

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
          its improved, unbiased version named SAGA. Both methods also use an
          iterative procedure, and are often faster than other solvers when
          both n_samples and n_features are large. Note that 'sag' and
          'saga' fast convergence is only guaranteed on features with
          approximately the same scale. You can preprocess the data with a
          scaler from sklearn.preprocessing.

        All last five solvers support both dense and sparse data. However, only
        'sparse_cg' supports sparse input when `fit_intercept` is True.

        .. versionadded:: 0.17
           Stochastic Average Gradient descent solver.
        .. versionadded:: 0.19
           SAGA solver.

    random_state : int, RandomState instance or None, optional, default None
        The seed of the pseudo random number generator to use when shuffling
        the data.  If int, random_state is the seed used by the random number
        generator; If RandomState instance, random_state is the random number
        generator; If None, the random number generator is the RandomState
        instance used by `np.random`. Used when ``solver`` == 'sag'.

        .. versionadded:: 0.17
           *random_state* to support Stochastic Average Gradient.

    Attributes
    ----------
    coef_ : array of shape (n_features,) or (n_targets, n_features)
        Weight vector(s).

    intercept_ : float or array of shape (n_targets,)
        Independent term in decision function. Set to 0.0 if
        ``fit_intercept = False``.

    n_iter_ : None or array of shape (n_targets,)
        Actual number of iterations for each target. Available only for
        sag and lsqr solvers. Other solvers will return None.

        .. versionadded:: 0.17

    See also
    --------
    RidgeClassifier : Ridge classifier
    RidgeCV : Ridge regression with built-in cross validation
    :class:`sklearn.kernel_ridge.KernelRidge` : Kernel ridge regression
        combines ridge regression with the kernel trick

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

    def __init__(self, alpha=1.0, fit_intercept=True, normalize=False,
                 copy_X=True, max_iter=None, tol=1e-3, solver="auto",
                 random_state=None):
        super().__init__(
            alpha=alpha, fit_intercept=fit_intercept,
            normalize=normalize, copy_X=copy_X,
            max_iter=max_iter, tol=tol, solver=solver,
            random_state=random_state)

    def fit(self, X, y, sample_weight=None):
        """Fit Ridge regression model.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            Training data

        y : array-like of shape (n_samples,) or (n_samples, n_targets)
            Target values

        sample_weight : float or array-like of shape (n_samples,), default=None
            Individual weights for each sample. If given a float, every sample
            will have the same weight.

        Returns
        -------
        self : returns an instance of self.
        """
        return super().fit(X, y, sample_weight=sample_weight)


class RidgeClassifier(LinearClassifierMixin, _BaseRidge):
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
        Alpha corresponds to ``C^-1`` in other linear models such as
        LogisticRegression or LinearSVC.

    fit_intercept : bool, default=True
        Whether to calculate the intercept for this model. If set to false, no
        intercept will be used in calculations (e.g. data is expected to be
        already centered).

    normalize : bool, default=False
        This parameter is ignored when ``fit_intercept`` is set to False.
        If True, the regressors X will be normalized before regression by
        subtracting the mean and dividing by the l2-norm.
        If you wish to standardize, please use
        :class:`sklearn.preprocessing.StandardScaler` before calling ``fit``
        on an estimator with ``normalize=False``.

    copy_X : bool, default=True
        If True, X will be copied; else, it may be overwritten.

    max_iter : int, optional
        Maximum number of iterations for conjugate gradient solver.
        The default value is determined by scipy.sparse.linalg.

    tol : float, default=1e-3
        Precision of the solution.

    class_weight : dict or 'balanced', optional
        Weights associated with classes in the form ``{class_label: weight}``.
        If not given, all classes are supposed to have weight one.

        The "balanced" mode uses the values of y to automatically adjust
        weights inversely proportional to class frequencies in the input data
        as ``n_samples / (n_classes * np.bincount(y))``.

    solver : {'auto', 'svd', 'cholesky', 'lsqr', 'sparse_cg', 'sag', 'saga'}
        Solver to use in the computational routines:

        - 'auto' chooses the solver automatically based on the type of data.

        - 'svd' uses a Singular Value Decomposition of X to compute the Ridge
          coefficients. More stable for singular matrices than
          'cholesky'.

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

    random_state : int, RandomState instance or None, default=None
        The seed of the pseudo random number generator to use when shuffling
        the data.  If int, random_state is the seed used by the random number
        generator; If RandomState instance, random_state is the random number
        generator; If None, the random number generator is the RandomState
        instance used by `np.random`. Used when ``solver`` == 'sag'.

    Attributes
    ----------
    coef_ : array of shape (1, n_features) or (n_classes, n_features)
        Coefficient of the features in the decision function.

        ``coef_`` is of shape (1, n_features) when the given problem is binary.

    intercept_ : float or array of shape (n_targets,)
        Independent term in decision function. Set to 0.0 if
        ``fit_intercept = False``.

    n_iter_ : None or array of shape (n_targets,)
        Actual number of iterations for each target. Available only for
        sag and lsqr solvers. Other solvers will return None.

    classes_ : array of shape (n_classes,)
        The classes labels.

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

    def __init__(self, alpha=1.0, fit_intercept=True, normalize=False,
                 copy_X=True, max_iter=None, tol=1e-3, class_weight=None,
                 solver="auto", random_state=None):
        super().__init__(
            alpha=alpha, fit_intercept=fit_intercept, normalize=normalize,
            copy_X=copy_X, max_iter=max_iter, tol=tol, solver=solver,
            random_state=random_state)
        self.class_weight = class_weight

    def fit(self, X, y, sample_weight=None):
        """Fit Ridge classifier model.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            Training data.

        y : array-like of shape (n_samples,)
            Target values.

        sample_weight : float or array-like of shape (n_samples,), default=None
            Individual weights for each sample. If given a float, every sample
            will have the same weight.

            .. versionadded:: 0.17
               *sample_weight* support to Classifier.

        Returns
        -------
        self : object
            Instance of the estimator.
        """
        _accept_sparse = _get_valid_accept_sparse(sparse.issparse(X),
                                                  self.solver)
        X, y = check_X_y(X, y, accept_sparse=_accept_sparse, multi_output=True,
                         y_numeric=False)
        sample_weight = _check_sample_weight(sample_weight, X, dtype=X.dtype)

        self._label_binarizer = LabelBinarizer(pos_label=1, neg_label=-1)
        Y = self._label_binarizer.fit_transform(y)
        if not self._label_binarizer.y_type_.startswith('multilabel'):
            y = column_or_1d(y, warn=True)
        else:
            # we don't (yet) support multi-label classification in Ridge
            raise ValueError(
                "%s doesn't support multi-label classification" % (
                    self.__class__.__name__))

        if self.class_weight:
            # modify the sample weights with the corresponding class weight
            sample_weight = (sample_weight *
                             compute_sample_weight(self.class_weight, y))

        super().fit(X, Y, sample_weight=sample_weight)
        return self

    @property
    def classes_(self):
        return self._label_binarizer.classes_


def _check_gcv_mode(X, gcv_mode):
    possible_gcv_modes = [None, 'auto', 'svd', 'eigen']
    if gcv_mode not in possible_gcv_modes:
        raise ValueError(
            "Unknown value for 'gcv_mode'. "
            "Got {} instead of one of {}" .format(
                gcv_mode, possible_gcv_modes))
    if gcv_mode in ['eigen', 'svd']:
        return gcv_mode
    # if X has more rows than columns, use decomposition of X^T.X,
    # otherwise X.X^T
    if X.shape[0] > X.shape[1]:
        return 'svd'
    return 'eigen'


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
    abs_cosine = np.abs(query.dot(vectors))
    index = np.argmax(abs_cosine)
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
        return safe_sparse_dot(
            self.X, v[:-1], dense_output=True
        ) - self.sqrt_sw * self.X_mean.dot(v[:-1]) + v[-1] * self.sqrt_sw

    def _matmat(self, v):
        return (
            safe_sparse_dot(self.X, v[:-1], dense_output=True) -
            self.sqrt_sw[:, None] * self.X_mean.dot(v[:-1]) + v[-1] *
            self.sqrt_sw[:, None])

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
        res[:-1] = (
            safe_sparse_dot(self.X.T, v, dense_output=True) -
            (self.X_mean * self.sqrt_sw.dot(v))
        )
        res[-1] = np.dot(v, self.sqrt_sw)
        return res

    def _matmat(self, v):
        n_features = self.shape[0]
        res = np.empty((n_features, v.shape[1]), dtype=self.X.dtype)
        res[:-1] = (
            safe_sparse_dot(self.X.T, v, dense_output=True) -
            self.X_mean[:, None] * self.sqrt_sw.dot(v)
        )
        res[-1] = np.dot(self.sqrt_sw, v)
        return res


class _RidgeGCV(LinearModel):
    """Ridge regression with built-in Generalized Cross-Validation

    It allows efficient Leave-One-Out cross-validation.

    This class is not intended to be used directly. Use RidgeCV instead.

    Notes
    -----

    We want to solve (K + alpha*Id)c = y,
    where K = X X^T is the kernel matrix.

    Let G = (K + alpha*Id).

    Dual solution: c = G^-1y
    Primal solution: w = X^T c

    Compute eigendecomposition K = Q V Q^T.
    Then G^-1 = Q (V + alpha*Id)^-1 Q^T,
    where (V + alpha*Id) is diagonal.
    It is thus inexpensive to inverse for many alphas.

    Let loov be the vector of prediction values for each example
    when the model was fitted with all examples but this example.

    loov = (KG^-1Y - diag(KG^-1)Y) / diag(I-KG^-1)

    Let looe be the vector of prediction errors for each example
    when the model was fitted with all examples but this example.

    looe = y - loov = c / diag(G^-1)

    References
    ----------
    http://cbcl.mit.edu/publications/ps/MIT-CSAIL-TR-2007-025.pdf
    https://www.mit.edu/~9.520/spring07/Classes/rlsslides.pdf
    """

    def __init__(self, alphas=(0.1, 1.0, 10.0),
                 fit_intercept=True, normalize=False,
                 scoring=None, copy_X=True,
                 gcv_mode=None, store_cv_values=False):
        self.alphas = np.asarray(alphas)
        self.fit_intercept = fit_intercept
        self.normalize = normalize
        self.scoring = scoring
        self.copy_X = copy_X
        self.gcv_mode = gcv_mode
        self.store_cv_values = store_cv_values

    def _decomp_diag(self, v_prime, Q):
        # compute diagonal of the matrix: dot(Q, dot(diag(v_prime), Q^T))
        return (v_prime * Q ** 2).sum(axis=-1)

    def _diag_dot(self, D, B):
        # compute dot(diag(D), B)
        if len(B.shape) > 1:
            # handle case where B is > 1-d
            D = D[(slice(None), ) + (np.newaxis, ) * (len(B.shape) - 1)]
        return D * B

    def _compute_gram(self, X, sqrt_sw):
        """Computes the Gram matrix XX^T with possible centering.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            The preprocessed design matrix.

        sqrt_sw : ndarray of shape (n_samples,)
            square roots of sample weights

        Returns
        -------
        gram : ndarray of shape (n_samples, n_samples)
            The Gram matrix.
        X_mean : ndarray of shape (n_feature,)
            The weighted mean of ``X`` for each feature.

        Notes
        -----
        When X is dense the centering has been done in preprocessing
        so the mean is 0 and we just compute XX^T.

        When X is sparse it has not been centered in preprocessing, but it has
        been scaled by sqrt(sample weights).

        When self.fit_intercept is False no centering is done.

        The centered X is never actually computed because centering would break
        the sparsity of X.
        """
        center = self.fit_intercept and sparse.issparse(X)
        if not center:
            # in this case centering has been done in preprocessing
            # or we are not fitting an intercept.
            X_mean = np.zeros(X.shape[1], dtype=X.dtype)
            return safe_sparse_dot(X, X.T, dense_output=True), X_mean
        # X is sparse
        n_samples = X.shape[0]
        sample_weight_matrix = sparse.dia_matrix(
            (sqrt_sw, 0), shape=(n_samples, n_samples))
        X_weighted = sample_weight_matrix.dot(X)
        X_mean, _ = mean_variance_axis(X_weighted, axis=0)
        X_mean *= n_samples / sqrt_sw.dot(sqrt_sw)
        X_mX = sqrt_sw[:, None] * safe_sparse_dot(
            X_mean, X.T, dense_output=True)
        X_mX_m = np.outer(sqrt_sw, sqrt_sw) * np.dot(X_mean, X_mean)
        return (safe_sparse_dot(X, X.T, dense_output=True) + X_mX_m
                - X_mX - X_mX.T, X_mean)

    def _compute_covariance(self, X, sqrt_sw):
        """Computes covariance matrix X^TX with possible centering.

        Parameters
        ----------
        X : sparse matrix of shape (n_samples, n_features)
            The preprocessed design matrix.

        sqrt_sw : ndarray of shape (n_samples,)
            square roots of sample weights

        Returns
        -------
        covariance : ndarray of shape (n_features, n_features)
            The covariance matrix.
        X_mean : ndarray of shape (n_feature,)
            The weighted mean of ``X`` for each feature.

        Notes
        -----
        Since X is sparse it has not been centered in preprocessing, but it has
        been scaled by sqrt(sample weights).

        When self.fit_intercept is False no centering is done.

        The centered X is never actually computed because centering would break
        the sparsity of X.
        """
        if not self.fit_intercept:
            # in this case centering has been done in preprocessing
            # or we are not fitting an intercept.
            X_mean = np.zeros(X.shape[1], dtype=X.dtype)
            return safe_sparse_dot(X.T, X, dense_output=True), X_mean
        # this function only gets called for sparse X
        n_samples = X.shape[0]
        sample_weight_matrix = sparse.dia_matrix(
            (sqrt_sw, 0), shape=(n_samples, n_samples))
        X_weighted = sample_weight_matrix.dot(X)
        X_mean, _ = mean_variance_axis(X_weighted, axis=0)
        X_mean = X_mean * n_samples / sqrt_sw.dot(sqrt_sw)
        weight_sum = sqrt_sw.dot(sqrt_sw)
        return (safe_sparse_dot(X.T, X, dense_output=True) -
                weight_sum * np.outer(X_mean, X_mean),
                X_mean)

    def _sparse_multidot_diag(self, X, A, X_mean, sqrt_sw):
        """Compute the diagonal of (X - X_mean).dot(A).dot((X - X_mean).T)
        without explicitely centering X nor computing X.dot(A)
        when X is sparse.

        Parameters
        ----------
        X : sparse matrix of shape (n_samples, n_features)

        A : ndarray of shape (n_features, n_features)

        X_mean : ndarray of shape (n_features,)

        sqrt_sw : ndarray of shape (n_features,)
            square roots of sample weights

        Returns
        -------
        diag : np.ndarray, shape (n_samples,)
            The computed diagonal.
        """
        intercept_col = scale = sqrt_sw
        batch_size = X.shape[1]
        diag = np.empty(X.shape[0], dtype=X.dtype)
        for start in range(0, X.shape[0], batch_size):
            batch = slice(start, min(X.shape[0], start + batch_size), 1)
            X_batch = np.empty(
                (X[batch].shape[0], X.shape[1] + self.fit_intercept),
                dtype=X.dtype
            )
            if self.fit_intercept:
                X_batch[:, :-1] = X[batch].A - X_mean * scale[batch][:, None]
                X_batch[:, -1] = intercept_col[batch]
            else:
                X_batch = X[batch].A
            diag[batch] = (X_batch.dot(A) * X_batch).sum(axis=1)
        return diag

    def _eigen_decompose_gram(self, X, y, sqrt_sw):
        """Eigendecomposition of X.X^T, used when n_samples <= n_features."""
        # if X is dense it has already been centered in preprocessing
        K, X_mean = self._compute_gram(X, sqrt_sw)
        if self.fit_intercept:
            # to emulate centering X with sample weights,
            # ie removing the weighted average, we add a column
            # containing the square roots of the sample weights.
            # by centering, it is orthogonal to the other columns
            K += np.outer(sqrt_sw, sqrt_sw)
        eigvals, Q = linalg.eigh(K)
        QT_y = np.dot(Q.T, y)
        return X_mean, eigvals, Q, QT_y

    def _solve_eigen_gram(self, alpha, y, sqrt_sw, X_mean, eigvals, Q, QT_y):
        """Compute dual coefficients and diagonal of G^-1.

        Used when we have a decomposition of X.X^T (n_samples <= n_features).
        """
        w = 1. / (eigvals + alpha)
        if self.fit_intercept:
            # the vector containing the square roots of the sample weights (1
            # when no sample weights) is the eigenvector of XX^T which
            # corresponds to the intercept; we cancel the regularization on
            # this dimension. the corresponding eigenvalue is
            # sum(sample_weight).
            normalized_sw = sqrt_sw / np.linalg.norm(sqrt_sw)
            intercept_dim = _find_smallest_angle(normalized_sw, Q)
            w[intercept_dim] = 0  # cancel regularization for the intercept

        c = np.dot(Q, self._diag_dot(w, QT_y))
        G_inverse_diag = self._decomp_diag(w, Q)
        # handle case where y is 2-d
        if len(y.shape) != 1:
            G_inverse_diag = G_inverse_diag[:, np.newaxis]
        return G_inverse_diag, c

    def _eigen_decompose_covariance(self, X, y, sqrt_sw):
        """Eigendecomposition of X^T.X, used when n_samples > n_features
        and X is sparse.
        """
        n_samples, n_features = X.shape
        cov = np.empty((n_features + 1, n_features + 1), dtype=X.dtype)
        cov[:-1, :-1], X_mean = self._compute_covariance(X, sqrt_sw)
        if not self.fit_intercept:
            cov = cov[:-1, :-1]
        # to emulate centering X with sample weights,
        # ie removing the weighted average, we add a column
        # containing the square roots of the sample weights.
        # by centering, it is orthogonal to the other columns
        # when all samples have the same weight we add a column of 1
        else:
            cov[-1] = 0
            cov[:, -1] = 0
            cov[-1, -1] = sqrt_sw.dot(sqrt_sw)
        nullspace_dim = max(0, X.shape[1] - X.shape[0])
        eigvals, V = linalg.eigh(cov)
        # remove eigenvalues and vectors in the null space of X^T.X
        eigvals = eigvals[nullspace_dim:]
        V = V[:, nullspace_dim:]
        return X_mean, eigvals, V, X

    def _solve_eigen_covariance_no_intercept(
            self, alpha, y, sqrt_sw, X_mean, eigvals, V, X):
        """Compute dual coefficients and diagonal of G^-1.

        Used when we have a decomposition of X^T.X
        (n_samples > n_features and X is sparse), and not fitting an intercept.
        """
        w = 1 / (eigvals + alpha)
        A = (V * w).dot(V.T)
        AXy = A.dot(safe_sparse_dot(X.T, y, dense_output=True))
        y_hat = safe_sparse_dot(X, AXy, dense_output=True)
        hat_diag = self._sparse_multidot_diag(X, A, X_mean, sqrt_sw)
        if len(y.shape) != 1:
            # handle case where y is 2-d
            hat_diag = hat_diag[:, np.newaxis]
        return (1 - hat_diag) / alpha, (y - y_hat) / alpha

    def _solve_eigen_covariance_intercept(
            self, alpha, y, sqrt_sw, X_mean, eigvals, V, X):
        """Compute dual coefficients and diagonal of G^-1.

        Used when we have a decomposition of X^T.X
        (n_samples > n_features and X is sparse),
        and we are fitting an intercept.
        """
        # the vector [0, 0, ..., 0, 1]
        # is the eigenvector of X^TX which
        # corresponds to the intercept; we cancel the regularization on
        # this dimension. the corresponding eigenvalue is
        # sum(sample_weight), e.g. n when uniform sample weights.
        intercept_sv = np.zeros(V.shape[0])
        intercept_sv[-1] = 1
        intercept_dim = _find_smallest_angle(intercept_sv, V)
        w = 1 / (eigvals + alpha)
        w[intercept_dim] = 1 / eigvals[intercept_dim]
        A = (V * w).dot(V.T)
        # add a column to X containing the square roots of sample weights
        X_op = _X_CenterStackOp(X, X_mean, sqrt_sw)
        AXy = A.dot(X_op.T.dot(y))
        y_hat = X_op.dot(AXy)
        hat_diag = self._sparse_multidot_diag(X, A, X_mean, sqrt_sw)
        # return (1 - hat_diag), (y - y_hat)
        if len(y.shape) != 1:
            # handle case where y is 2-d
            hat_diag = hat_diag[:, np.newaxis]
        return (1 - hat_diag) / alpha, (y - y_hat) / alpha

    def _solve_eigen_covariance(
            self, alpha, y, sqrt_sw, X_mean, eigvals, V, X):
        """Compute dual coefficients and diagonal of G^-1.

        Used when we have a decomposition of X^T.X
        (n_samples > n_features and X is sparse).
        """
        if self.fit_intercept:
            return self._solve_eigen_covariance_intercept(
                alpha, y, sqrt_sw, X_mean, eigvals, V, X)
        return self._solve_eigen_covariance_no_intercept(
            alpha, y, sqrt_sw, X_mean, eigvals, V, X)

    def _svd_decompose_design_matrix(self, X, y, sqrt_sw):
        # X already centered
        X_mean = np.zeros(X.shape[1], dtype=X.dtype)
        if self.fit_intercept:
            # to emulate fit_intercept=True situation, add a column
            # containing the square roots of the sample weights
            # by centering, the other columns are orthogonal to that one
            intercept_column = sqrt_sw[:, None]
            X = np.hstack((X, intercept_column))
        U, singvals, _ = linalg.svd(X, full_matrices=0)
        singvals_sq = singvals ** 2
        UT_y = np.dot(U.T, y)
        return X_mean, singvals_sq, U, UT_y

    def _solve_svd_design_matrix(
            self, alpha, y, sqrt_sw, X_mean, singvals_sq, U, UT_y):
        """Compute dual coefficients and diagonal of G^-1.

        Used when we have an SVD decomposition of X
        (n_samples > n_features and X is dense).
        """
        w = ((singvals_sq + alpha) ** -1) - (alpha ** -1)
        if self.fit_intercept:
            # detect intercept column
            normalized_sw = sqrt_sw / np.linalg.norm(sqrt_sw)
            intercept_dim = _find_smallest_angle(normalized_sw, U)
            # cancel the regularization for the intercept
            w[intercept_dim] = - (alpha ** -1)
        c = np.dot(U, self._diag_dot(w, UT_y)) + (alpha ** -1) * y
        G_inverse_diag = self._decomp_diag(w, U) + (alpha ** -1)
        if len(y.shape) != 1:
            # handle case where y is 2-d
            G_inverse_diag = G_inverse_diag[:, np.newaxis]
        return G_inverse_diag, c

    def fit(self, X, y, sample_weight=None):
        """Fit Ridge regression model with gcv.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            Training data. Will be cast to float64 if necessary.

        y : array-like of shape (n_samples,) or (n_samples, n_targets)
            Target values. Will be cast to float64 if necessary.

        sample_weight : float or array-like of shape (n_samples,), default=None
            Individual weights for each sample. If given a float, every sample
            will have the same weight.

        Returns
        -------
        self : object
        """
        X, y = check_X_y(X, y, ['csr', 'csc', 'coo'], dtype=[np.float64],
                         multi_output=True, y_numeric=True)

        if sample_weight is not None:
            sample_weight = _check_sample_weight(sample_weight, X,
                                                 dtype=X.dtype)

        if np.any(self.alphas <= 0):
            raise ValueError(
                "alphas must be positive. Got {} containing some "
                "negative or null value instead.".format(self.alphas))

        n_samples, n_features = X.shape

        X, y, X_offset, y_offset, X_scale = LinearModel._preprocess_data(
            X, y, self.fit_intercept, self.normalize, self.copy_X,
            sample_weight=sample_weight)

        gcv_mode = _check_gcv_mode(X, self.gcv_mode)

        if gcv_mode == 'eigen':
            decompose = self._eigen_decompose_gram
            solve = self._solve_eigen_gram
        elif gcv_mode == 'svd':
            if sparse.issparse(X):
                decompose = self._eigen_decompose_covariance
                solve = self._solve_eigen_covariance
            else:
                decompose = self._svd_decompose_design_matrix
                solve = self._solve_svd_design_matrix

        if sample_weight is not None:
            X, y = _rescale_data(X, y, sample_weight)
            sqrt_sw = np.sqrt(sample_weight)
        else:
            sqrt_sw = np.ones(X.shape[0], dtype=X.dtype)

        scorer = check_scoring(self, scoring=self.scoring, allow_none=True)
        error = scorer is None

        n_y = 1 if len(y.shape) == 1 else y.shape[1]
        cv_values = np.zeros((n_samples * n_y, len(self.alphas)),
                             dtype=X.dtype)
        C = []
        X_mean, *decomposition = decompose(X, y, sqrt_sw)
        for i, alpha in enumerate(self.alphas):
            G_inverse_diag, c = solve(
                float(alpha), y, sqrt_sw, X_mean, *decomposition)
            if error:
                squared_errors = (c / G_inverse_diag) ** 2
                cv_values[:, i] = squared_errors.ravel()
            else:
                predictions = y - (c / G_inverse_diag)
                cv_values[:, i] = predictions.ravel()
            C.append(c)

        if error:
            best = cv_values.mean(axis=0).argmin()
        else:
            # The scorer want an object that will make the predictions but
            # they are already computed efficiently by _RidgeGCV. This
            # identity_estimator will just return them
            def identity_estimator():
                pass
            identity_estimator.decision_function = lambda y_predict: y_predict
            identity_estimator.predict = lambda y_predict: y_predict

            # signature of scorer is (estimator, X, y)
            out = [scorer(identity_estimator, cv_values[:, i], y.ravel())
                   for i in range(len(self.alphas))]
            best = np.argmax(out)

        self.alpha_ = self.alphas[best]
        self.dual_coef_ = C[best]
        self.coef_ = safe_sparse_dot(self.dual_coef_.T, X)

        X_offset += X_mean * X_scale
        self._set_intercept(X_offset, y_offset, X_scale)

        if self.store_cv_values:
            if len(y.shape) == 1:
                cv_values_shape = n_samples, len(self.alphas)
            else:
                cv_values_shape = n_samples, n_y, len(self.alphas)
            self.cv_values_ = cv_values.reshape(cv_values_shape)

        return self


class _BaseRidgeCV(LinearModel):
    def __init__(self, alphas=(0.1, 1.0, 10.0),
                 fit_intercept=True, normalize=False, scoring=None,
                 cv=None, gcv_mode=None,
                 store_cv_values=False):
        self.alphas = np.asarray(alphas)
        self.fit_intercept = fit_intercept
        self.normalize = normalize
        self.scoring = scoring
        self.cv = cv
        self.gcv_mode = gcv_mode
        self.store_cv_values = store_cv_values

    def fit(self, X, y, sample_weight=None):
        """Fit Ridge regression model with cv.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            Training data. If using GCV, will be cast to float64
            if necessary.

        y : array-like of shape (n_samples,) or (n_samples, n_targets)
            Target values. Will be cast to X's dtype if necessary.

        sample_weight : float or array-like of shape (n_samples,), default=None
            Individual weights for each sample. If given a float, every sample
            will have the same weight.

        Returns
        -------
        self : object

        Notes
        -----
        When sample_weight is provided, the selected hyperparameter may depend
        on whether we use generalized cross-validation (cv=None or cv='auto')
        or another form of cross-validation, because only generalized
        cross-validation takes the sample weights into account when computing
        the validation score.
        """
        cv = self.cv
        if cv is None:
            estimator = _RidgeGCV(self.alphas,
                                  fit_intercept=self.fit_intercept,
                                  normalize=self.normalize,
                                  scoring=self.scoring,
                                  gcv_mode=self.gcv_mode,
                                  store_cv_values=self.store_cv_values)
            estimator.fit(X, y, sample_weight=sample_weight)
            self.alpha_ = estimator.alpha_
            if self.store_cv_values:
                self.cv_values_ = estimator.cv_values_
        else:
            if self.store_cv_values:
                raise ValueError("cv!=None and store_cv_values=True "
                                 " are incompatible")
            parameters = {'alpha': self.alphas}
            solver = 'sparse_cg' if sparse.issparse(X) else 'auto'
            gs = GridSearchCV(Ridge(fit_intercept=self.fit_intercept,
                                    normalize=self.normalize,
                                    solver=solver),
                              parameters, cv=cv, scoring=self.scoring)
            gs.fit(X, y, sample_weight=sample_weight)
            estimator = gs.best_estimator_
            self.alpha_ = gs.best_estimator_.alpha

        self.coef_ = estimator.coef_
        self.intercept_ = estimator.intercept_

        return self


class RidgeCV(MultiOutputMixin, RegressorMixin, _BaseRidgeCV):
    """Ridge regression with built-in cross-validation.

    See glossary entry for :term:`cross-validation estimator`.

    By default, it performs Generalized Cross-Validation, which is a form of
    efficient Leave-One-Out cross-validation.

    Read more in the :ref:`User Guide <ridge_regression>`.

    Parameters
    ----------
    alphas : ndarray of shape (n_alphas,), default=(0.1, 1.0, 10.0)
        Array of alpha values to try.
        Regularization strength; must be a positive float. Regularization
        improves the conditioning of the problem and reduces the variance of
        the estimates. Larger values specify stronger regularization.
        Alpha corresponds to ``C^-1`` in other linear models such as
        LogisticRegression or LinearSVC.
        If using generalized cross-validation, alphas must be positive.

    fit_intercept : bool, default=True
        Whether to calculate the intercept for this model. If set
        to false, no intercept will be used in calculations
        (i.e. data is expected to be centered).

    normalize : bool, default=False
        This parameter is ignored when ``fit_intercept`` is set to False.
        If True, the regressors X will be normalized before regression by
        subtracting the mean and dividing by the l2-norm.
        If you wish to standardize, please use
        :class:`sklearn.preprocessing.StandardScaler` before calling ``fit``
        on an estimator with ``normalize=False``.

    scoring : string, callable or None, default=None
        A string (see model evaluation documentation) or
        a scorer callable object / function with signature
        ``scorer(estimator, X, y)``.
        If None, the negative mean squared error if cv is 'auto' or None
        (i.e. when using generalized cross-validation), and r2 score otherwise.

    cv : int, cross-validation generator or an iterable, optional
        Determines the cross-validation splitting strategy.
        Possible inputs for cv are:

        - None, to use the efficient Leave-One-Out cross-validation
          (also known as Generalized Cross-Validation).
        - integer, to specify the number of folds.
        - :term:`CV splitter`,
        - An iterable yielding (train, test) splits as arrays of indices.

        For integer/None inputs, if ``y`` is binary or multiclass,
        :class:`sklearn.model_selection.StratifiedKFold` is used, else,
        :class:`sklearn.model_selection.KFold` is used.

        Refer :ref:`User Guide <cross_validation>` for the various
        cross-validation strategies that can be used here.

    gcv_mode : {None, 'auto', 'svd', eigen'}, optional
        Flag indicating which strategy to use when performing
        Generalized Cross-Validation. Options are::

            'auto' : use 'svd' if n_samples > n_features, otherwise use 'eigen'
            'svd' : force use of singular value decomposition of X when X is
                dense, eigenvalue decomposition of X^T.X when X is sparse.
            'eigen' : force computation via eigendecomposition of X.X^T

        The 'auto' mode is the default and is intended to pick the cheaper
        option of the two depending on the shape of the training data.

    store_cv_values : boolean, default=False
        Flag indicating if the cross-validation values corresponding to
        each alpha should be stored in the ``cv_values_`` attribute (see
        below). This flag is only compatible with ``cv=None`` (i.e. using
        Generalized Cross-Validation).

    Attributes
    ----------
    cv_values_ : array of shape (n_samples, n_alphas) or \
        shape (n_samples, n_targets, n_alphas), optional
        Cross-validation values for each alpha (if ``store_cv_values=True``\
        and ``cv=None``). After ``fit()`` has been called, this attribute \
        will contain the mean squared errors (by default) or the values \
        of the ``{loss,score}_func`` function (if provided in the constructor).

    coef_ : array of shape (n_features) or (n_targets, n_features)
        Weight vector(s).

    intercept_ : float or array of shape (n_targets,)
        Independent term in decision function. Set to 0.0 if
        ``fit_intercept = False``.

    alpha_ : float
        Estimated regularization parameter.

    Examples
    --------
    >>> from sklearn.datasets import load_diabetes
    >>> from sklearn.linear_model import RidgeCV
    >>> X, y = load_diabetes(return_X_y=True)
    >>> clf = RidgeCV(alphas=[1e-3, 1e-2, 1e-1, 1]).fit(X, y)
    >>> clf.score(X, y)
    0.5166...

    See also
    --------
    Ridge : Ridge regression
    RidgeClassifier : Ridge classifier
    RidgeClassifierCV : Ridge classifier with built-in cross validation
    """
    pass


class RidgeClassifierCV(LinearClassifierMixin, _BaseRidgeCV):
    """Ridge classifier with built-in cross-validation.

    See glossary entry for :term:`cross-validation estimator`.

    By default, it performs Generalized Cross-Validation, which is a form of
    efficient Leave-One-Out cross-validation. Currently, only the n_features >
    n_samples case is handled efficiently.

    Read more in the :ref:`User Guide <ridge_regression>`.

    Parameters
    ----------
    alphas : ndarray of shape (n_alphas,), default=(0.1, 1.0, 10.0)
        Array of alpha values to try.
        Regularization strength; must be a positive float. Regularization
        improves the conditioning of the problem and reduces the variance of
        the estimates. Larger values specify stronger regularization.
        Alpha corresponds to ``C^-1`` in other linear models such as
        LogisticRegression or LinearSVC.

    fit_intercept : bool, default=True
        Whether to calculate the intercept for this model. If set
        to false, no intercept will be used in calculations
        (i.e. data is expected to be centered).

    normalize : bool, default=False
        This parameter is ignored when ``fit_intercept`` is set to False.
        If True, the regressors X will be normalized before regression by
        subtracting the mean and dividing by the l2-norm.
        If you wish to standardize, please use
        :class:`sklearn.preprocessing.StandardScaler` before calling ``fit``
        on an estimator with ``normalize=False``.

    scoring : string, callable or None, default=None
        A string (see model evaluation documentation) or
        a scorer callable object / function with signature
        ``scorer(estimator, X, y)``.

    cv : int, cross-validation generator or an iterable, optional
        Determines the cross-validation splitting strategy.
        Possible inputs for cv are:

        - None, to use the efficient Leave-One-Out cross-validation
        - integer, to specify the number of folds.
        - :term:`CV splitter`,
        - An iterable yielding (train, test) splits as arrays of indices.

        Refer :ref:`User Guide <cross_validation>` for the various
        cross-validation strategies that can be used here.

    class_weight : dict or 'balanced', optional
        Weights associated with classes in the form ``{class_label: weight}``.
        If not given, all classes are supposed to have weight one.

        The "balanced" mode uses the values of y to automatically adjust
        weights inversely proportional to class frequencies in the input data
        as ``n_samples / (n_classes * np.bincount(y))``

    store_cv_values : boolean, default=False
        Flag indicating if the cross-validation values corresponding to
        each alpha should be stored in the ``cv_values_`` attribute (see
        below). This flag is only compatible with ``cv=None`` (i.e. using
        Generalized Cross-Validation).

    Attributes
    ----------
    cv_values_ : array of shape (n_samples, n_targets, n_alphas), optional
        Cross-validation values for each alpha (if ``store_cv_values=True`` and
        ``cv=None``). After ``fit()`` has been called, this attribute will
        contain the mean squared errors (by default) or the values of the
        ``{loss,score}_func`` function (if provided in the constructor). This
        attribute exists only when ``store_cv_values`` is True.

    coef_ : array of shape (1, n_features) or (n_targets, n_features)
        Coefficient of the features in the decision function.

        ``coef_`` is of shape (1, n_features) when the given problem is binary.

    intercept_ : float or array of shape (n_targets,)
        Independent term in decision function. Set to 0.0 if
        ``fit_intercept = False``.

    alpha_ : float
        Estimated regularization parameter

    classes_ : array of shape (n_classes,)
        The classes labels.

    Examples
    --------
    >>> from sklearn.datasets import load_breast_cancer
    >>> from sklearn.linear_model import RidgeClassifierCV
    >>> X, y = load_breast_cancer(return_X_y=True)
    >>> clf = RidgeClassifierCV(alphas=[1e-3, 1e-2, 1e-1, 1]).fit(X, y)
    >>> clf.score(X, y)
    0.9630...

    See also
    --------
    Ridge : Ridge regression
    RidgeClassifier : Ridge classifier
    RidgeCV : Ridge regression with built-in cross validation

    Notes
    -----
    For multi-class classification, n_class classifiers are trained in
    a one-versus-all approach. Concretely, this is implemented by taking
    advantage of the multi-variate response support in Ridge.
    """

    def __init__(self, alphas=(0.1, 1.0, 10.0), fit_intercept=True,
                 normalize=False, scoring=None, cv=None, class_weight=None,
                 store_cv_values=False):
        super().__init__(
            alphas=alphas, fit_intercept=fit_intercept, normalize=normalize,
            scoring=scoring, cv=cv, store_cv_values=store_cv_values)
        self.class_weight = class_weight

    def fit(self, X, y, sample_weight=None):
        """Fit Ridge classifier with cv.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            Training vectors, where n_samples is the number of samples
            and n_features is the number of features. When using GCV,
            will be cast to float64 if necessary.

        y : array-like of shape (n_samples,)
            Target values. Will be cast to X's dtype if necessary.

        sample_weight : float or array-like of shape (n_samples,), default=None
            Individual weights for each sample. If given a float, every sample
            will have the same weight.

        Returns
        -------
        self : object
        """
        X, y = check_X_y(X, y, accept_sparse=['csr', 'csc', 'coo'],
                         multi_output=True, y_numeric=False)
        sample_weight = _check_sample_weight(sample_weight, X, dtype=X.dtype)

        self._label_binarizer = LabelBinarizer(pos_label=1, neg_label=-1)
        Y = self._label_binarizer.fit_transform(y)
        if not self._label_binarizer.y_type_.startswith('multilabel'):
            y = column_or_1d(y, warn=True)

        if self.class_weight:
            # modify the sample weights with the corresponding class weight
            sample_weight = (sample_weight *
                             compute_sample_weight(self.class_weight, y))

        _BaseRidgeCV.fit(self, X, Y, sample_weight=sample_weight)
        return self

    @property
    def classes_(self):
        return self._label_binarizer.classes_
