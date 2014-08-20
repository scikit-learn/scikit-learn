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

from .base import LinearClassifierMixin, LinearModel
from ..base import RegressorMixin
from ..utils.extmath import safe_sparse_dot
from ..utils import check_X_y
from ..utils import compute_class_weight
from ..utils import column_or_1d
from ..preprocessing import LabelBinarizer
from ..grid_search import GridSearchCV
from ..externals import six
from ..metrics.scorer import check_scoring
from ..cross_validation import check_cv


def _ridge_path_svd(X_train, Y_train, alphas,
                    X_test=None, sg_val_thresh=1e-15):

    U, S, VT = linalg.svd(X_train, full_matrices=False)
    nnz_max = np.where(S > sg_val_thresh)[0].max() + 1
    UTY = U.T[:nnz_max].dot(Y_train)
    s_alpha = S[np.newaxis, :nnz_max, np.newaxis] / (
        S[np.newaxis, :nnz_max, np.newaxis] ** 2 + alphas[:, np.newaxis])
    s_alpha_UTY = s_alpha * UTY[np.newaxis]
    if X_test is not None:
        dekernelize = X_test.dot(VT[:nnz_max].T)
    else:
        dekernelize = VT[:nnz_max].T
    output = dekernelize.dot(s_alpha_UTY).transpose(1, 0, 2)
    return output


def _ridge_gcv_path_svd(X_train, Y_train, alphas,
                        sample_weight=None,
                        mode='looe',
                        sg_val_thresh=1e-15,
                        copy=True):

    if sample_weight is not None:
        has_sw = True
        sqrt_sw = np.sqrt(sample_weight).ravel()[:, np.newaxis]
        if copy:
            X_train = X_train * sqrt_sw
            Y_train = Y_train * sqrt_sw
        else:
            X_train *= sqrt_sw
            Y_train *= sqrt_sw
    else:
        has_sw = False

    U, S, VT = linalg.svd(X_train, full_matrices=False)
    del VT
    nnz_max = np.where(S > sg_val_thresh)[0].max() + 1
    U = U[:, :nnz_max]
    S = S[:nnz_max]

    inv_S_and_alpha = (1. / (S[np.newaxis, :, np.newaxis] ** 2 +
                             alphas[:, np.newaxis, :]) -
                       1. / alphas[:, np.newaxis, :])

    UTY = U.T.dot(Y_train)
    numerator = (
        U.dot(inv_S_and_alpha * UTY[np.newaxis]).transpose(1, 0, 2) +
        Y_train[np.newaxis] / alphas[:, np.newaxis, :])

    denominator = ((U ** 2).dot(inv_S_and_alpha).transpose(1, 0, 2)
                   + 1. / alphas[:, np.newaxis, :])

    looe = numerator / denominator

    if has_sw:
        looe /= sqrt_sw[np.newaxis]
        numerator *= sqrt_sw[np.newaxis]

    if sample_weight is not None:
        if copy:
            Y_train = Y_train / sqrt_sw
        else:
            Y_train /= sqrt_sw
            X_train /= sqrt_sw

    if mode == 'loov':
        return Y_train[np.newaxis] - looe, numerator
    elif mode == 'looe':
        return looe, numerator
    else:
        raise ValueError("mode must be in ['looe', 'loov']. Got %s" % mode)


def _precomp_kernel_ridge_path_eigen(gramX_train, Y_train, alphas,
                                     gramX_test=None,
                                     sample_weight=None,
                                     mode='normal',
                                     eig_val_thresh=1e-15,
                                     copy=True):
    """Precomputed kernel ridge without intercept.

    Parameters
    ----------

    gramX_train: ndarray, shape=(n_samples, n_samples)
        Kernel Gram matrix of the training data

    Y_train: ndarray, shape=(n_samples, n_targets)
        Target variables for regression

    alphas: ndarray, shape=(n_penalties, n_targets) or (n_penalties, 1)
        L2 penalties for ridge regression, several per target possible

    gramX_test: ndarray, shape=(n_test_samples, n_train_samples)
        Kernel similarities of test samples with the train samples

    sample_weight: ndarray or None, shape=(n_samples,), default None
        Weights for each sample

    mode: str, {'normal', 'looe', 'loov'}, default 'normal'
        Indicates the desired output.

        - 'normal' yields ridge coefficients or predictions
        - 'looe' yields a vector of leave one out residuals on the train
          set. A warning is raised if gramX_test is passed, because this is
          not used here.
        - 'loov' yields a vector of leave one out predictions

    eig_val_thresh: float, default 1e-15
        Threshold for eigenvalues of the Gram matrix. Any eigenvalue under
        eig_val_thresh will be set to zero. This is done to avoid explosion
        of the inverse at very low penalties

    """
    if gramX_test is not None and mode in ['looe', 'loov']:
        raise ValueError("gramX_test provided in loo mode. Remove X_test "
                         "or change mode.")

    if sample_weight is not None:
        has_sw = True
        sqrt_sw = np.sqrt(sample_weight.reshape(-1, 1))
        if copy:
            gramX_train = (sqrt_sw * gramX_train) * sqrt_sw.T
            Y_train = sqrt_sw * Y_train
        else:
            gramX_train *= sqrt_sw
            gramX_train *= sqrt_sw.T
            Y_train *= sqrt_sw
    else:
        has_sw = False

    eig_val_train, eig_vect_train = linalg.eigh(gramX_train)
    VTY = eig_vect_train.T.dot(Y_train)
    v_alpha = (eig_val_train[np.newaxis, :, np.newaxis] +
               alphas[:, np.newaxis])
    v_alpha_inv = np.zeros_like(v_alpha)
    passes_threshold = (v_alpha > eig_val_thresh)
    v_alpha_inv[passes_threshold] = 1. / v_alpha[passes_threshold]
    v_alpha_VTY = v_alpha_inv * VTY[np.newaxis]

    if gramX_test is None:
        dual_coef = eig_vect_train.dot(v_alpha_VTY).transpose(1, 0, 2)
        if has_sw:
            dual_coef *= sqrt_sw
        if mode == 'normal':
            output = dual_coef
        elif mode in ['looe', 'loov']:
            diag_rescale = (eig_vect_train ** 2).dot(
                v_alpha_inv).transpose(1, 0, 2)
            looe = dual_coef / diag_rescale
            if has_sw:
                looe /= sqrt_sw ** 2
            output = looe, dual_coef
            if mode == 'loov':
                loov = Y_train[np.newaxis] - looe
                output = loov, dual_coef
        else:
            raise ValueError("mode not understood")
    else:
        predict_mat = gramX_test.dot(eig_vect_train)
        if has_sw:
            predict_mat *= sqrt_sw
        predictions = predict_mat.dot(v_alpha_VTY).transpose(1, 0, 2)
        output = predictions

    return output


def _linear_kernel(X, Y):
    return safe_sparse_dot(X, Y.T, dense_output=True)


def _kernel_ridge_path_eigen(X_train, Y_train, alphas, X_test=None,
                             sample_weight=None,
                             kernel=_linear_kernel,
                             mode='normal'):
    gramX_train = kernel(X_train, X_train)
    gramX_test = None
    if X_test is not None:
        gramX_test = kernel(X_test, X_train)

    return _precomp_kernel_ridge_path_eigen(
        gramX_train, Y_train, alphas,
        gramX_test, sample_weight=sample_weight,
        mode=mode, copy=True)


def _feature_ridge_path_eigen(X_train, Y_train, alphas, X_test=None,
                              eig_val_thresh=1e-15):

    XTY = X_train.T.dot(Y_train)
    eig_val_train, eig_vect_train = linalg.eigh(X_train.T.dot(X_train))
    VTXTY = eig_vect_train.T.dot(XTY)  # under some circumstances it may be
                                       # better to do (V.T.dot(X.T)).dot(Y)

    v_alpha = (eig_val_train[np.newaxis, :, np.newaxis] +
                        alphas[:, np.newaxis])
    v_alpha_inv = np.zeros_like(v_alpha)
    passes_threshold = v_alpha > eig_val_thresh
    v_alpha_inv[passes_threshold] = 1. / v_alpha[passes_threshold]
    coef = eig_vect_train.dot(
        v_alpha_inv * VTXTY[np.newaxis]).transpose(1, 0, 2)

    if X_test is None:
        return coef
    else:
        return X_test.dot(coef).transpose(1, 0, 2)


def _multi_r2(y_true, y_pred):
    err = ((y_pred - y_true) ** 2).sum(-2)
    sq_norm_y_true = (y_true ** 2).sum(-2)
    return 1 - err / (sq_norm_y_true + 1e-18)


def ridge_path(X_train, Y_train, alphas, X_test=None, solver="eigen"):
    """Perform Ridge regression along a regularization path

    Parameters
    ----------

    X_train : ndarray, shape = (n_samples, n_features)
        Training data

    Y_train : ndarray, shape = (n_samples, n_targets) or (n_samples,)
        Training targets

    alphas : ndarray, shape=(n_penalties, n_targets) or (n_penalties, 1)
        Penalties on regularization path. Can be global or per target

    X_test : ndarray, shape=(n_test_samples, n_features), optional
        Test set. If specified, ridge_path returns predictions on X_test,
        otherwise it returns coefficients. Recommended if only predictions
        are of interest, especially in the n_samples << n_features case,
        since it bypasses the calculation of the feature coefficients,
        which can result in a memory and speed gain if n_features is large.

    solver : str, {'eigen', 'svd'}, (default 'eigen')
        The solver to use for ridge_path.


    Returns
    -------
    coef or predictions : ndarray
        shape = (n_penalties, n_features, n_targets) (coef) or
                (n_penalties, n_features) (coef if only one target) or
                (n_penalties, n_samples, n_targets) (predictions)
        If X_test is provided, predictions are returned, otherwise coef.

    Notes
    -----
    This solver does not fit the intercept.
"""

    y_raveled = Y_train.ndim == 1
    Y_train = np.atleast_2d(Y_train.T).T

    n_samples, n_features = X_train.shape
    n_samples_, n_targets = Y_train.shape

    if n_samples != n_samples_:
        raise ValueError(
            ("Number of samples in X_train (%d) and Y_train (%d) must"
            " correspond.") % (n_samples, n_samples_))

    alphas = np.atleast_1d(alphas)
    if alphas.ndim == 1:
        if len(alphas) == n_targets:
            raise ValueError(
                ("You have specified as many penalties as targets (%d). If"
                 " you want these to act as individual penalties, please"
                 " shape them as (1, n_targets) == (1, %d). If you wish"
                 " every penalty to be used on every target, please shape"
                 " alphas as (n_penalties, 1) == (%d, 1)") %
                (n_targets, n_targets, n_targets))
        else:
            alphas = alphas[:, np.newaxis]

    if solver not in ['eigen', 'svd']:
        raise NotImplementedError("Solver %s not implemented" % solver)

    if solver == 'eigen':
        if n_samples < n_features:
            # A kernelized implementation is more efficient, because it
            # work in the dual space (the sample axis).
            dual_path_or_predictions = _kernel_ridge_path_eigen(
                X_train, Y_train, alphas, X_test, kernel=_linear_kernel)
            if X_test is None:
                primal_path = X_train.T.dot(
                    dual_path_or_predictions).transpose(1, 0, 2)
                path = primal_path
            else:
                prediction_path = dual_path_or_predictions
                path = prediction_path
        else:
            path = _feature_ridge_path_eigen(X_train, Y_train, alphas, X_test)
    elif solver == 'svd':
        path = _ridge_path_svd(X_train, Y_train, alphas, X_test)

    if y_raveled:
        return path[:, :, 0]
    else:
        return path


def _solve_sparse_cg(X, y, alpha, max_iter=None, tol=1e-3):
    n_samples, n_features = X.shape
    X1 = sp_linalg.aslinearoperator(X)
    coefs = np.empty((y.shape[1], n_features))

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
            coef, info = sp_linalg.cg(C, y_column, tol=tol)
            coefs[i] = X1.rmatvec(coef)
        else:
            # linear ridge
            # w = inv(X^t X + alpha*Id) * X.T y
            y_column = X1.rmatvec(y_column)
            C = sp_linalg.LinearOperator(
                (n_features, n_features), matvec=mv, dtype=X.dtype)
            coefs[i], info = sp_linalg.cg(C, y_column, maxiter=max_iter,
                                          tol=tol)
        if info != 0:
            raise ValueError("Failed with error code %d" % info)

    return coefs


def _solve_lsqr(X, y, alpha, max_iter=None, tol=1e-3):
    n_samples, n_features = X.shape
    coefs = np.empty((y.shape[1], n_features))

    # According to the lsqr documentation, alpha = damp^2.
    sqrt_alpha = np.sqrt(alpha)

    for i in range(y.shape[1]):
        y_column = y[:, i]
        coefs[i] = sp_linalg.lsqr(X, y_column, damp=sqrt_alpha[i],
                                  atol=tol, btol=tol, iter_lim=max_iter)[0]

    return coefs


def _solve_cholesky(X, y, alpha, sample_weight=None):
    # w = inv(X^t X + alpha*Id) * X.T y
    n_samples, n_features = X.shape
    n_targets = y.shape[1]

    has_sw = sample_weight is not None

    if has_sw:
        sample_weight = sample_weight * np.ones(n_samples)
        sample_weight_matrix = sparse.dia_matrix((sample_weight, 0),
            shape=(n_samples, n_samples))
        weighted_X = safe_sparse_dot(sample_weight_matrix, X)
        A = safe_sparse_dot(weighted_X.T, X, dense_output=True)
        Xy = safe_sparse_dot(weighted_X.T, y, dense_output=True)
    else:
        A = safe_sparse_dot(X.T, X, dense_output=True)
        Xy = safe_sparse_dot(X.T, y, dense_output=True)

    one_alpha = np.array_equal(alpha, len(alpha) * [alpha[0]])

    if one_alpha:
        A.flat[::n_features + 1] += alpha[0]
        return linalg.solve(A, Xy, sym_pos=True,
                            overwrite_a=True).T
    else:
        coefs = np.empty([n_targets, n_features])
        for coef, target, current_alpha in zip(coefs, Xy.T, alpha):
            A.flat[::n_features + 1] += current_alpha
            coef[:] = linalg.solve(A, target, sym_pos=True,
                                   overwrite_a=False).ravel()
            A.flat[::n_features + 1] -= current_alpha
        return coefs


def _solve_cholesky_kernel(K, y, alpha, sample_weight=None):
    # dual_coef = inv(X X^t + alpha*Id) y
    n_samples = K.shape[0]
    n_targets = y.shape[1]

    one_alpha = np.array_equal(alpha, len(alpha) * [alpha[0]])

    has_sw = sample_weight is not None

    if has_sw:
        sw = np.sqrt(np.atleast_1d(sample_weight))
        y = y * sw[:, np.newaxis]
        K *= np.outer(sw, sw)

    if one_alpha:
        # Only one penalty, we can solve multi-target problems in one time.
        K.flat[::n_samples + 1] += alpha[0]

        dual_coef = linalg.solve(K, y, sym_pos=True, overwrite_a=True)

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


def _deprecate_dense_cholesky(solver):
    if solver == 'dense_cholesky':
        warnings.warn(DeprecationWarning("The name 'dense_cholesky' is "
                                         "deprecated. Using 'cholesky' "
                                         "instead. Changed in 0.15"))
        solver = 'cholesky'

    return solver


def ridge_regression(X, y, alpha, sample_weight=None, solver='auto',
                     max_iter=None, tol=1e-3):
    """Solve the ridge equation by the method of normal equations.

    Parameters
    ----------
    X : {array-like, sparse matrix, LinearOperator},
        shape = [n_samples, n_features]
        Training data

    y : array-like, shape = [n_samples] or [n_samples, n_targets]
        Target values

    alpha : {float, array-like},
        shape = [n_targets] if array-like
        The l_2 penalty to be used. If an array is passed, penalties are
        assumed to be specific to targets

    max_iter : int, optional
        Maximum number of iterations for conjugate gradient solver.
        The default value is determined by scipy.sparse.linalg.

    sample_weight : float or numpy array of shape [n_samples]
        Individual weights for each sample. If sample_weight is set, then
        the solver will automatically be set to 'cholesky'

    solver : {'auto', 'svd', 'cholesky', 'lsqr', 'sparse_cg', 'eigen'}
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
          scipy.sparse.linalg.lsqr. It is the fatest but may not be available
          in old scipy versions. It also uses an iterative procedure.

        - 'eigen' uses an eigenvalue decomposition of X.T.dot(X) or 
          X.dot(X.T) to compute the Ridge coefficients

        All three solvers support both dense and sparse data.

    tol: float
        Precision of the solution.

    Returns
    -------
    coef: array, shape = [n_features] or [n_targets, n_features]
        Weight vector(s).

    Notes
    -----
    This function won't compute the intercept.
    """

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

    has_sw = sample_weight is not None

    solver = _deprecate_dense_cholesky(solver)

    if solver == 'auto':
        # cholesky if it's a dense array and cg in
        # any other case
        if not sparse.issparse(X) or has_sw:
            solver = 'cholesky'
        else:
            solver = 'sparse_cg'

    elif solver == 'lsqr' and not hasattr(sp_linalg, 'lsqr'):
        warnings.warn("""lsqr not available on this machine, falling back
                      to sparse_cg.""")
        solver = 'sparse_cg'

    if has_sw:
        if np.atleast_1d(sample_weight).ndim > 1:
            raise ValueError("Sample weights must be 1D array or scalar")

        if solver != "cholesky":
            warnings.warn("sample_weight and class_weight not"
                          " supported in %s, fall back to "
                          "cholesky." % solver)
            solver = 'cholesky'

    # There should be either 1 or n_targets penalties
    alpha = np.asarray(alpha).ravel()
    if alpha.size not in [1, n_targets]:
        raise ValueError("Number of targets and number of penalties "
                         "do not correspond: %d != %d"
                         % (alpha.size, n_targets))

    if alpha.size == 1 and n_targets > 1:
        alpha = np.repeat(alpha, n_targets)

    if solver not in ('sparse_cg', 'cholesky', 'svd', 'lsqr', 'eigen'):
        raise ValueError('Solver %s not understood' % solver)

    if solver == 'sparse_cg':
        coef = _solve_sparse_cg(X, y, alpha, max_iter, tol)

    elif solver == "lsqr":
        coef = _solve_lsqr(X, y, alpha, max_iter, tol)

    elif solver == 'cholesky':
        if n_features > n_samples:
            K = safe_sparse_dot(X, X.T, dense_output=True)
            try:
                dual_coef = _solve_cholesky_kernel(K, y, alpha,
                                                         sample_weight)

                coef = safe_sparse_dot(X.T, dual_coef, dense_output=True).T
            except linalg.LinAlgError:
                # use SVD solver if matrix is singular
                solver = 'svd'

        else:
            try:
                coef = _solve_cholesky(X, y, alpha, sample_weight)
            except linalg.LinAlgError:
                # use SVD solver if matrix is singular
                solver = 'svd'

    if solver == 'svd':
        coef = ridge_path(X, y, np.atleast_2d(alpha), solver='svd')[0].T

    if solver == 'eigen':
        coef = ridge_path(X, y, np.atleast_2d(alpha), solver='eigen')[0].T
    if ravel:
        # When y was passed as a 1d-array, we flatten the coefficients.
        coef = coef.ravel()

    return coef


class _BaseRidge(six.with_metaclass(ABCMeta, LinearModel)):

    @abstractmethod
    def __init__(self, alpha=1.0, fit_intercept=True, normalize=False,
                 copy_X=True, max_iter=None, tol=1e-3, solver="auto"):
        self.alpha = alpha
        self.fit_intercept = fit_intercept
        self.normalize = normalize
        self.copy_X = copy_X
        self.max_iter = max_iter
        self.tol = tol
        self.solver = solver

    def fit(self, X, y, sample_weight=None):
        X, y = check_X_y(X, y, ['csr', 'csc', 'coo'], dtype=np.float, multi_output=True)

        if ((sample_weight is not None) and
                np.atleast_1d(sample_weight).ndim > 1):
            raise ValueError("Sample weights must be 1D array or scalar")

        X, y, X_mean, y_mean, X_std = self._center_data(
            X, y, self.fit_intercept, self.normalize, self.copy_X,
            sample_weight=sample_weight)

        solver = _deprecate_dense_cholesky(self.solver)

        self.coef_ = ridge_regression(X, y,
                                      alpha=self.alpha,
                                      sample_weight=sample_weight,
                                      max_iter=self.max_iter,
                                      tol=self.tol,
                                      solver=solver)
        self._set_intercept(X_mean, y_mean, X_std)
        return self


class Ridge(_BaseRidge, RegressorMixin):
    """Linear least squares with l2 regularization.

    This model solves a regression model where the loss function is
    the linear least squares function and regularization is given by
    the l2-norm. Also known as Ridge Regression or Tikhonov regularization.
    This estimator has built-in support for multi-variate regression
    (i.e., when y is a 2d-array of shape [n_samples, n_targets]).

    Parameters
    ----------
    alpha : {float, array-like}
        shape = [n_targets]
        Small positive values of alpha improve the conditioning of the problem
        and reduce the variance of the estimates.  Alpha corresponds to
        ``(2*C)^-1`` in other linear models such as LogisticRegression or
        LinearSVC. If an array is passed, penalties are assumed to be specific
        to the targets. Hence they must correspond in number.

    copy_X : boolean, optional, default True
        If True, X will be copied; else, it may be overwritten.

    fit_intercept : boolean
        Whether to calculate the intercept for this model. If set
        to false, no intercept will be used in calculations
        (e.g. data is expected to be already centered).

    max_iter : int, optional
        Maximum number of iterations for conjugate gradient solver.
        The default value is determined by scipy.sparse.linalg.

    normalize : boolean, optional, default False
        If True, the regressors X will be normalized before regression.

    solver : {'auto', 'svd', 'cholesky', 'lsqr', 'sparse_cg'}
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
          scipy.sparse.linalg.lsqr. It is the fatest but may not be available
          in old scipy versions. It also uses an iterative procedure.

        - 'eigen' uses an eigenvalue decomposition of X.T.dot(X) or 
          X.dot(X.T) to compute the Ridge coefficients

        All three solvers support both dense and sparse data.

    tol : float
        Precision of the solution.

    Attributes
    ----------
    coef_ : array, shape = [n_features] or [n_targets, n_features]
        Weight vector(s).

    See also
    --------
    RidgeClassifier, RidgeCV

    Examples
    --------
    >>> from sklearn.linear_model import Ridge
    >>> import numpy as np
    >>> n_samples, n_features = 10, 5
    >>> np.random.seed(0)
    >>> y = np.random.randn(n_samples)
    >>> X = np.random.randn(n_samples, n_features)
    >>> clf = Ridge(alpha=1.0)
    >>> clf.fit(X, y) # doctest: +NORMALIZE_WHITESPACE
    Ridge(alpha=1.0, copy_X=True, fit_intercept=True, max_iter=None,
          normalize=False, solver='auto', tol=0.001)
    """
    def __init__(self, alpha=1.0, fit_intercept=True, normalize=False,
                 copy_X=True, max_iter=None, tol=1e-3, solver="auto"):
        super(Ridge, self).__init__(alpha=alpha, fit_intercept=fit_intercept,
                                    normalize=normalize, copy_X=copy_X,
                                    max_iter=max_iter, tol=tol, solver=solver)

    def fit(self, X, y, sample_weight=None):
        """Fit Ridge regression model

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
        return super(Ridge, self).fit(X, y, sample_weight=sample_weight)


class RidgeClassifier(LinearClassifierMixin, _BaseRidge):
    """Classifier using Ridge regression.

    Parameters
    ----------
    alpha : float
        Small positive values of alpha improve the conditioning of the problem
        and reduce the variance of the estimates.  Alpha corresponds to
        ``(2*C)^-1`` in other linear models such as LogisticRegression or
        LinearSVC.

    class_weight : dict, optional
        Weights associated with classes in the form
        ``{class_label : weight}``. If not given, all classes are
        supposed to have weight one.

    copy_X : boolean, optional, default True
        If True, X will be copied; else, it may be overwritten.

    fit_intercept : boolean
        Whether to calculate the intercept for this model. If set to false, no
        intercept will be used in calculations (e.g. data is expected to be
        already centered).

    max_iter : int, optional
        Maximum number of iterations for conjugate gradient solver.
        The default value is determined by scipy.sparse.linalg.

    normalize : boolean, optional, default False
        If True, the regressors X will be normalized before regression.

    solver : {'auto', 'svd', 'cholesky', 'lsqr', 'sparse_cg'}
        Solver to use in the computational
        routines. 'svd' will use a Singular value decomposition to obtain
        the solution, 'cholesky' will use the standard
        scipy.linalg.solve function, 'sparse_cg' will use the
        conjugate gradient solver as found in
        scipy.sparse.linalg.cg while 'auto' will chose the most
        appropriate depending on the matrix X. 'lsqr' uses
        a direct regularized least-squares routine provided by scipy.

    tol : float
        Precision of the solution.

    Attributes
    ----------
    coef_ : array, shape = [n_features] or [n_classes, n_features]
        Weight vector(s).

    See also
    --------
    Ridge, RidgeClassifierCV

    Notes
    -----
    For multi-class classification, n_class classifiers are trained in
    a one-versus-all approach. Concretely, this is implemented by taking
    advantage of the multi-variate response support in Ridge.
    """
    def __init__(self, alpha=1.0, fit_intercept=True, normalize=False,
                 copy_X=True, max_iter=None, tol=1e-3, class_weight=None,
                 solver="auto"):
        super(RidgeClassifier, self).__init__(
            alpha=alpha, fit_intercept=fit_intercept, normalize=normalize,
            copy_X=copy_X, max_iter=max_iter, tol=tol, solver=solver)
        self.class_weight = class_weight

    def fit(self, X, y):
        """Fit Ridge regression model.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape = [n_samples,n_features]
            Training data

        y : array-like, shape = [n_samples]
            Target values

        Returns
        -------
        self : returns an instance of self.
        """
        self._label_binarizer = LabelBinarizer(pos_label=1, neg_label=-1)
        Y = self._label_binarizer.fit_transform(y)
        if not self._label_binarizer.y_type_.startswith('multilabel'):
            y = column_or_1d(y, warn=True)

        if self.class_weight:
            cw = compute_class_weight(self.class_weight,
                                      self.classes_, y)
            # get the class weight corresponding to each sample
            sample_weight = cw[np.searchsorted(self.classes_, y)]
        else:
            sample_weight = None

        super(RidgeClassifier, self).fit(X, Y, sample_weight=sample_weight)
        return self

    @property
    def classes_(self):
        return self._label_binarizer.classes_


class _RidgeGCV(LinearModel):
    """Ridge regression with built-in Generalized Cross-Validation

    It allows efficient Leave-One-Out cross-validation.

    This class is not intended to be used directly. Use RidgeCV instead.

    Notes
    -----

    We want to solve (K + alpha*Id)c = y,
    where K = X X^T is the kernel matrix.

    Let G = (K + alpha*Id)^-1.

    Dual solution: c = Gy
    Primal solution: w = X^T c

    Compute eigendecomposition K = Q V Q^T.
    Then G = Q (V + alpha*Id)^-1 Q^T,
    where (V + alpha*Id) is diagonal.
    It is thus inexpensive to inverse for many alphas.

    Let loov be the vector of prediction values for each example
    when the model was fitted with all examples but this example.

    loov = (KGY - diag(KG)Y) / diag(I-KG)

    Let looe be the vector of prediction errors for each example
    when the model was fitted with all examples but this example.

    looe = y - loov = c / diag(G)

    References
    ----------
    http://cbcl.mit.edu/projects/cbcl/publications/ps/MIT-CSAIL-TR-2007-025.pdf
    http://www.mit.edu/~9.520/spring07/Classes/rlsslides.pdf
    """

    def __init__(self, alphas=[0.1, 1.0, 10.0],
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

    def fit(self, X, y, sample_weight=None):
        """Fit Ridge regression model

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape = [n_samples, n_features]
            Training data

        y : array-like, shape = [n_samples] or [n_samples, n_targets]
            Target values

        sample_weight : float or array-like of shape [n_samples]
            Sample weight

        Returns
        -------
        self : Returns self.
        """
        X, y = check_X_y(X, y, ['csr', 'csc', 'coo'], dtype=np.float, multi_output=True)

        n_samples, n_features = X.shape

        X, y, X_mean, y_mean, X_std = LinearModel._center_data(
            X, y, self.fit_intercept, self.normalize, self.copy_X,
            sample_weight=sample_weight)

        gcv_mode = self.gcv_mode
        with_sw = len(np.shape(sample_weight))

        if gcv_mode is None or gcv_mode == 'auto':
            if sparse.issparse(X) or n_features > n_samples or with_sw:
                gcv_mode = 'eigen'
            else:
                gcv_mode = 'svd'

        if gcv_mode not in ['eigen', 'svd']:
            raise ValueError('bad gcv_mode "%s"' % gcv_mode)

        n_y = 1 if len(y.shape) == 1 else y.shape[1]
        cv_values = np.zeros((n_samples * n_y, len(self.alphas)))
        C = []

        scorer = check_scoring(self, scoring=self.scoring, allow_none=True,
                               score_overrides_loss=True)
        error = scorer is None

        alphas = np.atleast_2d(self.alphas.T).T
        mode = 'looe' if error else 'loov'
        y_is_raveled = y.ndim == 1
        y = np.atleast_2d(y.T).T
        if gcv_mode == 'eigen':
            out, C = _kernel_ridge_path_eigen(
                X, y, alphas, sample_weight=sample_weight, mode=mode)
        else:
            out, C = _ridge_gcv_path_svd(
                X, y, alphas, sample_weight=sample_weight, mode=mode)

        if mode == 'looe':
            out = out ** 2

        cv_values = out.transpose(1, 2, 0)

        if error:
            best = cv_values.mean(axis=0).argmin(axis=1)
        else:
            # The scorer want an object that will make the predictions but
            # they are already computed efficiently by _RidgeGCV. This
            # identity_estimator will just return them
            def identity_estimator():
                pass
            identity_estimator.decision_function = lambda y_predict: y_predict
            identity_estimator.predict = lambda y_predict: y_predict

            # XXX
            # This is the wrong application of the scorer: It is given all
            # of y and all cv values, although it should only be given one
            # value at a time, because the CV we are using here is LOO!
            out = np.array([scorer(identity_estimator, y, cv_values[..., i])
                   for i in range(len(self.alphas))])
            best = np.argmax(out)
            best = np.atleast_1d(best)

        self.alpha_ = self.alphas[best]
        self.dual_coef_ = C[best, :, np.arange(len(best))].T

        if y_is_raveled:
            C = C[:, :, 0]
            self.dual_coef_ = self.dual_coef_.ravel()

        self.coef_ = safe_sparse_dot(self.dual_coef_.T, X)
        self._set_intercept(X_mean, y_mean, X_std)

        if self.store_cv_values:
            if y_is_raveled:
                cv_values_shape = n_samples, len(self.alphas)
            else:
                cv_values_shape = n_samples, n_y, len(self.alphas)
            self.cv_values_ = cv_values.reshape(cv_values_shape)

        return self


class _BaseRidgeCV(LinearModel):
    def __init__(self, alphas=np.array([0.1, 1.0, 10.0]),
                 fit_intercept=True, normalize=False,
                 scoring=None, cv=None, gcv_mode=None,
                 store_cv_values=False, solver=None, copyX=True):
        self.alphas = alphas
        self.fit_intercept = fit_intercept
        self.normalize = normalize
        self.scoring = scoring
        self.cv = cv
        self.gcv_mode = gcv_mode
        self.store_cv_values = store_cv_values
        self.solver = solver
        self.copyX = copyX

    def fit(self, X, y, sample_weight=None):
        """Fit Ridge regression model

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]
            Training data

        y : array-like, shape = [n_samples] or [n_samples, n_targets]
            Target values

        sample_weight : float or array-like of shape [n_samples]
            Sample weight

        Returns
        -------
        self : Returns self.
        """
        if self.cv is None:
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
            self.coef_ = estimator.coef_
            self.intercept_ = estimator.intercept_
        else:
            if self.solver in ['eigen', 'svd']:
                cv = check_cv(self.cv)
                scorer = check_scoring(self,
                                       scoring=self.scoring,
                                       allow_none=True,
                                       score_overrides_loss=True)

                alphas = np.atleast_2d(np.array(self.alphas).T).T

                # need to make an object with methods that scorers like
                def identity_estimator():
                    pass
                identity_estimator.predict = lambda pred: pred
                identity_estimator.decision_function = lambda dec: dec
                identity_estimator.score = _multi_r2

                # Need to center data before regression paths ...
                # Means we probably need to write _center_and_ridge_path
                # Or even a center + ridge_path + predict + score
                # The following is rudimentary to ensure functionality
                all_scores = np.zeros([len(cv), len(alphas), y.shape[1]])
                for (train, test), scores in zip(cv, all_scores):
                    X_train, y_train = X[train], y[train]
                    (X_train, y_train,
                     X_train_mean, y_train_mean,
                     X_train_std) = self._center_data(
                        X_train, y_train, self.fit_intercept,
                        self.normalize, copy=True,
                        sample_weight=sample_weight)

                    predictions = ridge_path(
                        X_train, y_train, alphas,
                        (X[test] - X_train_mean) / X_train_std,
                        solver=self.solver) + y_train_mean[np.newaxis]
                    for score, prediction in zip(scores, predictions):
                        score[:] = scorer(
                            identity_estimator, y[test], prediction)

                arg_best_penalty = all_scores.mean(axis=0).argmax(axis=0)
                self.best_alphas_ = alphas[
                    arg_best_penalty,
                    np.maximum(
                        np.arange(len(arg_best_penalty)),
                        alphas.shape[1] - 1)]

                (X, y, X_mean, y_mean, X_std) = self._center_data(
                    X, y, self.fit_intercept, self.normalize,
                    copy=self.copyX, sample_weight=sample_weight)
                self.coef_ = ridge_path(X, y,
                                        np.atleast_2d(self.best_alphas_),
                                        solver=self.solver)[0].T
                self._set_intercept(X_mean, y_mean, X_std)
            else:
                # do the old grid search
                if self.store_cv_values:
                    raise ValueError("cv!=None and store_cv_values=True "
                                 " are incompatible")
                parameters = {'alpha': self.alphas}
                # FIXME: sample_weight must be split into training/validation data
                #        too!
                #fit_params = {'sample_weight' : sample_weight}
                fit_params = {}
                gs = GridSearchCV(Ridge(fit_intercept=self.fit_intercept),
                              parameters, fit_params=fit_params, cv=self.cv)
                gs.fit(X, y)
                estimator = gs.best_estimator_
                self.alpha_ = gs.best_estimator_.alpha

                self.coef_ = estimator.coef_
                self.intercept_ = estimator.intercept_

        return self


class RidgeCV(_BaseRidgeCV, RegressorMixin):
    """Ridge regression with built-in cross-validation.

    By default, it performs Generalized Cross-Validation, which is a form of
    efficient Leave-One-Out cross-validation.

    Parameters
    ----------
    alphas: numpy array of shape [n_alphas]
        Array of alpha values to try.
        Small positive values of alpha improve the conditioning of the
        problem and reduce the variance of the estimates.
        Alpha corresponds to ``(2*C)^-1`` in other linear models such as
        LogisticRegression or LinearSVC.

    fit_intercept : boolean
        Whether to calculate the intercept for this model. If set
        to false, no intercept will be used in calculations
        (e.g. data is expected to be already centered).

    normalize : boolean, optional, default False
        If True, the regressors X will be normalized before regression.

    scoring : string, callable or None, optional, default: None
        A string (see model evaluation documentation) or
        a scorer callable object / function with signature
        ``scorer(estimator, X, y)``.

    cv : cross-validation generator, optional
        If None, Generalized Cross-Validation (efficient Leave-One-Out)
        will be used.

    gcv_mode : {None, 'auto', 'svd', eigen'}, optional
        Flag indicating which strategy to use when performing
        Generalized Cross-Validation. Options are::

            'auto' : use svd if n_samples > n_features or when X is a sparse
                     matrix, otherwise use eigen
            'svd' : force computation via singular value decomposition of X
                    (does not work for sparse matrices)
            'eigen' : force computation via eigendecomposition of X^T X

        The 'auto' mode is the default and is intended to pick the cheaper
        option of the two depending upon the shape and format of the training
        data.

    store_cv_values : boolean, default=False
        Flag indicating if the cross-validation values corresponding to
        each alpha should be stored in the `cv_values_` attribute (see
        below). This flag is only compatible with `cv=None` (i.e. using
        Generalized Cross-Validation).

    Attributes
    ----------
    cv_values_ : array, shape = [n_samples, n_alphas] or \
        shape = [n_samples, n_targets, n_alphas], optional
        Cross-validation values for each alpha (if `store_cv_values=True` and \
        `cv=None`). After `fit()` has been called, this attribute will \
        contain the mean squared errors (by default) or the values of the \
        `{loss,score}_func` function (if provided in the constructor).

    coef_ : array, shape = [n_features] or [n_targets, n_features]
        Weight vector(s).

    alpha_ : float
        Estimated regularization parameter.

    intercept_ : float | array, shape = (n_targets,)
        Independent term in decision function. Set to 0.0 if
        ``fit_intercept = False``.

    See also
    --------
    Ridge: Ridge regression
    RidgeClassifier: Ridge classifier
    RidgeClassifierCV: Ridge classifier with built-in cross validation
    """
    pass


class RidgeClassifierCV(LinearClassifierMixin, _BaseRidgeCV):
    """Ridge classifier with built-in cross-validation.

    By default, it performs Generalized Cross-Validation, which is a form of
    efficient Leave-One-Out cross-validation. Currently, only the n_features >
    n_samples case is handled efficiently.

    Parameters
    ----------
    alphas: numpy array of shape [n_alphas]
        Array of alpha values to try.
        Small positive values of alpha improve the conditioning of the
        problem and reduce the variance of the estimates.
        Alpha corresponds to ``(2*C)^-1`` in other linear models such as
        LogisticRegression or LinearSVC.

    fit_intercept : boolean
        Whether to calculate the intercept for this model. If set
        to false, no intercept will be used in calculations
        (e.g. data is expected to be already centered).

    normalize : boolean, optional, default False
        If True, the regressors X will be normalized before regression.

    scoring : string, callable or None, optional, default: None
        A string (see model evaluation documentation) or
        a scorer callable object / function with signature
        ``scorer(estimator, X, y)``.

    cv : cross-validation generator, optional
        If None, Generalized Cross-Validation (efficient Leave-One-Out)
        will be used.

    class_weight : dict, optional
        Weights associated with classes in the form
        ``{class_label : weight}``. If not given, all classes are
        supposed to have weight one.

    Attributes
    ----------
    cv_values_ : array, shape = [n_samples, n_alphas] or \
    shape = [n_samples, n_responses, n_alphas], optional
        Cross-validation values for each alpha (if `store_cv_values=True` and
    `cv=None`). After `fit()` has been called, this attribute will contain \
    the mean squared errors (by default) or the values of the \
    `{loss,score}_func` function (if provided in the constructor).

    coef_ : array, shape = [n_features] or [n_targets, n_features]
        Weight vector(s).

    alpha_ : float
        Estimated regularization parameter

    See also
    --------
    Ridge: Ridge regression
    RidgeClassifier: Ridge classifier
    RidgeCV: Ridge regression with built-in cross validation

    Notes
    -----
    For multi-class classification, n_class classifiers are trained in
    a one-versus-all approach. Concretely, this is implemented by taking
    advantage of the multi-variate response support in Ridge.
    """
    def __init__(self, alphas=np.array([0.1, 1.0, 10.0]), fit_intercept=True,
                 normalize=False, scoring=None, cv=None, class_weight=None):
        super(RidgeClassifierCV, self).__init__(
            alphas=alphas, fit_intercept=fit_intercept, normalize=normalize,
            scoring=scoring, cv=cv)
        self.class_weight = class_weight

    def fit(self, X, y, sample_weight=None):
        """Fit the ridge classifier.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            Training vectors, where n_samples is the number of samples
            and n_features is the number of features.

        y : array-like, shape (n_samples,)
            Target values.

        sample_weight : float or numpy array of shape (n_samples,)
            Sample weight.

        Returns
        -------
        self : object
            Returns self.
        """
        if sample_weight is None:
            sample_weight = 1.

        self._label_binarizer = LabelBinarizer(pos_label=1, neg_label=-1)
        Y = self._label_binarizer.fit_transform(y)
        if not self._label_binarizer.y_type_.startswith('multilabel'):
            y = column_or_1d(y, warn=True)
        cw = compute_class_weight(self.class_weight,
                                  self.classes_, Y)
        # modify the sample weights with the corresponding class weight
        sample_weight *= cw[np.searchsorted(self.classes_, y)]
        _BaseRidgeCV.fit(self, X, Y, sample_weight=sample_weight)
        return self

    @property
    def classes_(self):
        return self._label_binarizer.classes_
