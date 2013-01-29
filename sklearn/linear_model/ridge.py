"""
Ridge regression
"""

# Author: Mathieu Blondel <mathieu@mblondel.org>
#         Reuben Fletcher-Costin <reuben.fletchercostin@gmail.com>
#         Fabian Pedregosa <fabian@fseoane.net>
# License: Simplified BSD


from abc import ABCMeta, abstractmethod
import warnings

import numpy as np
from scipy import linalg
from scipy import sparse
from scipy.sparse import linalg as sp_linalg

from .base import LinearClassifierMixin, LinearModel
from ..base import RegressorMixin
from ..utils.extmath import safe_sparse_dot
from ..utils import safe_asarray
from ..preprocessing import LabelBinarizer
from ..grid_search import GridSearchCV


def ridge_regression(X, y, alpha, sample_weight=1.0, solver='auto',
                     max_iter=None, tol=1e-3):
    """Solve the ridge equation by the method of normal equations.

    Parameters
    ----------
    X : {array-like, sparse matrix, LinearOperator},
        shape = [n_samples, n_features]
        Training data

    y : array-like, shape = [n_samples] or [n_samples, n_targets]
        Target values

    max_iter : int, optional
        Maximum number of iterations for conjugate gradient solver.
        The default value is determined by scipy.sparse.linalg.

    sample_weight : float or numpy array of shape [n_samples]
        Individual weights for each sample

    solver : {'auto', 'dense_cholesky', 'lsqr', 'sparse_cg'}
        Solver to use in the computational routines:

        - 'auto' chooses the solver automatically based on the type of data.

        - 'dense_cholesky' uses the standard scipy.linalg.solve function to
          obtain a closed-form solution.

        - 'sparse_cg' uses the conjugate gradient solver as found in
          scipy.sparse.linalg.cg. As an iterative algorithm, this solver is
          more appropriate than 'dense_cholesky' for large-scale data
          (possibility to set `tol` and `max_iter`).

        - 'lsqr' uses the dedicated regularized least-squares routine
          scipy.sparse.linalg.lsqr. It is the fatest but may not be available
          in old scipy versions. It also uses an iterative procedure.

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
    has_sw = isinstance(sample_weight, np.ndarray) or sample_weight != 1.0

    if solver == 'auto':
        # cholesky if it's a dense array and cg in
        # any other case
        if hasattr(X, '__array__'):
            solver = 'dense_cholesky'
        else:
            solver = 'sparse_cg'

    elif solver == 'lsqr' and not hasattr(sp_linalg, 'lsqr'):
        warnings.warn("""lsqr not available on this machine, falling back
                      to sparse_cg.""")
        solver = 'sparse_cg'

    if has_sw:
        solver = 'dense_cholesky'

    if solver == 'sparse_cg':
        # gradient descent
        X1 = sp_linalg.aslinearoperator(X)
        if y.ndim == 1:
            y1 = np.reshape(y, (-1, 1))
        else:
            y1 = y
        coefs = np.empty((y1.shape[1], n_features))

        if n_features > n_samples:
            def mv(x):
                return X1.matvec(X1.rmatvec(x)) + alpha * x
        else:
            def mv(x):
                return X1.rmatvec(X1.matvec(x)) + alpha * x

        for i in range(y1.shape[1]):
            y_column = y1[:, i]
            if n_features > n_samples:
                # kernel ridge
                # w = X.T * inv(X X^t + alpha*Id) y
                C = sp_linalg.LinearOperator(
                    (n_samples, n_samples), matvec=mv, dtype=X.dtype)
                coef, info = sp_linalg.cg(C, y_column, tol=tol)
                coefs[i] = X1.rmatvec(coef)
            else:
                # ridge
                # w = inv(X^t X + alpha*Id) * X.T y
                y_column = X1.rmatvec(y_column)
                C = sp_linalg.LinearOperator(
                    (n_features, n_features), matvec=mv, dtype=X.dtype)
                coefs[i], info = sp_linalg.cg(C, y_column, maxiter=max_iter,
                                              tol=tol)
            if info != 0:
                raise ValueError("Failed with error code %d" % info)

        if y.ndim == 1:
            coefs = np.ravel(coefs)

        return coefs
    elif solver == "lsqr":
        if y.ndim == 1:
            y1 = np.reshape(y, (-1, 1))
        else:
            y1 = y
        coefs = np.empty((y1.shape[1], n_features))

        # According to the lsqr documentation, alpha = damp^2.
        sqrt_alpha = np.sqrt(alpha)

        for i in range(y1.shape[1]):
            y_column = y1[:, i]
            coefs[i] = sp_linalg.lsqr(X, y_column, damp=sqrt_alpha,
                                      atol=tol, btol=tol, iter_lim=max_iter)[0]

        if y.ndim == 1:
            coefs = np.ravel(coefs)

        return coefs
    else:
        # normal equations (cholesky) method
        if n_features > n_samples or has_sw:
            # kernel ridge
            # w = X.T * inv(X X^t + alpha*Id) y
            K = safe_sparse_dot(X, X.T, dense_output=True)
            if has_sw:
                # We are doing a little danse with the sample weights to
                # avoid copying the original X, which could be big
                sw = np.sqrt(sample_weight)
                if y.ndim == 1:
                    y = y * sw
                else:
                    # Deal with multiple-output problems
                    y = y * sw[:, np.newaxis]
                K *= np.outer(sw, sw)
            K.flat[::n_samples + 1] += alpha
            dual_coef = linalg.solve(K, y,
                                     sym_pos=True, overwrite_a=True)
            if has_sw:
                if dual_coef.ndim == 1:
                    dual_coef *= sw
                else:
                    # Deal with multiple-output problems
                    dual_coef *= sw[:, np.newaxis]
            coef = safe_sparse_dot(X.T, dual_coef, dense_output=True)
        else:
            # ridge
            # w = inv(X^t X + alpha*Id) * X.T y
            A = safe_sparse_dot(X.T, X, dense_output=True)
            A.flat[::n_features + 1] += alpha
            Xy = safe_sparse_dot(X.T, y, dense_output=True)
            coef = linalg.solve(A, Xy, sym_pos=True, overwrite_a=True)

        return coef.T


class _BaseRidge(LinearModel):
    __metaclass__ = ABCMeta

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

    def fit(self, X, y, sample_weight=1.0, solver=None):
        X = safe_asarray(X, dtype=np.float)
        y = np.asarray(y, dtype=np.float)

        X, y, X_mean, y_mean, X_std = self._center_data(
            X, y, self.fit_intercept, self.normalize, self.copy_X,
            sample_weight=sample_weight)

        self.coef_ = ridge_regression(X, y,
                                      alpha=self.alpha,
                                      sample_weight=sample_weight,
                                      solver=solver,
                                      max_iter=self.max_iter,
                                      tol=self.tol)
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
    alpha : float
        Small positive values of alpha improve the conditioning of the problem
        and reduce the variance of the estimates.  Alpha corresponds to
        ``(2*C)^-1`` in other linear models such as LogisticRegression or
        LinearSVC.

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

    solver : {'auto', 'dense_cholesky', 'lsqr', 'sparse_cg'}
        Solver to use in the computational routines:

        - 'auto' chooses the solver automatically based on the type of data.

        - 'dense_cholesky' uses the standard scipy.linalg.solve function to
          obtain a closed-form solution.

        - 'sparse_cg' uses the conjugate gradient solver as found in
          scipy.sparse.linalg.cg. As an iterative algorithm, this solver is
          more appropriate than 'dense_cholesky' for large-scale data
          (possibility to set `tol` and `max_iter`).

        - 'lsqr' uses the dedicated regularized least-squares routine
          scipy.sparse.linalg.lsqr. It is the fatest but may not be available
          in old scipy versions. It also uses an iterative procedure.

        All three solvers support both dense and sparse data.

    tol : float
        Precision of the solution.

    Attributes
    ----------
    `coef_` : array, shape = [n_features] or [n_targets, n_features]
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

    def fit(self, X, y, sample_weight=1.0, solver=None):
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
        if solver is None:
            solver = self.solver
        else:
            # The fit method should be removed from Ridge when this warning is
            # removed
            warnings.warn("""solver option in fit is deprecated and will be
                          removed in v0.14.""")

        return _BaseRidge.fit(self, X, y, solver=solver,
                              sample_weight=sample_weight)


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
        {class_label : weight}. If not given, all classes are
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

    solver : {'auto', 'dense_cholesky', 'lsqr', 'sparse_cg'}
        Solver to use in the computational
        routines. 'dense_cholesky' will use the standard
        scipy.linalg.solve function, 'sparse_cg' will use the
        conjugate gradient solver as found in
        scipy.sparse.linalg.cg while 'auto' will chose the most
        appropriate depending on the matrix X. 'lsqr' uses
        a direct regularized least-squares routine provided by scipy.

    tol : float
        Precision of the solution.

    Attributes
    ----------
    `coef_` : array, shape = [n_features] or [n_classes, n_features]
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

    def fit(self, X, y, solver=None):
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
        if self.class_weight is None:
            class_weight = {}
        else:
            class_weight = self.class_weight

        if solver is None:
            solver = self.solver
        else:
            warnings.warn("""solver option in fit is deprecated and will be
                          removed in v0.14.""")

        sample_weight_classes = np.array([class_weight.get(k, 1.0) for k in y])
        self._label_binarizer = LabelBinarizer(pos_label=1, neg_label=-1)
        Y = self._label_binarizer.fit_transform(y)
        _BaseRidge.fit(self, X, Y, solver=solver,
                       sample_weight=sample_weight_classes)
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

    def __init__(self, alphas=[0.1, 1.0, 10.0], fit_intercept=True,
                 normalize=False, score_func=None, loss_func=None,
                 copy_X=True, gcv_mode=None, store_cv_values=False):
        self.alphas = np.asarray(alphas)
        self.fit_intercept = fit_intercept
        self.normalize = normalize
        self.score_func = score_func
        self.loss_func = loss_func
        self.copy_X = copy_X
        self.gcv_mode = gcv_mode
        self.store_cv_values = store_cv_values

    def _pre_compute(self, X, y):
        # even if X is very sparse, K is usually very dense
        K = safe_sparse_dot(X, X.T, dense_output=True)
        v, Q = linalg.eigh(K)
        QT_y = np.dot(Q.T, y)
        return v, Q, QT_y

    def _decomp_diag(self, v_prime, Q):
        # compute diagonal of the matrix: dot(Q, dot(diag(v_prime), Q^T))
        return (v_prime * Q ** 2).sum(axis=-1)

    def _diag_dot(self, D, B):
        # compute dot(diag(D), B)
        if len(B.shape) > 1:
            # handle case where B is > 1-d
            D = D[(slice(None), ) + (np.newaxis, ) * (len(B.shape) - 1)]
        return D * B

    def _errors(self, alpha, y, v, Q, QT_y):
        # don't construct matrix G, instead compute action on y & diagonal
        w = 1.0 / (v + alpha)
        c = np.dot(Q, self._diag_dot(w, QT_y))
        G_diag = self._decomp_diag(w, Q)
        # handle case where y is 2-d
        if len(y.shape) != 1:
            G_diag = G_diag[:, np.newaxis]
        return (c / G_diag) ** 2, c

    def _values(self, alpha, y, v, Q, QT_y):
        # don't construct matrix G, instead compute action on y & diagonal
        w = 1.0 / (v + alpha)
        c = np.dot(Q, self._diag_dot(w, QT_y))
        G_diag = self._decomp_diag(w, Q)
        # handle case where y is 2-d
        if len(y.shape) != 1:
            G_diag = G_diag[:, np.newaxis]
        return y - (c / G_diag), c

    def _pre_compute_svd(self, X, y):
        if sparse.issparse(X) and hasattr(X, 'toarray'):
            X = X.toarray()
        U, s, _ = np.linalg.svd(X, full_matrices=0)
        v = s ** 2
        UT_y = np.dot(U.T, y)
        return v, U, UT_y

    def _errors_svd(self, alpha, y, v, U, UT_y):
        w = ((v + alpha) ** -1) - (alpha ** -1)
        c = np.dot(U, self._diag_dot(w, UT_y)) + (alpha ** -1) * y
        G_diag = self._decomp_diag(w, U) + (alpha ** -1)
        if len(y.shape) != 1:
            # handle case where y is 2-d
            G_diag = G_diag[:, np.newaxis]
        return (c / G_diag) ** 2, c

    def _values_svd(self, alpha, y, v, U, UT_y):
        w = ((v + alpha) ** -1) - (alpha ** -1)
        c = np.dot(U, self._diag_dot(w, UT_y)) + (alpha ** -1) * y
        G_diag = self._decomp_diag(w, U) + (alpha ** -1)
        if len(y.shape) != 1:
            # handle case when y is 2-d
            G_diag = G_diag[:, np.newaxis]
        return y - (c / G_diag), c

    def fit(self, X, y, sample_weight=1.0):
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
        X = safe_asarray(X, dtype=np.float)
        y = np.asarray(y, dtype=np.float)

        n_samples, n_features = X.shape

        X, y, X_mean, y_mean, X_std = LinearModel._center_data(
            X, y, self.fit_intercept, self.normalize, self.copy_X,
            sample_weight=sample_weight)

        gcv_mode = self.gcv_mode
        with_sw = len(np.shape(sample_weight))

        if gcv_mode is None or gcv_mode == 'auto':
            if n_features > n_samples or with_sw:
                gcv_mode = 'eigen'
            else:
                gcv_mode = 'svd'
        elif gcv_mode == "svd" and with_sw:
            # FIXME non-uniform sample weights not yet supported
            warnings.warn("non-uniform sample weights unsupported for svd, "
                          "forcing usage of eigen")
            gcv_mode = 'eigen'

        if gcv_mode == 'eigen':
            _pre_compute = self._pre_compute
            _errors = self._errors
            _values = self._values
        elif gcv_mode == 'svd':
            # assert n_samples >= n_features
            _pre_compute = self._pre_compute_svd
            _errors = self._errors_svd
            _values = self._values_svd
        else:
            raise ValueError('bad gcv_mode "%s"' % gcv_mode)

        v, Q, QT_y = _pre_compute(X, y)
        n_y = 1 if len(y.shape) == 1 else y.shape[1]
        cv_values = np.zeros((n_samples * n_y, len(self.alphas)))
        C = []

        error = self.score_func is None and self.loss_func is None

        for i, alpha in enumerate(self.alphas):
            if error:
                out, c = _errors(sample_weight * alpha, y, v, Q, QT_y)
            else:
                out, c = _values(sample_weight * alpha, y, v, Q, QT_y)
            cv_values[:, i] = out.ravel()
            C.append(c)

        if error:
            best = cv_values.mean(axis=0).argmin()
        else:
            func = self.score_func if self.score_func else self.loss_func
            out = [func(y.ravel(), cv_values[:, i])
                   for i in range(len(self.alphas))]
            best = np.argmax(out) if self.score_func else np.argmin(out)

        self.alpha_ = self.alphas[best]
        self.dual_coef_ = C[best]
        self.coef_ = safe_sparse_dot(self.dual_coef_.T, X)

        self._set_intercept(X_mean, y_mean, X_std)

        if self.store_cv_values:
            if len(y.shape) == 1:
                cv_values_shape = n_samples, len(self.alphas)
            else:
                cv_values_shape = n_samples, n_y, len(self.alphas)
            self.cv_values_ = cv_values.reshape(cv_values_shape)

        return self

    @property
    def best_alpha(self):
        warnings.warn("Use alpha_. Using best_alpha is deprecated"
                      "since version 0.12, and backward compatibility "
                      "won't be maintained from version 0.14 onward. ",
                      DeprecationWarning, stacklevel=2)
        return self.alpha_


class _BaseRidgeCV(LinearModel):

    def __init__(self, alphas=np.array([0.1, 1.0, 10.0]),
                 fit_intercept=True, normalize=False, score_func=None,
                 loss_func=None, cv=None, gcv_mode=None,
                 store_cv_values=False):
        self.alphas = alphas
        self.fit_intercept = fit_intercept
        self.normalize = normalize
        self.score_func = score_func
        self.loss_func = loss_func
        self.cv = cv
        self.gcv_mode = gcv_mode
        self.store_cv_values = store_cv_values

    def fit(self, X, y, sample_weight=1.0):
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
                                  score_func=self.score_func,
                                  loss_func=self.loss_func,
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

    score_func: callable, optional
        function that takes 2 arguments and compares them in
        order to evaluate the performance of prediction (big is good)
        if None is passed, the score of the estimator is maximized

    loss_func: callable, optional
        function that takes 2 arguments and compares them in
        order to evaluate the performance of prediction (small is good)
        if None is passed, the score of the estimator is maximized

    cv : cross-validation generator, optional
        If None, Generalized Cross-Validation (efficient Leave-One-Out)
        will be used.

    gcv_mode : {None, 'auto', 'svd', eigen'}, optional
        Flag indicating which strategy to use when performing
        Generalized Cross-Validation. Options are::

            'auto' : use svd if n_samples > n_features, otherwise use eigen
            'svd' : force computation via singular value decomposition of X
            'eigen' : force computation via eigendecomposition of X^T X

        The 'auto' mode is the default and is intended to pick the cheaper \
        option of the two depending upon the shape of the training data.

    store_cv_values : boolean, default=False
        Flag indicating if the cross-validation values corresponding to
        each alpha should be stored in the `cv_values_` attribute (see
        below). This flag is only compatible with `cv=None` (i.e. using
        Generalized Cross-Validation).

    Attributes
    ----------
    `cv_values_` : array, shape = [n_samples, n_alphas] or \
        shape = [n_samples, n_targets, n_alphas], optional
        Cross-validation values for each alpha (if `store_cv_values=True` and \
        `cv=None`). After `fit()` has been called, this attribute will \
        contain the mean squared errors (by default) or the values of the \
        `{loss,score}_func` function (if provided in the constructor).

    `coef_` : array, shape = [n_features] or [n_targets, n_features]
        Weight vector(s).

    `alpha_` : float
        Estimated regularization parameter.

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
        Alpha corresponds to (2*C)^-1 in other linear models such as
        LogisticRegression or LinearSVC.

    fit_intercept : boolean
        Whether to calculate the intercept for this model. If set
        to false, no intercept will be used in calculations
        (e.g. data is expected to be already centered).

    normalize : boolean, optional, default False
        If True, the regressors X will be normalized before regression.

    score_func: callable, optional
        function that takes 2 arguments and compares them in
        order to evaluate the performance of prediction (big is good)
        if None is passed, the score of the estimator is maximized

    loss_func: callable, optional
        function that takes 2 arguments and compares them in
        order to evaluate the performance of prediction (small is good)
        if None is passed, the score of the estimator is maximized

    cv : cross-validation generator, optional
        If None, Generalized Cross-Validation (efficient Leave-One-Out)
        will be used.

    class_weight : dict, optional
        Weights associated with classes in the form
        {class_label : weight}. If not given, all classes are
        supposed to have weight one.

    Attributes
    ----------
    `cv_values_` : array, shape = [n_samples, n_alphas] or \
    shape = [n_samples, n_responses, n_alphas], optional
        Cross-validation values for each alpha (if `store_cv_values=True` and
    `cv=None`). After `fit()` has been called, this attribute will contain \
    the mean squared errors (by default) or the values of the \
    `{loss,score}_func` function (if provided in the constructor).

    `coef_` : array, shape = [n_features] or [n_targets, n_features]
        Weight vector(s).

    `alpha_` : float
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
                 normalize=False, score_func=None, loss_func=None, cv=None,
                 class_weight=None):
        super(RidgeClassifierCV, self).__init__(
            alphas=alphas, fit_intercept=fit_intercept, normalize=normalize,
            score_func=score_func, loss_func=loss_func, cv=cv)
        self.class_weight = class_weight

    def fit(self, X, y, sample_weight=1.0, class_weight=None):
        """Fit the ridge classifier.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]
            Training vectors, where n_samples is the number of samples
            and n_features is the number of features.

        y : array-like, shape = [n_samples]
            Target values.

        sample_weight : float or numpy array of shape [n_samples]
            Sample weight

        Returns
        -------
        self : object
            Returns self.
        """
        if self.class_weight is not None:
            get_cw = self.class_weight.get
            sample_weight = (sample_weight
                             * np.array([get_cw(k, 1.0) for k in y]))
        self._label_binarizer = LabelBinarizer(pos_label=1, neg_label=-1)
        Y = self._label_binarizer.fit_transform(y)
        _BaseRidgeCV.fit(self, X, Y, sample_weight=sample_weight)
        return self

    @property
    def classes_(self):
        return self._label_binarizer.classes_
