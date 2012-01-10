"""
Ridge regression
"""

# Author: Mathieu Blondel <mathieu@mblondel.org>
# License: Simplified BSD

import numpy as np

from .base import LinearModel
from ..utils.extmath import safe_sparse_dot
from ..utils import safe_asarray
from ..preprocessing import LabelBinarizer
from ..grid_search import GridSearchCV


def _solve(A, b, solver, tol):
    # helper method for ridge_regression, A is symmetric positive

    if solver == 'auto':
        if hasattr(A, 'todense'):
            solver = 'sparse_cg'
        else:
            solver = 'dense_cholesky'

    if solver == 'sparse_cg':
        if b.ndim < 2:
            from scipy.sparse import linalg as sp_linalg
            sol, error = sp_linalg.cg(A, b, tol=tol)
            if error:
                raise ValueError("Failed with error code %d" % error)
            return sol
        else:
            # sparse_cg cannot handle a 2-d b.
            sol = []
            for j in range(b.shape[1]):
                sol.append(_solve(A, b[:, j], solver="sparse_cg", tol=tol))
            return np.array(sol).T

    elif solver == 'dense_cholesky':
        from scipy import linalg
        if hasattr(A, 'todense'):
            A = A.todense()
        return linalg.solve(A, b, sym_pos=True, overwrite_a=True)
    else:
        raise NotImplementedError('Solver %s not implemented' % solver)


def ridge_regression(X, y, alpha, sample_weight=1.0, solver='auto', tol=1e-3):
    """Solve the ridge equation by the method of normal equations.

    Parameters
    ----------
    X : {array-like, sparse matrix}, shape = [n_samples, n_features]
        Training data

    y : array-like, shape = [n_samples] or [n_samples, n_responses]
        Target values

    sample_weight : float or numpy array of shape [n_samples]
        Individual weights for each sample

    solver : {'auto', 'dense_cholesky', 'sparse_cg'}, optional
        Solver to use in the computational routines. 'delse_cholesky'
        will use the standard scipy.linalg.solve function, 'sparse_cg'
        will use the a conjugate gradient solver as found in
        scipy.sparse.linalg.cg while 'auto' will chose the most
        appropiate depending on the matrix X.

    tol: float
        Precision of the solution.

    Returns
    -------
    coef: array, shape = [n_features] or [n_responses, n_features]
        Weight vector(s).

    Notes
    -----
    This function won't compute the intercept.
    """

    n_samples, n_features = X.shape
    is_sparse = False

    if hasattr(X, 'todense'):  # lazy import of scipy.sparse
        from scipy import sparse
        is_sparse = sparse.issparse(X)

    if is_sparse:
        if n_features > n_samples or \
           isinstance(sample_weight, np.ndarray) or \
           sample_weight != 1.0:

            I = sparse.lil_matrix((n_samples, n_samples))
            I.setdiag(np.ones(n_samples) * alpha * sample_weight)
            c = _solve(X * X.T + I, y, solver, tol)
            coef = X.T * c
        else:
            I = sparse.lil_matrix((n_features, n_features))
            I.setdiag(np.ones(n_features) * alpha)
            coef = _solve(X.T * X + I, X.T * y, solver, tol)
    else:
        if n_features > n_samples or \
           isinstance(sample_weight, np.ndarray) or \
           sample_weight != 1.0:

            # kernel ridge
            # w = X.T * inv(X X^t + alpha*Id) y
            A = np.dot(X, X.T)
            A.flat[::n_samples + 1] += alpha * sample_weight
            coef = np.dot(X.T, _solve(A, y, solver, tol))
        else:
            # ridge
            # w = inv(X^t X + alpha*Id) * X.T y
            A = np.dot(X.T, X)
            A.flat[::n_features + 1] += alpha
            coef = _solve(A, np.dot(X.T, y), solver, tol)

    return coef.T


class Ridge(LinearModel):
    """Linear least squares with l2 regularization.

    This model solves a regression model where the loss function is
    the linear least squares function and regularization is given by
    the l2-norm. Also known as Ridge Regression or Tikhonov regularization.

    Parameters
    ----------
    alpha : float
        Small positive values of alpha improve the conditioning of the
        problem and reduce the variance of the estimates.
        Alpha corresponds to (2*C)^-1 in other linear models such as
        LogisticRegression or LinearSVC.

    fit_intercept : boolean
        Whether to calculate the intercept for this model. If set
        to false, no intercept will be used in calculations
        (e.g. data is expected to be already centered).

    normalize : boolean, optional
        If True, the regressors X are normalized

    copy_X : boolean, optional, default True
        If True, X will be copied; else, it may be overwritten.

    tol: float
        Precision of the solution.

    Attributes
    ----------

    `coef_` : array, shape = [n_features] or [n_responses, n_features]
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
    Ridge(alpha=1.0, copy_X=True, fit_intercept=True, normalize=False,
       tol=0.001)
    """

    def __init__(self, alpha=1.0, fit_intercept=True, normalize=False,
            copy_X=True, tol=1e-3):
        self.alpha = alpha
        self.fit_intercept = fit_intercept
        self.normalize = normalize
        self.copy_X = copy_X
        self.tol = tol

    def fit(self, X, y, sample_weight=1.0, solver='auto'):
        """Fit Ridge regression model

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape = [n_samples, n_features]
            Training data

        y : array-like, shape = [n_samples] or [n_samples, n_responses]
            Target values

        sample_weight : float or numpy array of shape [n_samples]
            Individual weights for each sample

        solver : {'auto', 'dense_cholesky', 'sparse_cg'}
            Solver to use in the computational
            routines. 'delse_cholesky' will use the standard
            scipy.linalg.solve function, 'sparse_cg' will use the a
            conjugate gradient solver as found in
            scipy.sparse.linalg.cg while 'auto' will chose the most
            appropiate depending on the matrix X.

        Returns
        -------
        self : returns an instance of self.
        """
        X = safe_asarray(X, dtype=np.float)
        y = np.asarray(y, dtype=np.float)

        X, y, X_mean, y_mean, X_std = \
           self._center_data(X, y, self.fit_intercept,
                   self.normalize, self.copy_X)

        self.coef_ = ridge_regression(X, y, self.alpha, sample_weight,
                                      solver, self.tol)
        self._set_intercept(X_mean, y_mean, X_std)
        return self


class RidgeClassifier(Ridge):
    """Classifier using Ridge regression

    Parameters
    ----------
    alpha : float
        Small positive values of alpha improve the conditioning of the
        problem and reduce the variance of the estimates.
        Alpha corresponds to (2*C)^-1 in other linear models such as
        LogisticRegression or LinearSVC.

    fit_intercept : boolean
        Whether to calculate the intercept for this model. If set
        to false, no intercept will be used in calculations
        (e.g. data is expected to be already centered).

    normalize : boolean, optional
        If True, the regressors X are normalized

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
    a one-versus-all approach.
    """

    def fit(self, X, y, solver='auto'):
        """Fit Ridge regression model.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape = [n_samples,n_features]
            Training data

        y : array-like, shape = [n_samples]
            Target values

        solver : {'auto', 'dense_cholesky', 'sparse_cg'}
            Solver to use in the computational
            routines. 'delse_cholesky' will use the standard
            scipy.linalg.solve function, 'sparse_cg' will use the a
            conjugate gradient solver as found in
            scipy.sparse.linalg.cg while 'auto' will chose the most
            appropiate depending on the matrix X.

        Returns
        -------
        self : returns an instance of self.
        """
        self.label_binarizer = LabelBinarizer()
        Y = self.label_binarizer.fit_transform(y)
        Ridge.fit(self, X, Y, solver=solver)
        return self

    def decision_function(self, X):
        return Ridge.decision_function(self, X)

    def predict(self, X):
        """Predict target values according to the fitted model.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]

        Returns
        -------
        y : array, shape = [n_samples]
        """
        Y = self.decision_function(X)
        return self.label_binarizer.inverse_transform(Y)


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

    **References**:
    http://cbcl.mit.edu/projects/cbcl/publications/ps/MIT-CSAIL-TR-2007-025.pdf
    http://www.mit.edu/~9.520/spring07/Classes/rlsslides.pdf
    """

    def __init__(self, alphas=[0.1, 1.0, 10.0], fit_intercept=True,
            normalize=False, score_func=None, loss_func=None, copy_X=True):
        self.alphas = np.asarray(alphas)
        self.fit_intercept = fit_intercept
        self.normalize = normalize
        self.score_func = score_func
        self.loss_func = loss_func
        self.copy_X = copy_X

    def _pre_compute(self, X, y):
        # even if X is very sparse, K is usually very dense
        K = safe_sparse_dot(X, X.T, dense_output=True)
        from scipy import linalg
        v, Q = linalg.eigh(K)
        return K, v, Q

    def _errors(self, v, Q, y, alpha):
        G = np.dot(np.dot(Q, np.diag(1.0 / (v + alpha))), Q.T)
        c = np.dot(G, y)
        G_diag = np.diag(G)
        # handle case when y is 2-d
        G_diag = G_diag if len(y.shape) == 1 else G_diag[:, np.newaxis]
        return (c / G_diag) ** 2, c

    def _values(self, K, v, Q, y, alpha):
        n_samples = y.shape[0]

        G = np.dot(np.dot(Q, np.diag(1.0 / (v + alpha))), Q.T)
        c = np.dot(G, y)
        KG = np.dot(K, G)
        #KG = np.dot(np.dot(Q, np.diag(v / (v + alpha))), Q.T)
        KG_diag = np.diag(KG)

        denom = np.ones(n_samples) - KG_diag
        if len(y.shape) == 2:
            # handle case when y is 2-d
            KG_diag = KG_diag[:, np.newaxis]
            denom = denom[:, np.newaxis]

        num = np.dot(KG, y) - KG_diag * y

        return num / denom, c

    def fit(self, X, y, sample_weight=1.0):
        """Fit Ridge regression model

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape = [n_samples, n_features]
            Training data

        y : array-like, shape = [n_samples] or [n_samples, n_responses]
            Target values

        sample_weight : float or array-like of shape [n_samples]
            Sample weight

        Returns
        -------
        self : Returns self.
        """
        X = safe_asarray(X, dtype=np.float)
        y = np.asarray(y, dtype=np.float)

        n_samples = X.shape[0]

        X, y, X_mean, y_mean, X_std = LinearModel._center_data(X, y,
                self.fit_intercept, self.normalize, self.copy_X)

        K, v, Q = self._pre_compute(X, y)
        n_y = 1 if len(y.shape) == 1 else y.shape[1]
        M = np.zeros((n_samples * n_y, len(self.alphas)))
        C = []

        error = self.score_func is None and self.loss_func is None

        for i, alpha in enumerate(self.alphas):
            if error:
                out, c = self._errors(v, Q, y, sample_weight * alpha)
            else:
                out, c = self._values(K, v, Q, y, sample_weight * alpha)
            M[:, i] = out.ravel()
            C.append(c)

        if error:
            best = M.mean(axis=0).argmin()
        else:
            func = self.score_func if self.score_func else self.loss_func
            out = [func(y.ravel(), M[:, i]) for i in range(len(self.alphas))]
            best = np.argmax(out) if self.score_func else np.argmin(out)

        self.best_alpha = self.alphas[best]
        self.dual_coef_ = C[best]
        self.coef_ = safe_sparse_dot(self.dual_coef_.T, X)

        self._set_intercept(X_mean, y_mean, X_std)

        return self


class RidgeCV(LinearModel):
    """Ridge regression with built-in cross-validation.

    By default, it performs Generalized Cross-Validation, which is a form of
    efficient Leave-One-Out cross-validation. Currently, only the n_features >
    n_samples case is handled efficiently.

    Parameters
    ----------
    alphas: numpy array of shape [n_alpha]
        Array of alpha values to try.
        Small positive values of alpha improve the conditioning of the
        problem and reduce the variance of the estimates.
        Alpha corresponds to (2*C)^-1 in other linear models such as
        LogisticRegression or LinearSVC.

    fit_intercept : boolean
        Whether to calculate the intercept for this model. If set
        to false, no intercept will be used in calculations
        (e.g. data is expected to be already centered).

    normalize : boolean, optional
        If True, the regressors X are normalized

    loss_func: callable, optional
        function that takes 2 arguments and compares them in
        order to evaluate the performance of prediciton (small is good)
        if None is passed, the score of the estimator is maximized

    score_func: callable, optional
        function that takes 2 arguments and compares them in
        order to evaluate the performance of prediciton (big is good)
        if None is passed, the score of the estimator is maximized

    See also
    --------
    Ridge, RidgeClassifierCV
    """

    def __init__(self, alphas=np.array([0.1, 1.0, 10.0]), fit_intercept=True,
                   normalize=False, score_func=None, loss_func=None, cv=None):
        self.alphas = alphas
        self.fit_intercept = fit_intercept
        self.normalize = normalize
        self.score_func = score_func
        self.loss_func = loss_func
        self.cv = cv

    def fit(self, X, y, sample_weight=1.0):
        """Fit Ridge regression model

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]
            Training data

        y : array-like, shape = [n_samples] or [n_samples, n_responses]
            Target values

        sample_weight : float or array-like of shape [n_samples]
            Sample weight

        cv : cross-validation generator, optional
            If None, Generalized Cross-Validationn (efficient Leave-One-Out)
            will be used.

        Returns
        -------
        self : Returns self.
        """
        if self.cv is None:
            estimator = _RidgeGCV(self.alphas, self.fit_intercept,
                                  self.score_func, self.loss_func)
            estimator.fit(X, y, sample_weight=sample_weight)
            self.best_alpha = estimator.best_alpha
        else:
            parameters = {'alpha': self.alphas}
            # FIXME: sample_weight must be split into training/validation data
            #        too!
            #fit_params = {'sample_weight' : sample_weight}
            fit_params = {}
            gs = GridSearchCV(Ridge(fit_intercept=self.fit_intercept),
                              parameters, fit_params=fit_params, cv=self.cv)
            gs.fit(X, y)
            estimator = gs.best_estimator_
            self.best_alpha = gs.best_estimator_.alpha

        self.coef_ = estimator.coef_
        self.intercept_ = estimator.intercept_

        return self


class RidgeClassifierCV(RidgeCV):

    def fit(self, X, y, sample_weight=1.0, class_weight=None):
        """Fit the ridge classifier.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]
            Training vectors, where n_samples is the number of samples
            and n_features is the number of features.

        y : array-like, shape = [n_samples]
            Target values.

        class_weight : dict, optional
            Weights associated with classes in the form
            {class_label : weight}. If not given, all classes are
            supposed to have weight one.

        sample_weight : float or numpy array of shape [n_samples]
            Sample weight

        Returns
        -------
        self : object
            Returns self.
        """
        if class_weight is None:
            class_weight = {}
        sample_weight2 = np.array([class_weight.get(k, 1.0) for k in y])
        self.label_binarizer = LabelBinarizer()
        Y = self.label_binarizer.fit_transform(y)
        RidgeCV.fit(self, X, Y,
                    sample_weight=sample_weight * sample_weight2)
        return self

    def decision_function(self, X):
        return RidgeCV.decision_function(self, X)

    def predict(self, X):
        """Predict target values according to the fitted model.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]

        Returns
        -------
        y : array, shape = [n_samples]
        """
        Y = self.decision_function(X)
        return self.label_binarizer.inverse_transform(Y)
