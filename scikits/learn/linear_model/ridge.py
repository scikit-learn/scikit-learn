"""
Ridge regression
"""

# Author: Mathieu Blondel <mathieu@mblondel.org>
# License: Simplified BSD

import numpy as np
import scipy.sparse as sp
from scipy import linalg
from scipy.sparse import linalg as sp_linalg

from .base import LinearModel
from ..utils.extmath import safe_sparse_dot
from ..utils import safe_asanyarray
from ..preprocessing import LabelBinarizer
from ..grid_search import GridSearchCV


class Ridge(LinearModel):
    """
    Ridge regression.

    Parameters
    ----------
    alpha : float
        Small positive values of alpha improve the coditioning of the
        problem and reduce the variance of the estimates.

    fit_intercept : boolean
        wether to calculate the intercept for this model. If set
        to false, no intercept will be used in calculations
        (e.g. data is expected to be already centered).

    Examples
    --------
    >>> from scikits.learn.linear_model import Ridge
    >>> import numpy as np
    >>> n_samples, n_features = 10, 5
    >>> np.random.seed(0)
    >>> y = np.random.randn(n_samples)
    >>> X = np.random.randn(n_samples, n_features)
    >>> clf = Ridge(alpha=1.0)
    >>> clf.fit(X, y)
    Ridge(alpha=1.0, fit_intercept=True)
    """

    def __init__(self, alpha=1.0, fit_intercept=True):
        self.alpha = alpha
        self.fit_intercept = fit_intercept

    def fit(self, X, y, solver="default", **params):
        """
        Fit Ridge regression model

        Parameters
        ----------
        X : numpy array of shape [n_samples,n_features]
            Training data

        y : numpy array of shape [n_samples]
            Target values

        Returns
        -------
        self : returns an instance of self.
        """
        self._set_params(**params)
        self.solver = solver

        X = safe_asanyarray(X, dtype=np.float)
        y = np.asanyarray(y, dtype=np.float)

        X, y, Xmean, ymean = LinearModel._center_data(X, y, self.fit_intercept)

        if sp.issparse(X):
            self._solve_sparse(X, y)
        else:
            self._solve_dense(X, y)

        self._set_intercept(Xmean, ymean)

        return self

    def _solve_dense(self, X, y):
        n_samples, n_features = X.shape

        if n_samples > n_features:
            # w = inv(X^t X + alpha*Id) * X.T y
            A = np.dot(X.T, X)
            A.flat[::n_features + 1] += self.alpha
            self.coef_ = self._solve(A, np.dot(X.T, y))
        else:
            # w = X.T * inv(X X^t + alpha*Id) y
            A = np.dot(X, X.T)
            A.flat[::n_samples + 1] += self.alpha
            self.coef_ = np.dot(X.T, self._solve(A, y))

    def _solve_sparse(self, X, y):
        n_samples, n_features = X.shape

        if n_samples > n_features:
            I = sp.lil_matrix((n_features, n_features))
            I.setdiag(np.ones(n_features) * self.alpha)
            self.coef_ = self._solve(X.T * X + I, X.T * y)
        else:
            I = sp.lil_matrix((n_samples, n_samples))
            I.setdiag(np.ones(n_samples) * self.alpha)
            c = self._solve(X * X.T + I, y)
            self.coef_ = X.T * c

    def _solve(self, A, b):
        if self.solver == "cg":
            # this solver cannot handle a 2-d b.
            sol, error = sp_linalg.cg(A, b)
            if error:
                raise ValueError("Failed with error code %d" % error)
            return sol
        else:
            # we are working with dense symmetric positive A
            if sp.issparse(A):
                A = A.todense()
            return linalg.solve(A, b, sym_pos=True, overwrite_a=True)


class RidgeClassifier(Ridge):

    def fit(self, X, y):
        self.lb = LabelBinarizer()
        Y = self.lb.fit_transform(y)
        Ridge.fit(self, X, Y)
        return self

    def decision_function(self, X):
        return Ridge.predict(self, X)

    def predict(self, X):
        Y = self.decision_function(X)
        return self.lb.inverse_transform(Y)


class _RidgeLOO(LinearModel):
    """
    Ridge regression with built-in efficient Leave-One-Out cross-validation.

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

    Reference
    ---------

    http://www.mit.edu/~9.520/spring07/Classes/rlsslides.pdf
    """

    def __init__(self, alphas=np.array([0.1, 1.0, 10.0]), fit_intercept=True,
                       score_func=None, loss_func=None):
        self.alphas = alphas
        self.fit_intercept = fit_intercept
        self.score_func = score_func
        self.loss_func = loss_func

    def _pre_compute(self, X, y):
        # even if X is very sparse, K is usually very dense
        K = safe_sparse_dot(X, X.T, dense_output=True)
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
        X : numpy array of shape [n_samples, n_features]
            Training data

        y : numpy array of shape [n_samples] or [n_samples, n_responses]
            Target values

        sample_weight : float or numpy array of shape [n_samples]
            Sample weight

        Returns
        -------
        self : Returns self.
        """
        n_samples = X.shape[0]

        X, y, Xmean, ymean = LinearModel._center_data(X, y, self.fit_intercept)

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
        self.coef_ = safe_sparse_dot(X.T, self.dual_coef_)

        self._set_intercept(Xmean, ymean)

        return self


class RidgeCV(LinearModel):
    """
    Ridge regression with built-in cross-validation.

    Defaults to Leave-One-Out, which is handled efficiently.
    """

    def __init__(self, alphas=np.array([0.1, 1.0, 10.0]), fit_intercept=True,
                       score_func=None, loss_func=None):
        self.alphas = alphas
        self.fit_intercept = fit_intercept
        self.score_func = score_func
        self.loss_func = loss_func

    def fit(self, X, y, sample_weight=1.0, cv=None):
        """Fit Ridge regression model

        Parameters
        ----------
        X : numpy array of shape [n_samples, n_features]
            Training data

        y : numpy array of shape [n_samples] or [n_samples, n_responses]
            Target values

        sample_weight : float or numpy array of shape [n_samples]
            Sample weight

        cv : cross-validation generator, optional
            If None, Leave-One-Out will be used.

        Returns
        -------
        self : Returns self.
        """
        if cv is None:
            estimator = _RidgeLOO(self.alphas, self.fit_intercept,
                                  self.score_func, self.loss_func)
            estimator.fit(X, y, sample_weight)
        else:
            parameters = {'alpha': self.alphas}
            gs = GridSearchCV(Ridge(fit_intercept=self.fit_intercept),
                              parameters)
            # FIXME: need to implement sample_weight in Ridge
            gs.fit(X, y, cv=cv)
            estimator = gs.best_estimator

        self.coef_ = estimator.coef_
        self.intercept_ = estimator.intercept_

        return self


class RidgeClassifierCV(RidgeCV):

    def fit(self, X, y, sample_weight=1.0, class_weight={}, cv=None):
        """
        Fit the ridge classifier.

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
        sample_weight2 = np.array([class_weight.get(k, 1.0) for k in y])
        self.lb = LabelBinarizer()
        Y = self.lb.fit_transform(y)
        RidgeCV.fit(self, X, Y,
                    sample_weight=sample_weight * sample_weight2,
                    cv=cv)
        return self

    def decision_function(self, X):
        return RidgeCV.predict(self, X)

    def predict(self, X):
        Y = self.decision_function(X)
        return self.lb.inverse_transform(Y)
