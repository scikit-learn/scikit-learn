"""
Ridge regression
"""

import numpy as np
from scipy import linalg

from .base import LinearModel


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

    def fit(self, X, y, **params):
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

        X = np.asanyarray(X, dtype=np.float)
        y = np.asanyarray(y, dtype=np.float)

        n_samples, n_features = X.shape

        X, y, Xmean, ymean = LinearModel._center_data(X, y, self.fit_intercept)

        if n_samples > n_features:
            # w = inv(X^t X + alpha*Id) * X.T y
            A = np.dot(X.T, X)
            A.flat[::n_features+1] += self.alpha
            self.coef_ = linalg.solve(A, np.dot(X.T, y),
                                      overwrite_a=True, sym_pos=True)
        else:
            # w = X.T * inv(X X^t + alpha*Id) y
            A = np.dot(X, X.T)
            A.flat[::n_samples+1] += self.alpha
            self.coef_ = np.dot(X.T, linalg.solve(A, y, overwrite_a=True,
                                                  sym_pos=True))

        self._set_intercept(Xmean, ymean)
        return self

# Author: Mathieu Blondel <mathieu@mblondel.org>
# License: Simplified BSD

class RidgeLOO(LinearModel):
    """Ridge regression with built-in efficient
       Leave-One-Out cross-validation

    Explanations
    ------------

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
        K = np.dot(X, X.T)
        v, Q = linalg.eigh(K)
        return K, v, Q

    def _errors(self, v, Q, y, alpha):
        G = np.dot(np.dot(Q, np.diag(1.0 / (v + alpha))), Q.T)
        c = np.dot(G, y)
        G_diag = np.diag(G)
        # handle case when y is 2-d
        G_diag = G_diag if len(y.shape) == 1 else G_diag[:,np.newaxis]
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
            KG_diag = KG_diag[:,np.newaxis]
            denom = denom[:,np.newaxis]

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
            M[:,i] = out.ravel()
            C.append(c)

        if error:
            best = M.mean(axis=0).argmin()
        else:
            func = self.score_func if self.score_func else self.loss_func
            out = [func(y.ravel(), M[:,i]) for i in range(len(self.alphas))]
            best = np.argmax(out) if self.score_func else np.argmin(out)

        self.best_alpha = self.alphas[best]
        self.dual_coef_ = C[best]
        self._set_coef_(X)

        self._set_intercept(Xmean, ymean)

        return self

    def _set_coef_(self, X):
        self.coef_ = np.dot(X.T, self.dual_coef_)

