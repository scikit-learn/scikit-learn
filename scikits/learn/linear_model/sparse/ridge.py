
# Author: Mathieu Blondel <mathieu@mblondel.org>
# License: Simplified BSD

import numpy as np
import scipy.sparse as sp
from scipy import linalg
from scipy.sparse import linalg as sp_linalg

from ..base import LinearModel
from ..ridge import RidgeLOO as DenseRidgeLOO
from ...preprocessing import LabelBinarizer

class Ridge(LinearModel):

    def __init__(self, alpha=1.0, fit_intercept=False, solver="default"):
        self.alpha = alpha
        self.intercept_ = 0
        self.fit_intercept = fit_intercept
        self.solver = solver

    def _solve(self, A, b):
        if self.solver == "cg":
            # this solver cannot handle a 2-d b.
            sol, error = sp_linalg.cg(A, b)
            if error:
                raise ValueError, "Failed with error code %d" % error
            return sol
        else:
            # we are working with dense symmetric positive A
            return linalg.solve(A.todense(), b, sym_pos=True)

    def fit(self, X, y):
        X = sp.csr_matrix(X)
        n_samples, n_features = X.shape

        y = np.array(y)

        if n_samples > n_features:
            I = sp.lil_matrix((n_features, n_features))
            I.setdiag(np.ones(n_features) * self.alpha)
            self.coef_ = self._solve(X.T * X + I, X.T * y)
        else:
            I = sp.lil_matrix((n_samples, n_samples))
            I.setdiag(np.ones(n_samples) * self.alpha)
            c = self._solve(X * X.T + I, y)
            self.coef_ = X.T * c

        # FIXME: handle fit_intercept

        return self

    def predict(self, X):
        X = sp.csr_matrix(X)
        return X * self.coef_ + self.intercept_

class RidgeLOO(DenseRidgeLOO):

    def __init__(self, alphas=np.array([0.1, 1.0, 10.0]), fit_intercept=False,
                       score_func=None, loss_func=None):
        self.alphas = alphas
        self.fit_intercept = False # True not supported yet
        self.score_func = score_func
        self.loss_func = loss_func

    def _pre_compute(self, X, y):
        X = sp.csr_matrix(X)
        # the kernel matrix is dense
        K = (X * X.T).toarray()
        v, Q = linalg.eigh(K)
        return K, v, Q

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
        X = sp.csr_matrix(X)
        DenseRidgeLOO.fit(self, X, y, sample_weight)
        return self

    def _set_coef_(self, X):
        self.coef_ = X.T * self.dual_coef_

    def predict(self, X):
        X = sp.csr_matrix(X)
        return X * self.coef_ + self.intercept_

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

class RidgeClassifierLOO(RidgeLOO):

    def fit(self, X, y, class_weight={}):
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

        Returns
        -------
        self : object
            Returns self.
        """
        sample_weight = np.array([class_weight.get(k, 1.0) for k in y])
        self.lb = LabelBinarizer()
        Y = self.lb.fit_transform(y)
        RidgeLOO.fit(self, X, Y, sample_weight)
        return self

    def decision_function(self, X):
        return RidgeLOO.predict(self, X)

    def predict(self, X):
        Y = self.decision_function(X)
        return self.lb.inverse_transform(Y)

