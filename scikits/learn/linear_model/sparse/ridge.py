
# Author: Mathieu Blondel <mathieu@mblondel.org>
# License: Simplified BSD

import numpy as np
import scipy.sparse as sp
from scipy import linalg
from scipy.sparse import linalg as sp_linalg

from ..base import LinearModel
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

