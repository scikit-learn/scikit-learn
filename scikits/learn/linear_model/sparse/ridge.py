
import numpy as np
import scipy.sparse as sp
from scipy.sparse import linalg

from ..base import LinearModel

class Ridge(LinearModel):

    def __init__(self, alpha=1.0, fit_intercept=False):
        self.alpha = alpha
        self.intercept_ = 0

    def _check_error(self, error):
        if error:
            raise ValueError, "Failed with error code %d" % error

    def fit(self, X, y):
        n_samples, n_features = X.shape

        y = np.array(y, dtype=np.int32)

        if n_samples > n_features:
            I = sp.lil_matrix((n_features, n_features))
            I.setdiag(np.ones(n_features) * self.alpha)
            self.coef_, error = linalg.cg(X.T * X + I, X.T * y)
            self._check_error(error)
        else:
            I = sp.lil_matrix((n_samples, n_samples))
            I.setdiag(np.ones(n_samples) * self.alpha)
            sol, error = linalg.cg(X * X.T + I, y)
            self._check_error(error)
            self.coef_ = X.T * sol

        # FIXME: handle fit_intercept

        return self

    def predict(self, X):
        return X * self.coef_ + self.intercept_
