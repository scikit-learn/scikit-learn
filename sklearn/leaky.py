import scipy.sparse as sp
import resource
from .base import BaseEstimator


class LeakyEstimator(BaseEstimator):
    def fit(self, X, y):
        if sp.issparse(X):
            self.whoops = X.todense()

    def predict(self, X):
        return (X.todense() * 2)[:, 0]

    def _more_tags(self):
        return {'multioutput_only': True}
