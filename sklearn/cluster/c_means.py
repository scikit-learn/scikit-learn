"""C-means clustering"""

import numpy as np

from ..base import BaseEstimator, ClusterMixin, TransformerMixin
from ..utils import check_array

def c_means(X, n_clusters, n_init=10, max_iter=300, verbose=False, tol=1e-4,
            random_state=None, algorithm="auto", return_n_iter=False):
    if n_init <= 0:
        raise ValueError("Invalid number of initializations."
                         " n_init={:d}")


class CMeans(BaseEstimator, ClusterMixin, TransformerMixin):

    def __init__(self, n_clusters=8, n_init=10, max_iter=300, tol=1e-4,
                 random_state=None, algorithm='auto'):

        self.n_clusters = n_clusters
        self.max_iter = max_iter
        self.tol = tol
        self.n_init = n_init
        self.random_state = random_state
        self.algorithm = algorithm

    def fit(self, X, y=None):
        c_means(
            X, n_clusters=self.n_clusters, n_init=self.n_init,
            max_iter=self.max_iter, tol=self.tol,
            random_state=self.random_state, algorithm=self.algorithm)
        return self