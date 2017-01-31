"""C-means clustering"""

import numpy as np

from ..base import BaseEstimator, ClusterMixin, TransformerMixin
from ..utils import check_array
from ..utils import check_random_state
from ..utils import as_float_array
from ..externals.six import string_types


def _validate_center_shape(X, n_centers, centers):
    """Check if centers is compatible with X and n_centers"""
    if len(centers) != n_centers:
        raise ValueError('The shape of the initial centers (%s) '
                         'does not match the number of clusters %i'
                         % (centers.shape, n_centers))
    if centers.shape[1] != X.shape[1]:
        raise ValueError(
            "The number of features of the initial centers %s "
            "does not match the number of features of the data %s."
            % (centers.shape[1], X.shape[1]))


def c_means(X, n_clusters, n_init=10, max_iter=300, verbose=False, tol=1e-4,
            random_state=None, algorithm="auto", return_n_iter=False,
            copy_x=True):
    if n_init <= 0:
        raise ValueError('Number of initializations should be a positive'
                         ' number, got {:d} instead.'.format(n_init))

    random_state = check_random_state(random_state)

    if max_iter <= 0:
        raise ValueError('Number of iterations should be a positive number,'
                         ' got {:d} instead.'.format(max_iter))

    X = as_float_array(X, copy=copy_x)

    X_mean = X.mean(axis=0)
    X -= X_mean

    if not copy_x:
        X += X_mean

def _cmeans_single_probabilistic(X, n_clusters, max_iter=300, random_state=None,
                                 tol=1e-4):
    random_state = check_random_state(random_state)

    memberships_best, intertia_best, centers_best = None, None, None


def _init_centroids(X, k, init, random_state=None):
    random_state = check_random_state(random_state)
    n_samples = X.shape[0]

    if n_samples < k:
        raise ValueError(
            "`n_samples={:d} should be larger than n_clusters={:d}".format(
                n_samples, k
            )
        )

    if isinstance(init, string_types) and init == 'random':
        seeds = random_state.permutation(n_samples)[:k]
        centers = X[seeds]
    else:
        raise ValueError("`init` should be 'random'. '{}' (type '{}') "
                         "was passed.".format(init, type(init)))

    _validate_center_shape(X, k, centers)
    return centers


class CMeans(BaseEstimator, ClusterMixin, TransformerMixin):

    def __init__(self, n_clusters=8, n_init=10, max_iter=300, tol=1e-4,
                 random_state=None, algorithm='auto', copy_x=True):

        self.n_clusters = n_clusters
        self.max_iter = max_iter
        self.tol = tol
        self.n_init = n_init
        self.random_state = random_state
        self.algorithm = algorithm
        self.copy_x = copy_x

    def fit(self, X, y=None):
        c_means(
            X, n_clusters=self.n_clusters, n_init=self.n_init,
            max_iter=self.max_iter, tol=self.tol,
            random_state=self.random_state, algorithm=self.algorithm,
            copy_x=self.copy_x)
        return self
