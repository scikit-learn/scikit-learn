"""C-means clustering"""

import numpy as np

from sklearn.cluster._c_means import (
    _centers_probabilistic,
    _memberships_probabilistic,
    _centers_possibilistic,
    _memberships_possibilistic,
)
from ..base import BaseEstimator, ClusterMixin, TransformerMixin
from ..externals.six import string_types
from ..metrics.pairwise import euclidean_distances
from ..utils import as_float_array
from ..utils import check_array
from ..utils import check_random_state


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


def c_means(X, n_clusters, m=2, n_init=10, max_iter=300, init='random',
            verbose=False, tol=1e-4, random_state=None, algorithm="auto",
            return_n_iter=False, copy_x=True):
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

    membership_best, inertia_best, centers_best, n_iter_best = (
        None, None, None, None
    )

    if algorithm == "probabilistic":
        cmeans_single = _cmeans_single_probabilistic
    elif algorithm == "possibilistic":
        cmeans_single = _cmeans_single_possibilistic
    elif algorithm == "auto":
        cmeans_single = _cmeans_single_probabilistic
    else:
        cmeans_single = _cmeans_single_probabilistic

    for it in range(n_init):
        membership, inertia, centers, n_iter_ = cmeans_single(
            X, n_clusters, max_iter=max_iter, init=init, tol=tol,
            random_state=random_state)
        if inertia_best is None or inertia < inertia_best:
            membership_best = membership
            centers_best = centers
            inertia_best = inertia
            n_iter_best = n_iter_

    if not copy_x:
        X += X_mean

    return membership_best, centers_best, inertia_best, n_iter_best


def _cmeans_single_probabilistic(X, n_clusters, m=2., max_iter=300,
                                 init='random', random_state=None, tol=1e-4):
    random_state = check_random_state(random_state)
    memberships_best, inertia_best, centers_best = None, None, None

    centers = _init_centroids(X, n_clusters, init, random_state=random_state)
    inertia = np.infty

    for i in range(max_iter):
        inertia_old = inertia
        distances = euclidean_distances(X, centers)
        memberships = _memberships_probabilistic(distances, m)
        inertia = np.sum(memberships ** m * distances)
        centers = _centers_probabilistic(X, memberships, m)

        if inertia_best is None or inertia < inertia_best:
            memberships_best = memberships.copy()
            centers_best = centers.copy()
            inertia_best = inertia

        if inertia - inertia_old < tol:
            break

    return memberships_best, inertia_best, centers_best, i+1


def _cmeans_single_possibilistic(X, n_clusters, m=2., max_iter=300,
                                 init='probabilistic', random_state=None,
                                 tol=1e-4):
    random_state = check_random_state(random_state)
    memberships_best, inertia_best, centers_best = None, None, None

    # Initialize using a probabilistic run
    memberships, centers, inertia, _ = c_means(X,
                                               n_clusters=n_clusters,
                                               random_state=random_state)

    distances = euclidean_distances(X, centers)
    weights = np.sum(memberships ** m * distances,
                     axis=0) / np.sum(memberships ** m, axis=0)

    for i in range(max_iter):
        inertia_old = inertia
        distances = euclidean_distances(X, centers)
        memberships = _memberships_possibilistic(distances, weights, m)
        inertia = np.sum(memberships ** m * distances)
        centers = _centers_possibilistic(X, memberships, m)

        if inertia_best is None or inertia < inertia_best:
            memberships_best = memberships.copy()
            centers_best = centers.copy()
            inertia_best = inertia

        if inertia - inertia_old < tol:
            break

    return memberships_best, inertia_best, centers_best, i+1


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

    def __init__(self, n_clusters=8, m=2, n_init=10, max_iter=300, init='random',
                 tol=1e-4, random_state=None, algorithm='auto', copy_x=True):

        self.n_clusters = n_clusters
        self.m = m
        self.max_iter = max_iter
        self.init = init
        self.tol = tol
        self.n_init = n_init
        self.random_state = random_state
        self.algorithm = algorithm
        self.copy_x = copy_x

    def _check_fit_data(self, X):
        """Verify that the number of samples given is larger than k"""
        X = check_array(X, accept_sparse='csr', dtype=[np.float64, np.float32])
        if X.shape[0] < self.n_clusters:
            raise ValueError("n_samples=%d should be >= n_clusters=%d" % (
                X.shape[0], self.n_clusters))
        return X

    def fit(self, X, y=None):
        random_state = check_random_state(self.random_state)
        X = self._check_fit_data(X)

        self.memberships_, self.centers_, self.inertia_, self.n_iter_ = c_means(
            X, n_clusters=self.n_clusters, m=self.m, n_init=self.n_init,
            max_iter=self.max_iter, init=self.init, tol=self.tol,
            random_state=random_state, algorithm=self.algorithm,
            copy_x=self.copy_x)
        return self

    @property
    def labels_(self):
        return np.argmax(self.memberships_, axis=1)
