"""C-means clustering"""

import numpy as np

from sklearn.cluster._c_means import Probabilistic, Possibilistic
from ..base import BaseEstimator, ClusterMixin, TransformerMixin
from ..externals.six import string_types
from ..utils import as_float_array
from ..utils import check_array
from ..utils import check_random_state
from ..utils.validation import FLOAT_DTYPES
from ..utils.validation import check_is_fitted

from sklearn.metrics import euclidean_distances


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

    results_best = ()
    inertia_best = None

    if algorithm == "probabilistic":
        single = probabilistic_single
    elif algorithm == "possibilistic":
        single = possibilistic_single
    elif algorithm == "auto":
        single = probabilistic_single
    else:
        raise NotImplementedError(
            "{} is not an implemented algorithm.".format(algorithm))

    for it in range(n_init):
        results = single(
            X, n_clusters, max_iter=max_iter, init=init, tol=tol,
            random_state=random_state)
        inertia = results['inertia']
        if inertia_best is None or inertia < inertia_best:
            results_best = results
            inertia_best = results['inertia']
    if not copy_x:
        X += X_mean
    results_best['centers'] += X_mean

    return results_best


def probabilistic_single(X, n_clusters, m=2., max_iter=300,
                         init='random', random_state=None, tol=1e-4):
    random_state = check_random_state(random_state)
    memberships_best, inertia_best, centers_best = None, None, None

    centers = _init_centroids(X, n_clusters, init,
                              random_state=random_state)
    inertia = np.infty

    for i in range(max_iter):
        inertia_old = inertia
        distances = euclidean_distances(X, centers)
        memberships = Probabilistic.memberships(distances, m)
        inertia = np.sum(memberships ** m * distances)
        centers = Probabilistic.centers(X, memberships, m)

        if inertia_best is None or inertia < inertia_best:
            memberships_best = memberships.copy()
            centers_best = centers.copy()
            inertia_best = inertia

        if inertia - inertia_old < tol:
            break

    results = {
        'memberships': memberships_best,
        'inertia': inertia_best,
        'centers': centers_best,
        'n_iter': i + 1,
    }

    return results


def possibilistic_single(X, n_clusters, m=2., max_iter=300,
                         init='probabilistic', random_state=None,
                         tol=1e-4):
    random_state = check_random_state(random_state)
    memberships_best, inertia_best, centers_best = None, None, None

    # Initialize using a probabilistic run
    results_init = c_means(
        X, n_clusters=n_clusters, algorithm="probabilistic",
        random_state=random_state)

    memberships = results_init['memberships']
    centers = results_init['centers']
    inertia = results_init['inertia']

    distances = euclidean_distances(X, centers)
    weights = np.sum(memberships ** m * distances,
                     axis=0) / np.sum(memberships ** m, axis=0)

    for i in range(max_iter):
        inertia_old = inertia
        distances = euclidean_distances(X, centers)
        memberships = Possibilistic.memberships(distances, m)
        inertia = np.sum(memberships ** m * distances)
        centers = Possibilistic.centers(X, memberships, m)

        if inertia_best is None or inertia < inertia_best:
            memberships_best = memberships.copy()
            centers_best = centers.copy()
            inertia_best = inertia

        if inertia - inertia_old < tol:
            break

    results = {
        'memberships': memberships_best,
        'inertia': inertia_best,
        'centers': centers_best,
        'weights': weights,
        'n_iter': i + 1,
    }

    return results


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

    def _check_test_data(self, X):
        X = check_array(X, accept_sparse='csr', dtype=FLOAT_DTYPES)
        n_samples, n_features = X.shape
        expected_n_features = self.centers_.shape[1]
        if not n_features == expected_n_features:
            raise ValueError("Incorrect number of features. "
                             "Got %d features, expected %d" % (
                                 n_features, expected_n_features))
        return X

    @property
    def labels_(self):
        return np.argmax(self.memberships_, axis=1)


class ProbabilisticCMeans(CMeans):

    def fit(self, X, y=None):
        random_state = check_random_state(self.random_state)
        X = self._check_fit_data(X)

        results = c_means(
            X, n_clusters=self.n_clusters, m=self.m, n_init=self.n_init,
            max_iter=self.max_iter, init=self.init, tol=self.tol,
            random_state=random_state, algorithm="probabilistic",
            copy_x=self.copy_x)
        self.memberships_ = results['memberships']
        self.centers_ = results['centers']
        self.inertia_ = results['inertia']
        self.n_iter_ = results['n_iter']
        return self

    def predict(self, X):
        """Predict the memberships for each sample in X.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            Data to predict

        Returns
        -------
        memberships : ndarray, shape (n_samples, n_clusters)
            Membership of each sample to each cluster.

        """
        check_is_fitted(self, 'centers_')
        X = self._check_test_data(X)
        distances = euclidean_distances(X, self.centers_)
        return Probabilistic.memberships(distances, m=self.m)


class PossibilisticCMeans(CMeans):

    def fit(self, X, y=None):
        random_state = check_random_state(self.random_state)
        X = self._check_fit_data(X)

        results = c_means(
            X, n_clusters=self.n_clusters, m=self.m, n_init=self.n_init,
            max_iter=self.max_iter, init=self.init, tol=self.tol,
            random_state=random_state, algorithm="possibilistic",
            copy_x=self.copy_x)
        self.memberships_ = results['memberships']
        self.centers_ = results['centers']
        self.inertia_ = results['inertia']
        self.n_iter_ = results['n_iter']
        self.weights_ = results['weights']
        return self

    def predict(self, X):
        """Predict the memberships for each sample in X.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            Data to predict

        Returns
        -------
        memberships : ndarray, shape (n_samples, n_clusters)
            Membership of each sample to each cluster.

        """
        check_is_fitted(self, 'centers_')
        X = self._check_test_data(X)
        distances = euclidean_distances(X, self.centers_)
        return Possibilistic.memberships(distances, self.weights_, m=self.m)
