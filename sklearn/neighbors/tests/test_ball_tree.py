import itertools

import numpy as np
import pytest
from numpy.testing import assert_array_almost_equal
from sklearn.neighbors._ball_tree import BallTree
from sklearn.neighbors import DistanceMetric
from sklearn.utils import check_random_state

rng = np.random.RandomState(10)
V_mahalanobis = rng.rand(3, 3)
V_mahalanobis = np.dot(V_mahalanobis, V_mahalanobis.T)

DIMENSION = 3

METRICS = {'euclidean': {},
           'manhattan': {},
           'minkowski': dict(p=3),
           'chebyshev': {},
           'seuclidean': dict(V=rng.random_sample(DIMENSION)),
           'wminkowski': dict(p=3, w=rng.random_sample(DIMENSION)),
           'mahalanobis': dict(V=V_mahalanobis)}

DISCRETE_METRICS = ['hamming',
                    'canberra',
                    'braycurtis']

BOOLEAN_METRICS = ['matching', 'jaccard', 'dice', 'kulsinski',
                   'rogerstanimoto', 'russellrao', 'sokalmichener',
                   'sokalsneath']


def brute_force_neighbors(X, Y, k, metric, **kwargs):
    D = DistanceMetric.get_metric(metric, **kwargs).pairwise(Y, X)
    ind = np.argsort(D, axis=1)[:, :k]
    dist = D[np.arange(Y.shape[0])[:, None], ind]
    return dist, ind


@pytest.mark.parametrize('metric',
                         itertools.chain(BOOLEAN_METRICS, DISCRETE_METRICS))
def test_ball_tree_query_metrics(metric):
    rng = check_random_state(0)
    if metric in BOOLEAN_METRICS:
        X = rng.random_sample((40, 10)).round(0)
        Y = rng.random_sample((10, 10)).round(0)
    elif metric in DISCRETE_METRICS:
        X = (4 * rng.random_sample((40, 10))).round(0)
        Y = (4 * rng.random_sample((10, 10))).round(0)

    k = 5

    bt = BallTree(X, leaf_size=1, metric=metric)
    dist1, ind1 = bt.query(Y, k)
    dist2, ind2 = brute_force_neighbors(X, Y, k, metric)
    assert_array_almost_equal(dist1, dist2)


def test_query_haversine():
    rng = check_random_state(0)
    X = 2 * np.pi * rng.random_sample((40, 2))
    bt = BallTree(X, leaf_size=1, metric='haversine')
    dist1, ind1 = bt.query(X, k=5)
    dist2, ind2 = brute_force_neighbors(X, X, k=5, metric='haversine')

    assert_array_almost_equal(dist1, dist2)
    assert_array_almost_equal(ind1, ind2)
