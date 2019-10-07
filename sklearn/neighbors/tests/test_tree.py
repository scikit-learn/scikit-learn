import itertools

import numpy as np
import pytest

from sklearn.neighbors.dist_metrics import DistanceMetric
from sklearn.neighbors.ball_tree import BallTree
from sklearn.neighbors.kd_tree import KDTree

from sklearn.utils import check_random_state
from sklearn.utils.testing import assert_allclose
from numpy.testing import assert_array_almost_equal

rng = np.random.RandomState(42)
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

KD_TREE_METRICS = ['euclidean', 'manhattan', 'chebyshev', 'minkowski']
BALL_TREE_METRICS = list(METRICS.keys())


def brute_force_neighbors(X, Y, k, metric, **kwargs):
    D = DistanceMetric.get_metric(metric, **kwargs).pairwise(Y, X)
    ind = np.argsort(D, axis=1)[:, :k]
    dist = D[np.arange(Y.shape[0])[:, None], ind]
    return dist, ind


@pytest.mark.parametrize(
        'Cls, metric',
        itertools.chain(
            [(KDTree, metric) for metric in KD_TREE_METRICS],
            [(BallTree, metric) for metric in BALL_TREE_METRICS]))
@pytest.mark.parametrize('k', (1, 3, 5))
@pytest.mark.parametrize('dualtree', (True, False))
@pytest.mark.parametrize('breadth_first', (True, False))
def test_kd_tree_query(Cls, metric, k, dualtree, breadth_first):
    rng = check_random_state(0)
    X = rng.random_sample((40, DIMENSION))
    Y = rng.random_sample((10, DIMENSION))

    kwargs = METRICS[metric]

    kdt = Cls(X, leaf_size=1, metric=metric, **kwargs)
    dist1, ind1 = kdt.query(Y, k, dualtree=dualtree,
                            breadth_first=breadth_first)
    dist2, ind2 = brute_force_neighbors(X, Y, k, metric, **kwargs)

    # don't check indices here: if there are any duplicate distances,
    # the indices may not match.  Distances should not have this problem.
    assert_array_almost_equal(dist1, dist2)
