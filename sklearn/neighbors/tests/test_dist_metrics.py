import itertools
import pickle

import numpy as np
from numpy.testing import assert_array_almost_equal

import pytest

from distutils.version import LooseVersion
from scipy import __version__ as scipy_version
from scipy.spatial.distance import cdist
from sklearn.neighbors.dist_metrics import DistanceMetric
from sklearn.neighbors import BallTree
from sklearn.utils import check_random_state
from sklearn.utils.testing import assert_raises_regex


def dist_func(x1, x2, p):
    return np.sum((x1 - x2) ** p) ** (1. / p)


rng = check_random_state(0)
d = 4
n1 = 20
n2 = 25
X1 = rng.random_sample((n1, d)).astype('float64')
X2 = rng.random_sample((n2, d)).astype('float64')

# make boolean arrays: ones and zeros
X1_bool = X1.round(0)
X2_bool = X2.round(0)

V = rng.random_sample((d, d))
VI = np.dot(V, V.T)

BOOL_METRICS = ['matching', 'jaccard', 'dice',
                'kulsinski', 'rogerstanimoto', 'russellrao',
                'sokalmichener', 'sokalsneath']

METRICS_DEFAULT_PARAMS = {'euclidean': {},
                          'cityblock': {},
                          'minkowski': dict(p=(1, 1.5, 2, 3)),
                          'chebyshev': {},
                          'seuclidean': dict(V=(rng.random_sample(d),)),
                          'wminkowski': dict(p=(1, 1.5, 3),
                                             w=(rng.random_sample(d),)),
                          'mahalanobis': dict(VI=(VI,)),
                          'hamming': {},
                          'canberra': {},
                          'braycurtis': {}}


@pytest.mark.parametrize('metric', METRICS_DEFAULT_PARAMS)
def test_cdist(metric):
    argdict = METRICS_DEFAULT_PARAMS[metric]
    keys = argdict.keys()
    for vals in itertools.product(*argdict.values()):
        kwargs = dict(zip(keys, vals))
        D_true = cdist(X1, X2, metric, **kwargs)
        check_cdist(metric, kwargs, D_true)


@pytest.mark.parametrize('metric', BOOL_METRICS)
def test_cdist_bool_metric(metric):
    D_true = cdist(X1_bool, X2_bool, metric)
    check_cdist_bool(metric, D_true)


def check_cdist(metric, kwargs, D_true):
    dm = DistanceMetric.get_metric(metric, **kwargs)
    D12 = dm.pairwise(X1, X2)
    assert_array_almost_equal(D12, D_true)


def check_cdist_bool(metric, D_true):
    dm = DistanceMetric.get_metric(metric)
    D12 = dm.pairwise(X1_bool, X2_bool)
    assert_array_almost_equal(D12, D_true)


@pytest.mark.parametrize('metric', METRICS_DEFAULT_PARAMS)
def test_pdist(metric):
    argdict = METRICS_DEFAULT_PARAMS[metric]
    keys = argdict.keys()
    for vals in itertools.product(*argdict.values()):
        kwargs = dict(zip(keys, vals))
        D_true = cdist(X1, X1, metric, **kwargs)
        check_pdist(metric, kwargs, D_true)


@pytest.mark.parametrize('metric', BOOL_METRICS)
def test_pdist_bool_metrics(metric):
    D_true = cdist(X1_bool, X1_bool, metric)
    check_pdist_bool(metric, D_true)


def check_pdist(metric, kwargs, D_true):
    dm = DistanceMetric.get_metric(metric, **kwargs)
    D12 = dm.pairwise(X1)
    assert_array_almost_equal(D12, D_true)


def check_pdist_bool(metric, D_true):
    dm = DistanceMetric.get_metric(metric)
    D12 = dm.pairwise(X1_bool)
    # Based on https://github.com/scipy/scipy/pull/7373
    # When comparing two all-zero vectors, scipy>=1.2.0 jaccard metric
    # was changed to return 0, instead of nan.
    if metric == 'jaccard' and LooseVersion(scipy_version) < '1.2.0':
        D_true[np.isnan(D_true)] = 0
    assert_array_almost_equal(D12, D_true)


@pytest.mark.parametrize('metric', METRICS_DEFAULT_PARAMS)
def test_pickle(metric):
    argdict = METRICS_DEFAULT_PARAMS[metric]
    keys = argdict.keys()
    for vals in itertools.product(*argdict.values()):
        kwargs = dict(zip(keys, vals))
        check_pickle(metric, kwargs)


@pytest.mark.parametrize('metric', BOOL_METRICS)
def test_pickle_bool_metrics(metric):
    dm = DistanceMetric.get_metric(metric)
    D1 = dm.pairwise(X1_bool)
    dm2 = pickle.loads(pickle.dumps(dm))
    D2 = dm2.pairwise(X1_bool)
    assert_array_almost_equal(D1, D2)


def check_pickle(metric, kwargs):
    dm = DistanceMetric.get_metric(metric, **kwargs)
    D1 = dm.pairwise(X1)
    dm2 = pickle.loads(pickle.dumps(dm))
    D2 = dm2.pairwise(X1)
    assert_array_almost_equal(D1, D2)


def test_haversine_metric():
    def haversine_slow(x1, x2):
        return 2 * np.arcsin(np.sqrt(np.sin(0.5 * (x1[0] - x2[0])) ** 2
                                     + np.cos(x1[0]) * np.cos(x2[0]) *
                                     np.sin(0.5 * (x1[1] - x2[1])) ** 2))

    X = np.random.random((10, 2))

    haversine = DistanceMetric.get_metric("haversine")

    D1 = haversine.pairwise(X)
    D2 = np.zeros_like(D1)
    for i, x1 in enumerate(X):
        for j, x2 in enumerate(X):
            D2[i, j] = haversine_slow(x1, x2)

    assert_array_almost_equal(D1, D2)
    assert_array_almost_equal(haversine.dist_to_rdist(D1),
                              np.sin(0.5 * D2) ** 2)


def test_pyfunc_metric():
    X = np.random.random((10, 3))

    euclidean = DistanceMetric.get_metric("euclidean")
    pyfunc = DistanceMetric.get_metric("pyfunc", func=dist_func, p=2)

    # Check if both callable metric and predefined metric initialized
    # DistanceMetric object is picklable
    euclidean_pkl = pickle.loads(pickle.dumps(euclidean))
    pyfunc_pkl = pickle.loads(pickle.dumps(pyfunc))

    D1 = euclidean.pairwise(X)
    D2 = pyfunc.pairwise(X)

    D1_pkl = euclidean_pkl.pairwise(X)
    D2_pkl = pyfunc_pkl.pairwise(X)

    assert_array_almost_equal(D1, D2)
    assert_array_almost_equal(D1_pkl, D2_pkl)


def test_bad_pyfunc_metric():
    def wrong_distance(x, y):
        return "1"

    X = np.ones((5, 2))
    assert_raises_regex(TypeError,
                        "Custom distance function must accept two vectors",
                        BallTree, X, metric=wrong_distance)


def test_input_data_size():
    # Regression test for #6288
    # Previoulsly, a metric requiring a particular input dimension would fail
    def custom_metric(x, y):
        assert x.shape[0] == 3
        return np.sum((x - y) ** 2)

    rng = check_random_state(0)
    X = rng.rand(10, 3)

    pyfunc = DistanceMetric.get_metric("pyfunc", func=dist_func, p=2)
    eucl = DistanceMetric.get_metric("euclidean")
    assert_array_almost_equal(pyfunc.pairwise(X), eucl.pairwise(X))
