import itertools

import numpy as np
import pytest
from numpy.testing import assert_array_almost_equal, assert_allclose, assert_equal
from sklearn.neighbors._ball_tree import BallTree
from sklearn.utils import check_random_state
from sklearn.utils.validation import check_array
from sklearn.utils._testing import _convert_container

rng = np.random.RandomState(10)
V_mahalanobis = rng.rand(3, 3)
V_mahalanobis = np.dot(V_mahalanobis, V_mahalanobis.T)

DIMENSION = 3

METRICS = {
    "euclidean": {},
    "manhattan": {},
    "minkowski": dict(p=3),
    "chebyshev": {},
}

ADDITIONAL_METRICS = {
    "seuclidean": dict(V=rng.random_sample(DIMENSION)),
    "wminkowski": dict(p=3, w=rng.random_sample(DIMENSION)),
    "mahalanobis": dict(V=V_mahalanobis),
}

DISCRETE_METRICS = ["hamming", "canberra", "braycurtis"]

BOOLEAN_METRICS = [
    "jaccard",
    "dice",
    "rogerstanimoto",
    "russellrao",
    "sokalmichener",
    "sokalsneath",
]


def brute_force_neighbors(X, Y, k, metric, **kwargs):
    from sklearn.metrics import DistanceMetric

    X, Y = check_array(X), check_array(Y)
    D = DistanceMetric.get_metric(metric, **kwargs).pairwise(Y, X)
    ind = np.argsort(D, axis=1)[:, :k]
    dist = D[np.arange(Y.shape[0])[:, None], ind]
    return dist, ind


@pytest.mark.parametrize("metric", itertools.chain(BOOLEAN_METRICS, DISCRETE_METRICS))
@pytest.mark.parametrize("array_type", ["list", "array"])
def test_ball_tree_query_metrics(metric, array_type):
    rng = check_random_state(0)
    if metric in BOOLEAN_METRICS:
        X = rng.random_sample((40, 10)).round(0)
        Y = rng.random_sample((10, 10)).round(0)
    elif metric in DISCRETE_METRICS:
        X = (4 * rng.random_sample((40, 10))).round(0)
        Y = (4 * rng.random_sample((10, 10))).round(0)
    X = _convert_container(X, array_type)
    Y = _convert_container(Y, array_type)

    k = 5

    bt = BallTree(X, leaf_size=1, metric=metric)
    dist1, ind1 = bt.query(Y, k)
    dist2, ind2 = brute_force_neighbors(X, Y, k, metric)
    assert_array_almost_equal(dist1, dist2)


def test_query_haversine():
    rng = check_random_state(0)
    X = 2 * np.pi * rng.random_sample((40, 2))
    bt = BallTree(X, leaf_size=1, metric="haversine")
    dist1, ind1 = bt.query(X, k=5)
    dist2, ind2 = brute_force_neighbors(X, X, k=5, metric="haversine")

    assert_array_almost_equal(dist1, dist2)
    assert_array_almost_equal(ind1, ind2)


def test_array_object_type():
    """Check that we do not accept object dtype array."""
    X = np.array([(1, 2, 3), (2, 5), (5, 5, 1, 2)], dtype=object)
    with pytest.raises(ValueError, match="Unexpected dtype object provided"):
        BallTree(X)


def test_bad_pyfunc_metric():
    def wrong_returned_value(x, y):
        return "1"

    def one_arg_func(x):
        return 1.0  # pragma: no cover

    X = np.ones((5, 2))
    msg = "Custom distance function must accept two vectors and return a float."
    with pytest.raises(TypeError, match=msg):
        BallTree(X, metric=wrong_returned_value)

    msg = "takes 1 positional argument but 2 were given"
    with pytest.raises(TypeError, match=msg):
        BallTree(X, metric=one_arg_func)


@pytest.mark.parametrize("metric", itertools.chain(METRICS, BOOLEAN_METRICS))
def test_ball_tree_numerical_consistency(global_random_seed, metric):
    X_64, X_32, Y_64, Y_32 = get_dataset_for_query_methods(
        random_seed=global_random_seed
    )

    metric_params = METRICS.get(metric, {})
    bt_64 = BallTree(X_64, leaf_size=1, metric=metric, **metric_params)
    bt_32 = BallTree(X_32, leaf_size=1, metric=metric, **metric_params)

    # Test consistency with respect to the `query` method
    k = 5
    dist_64, ind_64 = bt_64.query(Y_64, k=k)
    dist_32, ind_32 = bt_32.query(Y_32, k=k)
    assert_allclose(dist_64, dist_32, rtol=1e-5)
    assert_equal(ind_64, ind_32)
    assert dist_64.dtype == np.float64
    assert dist_32.dtype == np.float32

    # Test consistency with respect to the `query_radius` method
    r = 2.38
    ind_64, neighbors_64 = bt_64.query_radius(Y_64[0:2, :], r=r)
    ind_32, neighbors_32 = bt_32.query_radius(Y_32[0:2, :], r=r)
    assert_equal(ind_64, ind_32)
    assert_allclose(
        neighbors_64,
        neighbors_32,
    )

    # Test consistency with respect to the `query_radius` method
    # with return distances being true
    ind_64, dist_64 = bt_64.query_radius(Y_64[4:5, :], r=r, return_distance=True)
    ind_32, dist_32 = bt_32.query_radius(Y_32[4:5, :], r=r, return_distance=True)
    assert_equal(ind_64[0], ind_32[0])
    assert_allclose(dist_64[0], dist_32[0], rtol=1e-5)
    assert dist_64[0].dtype == np.float64
    assert dist_32[0].dtype == np.float32


@pytest.mark.parametrize("metric", itertools.chain(METRICS, BOOLEAN_METRICS))
def test_kernel_density_numerical_consistency(global_random_seed, metric):
    # Test consistency with respect to the `kernel_density` method
    X_64, X_32, Y_64, Y_32 = get_dataset_for_kernel_density(
        random_seed=global_random_seed
    )

    metric_params = METRICS.get(metric, {})
    bt_64 = BallTree(X_64, leaf_size=1, metric=metric, **metric_params)
    bt_32 = BallTree(X_32, leaf_size=1, metric=metric, **metric_params)

    kernel = "gaussian"
    h = 0.1
    density64 = bt_64.kernel_density(Y_64, h=h, kernel=kernel, breadth_first=True)
    density32 = bt_32.kernel_density(Y_32, h=h, kernel=kernel, breadth_first=True)
    assert_allclose(density64, density32, rtol=1e-5)
    assert density64.dtype == np.float64
    assert density32.dtype == np.float32


def test_two_point_correlation_numerical_consistency(global_random_seed):
    # Test consistency with respect to the `two_point_correlation` method
    rng = np.random.RandomState(global_random_seed)
    _X = rng.random_sample((100, 3))
    _Y = rng.random_sample((5, 3))

    X_64 = _X.astype(dtype=np.float64, copy=False)
    Y_64 = _Y.astype(dtype=np.float64, copy=False)

    X_32 = _X.astype(dtype=np.float32, copy=False)
    Y_32 = _Y.astype(dtype=np.float32, copy=False)

    bt_64 = BallTree(X_64, leaf_size=10)
    bt_32 = BallTree(X_32, leaf_size=10)

    r = np.linspace(0, 1, 10)

    counts_64 = bt_64.two_point_correlation(Y_64, r=r, dualtree=True)
    counts_32 = bt_32.two_point_correlation(Y_32, r=r, dualtree=True)
    assert_allclose(counts_64, counts_32)


def get_dataset_for_query_methods(random_seed):
    rng = np.random.RandomState(random_seed)
    _X = rng.rand(100, 50)
    _Y = rng.rand(5, 50)

    X_64 = _X.astype(dtype=np.float64, copy=False)
    Y_64 = _Y.astype(dtype=np.float64, copy=False)

    X_32 = _X.astype(dtype=np.float32, copy=False)
    Y_32 = _Y.astype(dtype=np.float32, copy=False)

    return X_64, X_32, Y_64, Y_32


def get_dataset_for_kernel_density(random_seed):
    rng = np.random.RandomState(random_seed)
    _X = rng.random_sample((100, 3))
    _Y = rng.random_sample((5, 3))

    X_64 = _X.astype(dtype=np.float64, copy=False)
    Y_64 = _Y.astype(dtype=np.float64, copy=False)

    X_32 = _X.astype(dtype=np.float32, copy=False)
    Y_32 = _Y.astype(dtype=np.float32, copy=False)

    return X_64, X_32, Y_64, Y_32
