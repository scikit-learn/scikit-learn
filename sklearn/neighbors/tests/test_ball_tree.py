import pickle
import numpy as np
from numpy.testing import assert_array_almost_equal
from sklearn.neighbors.ball_tree import (BallTree, NeighborsHeap,
                                         simultaneous_sort, kernel_norm,
                                         nodeheap_sort, DTYPE, ITYPE)
from sklearn.neighbors.dist_metrics import DistanceMetric
from sklearn.utils import check_random_state
from sklearn.utils.testing import SkipTest, assert_allclose

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


def dist_func(x1, x2, p):
    return np.sum((x1 - x2) ** p) ** (1. / p)


def brute_force_neighbors(X, Y, k, metric, **kwargs):
    D = DistanceMetric.get_metric(metric, **kwargs).pairwise(Y, X)
    ind = np.argsort(D, axis=1)[:, :k]
    dist = D[np.arange(Y.shape[0])[:, None], ind]
    return dist, ind


def test_ball_tree_query():
    rng = check_random_state(0)
    X = rng.random_sample((40, DIMENSION))
    Y = rng.random_sample((10, DIMENSION))

    def check_neighbors(dualtree, breadth_first, k, metric, kwargs):
        bt = BallTree(X, leaf_size=1, metric=metric, **kwargs)
        dist1, ind1 = bt.query(Y, k, dualtree=dualtree,
                               breadth_first=breadth_first)
        dist2, ind2 = brute_force_neighbors(X, Y, k, metric, **kwargs)

        # don't check indices here: if there are any duplicate distances,
        # the indices may not match.  Distances should not have this problem.
        assert_array_almost_equal(dist1, dist2)

    for (metric, kwargs) in METRICS.items():
        for k in (1, 3, 5):
            for dualtree in (True, False):
                for breadth_first in (True, False):
                    yield (check_neighbors,
                           dualtree, breadth_first,
                           k, metric, kwargs)


def test_ball_tree_query_boolean_metrics():
    rng = check_random_state(0)
    X = rng.random_sample((40, 10)).round(0)
    Y = rng.random_sample((10, 10)).round(0)
    k = 5

    def check_neighbors(metric):
        bt = BallTree(X, leaf_size=1, metric=metric)
        dist1, ind1 = bt.query(Y, k)
        dist2, ind2 = brute_force_neighbors(X, Y, k, metric)
        assert_array_almost_equal(dist1, dist2)

    for metric in BOOLEAN_METRICS:
        yield check_neighbors, metric


def test_ball_tree_query_discrete_metrics():
    rng = check_random_state(0)
    X = (4 * rng.random_sample((40, 10))).round(0)
    Y = (4 * rng.random_sample((10, 10))).round(0)
    k = 5

    def check_neighbors(metric):
        bt = BallTree(X, leaf_size=1, metric=metric)
        dist1, ind1 = bt.query(Y, k)
        dist2, ind2 = brute_force_neighbors(X, Y, k, metric)
        assert_array_almost_equal(dist1, dist2)

    for metric in DISCRETE_METRICS:
        yield check_neighbors, metric


def test_ball_tree_query_radius(n_samples=100, n_features=10):
    rng = check_random_state(0)
    X = 2 * rng.random_sample(size=(n_samples, n_features)) - 1
    query_pt = np.zeros(n_features, dtype=float)

    eps = 1E-15  # roundoff error can cause test to fail
    bt = BallTree(X, leaf_size=5)
    rad = np.sqrt(((X - query_pt) ** 2).sum(1))

    for r in np.linspace(rad[0], rad[-1], 100):
        ind = bt.query_radius([query_pt], r + eps)[0]
        i = np.where(rad <= r + eps)[0]

        ind.sort()
        i.sort()

        assert_array_almost_equal(i, ind)


def test_ball_tree_query_radius_distance(n_samples=100, n_features=10):
    rng = check_random_state(0)
    X = 2 * rng.random_sample(size=(n_samples, n_features)) - 1
    query_pt = np.zeros(n_features, dtype=float)

    eps = 1E-15  # roundoff error can cause test to fail
    bt = BallTree(X, leaf_size=5)
    rad = np.sqrt(((X - query_pt) ** 2).sum(1))

    for r in np.linspace(rad[0], rad[-1], 100):
        ind, dist = bt.query_radius([query_pt], r + eps, return_distance=True)

        ind = ind[0]
        dist = dist[0]

        d = np.sqrt(((query_pt - X[ind]) ** 2).sum(1))

        assert_array_almost_equal(d, dist)


def compute_kernel_slow(Y, X, kernel, h):
    d = np.sqrt(((Y[:, None, :] - X) ** 2).sum(-1))
    norm = kernel_norm(h, X.shape[1], kernel)

    if kernel == 'gaussian':
        return norm * np.exp(-0.5 * (d * d) / (h * h)).sum(-1)
    elif kernel == 'tophat':
        return norm * (d < h).sum(-1)
    elif kernel == 'epanechnikov':
        return norm * ((1.0 - (d * d) / (h * h)) * (d < h)).sum(-1)
    elif kernel == 'exponential':
        return norm * (np.exp(-d / h)).sum(-1)
    elif kernel == 'linear':
        return norm * ((1 - d / h) * (d < h)).sum(-1)
    elif kernel == 'cosine':
        return norm * (np.cos(0.5 * np.pi * d / h) * (d < h)).sum(-1)
    else:
        raise ValueError('kernel not recognized')


def check_results(kernel, h, atol, rtol, breadth_first, bt, Y, dens_true):
    dens = bt.kernel_density(Y, h, atol=atol, rtol=rtol,
                             kernel=kernel,
                             breadth_first=breadth_first)
    assert_allclose(dens, dens_true,
                    atol=atol, rtol=max(rtol, 1e-7))


def test_ball_tree_kde(n_samples=100, n_features=3):
    rng = check_random_state(0)
    X = rng.random_sample((n_samples, n_features))
    Y = rng.random_sample((n_samples, n_features))
    bt = BallTree(X, leaf_size=10)

    for kernel in ['gaussian', 'tophat', 'epanechnikov',
                   'exponential', 'linear', 'cosine']:
        for h in [0.01, 0.1, 1]:
            dens_true = compute_kernel_slow(Y, X, kernel, h)

            for rtol in [0, 1E-5]:
                for atol in [1E-6, 1E-2]:
                    for breadth_first in (True, False):
                        yield (check_results, kernel, h, atol, rtol,
                               breadth_first, bt, Y, dens_true)


def test_gaussian_kde(n_samples=1000):
    # Compare gaussian KDE results to scipy.stats.gaussian_kde
    from scipy.stats import gaussian_kde
    rng = check_random_state(0)
    x_in = rng.normal(0, 1, n_samples)
    x_out = np.linspace(-5, 5, 30)

    for h in [0.01, 0.1, 1]:
        bt = BallTree(x_in[:, None])
        try:
            gkde = gaussian_kde(x_in, bw_method=h / np.std(x_in))
        except TypeError:
            raise SkipTest("Old version of scipy, doesn't accept "
                           "explicit bandwidth.")

        dens_bt = bt.kernel_density(x_out[:, None], h) / n_samples
        dens_gkde = gkde.evaluate(x_out)

        assert_array_almost_equal(dens_bt, dens_gkde, decimal=3)


def test_ball_tree_two_point(n_samples=100, n_features=3):
    rng = check_random_state(0)
    X = rng.random_sample((n_samples, n_features))
    Y = rng.random_sample((n_samples, n_features))
    r = np.linspace(0, 1, 10)
    bt = BallTree(X, leaf_size=10)

    D = DistanceMetric.get_metric("euclidean").pairwise(Y, X)
    counts_true = [(D <= ri).sum() for ri in r]

    def check_two_point(r, dualtree):
        counts = bt.two_point_correlation(Y, r=r, dualtree=dualtree)
        assert_array_almost_equal(counts, counts_true)

    for dualtree in (True, False):
        yield check_two_point, r, dualtree


def test_ball_tree_pickle():
    rng = check_random_state(0)
    X = rng.random_sample((10, 3))

    bt1 = BallTree(X, leaf_size=1)
    # Test if BallTree with callable metric is picklable
    bt1_pyfunc = BallTree(X, metric=dist_func, leaf_size=1, p=2)

    ind1, dist1 = bt1.query(X)
    ind1_pyfunc, dist1_pyfunc = bt1_pyfunc.query(X)

    def check_pickle_protocol(protocol):
        s = pickle.dumps(bt1, protocol=protocol)
        bt2 = pickle.loads(s)

        s_pyfunc = pickle.dumps(bt1_pyfunc, protocol=protocol)
        bt2_pyfunc = pickle.loads(s_pyfunc)

        ind2, dist2 = bt2.query(X)
        ind2_pyfunc, dist2_pyfunc = bt2_pyfunc.query(X)

        assert_array_almost_equal(ind1, ind2)
        assert_array_almost_equal(dist1, dist2)

        assert_array_almost_equal(ind1_pyfunc, ind2_pyfunc)
        assert_array_almost_equal(dist1_pyfunc, dist2_pyfunc)

    for protocol in (0, 1, 2):
        yield check_pickle_protocol, protocol


def test_neighbors_heap(n_pts=5, n_nbrs=10):
    heap = NeighborsHeap(n_pts, n_nbrs)

    for row in range(n_pts):
        d_in = rng.random_sample(2 * n_nbrs).astype(DTYPE)
        i_in = np.arange(2 * n_nbrs, dtype=ITYPE)
        for d, i in zip(d_in, i_in):
            heap.push(row, d, i)

        ind = np.argsort(d_in)
        d_in = d_in[ind]
        i_in = i_in[ind]

        d_heap, i_heap = heap.get_arrays(sort=True)

        assert_array_almost_equal(d_in[:n_nbrs], d_heap[row])
        assert_array_almost_equal(i_in[:n_nbrs], i_heap[row])


def test_node_heap(n_nodes=50):
    vals = rng.random_sample(n_nodes).astype(DTYPE)

    i1 = np.argsort(vals)
    vals2, i2 = nodeheap_sort(vals)

    assert_array_almost_equal(i1, i2)
    assert_array_almost_equal(vals[i1], vals2)


def test_simultaneous_sort(n_rows=10, n_pts=201):
    dist = rng.random_sample((n_rows, n_pts)).astype(DTYPE)
    ind = (np.arange(n_pts) + np.zeros((n_rows, 1))).astype(ITYPE)

    dist2 = dist.copy()
    ind2 = ind.copy()

    # simultaneous sort rows using function
    simultaneous_sort(dist, ind)

    # simultaneous sort rows using numpy
    i = np.argsort(dist2, axis=1)
    row_ind = np.arange(n_rows)[:, None]
    dist2 = dist2[row_ind, i]
    ind2 = ind2[row_ind, i]

    assert_array_almost_equal(dist, dist2)
    assert_array_almost_equal(ind, ind2)


def test_query_haversine():
    rng = check_random_state(0)
    X = 2 * np.pi * rng.random_sample((40, 2))
    bt = BallTree(X, leaf_size=1, metric='haversine')
    dist1, ind1 = bt.query(X, k=5)
    dist2, ind2 = brute_force_neighbors(X, X, k=5, metric='haversine')

    assert_array_almost_equal(dist1, dist2)
    assert_array_almost_equal(ind1, ind2)
