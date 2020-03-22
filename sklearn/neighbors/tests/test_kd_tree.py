import numpy as np

from sklearn.neighbors import KDTree

DIMENSION = 3

METRICS = {
    "euclidean": {},
    "manhattan": {},
    "chebyshev": {},
    "minkowski": dict(p=3),
}


def test_kd_tree_filter_query():

    X = np.array(
        [
            [0, 2, 1],
            [2, 2, 3],
            [10, 3, 1],
            [4, 5, 6],
            [0, 0, 0],
            [-10, -3, 2],
            [1, -5, 2],
        ]
    )
    tree = KDTree(X, leaf_size=1)
    queries = np.array([[1, 1, 1]])

    dist, ind = tree.query_filtered(queries, [[1, -1, -1]], k=2)

    assert ind[0][0] == 0
    assert ind[0][1] == 4
    assert dist[0][0] == 1

    rand = np.random.RandomState(18)

    n_features = 6
    data_size = 2000
    query_num = 20

    X = rand.rand(data_size, n_features)

    tree2 = KDTree(X, leaf_size=4)

    queries = rand.rand(query_num, n_features)
    filter_radiuses = rand.rand(query_num, n_features) - 0.5

    bruteforce_ind = []
    bruteforce_dist = []
    for q, rads in zip(queries, filter_radiuses):
        effective_rads = np.array([r if r >= 0 else np.inf for r in rads])
        _dist_inds = rads < 0
        min_rdist = np.inf
        min_ind = 0
        for i, x_i in enumerate(X):
            if (np.abs(x_i - q) <= effective_rads).all():
                _d = ((x_i[_dist_inds] - q[_dist_inds]) ** 2).sum()
                if _d < min_rdist:
                    min_rdist = _d
                    min_ind = i
        bruteforce_ind.append([min_ind])
        bruteforce_dist.append([np.sqrt(min_rdist)])

    bf_dist = np.array(bruteforce_dist)
    bf_ind = np.array(bruteforce_ind)

    dist, ind = tree2.query_filtered(queries, filter_radiuses, k=1)

    for q, rad, i, bi, d, bd in zip(
        queries,
        filter_radiuses,
        ind[:, 0],
        bf_ind[:, 0],
        dist[:, 0],
        bf_dist[:, 0],
    ):
        if (rad < 0).any():
            assert i == bi
        assert d == bd
