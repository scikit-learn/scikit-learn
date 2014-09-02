from collections import defaultdict

import numpy as np
from numpy.testing import assert_array_almost_equal
from sklearn.utils.graph import (graph_shortest_path,
                                 modified_graph_shortest_path,
                                 single_source_shortest_path_length)


def floyd_warshall_slow(graph, directed=False):
    N = graph.shape[0]

    #set nonzero entries to infinity
    graph[np.where(graph == 0)] = np.inf

    #set diagonal to zero
    graph.flat[::N + 1] = 0

    if not directed:
        graph = np.minimum(graph, graph.T)

    for k in range(N):
        for i in range(N):
            for j in range(N):
                graph[i, j] = min(graph[i, j], graph[i, k] + graph[k, j])

    graph[np.where(np.isinf(graph))] = 0

    return graph


def generate_graph(N=20):
    #sparse grid of distances
    rng = np.random.RandomState(0)
    dist_matrix = rng.random_sample((N, N))

    #make symmetric: distances are not direction-dependent
    dist_matrix += dist_matrix.T

    #make graph sparse
    i = (rng.randint(N, size=N * N // 2), rng.randint(N, size=N * N // 2))
    dist_matrix[i] = 0

    #set diagonal to zero
    dist_matrix.flat[::N + 1] = 0

    return dist_matrix


def test_floyd_warshall():
    dist_matrix = generate_graph(20)

    for directed in (True, False):
        graph_FW = graph_shortest_path(dist_matrix, directed, 'FW')[0]
        graph_py = floyd_warshall_slow(dist_matrix.copy(), directed)

        assert_array_almost_equal(graph_FW, graph_py)


def test_dijkstra():
    dist_matrix = generate_graph(20)

    for directed in (True, False):
        graph_D = graph_shortest_path(dist_matrix, directed, 'D')[0]
        graph_py = floyd_warshall_slow(dist_matrix.copy(), directed)

        assert_array_almost_equal(graph_D, graph_py)


def test_shortest_path():
    dist_matrix = generate_graph(20)
    # We compare path length and not costs (-> set distances to 0 or 1)
    dist_matrix[dist_matrix != 0] = 1

    for directed in (True, False):
        if not directed:
            dist_matrix = np.minimum(dist_matrix, dist_matrix.T)

        graph_py = floyd_warshall_slow(dist_matrix.copy(), directed)
        for i in range(dist_matrix.shape[0]):
            # Non-reachable nodes have distance 0 in graph_py
            dist_dict = defaultdict(int)
            dist_dict.update(single_source_shortest_path_length(dist_matrix,
                                                                i))

            for j in range(graph_py[i].shape[0]):
                assert_array_almost_equal(dist_dict[j], graph_py[i, j])


def test_dijkstra_bug_fix():
    X = np.array([[0., 0., 4.],
                  [1., 0., 2.],
                  [0., 5., 0.]])
    dist_FW, pred_FW = graph_shortest_path(X, directed=False, method='FW')
    dist_D, pred_D = graph_shortest_path(X, directed=False, method='D')
    assert_array_almost_equal(dist_D, dist_FW)
    assert_array_almost_equal(pred_D, pred_FW)


def test_predecessors():
    kng = np.array([[0., 3., 0., 0., 1.],
                    [3., 0., 0., 1., 2.],
                    [0., 0., 0., 2., 0.],
                    [0., 1., 2., 0., 5.],
                    [1., 2., 0., 5., 0.]])

    pred_true = np.array([[-1,  0,  3,  1,  0],
                          [ 1, -1,  3,  1,  1],
                          [ 1,  3, -1,  2,  1],
                          [ 1,  3,  3, -1,  1],
                          [ 4,  4,  3,  1, -1]])
    _, pred_FW = graph_shortest_path(kng, directed=False, method='FW')
    _, pred_D = graph_shortest_path(kng, directed=False, method='D')
    assert_array_almost_equal(pred_true, pred_FW)
    assert_array_almost_equal(pred_true, pred_D)


def test_modified_shortest_path():
    kng = np.array([[0, 3,   0,   0, 1, 0, 0, 0, 0,   2],
                    [3, 0,   0,   1, 2, 0, 0, 2, 0,   1.5],
                    [0, 0,   0,   2, 0, 0, 4, 3, 4.5, 0],
                    [0, 1,   2,   0, 5, 0, 1, 0, 0,   0],
                    [1, 2,   0,   5, 0, 1, 0, 0, 0,   0],
                    [0, 0,   0,   0, 1, 0, 2, 0, 0,   0],
                    [0, 0,   4,   1, 0, 2, 0, 0, 0,   0],
                    [0, 2,   3,   0, 0, 0, 0, 0, 4,   0],
                    [0, 0,   4.5, 0, 0, 0, 0, 4, 0,   0],
                    [2, 1.5, 0,   0, 0, 0, 0, 0, 0,   0]])
    directed = False
    dist_FW, pred_FW = graph_shortest_path(kng, directed=directed, method='FW')
    dist_mod = np.zeros_like(dist_FW)
    pred_mod = np.ones_like(pred_FW) * -1
    vertices = range(dist_mod.shape[0])
    for v in range(dist_mod.shape[0]):
        if directed:
            C_u = np.array(vertices, dtype=np.int)
            C_u = np.setdiff1d(C_u, [v], assume_unique=True)
        else:
            vertices.remove(v)
            C_u = np.array(vertices, dtype=np.int)
        dist_mod, pred_mod = modified_graph_shortest_path(v, dist_mod.shape[0],
                                    C_u, kng, dist_mod, pred_mod,
                                    directed=directed)

    assert_array_almost_equal(dist_mod, dist_FW)
    assert_array_almost_equal(pred_mod, pred_FW)


if __name__ == '__main__':
    import nose
    nose.runmodule()
