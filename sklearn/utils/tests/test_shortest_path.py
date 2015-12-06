from collections import defaultdict

import numpy as np
from numpy.testing import (assert_equal,
                           assert_array_equal,
                           assert_array_almost_equal,
                           assert_array_less)

from sklearn.utils.graph import (graph_shortest_path,
                                 single_source_shortest_path_length)


def floyd_warshall_slow(graph, directed=False):
    N = graph.shape[0]

    # set nonzero entries to infinity
    graph[np.where(graph == 0)] = np.inf

    # set diagonal to zero
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
    # sparse grid of distances
    rng = np.random.RandomState(0)
    dist_matrix = rng.random_sample((N, N))

    # make symmetric: distances are not direction-dependent
    dist_matrix = dist_matrix + dist_matrix.T

    # make graph sparse
    i = (rng.randint(N, size=N * N // 2), rng.randint(N, size=N * N // 2))
    dist_matrix[i] = 0

    # set diagonal to zero
    dist_matrix.flat[::N + 1] = 0

    return dist_matrix


def get_random_vertices_subset(n_vertices, n_subset_size):
    vertices = np.arange(n_vertices)
    np.random.shuffle(vertices)
    vertices = vertices[:n_subset_size]
    vertices.sort()

    return vertices


def test_floyd_warshall():
    dist_matrix = generate_graph(20)

    for directed in (True, False):
        graph_FW = graph_shortest_path(dist_matrix, directed, 'FW')
        graph_py = floyd_warshall_slow(dist_matrix.copy(), directed)

        assert_array_almost_equal(graph_FW, graph_py)


def test_dijkstra():
    dist_matrix = generate_graph(20)

    for directed in (True, False):
        graph_D = graph_shortest_path(dist_matrix, directed, 'D')
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
    dist_FW = graph_shortest_path(X, directed=False, method='FW')
    dist_D = graph_shortest_path(X, directed=False, method='D')
    assert_array_almost_equal(dist_D, dist_FW)


def test_dijkstra_compute_vertices_subset():
    n = 10
    d = generate_graph(n)
    vertices = get_random_vertices_subset(n, 5)

    g = graph_shortest_path(d, directed=False, method='D',
                            only_vertices=vertices)

    assert_array_equal(g.shape, (len(vertices), n))

    for pos_in_graph, vertex in enumerate(vertices):
        # Assert dissimilarity with itself is zero.
        assert_equal(g[pos_in_graph, vertex], 0)

        # Assert all other dissimilarities are not zero.
        mask = np.ones(n, dtype=bool)
        mask[vertex] = 0
        assert_array_less(-1 * np.ones(n - 1), g[pos_in_graph][mask])


def test_dijkstra_equal_floyd_warshall():
    n = 10
    dist_matrix = generate_graph(n)
    vertices = get_random_vertices_subset(n, 5)

    g1 = graph_shortest_path(dist_matrix, directed=False, method='D')
    g2 = graph_shortest_path(dist_matrix, directed=False, method='FW')

    assert_array_equal(g1.shape, (n, n))
    assert_array_equal(g2.shape, (n, n))
    assert_array_almost_equal(g1, g2)

    g1 = graph_shortest_path(dist_matrix, directed=True, method='D')
    g2 = graph_shortest_path(dist_matrix, directed=True, method='FW')

    assert_array_equal(g1.shape, (n, n))
    assert_array_equal(g2.shape, (n, n))
    assert_array_almost_equal(g1, g2)

    g1 = graph_shortest_path(dist_matrix, directed=False, method='D',
                             only_vertices=vertices)
    g2 = graph_shortest_path(dist_matrix, directed=False, method='FW',
                             only_vertices=vertices)

    assert_array_equal(g1.shape, (len(vertices), n))
    assert_array_equal(g2.shape, (len(vertices), n))

    assert_array_almost_equal(g1, g2)

    g1 = graph_shortest_path(dist_matrix, directed=True, method='D',
                             only_vertices=vertices)
    g2 = graph_shortest_path(dist_matrix, directed=True, method='FW',
                             only_vertices=vertices)

    assert_array_equal(g1.shape, (len(vertices), n))
    assert_array_equal(g2.shape, (len(vertices), n))

    assert_array_almost_equal(g1, g2)
