import numpy as np
from numpy.testing import assert_array_almost_equal
from sklearn.utils.graph import graph_shortest_path


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
    dist_matrix = np.random.random((N, N))

    #make symmetric: distances are not direction-dependent
    dist_matrix += dist_matrix.T

    #make graph sparse
    i = (np.random.randint(N, size=N * N / 2),
         np.random.randint(N, size=N * N / 2))
    dist_matrix[i] = 0

    #set diagonal to zero
    dist_matrix.flat[::N + 1] = 0

    return dist_matrix


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


def test_dijkstra_bug_fix():
    X = np.array([[0., 0., 4.],
                  [1., 0., 2.],
                  [0., 5., 0.]])
    dist_FW = graph_shortest_path(X, directed=False, method='FW')
    dist_D = graph_shortest_path(X, directed=False, method='D')
    assert_array_almost_equal(dist_D, dist_FW)


if __name__ == '__main__':
    import nose
    nose.runmodule()
