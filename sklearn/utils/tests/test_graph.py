import numpy as np
from scipy.sparse.csgraph import connected_components

from sklearn.neighbors import kneighbors_graph
from sklearn.utils.graph import _fix_connected_components


def test_fix_connected_components():
    X = np.array([0, 1, 2, 5, 6, 7])[:, None]
    graph = kneighbors_graph(X, n_neighbors=2, mode="distance")

    n_connected_components, labels = connected_components(graph)
    assert n_connected_components > 1

    graph = _fix_connected_components(
        X, graph, n_connected_components, labels, metric="euclidean"
    )

    n_connected_components, labels = connected_components(graph)
    assert n_connected_components == 1
