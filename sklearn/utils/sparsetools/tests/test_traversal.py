import numpy as np
from scipy.sparse import csr_matrix
from numpy.testing import assert_array_equal, assert_equal
from .. import connected_components


def test_connected_components():
    """Tests connected components with a basic graph."""
    # Test with a symmetric matrix
    # Create a dense graph as a numpy array
    graph = [[0, 1, 0, 0, 0],
             [1, 0, 0, 0, 0],
             [0, 0, 0, 8, 5],
             [0, 0, 8, 0, 1],
             [0, 0, 5, 1, 0]]
    graph = np.asarray(graph)
    _ensure_correct_components(graph)
    # Try as a sparse matrix
    sparse_graph = csr_matrix(graph)
    _ensure_correct_components(graph)
    # Try as a sparse matrix, upper triangle only
    graph = [[0, 1, 0, 0, 0],
             [0, 0, 0, 0, 0],
             [0, 0, 0, 8, 5],
             [0, 0, 0, 0, 1],
             [0, 0, 0, 0, 0]]
    sparse_graph = csr_matrix(graph)
    _ensure_correct_components(graph)


def _ensure_correct_components(graph):
    # Ensures the graph components are correctly found.
    # This works only on the `graph` found in `test_connected_components`.
    n_components, labels = connected_components(graph)
    assert_equal(n_components, 2)
    expected_labels = np.array([0, 0, 1, 1, 1])
    assert_array_equal(labels, expected_labels)


def test_components_singleton():
    """Tests connected components with singletons in the graph."""
    # Test with a symmetric matrix
    # Create a dense graph as a numpy array
    graph = [[0, 1, 0, 0, 0],
             [1, 0, 0, 0, 0],
             [0, 0, 0, 0, 0],
             [0, 0, 0, 0, 0],
             [0, 0, 0, 0, 0]]
    graph = np.asarray(graph)
    _ensure_correct_components_singleton(graph)
    # Try as a sparse matrix
    sparse_graph = csr_matrix(graph)
    _ensure_correct_components_singleton(graph)
    # Try as a sparse matrix, upper triangle only
    graph = [[0, 1, 0, 0, 0],
             [0, 0, 0, 0, 0],
             [0, 0, 0, 0, 0],
             [0, 0, 0, 0, 0],
             [0, 0, 0, 0, 0]]
    sparse_graph = csr_matrix(graph)
    _ensure_correct_components_singleton(graph)


def _ensure_correct_components_singleton(graph):
    # Ensures the graph components are correctly found.
    # This works only on the `graph` found in `test_connected_components`.
    n_components, labels = connected_components(graph)
    assert_equal(n_components, 4)
    expected_labels = np.array([0, 0, 1, 2, 3])
    assert_array_equal(labels, expected_labels)
