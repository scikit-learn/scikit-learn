import numpy as np
from numpy.testing import assert_array_equal, assert_equal
import pytest

from scipy.sparse import csr_matrix, coo_matrix, diags
from scipy.sparse.csgraph import maximum_bipartite_matching


def test_raises_on_dense_input():
    with pytest.raises(TypeError):
        graph = np.array([[0, 1], [0, 0]])
        maximum_bipartite_matching(graph)


def test_empty_graph():
    graph = csr_matrix((0, 0))
    x = maximum_bipartite_matching(graph, perm_type='row')
    y = maximum_bipartite_matching(graph, perm_type='column')
    expected_matching = np.array([])
    assert_array_equal(expected_matching, x)
    assert_array_equal(expected_matching, y)


def test_empty_left_partition():
    graph = csr_matrix((2, 0))
    x = maximum_bipartite_matching(graph, perm_type='row')
    y = maximum_bipartite_matching(graph, perm_type='column')
    assert_array_equal(np.array([]), x)
    assert_array_equal(np.array([-1, -1]), y)


def test_empty_right_partition():
    graph = csr_matrix((0, 3))
    x = maximum_bipartite_matching(graph, perm_type='row')
    y = maximum_bipartite_matching(graph, perm_type='column')
    assert_array_equal(np.array([-1, -1, -1]), x)
    assert_array_equal(np.array([]), y)


def test_graph_with_no_edges():
    graph = csr_matrix((2, 2))
    x = maximum_bipartite_matching(graph, perm_type='row')
    y = maximum_bipartite_matching(graph, perm_type='column')
    assert_array_equal(np.array([-1, -1]), x)
    assert_array_equal(np.array([-1, -1]), y)


def test_graph_that_causes_augmentation():
    # In this graph, column 1 is initially assigned to row 1, but it should be
    # reassigned to make room for row 2.
    graph = csr_matrix([[1, 1], [1, 0]])
    x = maximum_bipartite_matching(graph, perm_type='column')
    y = maximum_bipartite_matching(graph, perm_type='row')
    expected_matching = np.array([1, 0])
    assert_array_equal(expected_matching, x)
    assert_array_equal(expected_matching, y)


def test_graph_with_more_rows_than_columns():
    graph = csr_matrix([[1, 1], [1, 0], [0, 1]])
    x = maximum_bipartite_matching(graph, perm_type='column')
    y = maximum_bipartite_matching(graph, perm_type='row')
    assert_array_equal(np.array([0, -1, 1]), x)
    assert_array_equal(np.array([0, 2]), y)


def test_graph_with_more_columns_than_rows():
    graph = csr_matrix([[1, 1, 0], [0, 0, 1]])
    x = maximum_bipartite_matching(graph, perm_type='column')
    y = maximum_bipartite_matching(graph, perm_type='row')
    assert_array_equal(np.array([0, 2]), x)
    assert_array_equal(np.array([0, -1, 1]), y)


def test_explicit_zeros_count_as_edges():
    data = [0, 0]
    indices = [1, 0]
    indptr = [0, 1, 2]
    graph = csr_matrix((data, indices, indptr), shape=(2, 2))
    x = maximum_bipartite_matching(graph, perm_type='row')
    y = maximum_bipartite_matching(graph, perm_type='column')
    expected_matching = np.array([1, 0])
    assert_array_equal(expected_matching, x)
    assert_array_equal(expected_matching, y)


def test_large_random_graph_with_one_edge_incident_to_each_vertex():
    np.random.seed(42)
    A = diags(np.ones(25), offsets=0, format='csr')
    rand_perm = np.random.permutation(25)
    rand_perm2 = np.random.permutation(25)

    Rrow = np.arange(25)
    Rcol = rand_perm
    Rdata = np.ones(25, dtype=int)
    Rmat = coo_matrix((Rdata, (Rrow, Rcol))).tocsr()

    Crow = rand_perm2
    Ccol = np.arange(25)
    Cdata = np.ones(25, dtype=int)
    Cmat = coo_matrix((Cdata, (Crow, Ccol))).tocsr()
    # Randomly permute identity matrix
    B = Rmat * A * Cmat

    # Row permute
    perm = maximum_bipartite_matching(B, perm_type='row')
    Rrow = np.arange(25)
    Rcol = perm
    Rdata = np.ones(25, dtype=int)
    Rmat = coo_matrix((Rdata, (Rrow, Rcol))).tocsr()
    C1 = Rmat * B

    # Column permute
    perm2 = maximum_bipartite_matching(B, perm_type='column')
    Crow = perm2
    Ccol = np.arange(25)
    Cdata = np.ones(25, dtype=int)
    Cmat = coo_matrix((Cdata, (Crow, Ccol))).tocsr()
    C2 = B * Cmat

    # Should get identity matrix back
    assert_equal(any(C1.diagonal() == 0), False)
    assert_equal(any(C2.diagonal() == 0), False)
