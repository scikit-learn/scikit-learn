#!/usr/bin/env python
from nose.tools import *
from nose import SkipTest
import networkx as nx
from networkx.algorithms import bipartite
from networkx.testing.utils import assert_edges_equal


class TestBiadjacencyMatrix:
    @classmethod
    def setupClass(cls):
        global np, sp, sparse, np_assert_equal
        try:
            import numpy as np
            import scipy as sp
            import scipy.sparse as sparse
            np_assert_equal = np.testing.assert_equal
        except ImportError:
            raise SkipTest('SciPy sparse library not available.')

    def test_biadjacency_matrix_weight(self):
        G = nx.path_graph(5)
        G.add_edge(0, 1, weight=2, other=4)
        X = [1, 3]
        Y = [0, 2, 4]
        M = bipartite.biadjacency_matrix(G, X, weight='weight')
        assert_equal(M[0, 0], 2)
        M = bipartite.biadjacency_matrix(G, X, weight='other')
        assert_equal(M[0, 0], 4)

    def test_biadjacency_matrix(self):
        tops = [2, 5, 10]
        bots = [5, 10, 15]
        for i in range(len(tops)):
            G = bipartite.random_graph(tops[i], bots[i], 0.2)
            top = [n for n, d in G.nodes(data=True) if d['bipartite'] == 0]
            M = bipartite.biadjacency_matrix(G, top)
            assert_equal(M.shape[0], tops[i])
            assert_equal(M.shape[1], bots[i])

    def test_biadjacency_matrix_order(self):
        G = nx.path_graph(5)
        G.add_edge(0, 1, weight=2)
        X = [3, 1]
        Y = [4, 2, 0]
        M = bipartite.biadjacency_matrix(G, X, Y, weight='weight')
        assert_equal(M[1, 2], 2)

    @raises(nx.NetworkXError)
    def test_null_graph(self):
        bipartite.biadjacency_matrix(nx.Graph(), [])

    @raises(nx.NetworkXError)
    def test_empty_graph(self):
        bipartite.biadjacency_matrix(nx.Graph([(1, 0)]), [])

    @raises(nx.NetworkXError)
    def test_duplicate_row(self):
        bipartite.biadjacency_matrix(nx.Graph([(1, 0)]), [1, 1])

    @raises(nx.NetworkXError)
    def test_duplicate_col(self):
        bipartite.biadjacency_matrix(nx.Graph([(1, 0)]), [0], [1, 1])

    @raises(nx.NetworkXError)
    def test_duplicate_col(self):
        bipartite.biadjacency_matrix(nx.Graph([(1, 0)]), [0], [1, 1])

    @raises(nx.NetworkXError)
    def test_format_keyword(self):
        bipartite.biadjacency_matrix(nx.Graph([(1, 0)]), [0], format='foo')

    def test_from_biadjacency_roundtrip(self):
        B1 = nx.path_graph(5)
        M = bipartite.biadjacency_matrix(B1, [0, 2, 4])
        B2 = bipartite.from_biadjacency_matrix(M)
        assert_true(nx.is_isomorphic(B1, B2))

    def test_from_biadjacency_weight(self):
        M = sparse.csc_matrix([[1, 2], [0, 3]])
        B = bipartite.from_biadjacency_matrix(M)
        assert_edges_equal(B.edges(), [(0, 2), (0, 3), (1, 3)])
        B = bipartite.from_biadjacency_matrix(M, edge_attribute='weight')
        e = [(0, 2, {'weight': 1}), (0, 3, {'weight': 2}), (1, 3, {'weight': 3})]
        assert_edges_equal(B.edges(data=True), e)

    def test_from_biadjacency_multigraph(self):
        M = sparse.csc_matrix([[1, 2], [0, 3]])
        B = bipartite.from_biadjacency_matrix(M, create_using=nx.MultiGraph())
        assert_edges_equal(B.edges(), [(0, 2), (0, 3), (0, 3), (1, 3), (1, 3), (1, 3)])
