#!/usr/bin/env python
from nose.tools import *
from nose import SkipTest
from nose.plugins.attrib import attr
import networkx as nx
from networkx.algorithms import bipartite


class TestBipartiteBasic:

    def test_is_bipartite(self):
        assert_true(bipartite.is_bipartite(nx.path_graph(4)))
        assert_true(bipartite.is_bipartite(nx.DiGraph([(1, 0)])))
        assert_false(bipartite.is_bipartite(nx.complete_graph(3)))

    def test_bipartite_color(self):
        G = nx.path_graph(4)
        c = bipartite.color(G)
        assert_equal(c, {0: 1, 1: 0, 2: 1, 3: 0})

    @raises(nx.NetworkXError)
    def test_not_bipartite_color(self):
        c = bipartite.color(nx.complete_graph(4))

    def test_bipartite_directed(self):
        G = bipartite.random_graph(10, 10, 0.1, directed=True)
        assert_true(bipartite.is_bipartite(G))

    def test_bipartite_sets(self):
        G = nx.path_graph(4)
        X, Y = bipartite.sets(G)
        assert_equal(X, {0, 2})
        assert_equal(Y, {1, 3})

    def test_bipartite_sets_directed(self):
        G = nx.path_graph(4)
        D = G.to_directed()
        X, Y = bipartite.sets(D)
        assert_equal(X, {0, 2})
        assert_equal(Y, {1, 3})

    def test_bipartite_sets_given_top_nodes(self):
        G = nx.path_graph(4)
        top_nodes = [0, 2]
        X, Y = bipartite.sets(G, top_nodes)
        assert_equal(X, {0, 2})
        assert_equal(Y, {1, 3})

    @raises(nx.AmbiguousSolution)
    def test_bipartite_sets_disconnected(self):
        G = nx.path_graph(4)
        G.add_edges_from([(5, 6), (6, 7)])
        X, Y = bipartite.sets(G)

    def test_is_bipartite_node_set(self):
        G = nx.path_graph(4)
        assert_true(bipartite.is_bipartite_node_set(G, [0, 2]))
        assert_true(bipartite.is_bipartite_node_set(G, [1, 3]))
        assert_false(bipartite.is_bipartite_node_set(G, [1, 2]))
        G.add_edge(10, 20)
        assert_true(bipartite.is_bipartite_node_set(G, [0, 2, 10]))
        assert_true(bipartite.is_bipartite_node_set(G, [0, 2, 20]))
        assert_true(bipartite.is_bipartite_node_set(G, [1, 3, 10]))
        assert_true(bipartite.is_bipartite_node_set(G, [1, 3, 20]))

    def test_bipartite_density(self):
        G = nx.path_graph(5)
        X, Y = bipartite.sets(G)
        density = float(len(list(G.edges()))) / (len(X) * len(Y))
        assert_equal(bipartite.density(G, X), density)
        D = nx.DiGraph(G.edges())
        assert_equal(bipartite.density(D, X), density / 2.0)
        assert_equal(bipartite.density(nx.Graph(), {}), 0.0)

    def test_bipartite_degrees(self):
        G = nx.path_graph(5)
        X = set([1, 3])
        Y = set([0, 2, 4])
        u, d = bipartite.degrees(G, Y)
        assert_equal(dict(u), {1: 2, 3: 2})
        assert_equal(dict(d), {0: 1, 2: 2, 4: 1})

    def test_bipartite_weighted_degrees(self):
        G = nx.path_graph(5)
        G.add_edge(0, 1, weight=0.1, other=0.2)
        X = set([1, 3])
        Y = set([0, 2, 4])
        u, d = bipartite.degrees(G, Y, weight='weight')
        assert_equal(dict(u), {1: 1.1, 3: 2})
        assert_equal(dict(d), {0: 0.1, 2: 2, 4: 1})
        u, d = bipartite.degrees(G, Y, weight='other')
        assert_equal(dict(u), {1: 1.2, 3: 2})
        assert_equal(dict(d), {0: 0.2, 2: 2, 4: 1})

    @attr('numpy')
    def test_biadjacency_matrix_weight(self):
        try:
            import scipy
        except ImportError:
            raise SkipTest('SciPy not available.')
        G = nx.path_graph(5)
        G.add_edge(0, 1, weight=2, other=4)
        X = [1, 3]
        Y = [0, 2, 4]
        M = bipartite.biadjacency_matrix(G, X, weight='weight')
        assert_equal(M[0, 0], 2)
        M = bipartite.biadjacency_matrix(G, X, weight='other')
        assert_equal(M[0, 0], 4)

    @attr('numpy')
    def test_biadjacency_matrix(self):
        try:
            import scipy
        except ImportError:
            raise SkipTest('SciPy not available.')
        tops = [2, 5, 10]
        bots = [5, 10, 15]
        for i in range(len(tops)):
            G = bipartite.random_graph(tops[i], bots[i], 0.2)
            top = [n for n, d in G.nodes(data=True) if d['bipartite'] == 0]
            M = bipartite.biadjacency_matrix(G, top)
            assert_equal(M.shape[0], tops[i])
            assert_equal(M.shape[1], bots[i])

    @attr('numpy')
    def test_biadjacency_matrix_order(self):
        try:
            import scipy
        except ImportError:
            raise SkipTest('SciPy not available.')
        G = nx.path_graph(5)
        G.add_edge(0, 1, weight=2)
        X = [3, 1]
        Y = [4, 2, 0]
        M = bipartite.biadjacency_matrix(G, X, Y, weight='weight')
        assert_equal(M[1, 2], 2)
