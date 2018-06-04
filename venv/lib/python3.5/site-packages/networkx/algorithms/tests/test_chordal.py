#!/usr/bin/env python
from nose.tools import *
import networkx as nx


class TestMCS:

    def setUp(self):
        # simple graph
        connected_chordal_G = nx.Graph()
        connected_chordal_G.add_edges_from([(1, 2), (1, 3), (2, 3), (2, 4), (3, 4),
                                            (3, 5), (3, 6), (4, 5), (4, 6), (5, 6)])
        self.connected_chordal_G = connected_chordal_G

        chordal_G = nx.Graph()
        chordal_G.add_edges_from([(1, 2), (1, 3), (2, 3), (2, 4), (3, 4),
                                  (3, 5), (3, 6), (4, 5), (4, 6), (5, 6), (7, 8)])
        chordal_G.add_node(9)
        self.chordal_G = chordal_G

        non_chordal_G = nx.Graph()
        non_chordal_G.add_edges_from([(1, 2), (1, 3), (2, 4), (2, 5), (3, 4), (3, 5)])
        self.non_chordal_G = non_chordal_G

    def test_is_chordal(self):
        assert_false(nx.is_chordal(self.non_chordal_G))
        assert_true(nx.is_chordal(self.chordal_G))
        assert_true(nx.is_chordal(self.connected_chordal_G))
        assert_true(nx.is_chordal(nx.complete_graph(3)))
        assert_true(nx.is_chordal(nx.cycle_graph(3)))
        assert_false(nx.is_chordal(nx.cycle_graph(5)))

    def test_induced_nodes(self):
        G = nx.generators.classic.path_graph(10)
        I = nx.find_induced_nodes(G, 1, 9, 2)
        assert_equal(I, set([1, 2, 3, 4, 5, 6, 7, 8, 9]))
        assert_raises(nx.NetworkXTreewidthBoundExceeded,
                      nx.find_induced_nodes, G, 1, 9, 1)
        I = nx.find_induced_nodes(self.chordal_G, 1, 6)
        assert_equal(I, set([1, 2, 4, 6]))
        assert_raises(nx.NetworkXError,
                      nx.find_induced_nodes, self.non_chordal_G, 1, 5)

    def test_chordal_find_cliques(self):
        cliques = set([frozenset([9]), frozenset([7, 8]), frozenset([1, 2, 3]),
                       frozenset([2, 3, 4]), frozenset([3, 4, 5, 6])])
        assert_equal(nx.chordal_graph_cliques(self.chordal_G), cliques)

    def test_chordal_find_cliques_path(self):
        G = nx.path_graph(10)
        cliqueset = nx.chordal_graph_cliques(G)
        for (u, v) in G.edges():
            assert_true(frozenset([u, v]) in cliqueset
                        or frozenset([v, u]) in cliqueset)

    def test_chordal_find_cliquesCC(self):
        cliques = set([frozenset([1, 2, 3]), frozenset([2, 3, 4]),
                       frozenset([3, 4, 5, 6])])
        assert_equal(nx.chordal_graph_cliques(self.connected_chordal_G), cliques)
