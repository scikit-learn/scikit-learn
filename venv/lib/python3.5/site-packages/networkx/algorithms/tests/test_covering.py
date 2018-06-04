#    Copyright 2016-2018 NetworkX developers.
#    Copyright (C) 2016 by
#    Nishant Nikhil <nishantiam@gmail.com>
#    All rights reserved.
#    BSD license.

from nose.tools import assert_equal, assert_true, assert_false
import networkx as nx


class TestMinEdgeCover:
    """Tests for :func:`networkx.algorithms.min_edge_cover`"""

    def test_empty_graph(self):
        G = nx.Graph()
        assert_equal(nx.min_edge_cover(G), set())

    def test_graph_with_loop(self):
        G = nx.Graph()
        G.add_edge(0, 0)
        assert_equal(nx.min_edge_cover(G), {(0, 0)})

    def test_graph_single_edge(self):
        G = nx.Graph()
        G.add_edge(0, 1)
        assert_equal(nx.min_edge_cover(G),
                     {(0, 1)})

    def test_bipartite_explicit(self):
        G = nx.Graph()
        G.add_nodes_from([1, 2, 3, 4], bipartite=0)
        G.add_nodes_from(['a', 'b', 'c'], bipartite=1)
        G.add_edges_from([(1, 'a'), (1, 'b'), (2, 'b'),
                          (2, 'c'), (3, 'c'), (4, 'a')])
        min_cover = nx.min_edge_cover(G, nx.algorithms.bipartite.matching.
                                      eppstein_matching)
        min_cover2 = nx.min_edge_cover(G)
        assert_true(nx.is_edge_cover(G, min_cover))
        assert_equal(len(min_cover), 8)

    def test_complete_graph(self):
        G = nx.complete_graph(10)
        min_cover = nx.min_edge_cover(G)
        assert_true(nx.is_edge_cover(G, min_cover))
        assert_equal(len(min_cover), 5)


class TestIsEdgeCover:
    """Tests for :func:`networkx.algorithms.is_edge_cover`"""

    def test_empty_graph(self):
        G = nx.Graph()
        assert_true(nx.is_edge_cover(G, set()))

    def test_graph_with_loop(self):
        G = nx.Graph()
        G.add_edge(1, 1)
        assert_true(nx.is_edge_cover(G, {(1, 1)}))

    def test_graph_single_edge(self):
        G = nx.Graph()
        G.add_edge(0, 1)
        assert_true(nx.is_edge_cover(G, {(0, 0), (1, 1)}))
        assert_true(nx.is_edge_cover(G, {(0, 1), (1, 0)}))
        assert_true(nx.is_edge_cover(G, {(0, 1)}))
        assert_false(nx.is_edge_cover(G, {(0, 0)}))
