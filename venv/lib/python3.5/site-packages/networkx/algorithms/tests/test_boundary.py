# test_boundary.py - unit tests for the boundary module
#
# Copyright 2015 NetworkX developers.
#
# This file is part of NetworkX.
#
# NetworkX is distributed under a BSD license; see LICENSE.txt for more
# information.
"""Unit tests for the :mod:`networkx.algorithms.boundary` module."""
from __future__ import division

from itertools import combinations

from nose.tools import assert_almost_equals
from nose.tools import assert_equal

import networkx as nx
from networkx.testing.utils import *
from networkx import convert_node_labels_to_integers as cnlti


class TestNodeBoundary(object):
    """Unit tests for the :func:`~networkx.node_boundary` function."""

    def test_null_graph(self):
        """Tests that the null graph has empty node boundaries."""
        null = nx.null_graph()
        assert_equal(nx.node_boundary(null, []), set())
        assert_equal(nx.node_boundary(null, [], []), set())
        assert_equal(nx.node_boundary(null, [1, 2, 3]), set())
        assert_equal(nx.node_boundary(null, [1, 2, 3], [4, 5, 6]), set())
        assert_equal(nx.node_boundary(null, [1, 2, 3], [3, 4, 5]), set())

    def test_path_graph(self):
        P10 = cnlti(nx.path_graph(10), first_label=1)
        assert_equal(nx.node_boundary(P10, []), set())
        assert_equal(nx.node_boundary(P10, [], []), set())
        assert_equal(nx.node_boundary(P10, [1, 2, 3]), {4})
        assert_equal(nx.node_boundary(P10, [4, 5, 6]), {3, 7})
        assert_equal(nx.node_boundary(P10, [3, 4, 5, 6, 7]), {2, 8})
        assert_equal(nx.node_boundary(P10, [8, 9, 10]), {7})
        assert_equal(nx.node_boundary(P10, [4, 5, 6], [9, 10]), set())

    def test_complete_graph(self):
        K10 = cnlti(nx.complete_graph(10), first_label=1)
        assert_equal(nx.node_boundary(K10, []), set())
        assert_equal(nx.node_boundary(K10, [], []), set())
        assert_equal(nx.node_boundary(K10, [1, 2, 3]), {4, 5, 6, 7, 8, 9, 10})
        assert_equal(nx.node_boundary(K10, [4, 5, 6]), {1, 2, 3, 7, 8, 9, 10})
        assert_equal(nx.node_boundary(K10, [3, 4, 5, 6, 7]), {1, 2, 8, 9, 10})
        assert_equal(nx.node_boundary(K10, [4, 5, 6], []), set())
        assert_equal(nx.node_boundary(K10, K10), set())
        assert_equal(nx.node_boundary(K10, [1, 2, 3], [3, 4, 5]), {4, 5})

    def test_petersen(self):
        """Check boundaries in the petersen graph

        cheeger(G,k)=min(|bdy(S)|/|S| for |S|=k, 0<k<=|V(G)|/2)

        """

        def cheeger(G, k):
            return min(len(nx.node_boundary(G, nn)) / k
                       for nn in combinations(G, k))

        P = nx.petersen_graph()
        assert_almost_equals(cheeger(P, 1), 3.00, places=2)
        assert_almost_equals(cheeger(P, 2), 2.00, places=2)
        assert_almost_equals(cheeger(P, 3), 1.67, places=2)
        assert_almost_equals(cheeger(P, 4), 1.00, places=2)
        assert_almost_equals(cheeger(P, 5), 0.80, places=2)

    def test_directed(self):
        """Tests the node boundary of a directed graph."""
        G = nx.DiGraph([(0, 1), (1, 2), (2, 3), (3, 4), (4, 0)])
        S = {0, 1}
        boundary = nx.node_boundary(G, S)
        expected = {2}
        assert_equal(boundary, expected)

    def test_multigraph(self):
        """Tests the node boundary of a multigraph."""
        G = nx.MultiGraph(list(nx.cycle_graph(5).edges()) * 2)
        S = {0, 1}
        boundary = nx.node_boundary(G, S)
        expected = {2, 4}
        assert_equal(boundary, expected)

    def test_multidigraph(self):
        """Tests the edge boundary of a multdiigraph."""
        edges = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 0)]
        G = nx.MultiDiGraph(edges * 2)
        S = {0, 1}
        boundary = nx.node_boundary(G, S)
        expected = {2}
        assert_equal(boundary, expected)


class TestEdgeBoundary(object):
    """Unit tests for the :func:`~networkx.edge_boundary` function."""

    def test_null_graph(self):
        null = nx.null_graph()
        assert_equal(list(nx.edge_boundary(null, [])), [])
        assert_equal(list(nx.edge_boundary(null, [], [])), [])
        assert_equal(list(nx.edge_boundary(null, [1, 2, 3])), [])
        assert_equal(list(nx.edge_boundary(null, [1, 2, 3], [4, 5, 6])), [])
        assert_equal(list(nx.edge_boundary(null, [1, 2, 3], [3, 4, 5])), [])

    def test_path_graph(self):
        P10 = cnlti(nx.path_graph(10), first_label=1)
        assert_equal(list(nx.edge_boundary(P10, [])), [])
        assert_equal(list(nx.edge_boundary(P10, [], [])), [])
        assert_equal(list(nx.edge_boundary(P10, [1, 2, 3])), [(3, 4)])
        assert_equal(sorted(nx.edge_boundary(P10, [4, 5, 6])),
                     [(4, 3), (6, 7)])
        assert_equal(sorted(nx.edge_boundary(P10, [3, 4, 5, 6, 7])),
                     [(3, 2), (7, 8)])
        assert_equal(list(nx.edge_boundary(P10, [8, 9, 10])), [(8, 7)])
        assert_equal(sorted(nx.edge_boundary(P10, [4, 5, 6], [9, 10])), [])
        assert_equal(list(nx.edge_boundary(P10, [1, 2, 3], [3, 4, 5])),
                     [(2, 3), (3, 4)])

    def test_complete_graph(self):
        K10 = cnlti(nx.complete_graph(10), first_label=1)

        def ilen(iterable): return sum(1 for i in iterable)
        assert_equal(list(nx.edge_boundary(K10, [])), [])
        assert_equal(list(nx.edge_boundary(K10, [], [])), [])
        assert_equal(ilen(nx.edge_boundary(K10, [1, 2, 3])), 21)
        assert_equal(ilen(nx.edge_boundary(K10, [4, 5, 6, 7])), 24)
        assert_equal(ilen(nx.edge_boundary(K10, [3, 4, 5, 6, 7])), 25)
        assert_equal(ilen(nx.edge_boundary(K10, [8, 9, 10])), 21)
        assert_edges_equal(nx.edge_boundary(K10, [4, 5, 6], [9, 10]),
                           [(4, 9), (4, 10), (5, 9), (5, 10), (6, 9), (6, 10)])
        assert_edges_equal(nx.edge_boundary(K10, [1, 2, 3], [3, 4, 5]),
                           [(1, 3), (1, 4), (1, 5), (2, 3), (2, 4),
                            (2, 5), (3, 4), (3, 5)])

    def test_directed(self):
        """Tests the edge boundary of a directed graph."""
        G = nx.DiGraph([(0, 1), (1, 2), (2, 3), (3, 4), (4, 0)])
        S = {0, 1}
        boundary = list(nx.edge_boundary(G, S))
        expected = [(1, 2)]
        assert_equal(boundary, expected)

    def test_multigraph(self):
        """Tests the edge boundary of a multigraph."""
        G = nx.MultiGraph(list(nx.cycle_graph(5).edges()) * 2)
        S = {0, 1}
        boundary = list(nx.edge_boundary(G, S))
        expected = [(0, 4), (0, 4), (1, 2), (1, 2)]
        assert_equal(boundary, expected)

    def test_multidigraph(self):
        """Tests the edge boundary of a multdiigraph."""
        edges = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 0)]
        G = nx.MultiDiGraph(edges * 2)
        S = {0, 1}
        boundary = list(nx.edge_boundary(G, S))
        expected = [(1, 2), (1, 2)]
        assert_equal(boundary, expected)
