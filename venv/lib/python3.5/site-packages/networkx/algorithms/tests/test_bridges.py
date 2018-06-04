# test_bridges.py - unit tests for bridge-finding algorithms
#
# Copyright 2004-2018 NetworkX developers.
#
# This file is part of NetworkX.
#
# NetworkX is distributed under a BSD license; see LICENSE.txt for more
# information.
"""Unit tests for bridge-finding algorithms."""
from unittest import TestCase
from nose.tools import assert_equal, assert_in

import networkx as nx


class TestBridges(TestCase):
    """Unit tests for the bridge-finding function."""

    def test_single_bridge(self):
        edges = [
            # DFS tree edges.
            (1, 2), (2, 3), (3, 4), (3, 5), (5, 6), (6, 7), (7, 8), (5, 9),
            (9, 10),
            # Nontree edges.
            (1, 3), (1, 4), (2, 5), (5, 10), (6, 8)
        ]
        G = nx.Graph(edges)
        source = 1
        bridges = list(nx.bridges(G, source))
        self.assertEqual(bridges, [(5, 6)])

    def test_barbell_graph(self):
        # The (3, 0) barbell graph has two triangles joined by a single edge.
        G = nx.barbell_graph(3, 0)
        source = 0
        bridges = list(nx.bridges(G, source))
        self.assertEqual(bridges, [(2, 3)])


class TestLocalBridges(TestCase):
    """Unit tests for the local_bridge function."""

    def setUp(self):
        self.BB = nx.barbell_graph(4, 0)
        self.square = nx.cycle_graph(4)
        self.tri = nx.cycle_graph(3)

    def test_nospan(self):
        expected = {(3, 4), (4, 3)}
        assert_in(next(nx.local_bridges(self.BB, with_span=False)), expected)
        assert_equal(set(nx.local_bridges(self.square, with_span=False)), self.square.edges)
        assert_equal(list(nx.local_bridges(self.tri, with_span=False)), [])

    def test_no_weight(self):
        inf = float('inf')
        expected = {(3, 4, inf), (4, 3, inf)}
        assert_in(next(nx.local_bridges(self.BB)), expected)
        expected = {(u, v, 3) for u, v, in self.square.edges}
        assert_equal(set(nx.local_bridges(self.square)), expected)
        assert_equal(list(nx.local_bridges(self.tri)), [])

    def test_weight(self):
        inf = float('inf')
        G = self.square.copy()

        G.edges[1, 2]['weight'] = 2
        expected = {(u, v, 5 - wt) for u, v, wt in G.edges(data='weight', default=1)}
        assert_equal(set(nx.local_bridges(G, weight='weight')), expected)

        expected = {(u, v, 6) for u, v in G.edges}
        lb = nx.local_bridges(G, weight=lambda u, v, d: 2)
        assert_equal(set(lb), expected)
