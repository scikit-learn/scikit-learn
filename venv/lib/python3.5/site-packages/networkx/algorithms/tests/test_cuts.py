# test_cuts.py - unit tests for the cuts module
#
# Copyright 2015 NetworkX developers.
#
# This file is part of NetworkX.
#
# NetworkX is distributed under a BSD license; see LICENSE.txt for more
# information.
"""Unit tests for the :mod:`networkx.algorithms.cuts` module."""
from __future__ import division

from nose.tools import assert_equal

import networkx as nx


class TestCutSize(object):
    """Unit tests for the :func:`~networkx.cut_size` function."""

    def test_symmetric(self):
        """Tests that the cut size is symmetric."""
        G = nx.barbell_graph(3, 0)
        S = {0, 1, 4}
        T = {2, 3, 5}
        assert_equal(nx.cut_size(G, S, T), 4)
        assert_equal(nx.cut_size(G, T, S), 4)

    def test_single_edge(self):
        """Tests for a cut of a single edge."""
        G = nx.barbell_graph(3, 0)
        S = {0, 1, 2}
        T = {3, 4, 5}
        assert_equal(nx.cut_size(G, S, T), 1)
        assert_equal(nx.cut_size(G, T, S), 1)

    def test_directed(self):
        """Tests that each directed edge is counted once in the cut."""
        G = nx.barbell_graph(3, 0).to_directed()
        S = {0, 1, 2}
        T = {3, 4, 5}
        assert_equal(nx.cut_size(G, S, T), 2)
        assert_equal(nx.cut_size(G, T, S), 2)

    def test_directed_symmetric(self):
        """Tests that a cut in a directed graph is symmetric."""
        G = nx.barbell_graph(3, 0).to_directed()
        S = {0, 1, 4}
        T = {2, 3, 5}
        assert_equal(nx.cut_size(G, S, T), 8)
        assert_equal(nx.cut_size(G, T, S), 8)

    def test_multigraph(self):
        """Tests that parallel edges are each counted for a cut."""
        G = nx.MultiGraph(['ab', 'ab'])
        assert_equal(nx.cut_size(G, {'a'}, {'b'}), 2)


class TestVolume(object):
    """Unit tests for the :func:`~networkx.volume` function."""

    def test_graph(self):
        G = nx.cycle_graph(4)
        assert_equal(nx.volume(G, {0, 1}), 4)

    def test_digraph(self):
        G = nx.DiGraph([(0, 1), (1, 2), (2, 3), (3, 0)])
        assert_equal(nx.volume(G, {0, 1}), 2)

    def test_multigraph(self):
        edges = list(nx.cycle_graph(4).edges())
        G = nx.MultiGraph(edges * 2)
        assert_equal(nx.volume(G, {0, 1}), 8)

    def test_multidigraph(self):
        edges = [(0, 1), (1, 2), (2, 3), (3, 0)]
        G = nx.MultiDiGraph(edges * 2)
        assert_equal(nx.volume(G, {0, 1}), 4)


class TestNormalizedCutSize(object):
    """Unit tests for the :func:`~networkx.normalized_cut_size`
    function.

    """

    def test_graph(self):
        G = nx.path_graph(4)
        S = {1, 2}
        T = set(G) - S
        size = nx.normalized_cut_size(G, S, T)
        # The cut looks like this: o-{-o--o-}-o
        expected = 2 * ((1 / 4) + (1 / 2))
        assert_equal(expected, size)

    def test_directed(self):
        G = nx.DiGraph([(0, 1), (1, 2), (2, 3)])
        S = {1, 2}
        T = set(G) - S
        size = nx.normalized_cut_size(G, S, T)
        # The cut looks like this: o-{->o-->o-}->o
        expected = 2 * ((1 / 2) + (1 / 1))
        assert_equal(expected, size)


class TestConductance(object):
    """Unit tests for the :func:`~networkx.conductance` function."""

    def test_graph(self):
        G = nx.barbell_graph(5, 0)
        # Consider the singleton sets containing the "bridge" nodes.
        # There is only one cut edge, and each set has volume five.
        S = {4}
        T = {5}
        conductance = nx.conductance(G, S, T)
        expected = 1 / 5
        assert_equal(expected, conductance)


class TestEdgeExpansion(object):
    """Unit tests for the :func:`~networkx.edge_expansion` function."""

    def test_graph(self):
        G = nx.barbell_graph(5, 0)
        S = set(range(5))
        T = set(G) - S
        expansion = nx.edge_expansion(G, S, T)
        expected = 1 / 5
        assert_equal(expected, expansion)


class TestNodeExpansion(object):
    """Unit tests for the :func:`~networkx.node_expansion` function.

    """

    def test_graph(self):
        G = nx.path_graph(8)
        S = {3, 4, 5}
        expansion = nx.node_expansion(G, S)
        # The neighborhood of S has cardinality five, and S has
        # cardinality three.
        expected = 5 / 3
        assert_equal(expected, expansion)


class TestBoundaryExpansion(object):
    """Unit tests for the :func:`~networkx.boundary_expansion` function.

    """

    def test_graph(self):
        G = nx.complete_graph(10)
        S = set(range(4))
        expansion = nx.boundary_expansion(G, S)
        # The node boundary of S has cardinality six, and S has
        # cardinality three.
        expected = 6 / 4
        assert_equal(expected, expansion)


class TestMixingExpansion(object):
    """Unit tests for the :func:`~networkx.mixing_expansion` function.

    """

    def test_graph(self):
        G = nx.barbell_graph(5, 0)
        S = set(range(5))
        T = set(G) - S
        expansion = nx.mixing_expansion(G, S, T)
        # There is one cut edge, and the total number of edges in the
        # graph is twice the total number of edges in a clique of size
        # five, plus one more for the bridge.
        expected = 1 / (2 * (5 * 4 + 1))
        assert_equal(expected, expansion)
