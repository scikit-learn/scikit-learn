# test_tournament.py - unit tests for the tournament module
#
# Copyright 2015 NetworkX developers.
#
# This file is part of NetworkX.
#
# NetworkX is distributed under a BSD license; see LICENSE.txt for more
# information.
"""Unit tests for the :mod:`networkx.algorithms.tournament` module."""
from itertools import combinations

from nose.tools import assert_equal
from nose.tools import assert_false
from nose.tools import assert_true

from networkx import DiGraph
from networkx.algorithms.tournament import is_reachable
from networkx.algorithms.tournament import is_strongly_connected
from networkx.algorithms.tournament import is_tournament
from networkx.algorithms.tournament import random_tournament
from networkx.algorithms.tournament import hamiltonian_path


class TestIsTournament(object):
    """Unit tests for the :func:`networkx.tournament.is_tournament`
    function.

    """

    def test_is_tournament(self):
        G = DiGraph()
        G.add_edges_from([(0, 1), (1, 2), (2, 3), (3, 0), (1, 3), (0, 2)])
        assert_true(is_tournament(G))

    def test_self_loops(self):
        """A tournament must have no self-loops."""
        G = DiGraph()
        G.add_edges_from([(0, 1), (1, 2), (2, 3), (3, 0), (1, 3), (0, 2)])
        G.add_edge(0, 0)
        assert_false(is_tournament(G))

    def test_missing_edges(self):
        """A tournament must not have any pair of nodes without at least
        one edge joining the pair.

        """
        G = DiGraph()
        G.add_edges_from([(0, 1), (1, 2), (2, 3), (3, 0), (1, 3)])
        assert_false(is_tournament(G))

    def test_bidirectional_edges(self):
        """A tournament must not have any pair of nodes with greater
        than one edge joining the pair.

        """
        G = DiGraph()
        G.add_edges_from([(0, 1), (1, 2), (2, 3), (3, 0), (1, 3), (0, 2)])
        G.add_edge(1, 0)
        assert_false(is_tournament(G))


class TestRandomTournament(object):
    """Unit tests for the :func:`networkx.tournament.random_tournament`
    function.

    """

    def test_graph_is_tournament(self):
        for n in range(10):
            G = random_tournament(5)
            assert_true(is_tournament(G))


class TestHamiltonianPath(object):
    """Unit tests for the :func:`networkx.tournament.hamiltonian_path`
    function.

    """

    def test_path_is_hamiltonian(self):
        G = DiGraph()
        G.add_edges_from([(0, 1), (1, 2), (2, 3), (3, 0), (1, 3), (0, 2)])
        path = hamiltonian_path(G)
        assert_equal(len(path), 4)
        assert_true(all(v in G[u] for u, v in zip(path, path[1:])))

    def test_hamiltonian_cycle(self):
        """Tests that :func:`networkx.tournament.hamiltonian_path`
        returns a Hamiltonian cycle when provided a strongly connected
        tournament.

        """
        G = DiGraph()
        G.add_edges_from([(0, 1), (1, 2), (2, 3), (3, 0), (1, 3), (0, 2)])
        path = hamiltonian_path(G)
        assert_equal(len(path), 4)
        assert_true(all(v in G[u] for u, v in zip(path, path[1:])))
        assert_true(path[0] in G[path[-1]])


class TestReachability(object):
    """Unit tests for the :func:`networkx.tournament.is_reachable`
    function.

    """

    def test_reachable_pair(self):
        """Tests for a reachable pair of nodes."""
        G = DiGraph([(0, 1), (1, 2), (2, 0)])
        assert_true(is_reachable(G, 0, 2))

    def test_same_node_is_reachable(self):
        """Tests that a node is always reachable from itself."""
        # G is an arbitrary tournament on ten nodes.
        G = DiGraph(sorted(p) for p in combinations(range(10), 2))
        assert_true(all(is_reachable(G, v, v) for v in G))

    def test_unreachable_pair(self):
        """Tests for an unreachable pair of nodes."""
        G = DiGraph([(0, 1), (0, 2), (1, 2)])
        assert_false(is_reachable(G, 1, 0))


class TestStronglyConnected(object):
    """Unit tests for the
    :func:`networkx.tournament.is_strongly_connected` function.

    """

    def test_is_strongly_connected(self):
        """Tests for a strongly connected tournament."""
        G = DiGraph([(0, 1), (1, 2), (2, 0)])
        assert_true(is_strongly_connected(G))

    def test_not_strongly_connected(self):
        """Tests for a tournament that is not strongly connected."""
        G = DiGraph([(0, 1), (0, 2), (1, 2)])
        assert_false(is_strongly_connected(G))
