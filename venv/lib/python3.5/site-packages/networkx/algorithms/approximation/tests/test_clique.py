# test_clique.py - unit tests for the approximation.clique module
#
# Copyright 2015 NetworkX developers.
#
# This file is part of NetworkX.
#
# NetworkX is distributed under a BSD license; see LICENSE.txt for more
# information.
"""Unit tests for the :mod:`networkx.algorithms.approximation.clique`
module.

"""
from __future__ import division

from nose.tools import assert_greater
from nose.tools import assert_true
from nose.tools import assert_equal

import networkx as nx
from networkx.algorithms.approximation import max_clique
from networkx.algorithms.approximation import clique_removal
from networkx.algorithms.approximation import large_clique_size


def is_independent_set(G, nodes):
    """Returns True if and only if `nodes` is a clique in `G`.

    `G` is a NetworkX graph. `nodes` is an iterable of nodes in
    `G`.

    """
    return G.subgraph(nodes).number_of_edges() == 0


def is_clique(G, nodes):
    """Returns True if and only if `nodes` is an independent set
    in `G`.

    `G` is an undirected simple graph. `nodes` is an iterable of
    nodes in `G`.

    """
    H = G.subgraph(nodes)
    n = len(H)
    return H.number_of_edges() == n * (n - 1) // 2


class TestCliqueRemoval(object):
    """Unit tests for the
    :func:`~networkx.algorithms.approximation.clique_removal` function.

    """

    def test_trivial_graph(self):
        G = nx.trivial_graph()
        independent_set, cliques = clique_removal(G)
        assert_true(is_independent_set(G, independent_set))
        assert_true(all(is_clique(G, clique) for clique in cliques))
        # In fact, we should only have 1-cliques, that is, singleton nodes.
        assert_true(all(len(clique) == 1 for clique in cliques))

    def test_complete_graph(self):
        G = nx.complete_graph(10)
        independent_set, cliques = clique_removal(G)
        assert_true(is_independent_set(G, independent_set))
        assert_true(all(is_clique(G, clique) for clique in cliques))

    def test_barbell_graph(self):
        G = nx.barbell_graph(10, 5)
        independent_set, cliques = clique_removal(G)
        assert_true(is_independent_set(G, independent_set))
        assert_true(all(is_clique(G, clique) for clique in cliques))


class TestMaxClique(object):
    """Unit tests for the :func:`networkx.algorithms.approximation.max_clique`
    function.

    """

    def test_null_graph(self):
        G = nx.null_graph()
        assert_equal(len(max_clique(G)), 0)

    def test_complete_graph(self):
        graph = nx.complete_graph(30)
        # this should return the entire graph
        mc = max_clique(graph)
        assert_equal(30, len(mc))

    def test_maximal_by_cardinality(self):
        """Tests that the maximal clique is computed according to maximum
        cardinality of the sets.

        For more information, see pull request #1531.

        """
        G = nx.complete_graph(5)
        G.add_edge(4, 5)
        clique = max_clique(G)
        assert_greater(len(clique), 1)

        G = nx.lollipop_graph(30, 2)
        clique = max_clique(G)
        assert_greater(len(clique), 2)


def test_large_clique_size():
    G = nx.complete_graph(9)
    nx.add_cycle(G, [9, 10, 11])
    G.add_edge(8, 9)
    G.add_edge(1, 12)
    G.add_node(13)

    assert_equal(large_clique_size(G), 9)
    G.remove_node(5)
    assert_equal(large_clique_size(G), 8)
    G.remove_edge(2, 3)
    assert_equal(large_clique_size(G), 7)
