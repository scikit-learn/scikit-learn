# test_quality.py - unit tests for the quality module
#
# Copyright 2015 NetworkX developers.
#
# This file is part of NetworkX.
#
# NetworkX is distributed under a BSD license; see LICENSE.txt for more
# information.
"""Unit tests for the :mod:`networkx.algorithms.community.quality`
module.

"""
from __future__ import division

from nose.tools import assert_equal
from nose.tools import assert_almost_equal

import networkx as nx
from networkx import barbell_graph
from networkx.algorithms.community import coverage
from networkx.algorithms.community import modularity
from networkx.algorithms.community import performance
from networkx.algorithms.community.quality import inter_community_edges


class TestPerformance(object):
    """Unit tests for the :func:`performance` function."""

    def test_bad_partition(self):
        """Tests that a poor partition has a low performance measure."""
        G = barbell_graph(3, 0)
        partition = [{0, 1, 4}, {2, 3, 5}]
        assert_almost_equal(8 / 15, performance(G, partition))

    def test_good_partition(self):
        """Tests that a good partition has a high performance measure.

        """
        G = barbell_graph(3, 0)
        partition = [{0, 1, 2}, {3, 4, 5}]
        assert_almost_equal(14 / 15, performance(G, partition))


class TestCoverage(object):
    """Unit tests for the :func:`coverage` function."""

    def test_bad_partition(self):
        """Tests that a poor partition has a low coverage measure."""
        G = barbell_graph(3, 0)
        partition = [{0, 1, 4}, {2, 3, 5}]
        assert_almost_equal(3 / 7, coverage(G, partition))

    def test_good_partition(self):
        """Tests that a good partition has a high coverage measure."""
        G = barbell_graph(3, 0)
        partition = [{0, 1, 2}, {3, 4, 5}]
        assert_almost_equal(6 / 7, coverage(G, partition))


def test_modularity():
    G = nx.barbell_graph(3, 0)
    C = [{0, 1, 4}, {2, 3, 5}]
    assert_almost_equal(-16 / (14 ** 2), modularity(G, C))
    C = [{0, 1, 2}, {3, 4, 5}]
    assert_almost_equal((35 * 2) / (14 ** 2), modularity(G, C))


def test_inter_community_edges_with_digraphs():
    G = nx.complete_graph(2, create_using=nx.DiGraph())
    partition = [{0}, {1}]
    assert_equal(inter_community_edges(G, partition), 2)

    G = nx.complete_graph(10, create_using=nx.DiGraph())
    partition = [{0}, {1, 2}, {3, 4, 5}, {6, 7, 8, 9}]
    assert_equal(inter_community_edges(G, partition), 70)

    G = nx.cycle_graph(4, create_using=nx.DiGraph())
    partition = [{0, 1}, {2, 3}]
    assert_equal(inter_community_edges(G, partition), 2)
