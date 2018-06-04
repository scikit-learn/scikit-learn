# test_stochastic.py - unit tests for the stochastic module
#
# Copyright 2010, 2011, 2012, 2013, 2014, 2015 NetworkX developers.
#
# This file is part of NetworkX.
#
# NetworkX is distributed under a BSD license; see LICENSE.txt for more
# information.
"""Unit tests for the :mod:`networkx.generators.stochastic` module."""
from nose.tools import assert_true, assert_equal, raises
import networkx as nx


class TestStochasticGraph(object):
    """Unit tests for the :func:`~networkx.stochastic_graph` function.

    """

    def test_default_weights(self):
        G = nx.DiGraph()
        G.add_edge(0, 1)
        G.add_edge(0, 2)
        S = nx.stochastic_graph(G)
        assert_true(nx.is_isomorphic(G, S))
        assert_equal(sorted(S.edges(data=True)),
                     [(0, 1, {'weight': 0.5}), (0, 2, {'weight': 0.5})])

    def test_in_place(self):
        """Tests for an in-place reweighting of the edges of the graph.

        """
        G = nx.DiGraph()
        G.add_edge(0, 1, weight=1)
        G.add_edge(0, 2, weight=1)
        nx.stochastic_graph(G, copy=False)
        assert_equal(sorted(G.edges(data=True)),
                     [(0, 1, {'weight': 0.5}), (0, 2, {'weight': 0.5})])

    def test_arbitrary_weights(self):
        G = nx.DiGraph()
        G.add_edge(0, 1, weight=1)
        G.add_edge(0, 2, weight=1)
        S = nx.stochastic_graph(G)
        assert_equal(sorted(S.edges(data=True)),
                     [(0, 1, {'weight': 0.5}), (0, 2, {'weight': 0.5})])

    def test_multidigraph(self):
        G = nx.MultiDiGraph()
        G.add_edges_from([(0, 1), (0, 1), (0, 2), (0, 2)])
        S = nx.stochastic_graph(G)
        d = dict(weight=0.25)
        assert_equal(sorted(S.edges(data=True)),
                     [(0, 1, d), (0, 1, d), (0, 2, d), (0, 2, d)])

    @raises(nx.NetworkXNotImplemented)
    def test_graph_disallowed(self):
        nx.stochastic_graph(nx.Graph())

    @raises(nx.NetworkXNotImplemented)
    def test_multigraph_disallowed(self):
        nx.stochastic_graph(nx.MultiGraph())
