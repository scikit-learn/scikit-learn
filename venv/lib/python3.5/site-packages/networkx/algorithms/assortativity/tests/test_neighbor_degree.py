#!/usr/bin/env python
from nose.tools import *
import networkx as nx


class TestAverageNeighbor(object):

    def test_degree_p4(self):
        G = nx.path_graph(4)
        answer = {0: 2, 1: 1.5, 2: 1.5, 3: 2}
        nd = nx.average_neighbor_degree(G)
        assert_equal(nd, answer)

        D = G.to_directed()
        nd = nx.average_neighbor_degree(D)
        assert_equal(nd, answer)

        D = G.to_directed()
        nd = nx.average_neighbor_degree(D)
        assert_equal(nd, answer)

        D = G.to_directed()
        nd = nx.average_neighbor_degree(D, source='in', target='in')
        assert_equal(nd, answer)

    def test_degree_p4_weighted(self):
        G = nx.path_graph(4)
        G[1][2]['weight'] = 4
        answer = {0: 2, 1: 1.8, 2: 1.8, 3: 2}
        nd = nx.average_neighbor_degree(G, weight='weight')
        assert_equal(nd, answer)

        D = G.to_directed()
        nd = nx.average_neighbor_degree(D, weight='weight')
        assert_equal(nd, answer)

        D = G.to_directed()
        nd = nx.average_neighbor_degree(D, weight='weight')
        assert_equal(nd, answer)
        nd = nx.average_neighbor_degree(D, source='out', target='out',
                                        weight='weight')
        assert_equal(nd, answer)

        D = G.to_directed()
        nd = nx.average_neighbor_degree(D, source='in', target='in',
                                        weight='weight')
        assert_equal(nd, answer)

    def test_degree_k4(self):
        G = nx.complete_graph(4)
        answer = {0: 3, 1: 3, 2: 3, 3: 3}
        nd = nx.average_neighbor_degree(G)
        assert_equal(nd, answer)

        D = G.to_directed()
        nd = nx.average_neighbor_degree(D)
        assert_equal(nd, answer)

        D = G.to_directed()
        nd = nx.average_neighbor_degree(D)
        assert_equal(nd, answer)

        D = G.to_directed()
        nd = nx.average_neighbor_degree(D, source='in', target='in')
        assert_equal(nd, answer)

    def test_degree_k4_nodes(self):
        G = nx.complete_graph(4)
        answer = {1: 3.0, 2: 3.0}
        nd = nx.average_neighbor_degree(G, nodes=[1, 2])
        assert_equal(nd, answer)

    def test_degree_barrat(self):
        G = nx.star_graph(5)
        G.add_edges_from([(5, 6), (5, 7), (5, 8), (5, 9)])
        G[0][5]['weight'] = 5
        nd = nx.average_neighbor_degree(G)[5]
        assert_equal(nd, 1.8)
        nd = nx.average_neighbor_degree(G, weight='weight')[5]
        assert_almost_equal(nd, 3.222222, places=5)
