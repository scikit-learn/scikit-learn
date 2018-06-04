from itertools import permutations

from nose.tools import assert_almost_equal
from nose.tools import assert_equal
from nose.tools import raises

import networkx as nx


class TestNeighborConnectivity(object):

    def test_degree_p4(self):
        G = nx.path_graph(4)
        answer = {1: 2.0, 2: 1.5}
        nd = nx.average_degree_connectivity(G)
        assert_equal(nd, answer)

        D = G.to_directed()
        answer = {2: 2.0, 4: 1.5}
        nd = nx.average_degree_connectivity(D)
        assert_equal(nd, answer)

        answer = {1: 2.0, 2: 1.5}
        D = G.to_directed()
        nd = nx.average_degree_connectivity(D, source='in', target='in')
        assert_equal(nd, answer)

        D = G.to_directed()
        nd = nx.average_degree_connectivity(D, source='in', target='in')
        assert_equal(nd, answer)

    def test_degree_p4_weighted(self):
        G = nx.path_graph(4)
        G[1][2]['weight'] = 4
        answer = {1: 2.0, 2: 1.8}
        nd = nx.average_degree_connectivity(G, weight='weight')
        assert_equal(nd, answer)
        answer = {1: 2.0, 2: 1.5}
        nd = nx.average_degree_connectivity(G)
        assert_equal(nd, answer)

        D = G.to_directed()
        answer = {2: 2.0, 4: 1.8}
        nd = nx.average_degree_connectivity(D, weight='weight')
        assert_equal(nd, answer)

        answer = {1: 2.0, 2: 1.8}
        D = G.to_directed()
        nd = nx.average_degree_connectivity(D, weight='weight', source='in',
                                            target='in')
        assert_equal(nd, answer)

        D = G.to_directed()
        nd = nx.average_degree_connectivity(D, source='in', target='out',
                                            weight='weight')
        assert_equal(nd, answer)

    def test_weight_keyword(self):
        G = nx.path_graph(4)
        G[1][2]['other'] = 4
        answer = {1: 2.0, 2: 1.8}
        nd = nx.average_degree_connectivity(G, weight='other')
        assert_equal(nd, answer)
        answer = {1: 2.0, 2: 1.5}
        nd = nx.average_degree_connectivity(G, weight=None)
        assert_equal(nd, answer)

        D = G.to_directed()
        answer = {2: 2.0, 4: 1.8}
        nd = nx.average_degree_connectivity(D, weight='other')
        assert_equal(nd, answer)

        answer = {1: 2.0, 2: 1.8}
        D = G.to_directed()
        nd = nx.average_degree_connectivity(D, weight='other', source='in',
                                            target='in')
        assert_equal(nd, answer)

        D = G.to_directed()
        nd = nx.average_degree_connectivity(D, weight='other', source='in',
                                            target='in')
        assert_equal(nd, answer)

    def test_degree_barrat(self):
        G = nx.star_graph(5)
        G.add_edges_from([(5, 6), (5, 7), (5, 8), (5, 9)])
        G[0][5]['weight'] = 5
        nd = nx.average_degree_connectivity(G)[5]
        assert_equal(nd, 1.8)
        nd = nx.average_degree_connectivity(G, weight='weight')[5]
        assert_almost_equal(nd, 3.222222, places=5)
        nd = nx.k_nearest_neighbors(G, weight='weight')[5]
        assert_almost_equal(nd, 3.222222, places=5)

    def test_zero_deg(self):
        G = nx.DiGraph()
        G.add_edge(1, 2)
        G.add_edge(1, 3)
        G.add_edge(1, 4)
        c = nx.average_degree_connectivity(G)
        assert_equal(c, {1: 0, 3: 1})
        c = nx.average_degree_connectivity(G, source='in', target='in')
        assert_equal(c, {0: 0, 1: 0})
        c = nx.average_degree_connectivity(G, source='in', target='out')
        assert_equal(c, {0: 0, 1: 3})
        c = nx.average_degree_connectivity(G, source='in', target='in+out')
        assert_equal(c, {0: 0, 1: 3})
        c = nx.average_degree_connectivity(G, source='out', target='out')
        assert_equal(c, {0: 0, 3: 0})
        c = nx.average_degree_connectivity(G, source='out', target='in')
        assert_equal(c, {0: 0, 3: 1})
        c = nx.average_degree_connectivity(G, source='out', target='in+out')
        assert_equal(c, {0: 0, 3: 1})

    def test_in_out_weight(self):
        G = nx.DiGraph()
        G.add_edge(1, 2, weight=1)
        G.add_edge(1, 3, weight=1)
        G.add_edge(3, 1, weight=1)
        for s, t in permutations(['in', 'out', 'in+out'], 2):
            c = nx.average_degree_connectivity(G, source=s, target=t)
            cw = nx.average_degree_connectivity(G, source=s, target=t,
                                                weight='weight')
            assert_equal(c, cw)

    @raises(ValueError)
    def test_invalid_source(self):
        G = nx.DiGraph()
        nx.average_degree_connectivity(G, source='bogus')

    @raises(ValueError)
    def test_invalid_target(self):
        G = nx.DiGraph()
        nx.average_degree_connectivity(G, target='bogus')

    def test_single_node(self):
        # TODO Is this really the intended behavior for providing a
        # single node as the argument `nodes`? Shouldn't the function
        # just return the connectivity value itself?
        G = nx.trivial_graph()
        conn = nx.average_degree_connectivity(G, nodes=0)
        assert_equal(conn, {0: 0})
