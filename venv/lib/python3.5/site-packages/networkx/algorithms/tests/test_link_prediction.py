import math

from functools import partial
from nose.tools import *

import networkx as nx


def _test_func(G, ebunch, expected, predict_func, **kwargs):
    result = predict_func(G, ebunch, **kwargs)
    exp_dict = dict((tuple(sorted([u, v])), score) for u, v, score in expected)
    res_dict = dict((tuple(sorted([u, v])), score) for u, v, score in result)

    assert_equal(len(exp_dict), len(res_dict))
    for p in exp_dict:
        assert_almost_equal(exp_dict[p], res_dict[p])


class TestResourceAllocationIndex():
    def setUp(self):
        self.func = nx.resource_allocation_index
        self.test = partial(_test_func, predict_func=self.func)

    def test_K5(self):
        G = nx.complete_graph(5)
        self.test(G, [(0, 1)], [(0, 1, 0.75)])

    def test_P3(self):
        G = nx.path_graph(3)
        self.test(G, [(0, 2)], [(0, 2, 0.5)])

    def test_S4(self):
        G = nx.star_graph(4)
        self.test(G, [(1, 2)], [(1, 2, 0.25)])

    @raises(nx.NetworkXNotImplemented)
    def test_digraph(self):
        G = nx.DiGraph()
        G.add_edges_from([(0, 1), (1, 2)])
        self.func(G, [(0, 2)])

    @raises(nx.NetworkXNotImplemented)
    def test_multigraph(self):
        G = nx.MultiGraph()
        G.add_edges_from([(0, 1), (1, 2)])
        self.func(G, [(0, 2)])

    @raises(nx.NetworkXNotImplemented)
    def test_multidigraph(self):
        G = nx.MultiDiGraph()
        G.add_edges_from([(0, 1), (1, 2)])
        self.func(G, [(0, 2)])

    def test_no_common_neighbor(self):
        G = nx.Graph()
        G.add_nodes_from([0, 1])
        self.test(G, [(0, 1)], [(0, 1, 0)])

    def test_equal_nodes(self):
        G = nx.complete_graph(4)
        self.test(G, [(0, 0)], [(0, 0, 1)])

    def test_all_nonexistent_edges(self):
        G = nx.Graph()
        G.add_edges_from([(0, 1), (0, 2), (2, 3)])
        self.test(G, None, [(0, 3, 0.5), (1, 2, 0.5), (1, 3, 0)])


class TestJaccardCoefficient():
    def setUp(self):
        self.func = nx.jaccard_coefficient
        self.test = partial(_test_func, predict_func=self.func)

    def test_K5(self):
        G = nx.complete_graph(5)
        self.test(G, [(0, 1)], [(0, 1, 0.6)])

    def test_P4(self):
        G = nx.path_graph(4)
        self.test(G, [(0, 2)], [(0, 2, 0.5)])

    @raises(nx.NetworkXNotImplemented)
    def test_digraph(self):
        G = nx.DiGraph()
        G.add_edges_from([(0, 1), (1, 2)])
        self.func(G, [(0, 2)])

    @raises(nx.NetworkXNotImplemented)
    def test_multigraph(self):
        G = nx.MultiGraph()
        G.add_edges_from([(0, 1), (1, 2)])
        self.func(G, [(0, 2)])

    @raises(nx.NetworkXNotImplemented)
    def test_multidigraph(self):
        G = nx.MultiDiGraph()
        G.add_edges_from([(0, 1), (1, 2)])
        self.func(G, [(0, 2)])

    def test_no_common_neighbor(self):
        G = nx.Graph()
        G.add_edges_from([(0, 1), (2, 3)])
        self.test(G, [(0, 2)], [(0, 2, 0)])

    def test_isolated_nodes(self):
        G = nx.Graph()
        G.add_nodes_from([0, 1])
        self.test(G, [(0, 1)], [(0, 1, 0)])

    def test_all_nonexistent_edges(self):
        G = nx.Graph()
        G.add_edges_from([(0, 1), (0, 2), (2, 3)])
        self.test(G, None, [(0, 3, 0.5), (1, 2, 0.5), (1, 3, 0)])


class TestAdamicAdarIndex():
    def setUp(self):
        self.func = nx.adamic_adar_index
        self.test = partial(_test_func, predict_func=self.func)

    def test_K5(self):
        G = nx.complete_graph(5)
        self.test(G, [(0, 1)], [(0, 1, 3 / math.log(4))])

    def test_P3(self):
        G = nx.path_graph(3)
        self.test(G, [(0, 2)], [(0, 2, 1 / math.log(2))])

    def test_S4(self):
        G = nx.star_graph(4)
        self.test(G, [(1, 2)], [(1, 2, 1 / math.log(4))])

    @raises(nx.NetworkXNotImplemented)
    def test_digraph(self):
        G = nx.DiGraph()
        G.add_edges_from([(0, 1), (1, 2)])
        self.func(G, [(0, 2)])

    @raises(nx.NetworkXNotImplemented)
    def test_multigraph(self):
        G = nx.MultiGraph()
        G.add_edges_from([(0, 1), (1, 2)])
        self.func(G, [(0, 2)])

    @raises(nx.NetworkXNotImplemented)
    def test_multidigraph(self):
        G = nx.MultiDiGraph()
        G.add_edges_from([(0, 1), (1, 2)])
        self.func(G, [(0, 2)])

    def test_no_common_neighbor(self):
        G = nx.Graph()
        G.add_nodes_from([0, 1])
        self.test(G, [(0, 1)], [(0, 1, 0)])

    def test_equal_nodes(self):
        G = nx.complete_graph(4)
        self.test(G, [(0, 0)], [(0, 0, 3 / math.log(3))])

    def test_all_nonexistent_edges(self):
        G = nx.Graph()
        G.add_edges_from([(0, 1), (0, 2), (2, 3)])
        self.test(G, None, [(0, 3, 1 / math.log(2)), (1, 2, 1 / math.log(2)),
                            (1, 3, 0)])


class TestPreferentialAttachment():
    def setUp(self):
        self.func = nx.preferential_attachment
        self.test = partial(_test_func, predict_func=self.func)

    def test_K5(self):
        G = nx.complete_graph(5)
        self.test(G, [(0, 1)], [(0, 1, 16)])

    def test_P3(self):
        G = nx.path_graph(3)
        self.test(G, [(0, 1)], [(0, 1, 2)])

    def test_S4(self):
        G = nx.star_graph(4)
        self.test(G, [(0, 2)], [(0, 2, 4)])

    @raises(nx.NetworkXNotImplemented)
    def test_digraph(self):
        G = nx.DiGraph()
        G.add_edges_from([(0, 1), (1, 2)])
        self.func(G, [(0, 2)])

    @raises(nx.NetworkXNotImplemented)
    def test_multigraph(self):
        G = nx.MultiGraph()
        G.add_edges_from([(0, 1), (1, 2)])
        self.func(G, [(0, 2)])

    @raises(nx.NetworkXNotImplemented)
    def test_multidigraph(self):
        G = nx.MultiDiGraph()
        G.add_edges_from([(0, 1), (1, 2)])
        self.func(G, [(0, 2)])

    def test_zero_degrees(self):
        G = nx.Graph()
        G.add_nodes_from([0, 1])
        self.test(G, [(0, 1)], [(0, 1, 0)])

    def test_all_nonexistent_edges(self):
        G = nx.Graph()
        G.add_edges_from([(0, 1), (0, 2), (2, 3)])
        self.test(G, None, [(0, 3, 2), (1, 2, 2), (1, 3, 1)])


class TestCNSoundarajanHopcroft():
    def setUp(self):
        self.func = nx.cn_soundarajan_hopcroft
        self.test = partial(_test_func, predict_func=self.func,
                            community='community')

    def test_K5(self):
        G = nx.complete_graph(5)
        G.nodes[0]['community'] = 0
        G.nodes[1]['community'] = 0
        G.nodes[2]['community'] = 0
        G.nodes[3]['community'] = 0
        G.nodes[4]['community'] = 1
        self.test(G, [(0, 1)], [(0, 1, 5)])

    def test_P3(self):
        G = nx.path_graph(3)
        G.nodes[0]['community'] = 0
        G.nodes[1]['community'] = 1
        G.nodes[2]['community'] = 0
        self.test(G, [(0, 2)], [(0, 2, 1)])

    def test_S4(self):
        G = nx.star_graph(4)
        G.nodes[0]['community'] = 1
        G.nodes[1]['community'] = 1
        G.nodes[2]['community'] = 1
        G.nodes[3]['community'] = 0
        G.nodes[4]['community'] = 0
        self.test(G, [(1, 2)], [(1, 2, 2)])

    @raises(nx.NetworkXNotImplemented)
    def test_digraph(self):
        G = nx.DiGraph()
        G.add_edges_from([(0, 1), (1, 2)])
        G.nodes[0]['community'] = 0
        G.nodes[1]['community'] = 0
        G.nodes[2]['community'] = 0
        self.func(G, [(0, 2)])

    @raises(nx.NetworkXNotImplemented)
    def test_multigraph(self):
        G = nx.MultiGraph()
        G.add_edges_from([(0, 1), (1, 2)])
        G.nodes[0]['community'] = 0
        G.nodes[1]['community'] = 0
        G.nodes[2]['community'] = 0
        self.func(G, [(0, 2)])

    @raises(nx.NetworkXNotImplemented)
    def test_multidigraph(self):
        G = nx.MultiDiGraph()
        G.add_edges_from([(0, 1), (1, 2)])
        G.nodes[0]['community'] = 0
        G.nodes[1]['community'] = 0
        G.nodes[2]['community'] = 0
        self.func(G, [(0, 2)])

    def test_no_common_neighbor(self):
        G = nx.Graph()
        G.add_nodes_from([0, 1])
        G.nodes[0]['community'] = 0
        G.nodes[1]['community'] = 0
        self.test(G, [(0, 1)], [(0, 1, 0)])

    def test_equal_nodes(self):
        G = nx.complete_graph(3)
        G.nodes[0]['community'] = 0
        G.nodes[1]['community'] = 0
        G.nodes[2]['community'] = 0
        self.test(G, [(0, 0)], [(0, 0, 4)])

    def test_different_community(self):
        G = nx.Graph()
        G.add_edges_from([(0, 1), (0, 2), (1, 3), (2, 3)])
        G.nodes[0]['community'] = 0
        G.nodes[1]['community'] = 0
        G.nodes[2]['community'] = 0
        G.nodes[3]['community'] = 1
        self.test(G, [(0, 3)], [(0, 3, 2)])

    @raises(nx.NetworkXAlgorithmError)
    def test_no_community_information(self):
        G = nx.complete_graph(5)
        list(self.func(G, [(0, 1)]))

    @raises(nx.NetworkXAlgorithmError)
    def test_insufficient_community_information(self):
        G = nx.Graph()
        G.add_edges_from([(0, 1), (0, 2), (1, 3), (2, 3)])
        G.nodes[0]['community'] = 0
        G.nodes[1]['community'] = 0
        G.nodes[3]['community'] = 0
        list(self.func(G, [(0, 3)]))

    def test_sufficient_community_information(self):
        G = nx.Graph()
        G.add_edges_from([(0, 1), (1, 2), (1, 3), (2, 4), (3, 4), (4, 5)])
        G.nodes[1]['community'] = 0
        G.nodes[2]['community'] = 0
        G.nodes[3]['community'] = 0
        G.nodes[4]['community'] = 0
        self.test(G, [(1, 4)], [(1, 4, 4)])

    def test_custom_community_attribute_name(self):
        G = nx.Graph()
        G.add_edges_from([(0, 1), (0, 2), (1, 3), (2, 3)])
        G.nodes[0]['cmty'] = 0
        G.nodes[1]['cmty'] = 0
        G.nodes[2]['cmty'] = 0
        G.nodes[3]['cmty'] = 1
        self.test(G, [(0, 3)], [(0, 3, 2)], community='cmty')

    def test_all_nonexistent_edges(self):
        G = nx.Graph()
        G.add_edges_from([(0, 1), (0, 2), (2, 3)])
        G.nodes[0]['community'] = 0
        G.nodes[1]['community'] = 1
        G.nodes[2]['community'] = 0
        G.nodes[3]['community'] = 0
        self.test(G, None, [(0, 3, 2), (1, 2, 1), (1, 3, 0)])


class TestRAIndexSoundarajanHopcroft():
    def setUp(self):
        self.func = nx.ra_index_soundarajan_hopcroft
        self.test = partial(_test_func, predict_func=self.func,
                            community='community')

    def test_K5(self):
        G = nx.complete_graph(5)
        G.nodes[0]['community'] = 0
        G.nodes[1]['community'] = 0
        G.nodes[2]['community'] = 0
        G.nodes[3]['community'] = 0
        G.nodes[4]['community'] = 1
        self.test(G, [(0, 1)], [(0, 1, 0.5)])

    def test_P3(self):
        G = nx.path_graph(3)
        G.nodes[0]['community'] = 0
        G.nodes[1]['community'] = 1
        G.nodes[2]['community'] = 0
        self.test(G, [(0, 2)], [(0, 2, 0)])

    def test_S4(self):
        G = nx.star_graph(4)
        G.nodes[0]['community'] = 1
        G.nodes[1]['community'] = 1
        G.nodes[2]['community'] = 1
        G.nodes[3]['community'] = 0
        G.nodes[4]['community'] = 0
        self.test(G, [(1, 2)], [(1, 2, 0.25)])

    @raises(nx.NetworkXNotImplemented)
    def test_digraph(self):
        G = nx.DiGraph()
        G.add_edges_from([(0, 1), (1, 2)])
        G.nodes[0]['community'] = 0
        G.nodes[1]['community'] = 0
        G.nodes[2]['community'] = 0
        self.func(G, [(0, 2)])

    @raises(nx.NetworkXNotImplemented)
    def test_multigraph(self):
        G = nx.MultiGraph()
        G.add_edges_from([(0, 1), (1, 2)])
        G.nodes[0]['community'] = 0
        G.nodes[1]['community'] = 0
        G.nodes[2]['community'] = 0
        self.func(G, [(0, 2)])

    @raises(nx.NetworkXNotImplemented)
    def test_multidigraph(self):
        G = nx.MultiDiGraph()
        G.add_edges_from([(0, 1), (1, 2)])
        G.nodes[0]['community'] = 0
        G.nodes[1]['community'] = 0
        G.nodes[2]['community'] = 0
        self.func(G, [(0, 2)])

    def test_no_common_neighbor(self):
        G = nx.Graph()
        G.add_nodes_from([0, 1])
        G.nodes[0]['community'] = 0
        G.nodes[1]['community'] = 0
        self.test(G, [(0, 1)], [(0, 1, 0)])

    def test_equal_nodes(self):
        G = nx.complete_graph(3)
        G.nodes[0]['community'] = 0
        G.nodes[1]['community'] = 0
        G.nodes[2]['community'] = 0
        self.test(G, [(0, 0)], [(0, 0, 1)])

    def test_different_community(self):
        G = nx.Graph()
        G.add_edges_from([(0, 1), (0, 2), (1, 3), (2, 3)])
        G.nodes[0]['community'] = 0
        G.nodes[1]['community'] = 0
        G.nodes[2]['community'] = 0
        G.nodes[3]['community'] = 1
        self.test(G, [(0, 3)], [(0, 3, 0)])

    @raises(nx.NetworkXAlgorithmError)
    def test_no_community_information(self):
        G = nx.complete_graph(5)
        list(self.func(G, [(0, 1)]))

    @raises(nx.NetworkXAlgorithmError)
    def test_insufficient_community_information(self):
        G = nx.Graph()
        G.add_edges_from([(0, 1), (0, 2), (1, 3), (2, 3)])
        G.nodes[0]['community'] = 0
        G.nodes[1]['community'] = 0
        G.nodes[3]['community'] = 0
        list(self.func(G, [(0, 3)]))

    def test_sufficient_community_information(self):
        G = nx.Graph()
        G.add_edges_from([(0, 1), (1, 2), (1, 3), (2, 4), (3, 4), (4, 5)])
        G.nodes[1]['community'] = 0
        G.nodes[2]['community'] = 0
        G.nodes[3]['community'] = 0
        G.nodes[4]['community'] = 0
        self.test(G, [(1, 4)], [(1, 4, 1)])

    def test_custom_community_attribute_name(self):
        G = nx.Graph()
        G.add_edges_from([(0, 1), (0, 2), (1, 3), (2, 3)])
        G.nodes[0]['cmty'] = 0
        G.nodes[1]['cmty'] = 0
        G.nodes[2]['cmty'] = 0
        G.nodes[3]['cmty'] = 1
        self.test(G, [(0, 3)], [(0, 3, 0)], community='cmty')

    def test_all_nonexistent_edges(self):
        G = nx.Graph()
        G.add_edges_from([(0, 1), (0, 2), (2, 3)])
        G.nodes[0]['community'] = 0
        G.nodes[1]['community'] = 1
        G.nodes[2]['community'] = 0
        G.nodes[3]['community'] = 0
        self.test(G, None, [(0, 3, 0.5), (1, 2, 0), (1, 3, 0)])


class TestWithinInterCluster():
    def setUp(self):
        self.delta = 0.001
        self.func = nx.within_inter_cluster
        self.test = partial(_test_func, predict_func=self.func,
                            delta=self.delta, community='community')

    def test_K5(self):
        G = nx.complete_graph(5)
        G.nodes[0]['community'] = 0
        G.nodes[1]['community'] = 0
        G.nodes[2]['community'] = 0
        G.nodes[3]['community'] = 0
        G.nodes[4]['community'] = 1
        self.test(G, [(0, 1)], [(0, 1, 2 / (1 + self.delta))])

    def test_P3(self):
        G = nx.path_graph(3)
        G.nodes[0]['community'] = 0
        G.nodes[1]['community'] = 1
        G.nodes[2]['community'] = 0
        self.test(G, [(0, 2)], [(0, 2, 0)])

    def test_S4(self):
        G = nx.star_graph(4)
        G.nodes[0]['community'] = 1
        G.nodes[1]['community'] = 1
        G.nodes[2]['community'] = 1
        G.nodes[3]['community'] = 0
        G.nodes[4]['community'] = 0
        self.test(G, [(1, 2)], [(1, 2, 1 / self.delta)])

    @raises(nx.NetworkXNotImplemented)
    def test_digraph(self):
        G = nx.DiGraph()
        G.add_edges_from([(0, 1), (1, 2)])
        G.nodes[0]['community'] = 0
        G.nodes[1]['community'] = 0
        G.nodes[2]['community'] = 0
        self.func(G, [(0, 2)])

    @raises(nx.NetworkXNotImplemented)
    def test_multigraph(self):
        G = nx.MultiGraph()
        G.add_edges_from([(0, 1), (1, 2)])
        G.nodes[0]['community'] = 0
        G.nodes[1]['community'] = 0
        G.nodes[2]['community'] = 0
        self.func(G, [(0, 2)])

    @raises(nx.NetworkXNotImplemented)
    def test_multidigraph(self):
        G = nx.MultiDiGraph()
        G.add_edges_from([(0, 1), (1, 2)])
        G.nodes[0]['community'] = 0
        G.nodes[1]['community'] = 0
        G.nodes[2]['community'] = 0
        self.func(G, [(0, 2)])

    def test_no_common_neighbor(self):
        G = nx.Graph()
        G.add_nodes_from([0, 1])
        G.nodes[0]['community'] = 0
        G.nodes[1]['community'] = 0
        self.test(G, [(0, 1)], [(0, 1, 0)])

    def test_equal_nodes(self):
        G = nx.complete_graph(3)
        G.nodes[0]['community'] = 0
        G.nodes[1]['community'] = 0
        G.nodes[2]['community'] = 0
        self.test(G, [(0, 0)], [(0, 0, 2 / self.delta)])

    def test_different_community(self):
        G = nx.Graph()
        G.add_edges_from([(0, 1), (0, 2), (1, 3), (2, 3)])
        G.nodes[0]['community'] = 0
        G.nodes[1]['community'] = 0
        G.nodes[2]['community'] = 0
        G.nodes[3]['community'] = 1
        self.test(G, [(0, 3)], [(0, 3, 0)])

    def test_no_inter_cluster_common_neighbor(self):
        G = nx.complete_graph(4)
        G.nodes[0]['community'] = 0
        G.nodes[1]['community'] = 0
        G.nodes[2]['community'] = 0
        G.nodes[3]['community'] = 0
        self.test(G, [(0, 3)], [(0, 3, 2 / self.delta)])

    @raises(nx.NetworkXAlgorithmError)
    def test_no_community_information(self):
        G = nx.complete_graph(5)
        list(self.func(G, [(0, 1)]))

    @raises(nx.NetworkXAlgorithmError)
    def test_insufficient_community_information(self):
        G = nx.Graph()
        G.add_edges_from([(0, 1), (0, 2), (1, 3), (2, 3)])
        G.nodes[0]['community'] = 0
        G.nodes[1]['community'] = 0
        G.nodes[3]['community'] = 0
        list(self.func(G, [(0, 3)]))

    def test_sufficient_community_information(self):
        G = nx.Graph()
        G.add_edges_from([(0, 1), (1, 2), (1, 3), (2, 4), (3, 4), (4, 5)])
        G.nodes[1]['community'] = 0
        G.nodes[2]['community'] = 0
        G.nodes[3]['community'] = 0
        G.nodes[4]['community'] = 0
        self.test(G, [(1, 4)], [(1, 4, 2 / self.delta)])

    @raises(nx.NetworkXAlgorithmError)
    def test_zero_delta(self):
        G = nx.complete_graph(3)
        G.nodes[0]['community'] = 0
        G.nodes[1]['community'] = 0
        G.nodes[2]['community'] = 0
        list(self.func(G, [(0, 1)], 0))

    @raises(nx.NetworkXAlgorithmError)
    def test_negative_delta(self):
        G = nx.complete_graph(3)
        G.nodes[0]['community'] = 0
        G.nodes[1]['community'] = 0
        G.nodes[2]['community'] = 0
        list(self.func(G, [(0, 1)], -0.5))

    def test_custom_community_attribute_name(self):
        G = nx.complete_graph(4)
        G.nodes[0]['cmty'] = 0
        G.nodes[1]['cmty'] = 0
        G.nodes[2]['cmty'] = 0
        G.nodes[3]['cmty'] = 0
        self.test(G, [(0, 3)], [(0, 3, 2 / self.delta)], community='cmty')

    def test_all_nonexistent_edges(self):
        G = nx.Graph()
        G.add_edges_from([(0, 1), (0, 2), (2, 3)])
        G.nodes[0]['community'] = 0
        G.nodes[1]['community'] = 1
        G.nodes[2]['community'] = 0
        G.nodes[3]['community'] = 0
        self.test(G, None, [(0, 3, 1 / self.delta), (1, 2, 0), (1, 3, 0)])
