from itertools import combinations
from nose.tools import assert_equal, assert_true, raises

import networkx as nx
from networkx.algorithms.flow import boykov_kolmogorov
from networkx.algorithms.flow import edmonds_karp
from networkx.algorithms.flow import preflow_push
from networkx.algorithms.flow import shortest_augmenting_path
from networkx.algorithms.flow import dinitz

flow_funcs = [
    boykov_kolmogorov,
    dinitz,
    edmonds_karp,
    preflow_push,
    shortest_augmenting_path,
]


class TestGomoryHuTree:

    def minimum_edge_weight(self, T, u, v):
        path = nx.shortest_path(T, u, v, weight='weight')
        return min((T[u][v]['weight'], (u, v)) for (u, v) in zip(path, path[1:]))

    def compute_cutset(self, G, T_orig, edge):
        T = T_orig.copy()
        T.remove_edge(*edge)
        U, V = list(nx.connected_components(T))
        cutset = set()
        for x, nbrs in ((n, G[n]) for n in U):
            cutset.update((x, y) for y in nbrs if y in V)
        return cutset

    def test_default_flow_function_karate_club_graph(self):
        G = nx.karate_club_graph()
        nx.set_edge_attributes(G, 1, 'capacity')
        T = nx.gomory_hu_tree(G)
        assert_true(nx.is_tree(T))
        for u, v in combinations(G, 2):
            cut_value, edge = self.minimum_edge_weight(T, u, v)
            assert_equal(nx.minimum_cut_value(G, u, v),
                         cut_value)

    def test_karate_club_graph(self):
        G = nx.karate_club_graph()
        nx.set_edge_attributes(G, 1, 'capacity')
        for flow_func in flow_funcs:
            T = nx.gomory_hu_tree(G, flow_func=flow_func)
            assert_true(nx.is_tree(T))
            for u, v in combinations(G, 2):
                cut_value, edge = self.minimum_edge_weight(T, u, v)
                assert_equal(nx.minimum_cut_value(G, u, v),
                             cut_value)

    def test_davis_southern_women_graph(self):
        G = nx.davis_southern_women_graph()
        nx.set_edge_attributes(G, 1, 'capacity')
        for flow_func in flow_funcs:
            T = nx.gomory_hu_tree(G, flow_func=flow_func)
            assert_true(nx.is_tree(T))
            for u, v in combinations(G, 2):
                cut_value, edge = self.minimum_edge_weight(T, u, v)
                assert_equal(nx.minimum_cut_value(G, u, v),
                             cut_value)

    def test_florentine_families_graph(self):
        G = nx.florentine_families_graph()
        nx.set_edge_attributes(G, 1, 'capacity')
        for flow_func in flow_funcs:
            T = nx.gomory_hu_tree(G, flow_func=flow_func)
            assert_true(nx.is_tree(T))
            for u, v in combinations(G, 2):
                cut_value, edge = self.minimum_edge_weight(T, u, v)
                assert_equal(nx.minimum_cut_value(G, u, v),
                             cut_value)

    def test_karate_club_graph_cutset(self):
        G = nx.karate_club_graph()
        nx.set_edge_attributes(G, 1, 'capacity')
        T = nx.gomory_hu_tree(G)
        assert_true(nx.is_tree(T))
        u, v = 0, 33
        cut_value, edge = self.minimum_edge_weight(T, u, v)
        cutset = self.compute_cutset(G, T, edge)
        assert_equal(cut_value, len(cutset))

    def test_wikipedia_example(self):
        # Example from https://en.wikipedia.org/wiki/Gomory%E2%80%93Hu_tree
        G = nx.Graph()
        G.add_weighted_edges_from((
            (0, 1, 1), (0, 2, 7), (1, 2, 1),
            (1, 3, 3), (1, 4, 2), (2, 4, 4),
            (3, 4, 1), (3, 5, 6), (4, 5, 2),
        ))
        for flow_func in flow_funcs:
            T = nx.gomory_hu_tree(G, capacity='weight', flow_func=flow_func)
            assert_true(nx.is_tree(T))
            for u, v in combinations(G, 2):
                cut_value, edge = self.minimum_edge_weight(T, u, v)
                assert_equal(nx.minimum_cut_value(G, u, v, capacity='weight'),
                             cut_value)

    @raises(nx.NetworkXNotImplemented)
    def test_directed_raises(self):
        G = nx.DiGraph()
        T = nx.gomory_hu_tree(G)

    @raises(nx.NetworkXError)
    def test_empty_raises(self):
        G = nx.empty_graph()
        T = nx.gomory_hu_tree(G)
