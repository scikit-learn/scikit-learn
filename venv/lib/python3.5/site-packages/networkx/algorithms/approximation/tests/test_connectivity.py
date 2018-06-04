import itertools
from nose.tools import assert_true, assert_equal, assert_raises

import networkx as nx
from networkx.algorithms import approximation as approx


def test_global_node_connectivity():
    # Figure 1 chapter on Connectivity
    G = nx.Graph()
    G.add_edges_from([(1, 2), (1, 3), (1, 4), (1, 5), (2, 3), (2, 6), (3, 4),
                      (3, 6), (4, 6), (4, 7), (5, 7), (6, 8), (6, 9), (7, 8),
                      (7, 10), (8, 11), (9, 10), (9, 11), (10, 11)])
    assert_equal(2, approx.local_node_connectivity(G, 1, 11))
    assert_equal(2, approx.node_connectivity(G))
    assert_equal(2, approx.node_connectivity(G, 1, 11))


def test_white_harary1():
    # Figure 1b white and harary (2001)
    # A graph with high adhesion (edge connectivity) and low cohesion
    # (node connectivity)
    G = nx.disjoint_union(nx.complete_graph(4), nx.complete_graph(4))
    G.remove_node(7)
    for i in range(4, 7):
        G.add_edge(0, i)
    G = nx.disjoint_union(G, nx.complete_graph(4))
    G.remove_node(G.order() - 1)
    for i in range(7, 10):
        G.add_edge(0, i)
    assert_equal(1, approx.node_connectivity(G))


def test_complete_graphs():
    for n in range(5, 25, 5):
        G = nx.complete_graph(n)
        assert_equal(n - 1, approx.node_connectivity(G))
        assert_equal(n - 1, approx.node_connectivity(G, 0, 3))


def test_empty_graphs():
    for k in range(5, 25, 5):
        G = nx.empty_graph(k)
        assert_equal(0, approx.node_connectivity(G))
        assert_equal(0, approx.node_connectivity(G, 0, 3))


def test_petersen():
    G = nx.petersen_graph()
    assert_equal(3, approx.node_connectivity(G))
    assert_equal(3, approx.node_connectivity(G, 0, 5))

# Approximation fails with tutte graph
# def test_tutte():
#    G = nx.tutte_graph()
#    assert_equal(3, approx.node_connectivity(G))


def test_dodecahedral():
    G = nx.dodecahedral_graph()
    assert_equal(3, approx.node_connectivity(G))
    assert_equal(3, approx.node_connectivity(G, 0, 5))


def test_octahedral():
    G = nx.octahedral_graph()
    assert_equal(4, approx.node_connectivity(G))
    assert_equal(4, approx.node_connectivity(G, 0, 5))

# Approximation can fail with icosahedral graph depending
# on iteration order.
# def test_icosahedral():
#    G=nx.icosahedral_graph()
#    assert_equal(5, approx.node_connectivity(G))
#    assert_equal(5, approx.node_connectivity(G, 0, 5))


def test_only_source():
    G = nx.complete_graph(5)
    assert_raises(nx.NetworkXError, approx.node_connectivity, G, s=0)


def test_only_target():
    G = nx.complete_graph(5)
    assert_raises(nx.NetworkXError, approx.node_connectivity, G, t=0)


def test_missing_source():
    G = nx.path_graph(4)
    assert_raises(nx.NetworkXError, approx.node_connectivity, G, 10, 1)


def test_missing_target():
    G = nx.path_graph(4)
    assert_raises(nx.NetworkXError, approx.node_connectivity, G, 1, 10)


def test_source_equals_target():
    G = nx.complete_graph(5)
    assert_raises(nx.NetworkXError, approx.local_node_connectivity, G, 0, 0)


def test_directed_node_connectivity():
    G = nx.cycle_graph(10, create_using=nx.DiGraph())  # only one direction
    D = nx.cycle_graph(10).to_directed()  # 2 reciprocal edges
    assert_equal(1, approx.node_connectivity(G))
    assert_equal(1, approx.node_connectivity(G, 1, 4))
    assert_equal(2,  approx.node_connectivity(D))
    assert_equal(2,  approx.node_connectivity(D, 1, 4))


class TestAllPairsNodeConnectivityApprox:

    def setUp(self):
        self.path = nx.path_graph(7)
        self.directed_path = nx.path_graph(7, create_using=nx.DiGraph())
        self.cycle = nx.cycle_graph(7)
        self.directed_cycle = nx.cycle_graph(7, create_using=nx.DiGraph())
        self.gnp = nx.gnp_random_graph(30, 0.1)
        self.directed_gnp = nx.gnp_random_graph(30, 0.1, directed=True)
        self.K20 = nx.complete_graph(20)
        self.K10 = nx.complete_graph(10)
        self.K5 = nx.complete_graph(5)
        self.G_list = [self.path, self.directed_path, self.cycle,
                       self.directed_cycle, self.gnp, self.directed_gnp, self.K10,
                       self.K5, self.K20]

    def test_cycles(self):
        K_undir = approx.all_pairs_node_connectivity(self.cycle)
        for source in K_undir:
            for target, k in K_undir[source].items():
                assert_true(k == 2)
        K_dir = approx.all_pairs_node_connectivity(self.directed_cycle)
        for source in K_dir:
            for target, k in K_dir[source].items():
                assert_true(k == 1)

    def test_complete(self):
        for G in [self.K10, self.K5, self.K20]:
            K = approx.all_pairs_node_connectivity(G)
            for source in K:
                for target, k in K[source].items():
                    assert_true(k == len(G) - 1)

    def test_paths(self):
        K_undir = approx.all_pairs_node_connectivity(self.path)
        for source in K_undir:
            for target, k in K_undir[source].items():
                assert_true(k == 1)
        K_dir = approx.all_pairs_node_connectivity(self.directed_path)
        for source in K_dir:
            for target, k in K_dir[source].items():
                if source < target:
                    assert_true(k == 1)
                else:
                    assert_true(k == 0)

    def test_cutoff(self):
        for G in [self.K10, self.K5, self.K20]:
            for mp in [2, 3, 4]:
                paths = approx.all_pairs_node_connectivity(G, cutoff=mp)
                for source in paths:
                    for target, K in paths[source].items():
                        assert_true(K == mp)

    def test_all_pairs_connectivity_nbunch(self):
        G = nx.complete_graph(5)
        nbunch = [0, 2, 3]
        C = approx.all_pairs_node_connectivity(G, nbunch=nbunch)
        assert_equal(len(C), len(nbunch))
