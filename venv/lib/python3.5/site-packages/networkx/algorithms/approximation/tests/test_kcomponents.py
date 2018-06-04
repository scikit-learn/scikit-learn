# Test for approximation to k-components algorithm
from nose.tools import assert_equal, assert_true, assert_false, assert_in
from nose.tools import assert_raises, raises, assert_greater_equal
import networkx as nx
from networkx.algorithms.approximation import k_components
from networkx.algorithms.approximation.kcomponents import _AntiGraph, _same


def build_k_number_dict(k_components):
    k_num = {}
    for k, comps in sorted(k_components.items()):
        for comp in comps:
            for node in comp:
                k_num[node] = k
    return k_num

##
# Some nice synthetic graphs
##


def graph_example_1():
    G = nx.convert_node_labels_to_integers(nx.grid_graph([5, 5]),
                                           label_attribute='labels')
    rlabels = nx.get_node_attributes(G, 'labels')
    labels = {v: k for k, v in rlabels.items()}

    for nodes in [(labels[(0, 0)], labels[(1, 0)]),
                  (labels[(0, 4)], labels[(1, 4)]),
                  (labels[(3, 0)], labels[(4, 0)]),
                  (labels[(3, 4)], labels[(4, 4)])]:
        new_node = G.order() + 1
        # Petersen graph is triconnected
        P = nx.petersen_graph()
        G = nx.disjoint_union(G, P)
        # Add two edges between the grid and P
        G.add_edge(new_node + 1, nodes[0])
        G.add_edge(new_node, nodes[1])
        # K5 is 4-connected
        K = nx.complete_graph(5)
        G = nx.disjoint_union(G, K)
        # Add three edges between P and K5
        G.add_edge(new_node + 2, new_node + 11)
        G.add_edge(new_node + 3, new_node + 12)
        G.add_edge(new_node + 4, new_node + 13)
        # Add another K5 sharing a node
        G = nx.disjoint_union(G, K)
        nbrs = G[new_node + 10]
        G.remove_node(new_node + 10)
        for nbr in nbrs:
            G.add_edge(new_node + 17, nbr)
        G.add_edge(new_node + 16, new_node + 5)
    return G


def torrents_and_ferraro_graph():
    G = nx.convert_node_labels_to_integers(nx.grid_graph([5, 5]),
                                           label_attribute='labels')
    rlabels = nx.get_node_attributes(G, 'labels')
    labels = {v: k for k, v in rlabels.items()}

    for nodes in [(labels[(0, 4)], labels[(1, 4)]),
                  (labels[(3, 4)], labels[(4, 4)])]:
        new_node = G.order() + 1
        # Petersen graph is triconnected
        P = nx.petersen_graph()
        G = nx.disjoint_union(G, P)
        # Add two edges between the grid and P
        G.add_edge(new_node + 1, nodes[0])
        G.add_edge(new_node, nodes[1])
        # K5 is 4-connected
        K = nx.complete_graph(5)
        G = nx.disjoint_union(G, K)
        # Add three edges between P and K5
        G.add_edge(new_node + 2, new_node + 11)
        G.add_edge(new_node + 3, new_node + 12)
        G.add_edge(new_node + 4, new_node + 13)
        # Add another K5 sharing a node
        G = nx.disjoint_union(G, K)
        nbrs = G[new_node + 10]
        G.remove_node(new_node + 10)
        for nbr in nbrs:
            G.add_edge(new_node + 17, nbr)
        # Commenting this makes the graph not biconnected !!
        # This stupid mistake make one reviewer very angry :P
        G.add_edge(new_node + 16, new_node + 8)

    for nodes in [(labels[(0, 0)], labels[(1, 0)]),
                  (labels[(3, 0)], labels[(4, 0)])]:
        new_node = G.order() + 1
        # Petersen graph is triconnected
        P = nx.petersen_graph()
        G = nx.disjoint_union(G, P)
        # Add two edges between the grid and P
        G.add_edge(new_node + 1, nodes[0])
        G.add_edge(new_node, nodes[1])
        # K5 is 4-connected
        K = nx.complete_graph(5)
        G = nx.disjoint_union(G, K)
        # Add three edges between P and K5
        G.add_edge(new_node + 2, new_node + 11)
        G.add_edge(new_node + 3, new_node + 12)
        G.add_edge(new_node + 4, new_node + 13)
        # Add another K5 sharing two nodes
        G = nx.disjoint_union(G, K)
        nbrs = G[new_node + 10]
        G.remove_node(new_node + 10)
        for nbr in nbrs:
            G.add_edge(new_node + 17, nbr)
        nbrs2 = G[new_node + 9]
        G.remove_node(new_node + 9)
        for nbr in nbrs2:
            G.add_edge(new_node + 18, nbr)
    return G

# Helper function


def _check_connectivity(G):
    result = k_components(G)
    for k, components in result.items():
        if k < 3:
            continue
        for component in components:
            C = G.subgraph(component)
            K = nx.node_connectivity(C)
            assert_greater_equal(K, k)


def test_torrents_and_ferraro_graph():
    G = torrents_and_ferraro_graph()
    _check_connectivity(G)


def test_example_1():
    G = graph_example_1()
    _check_connectivity(G)


def test_karate_0():
    G = nx.karate_club_graph()
    _check_connectivity(G)


def test_karate_1():
    karate_k_num = {0: 4, 1: 4, 2: 4, 3: 4, 4: 3, 5: 3, 6: 3, 7: 4, 8: 4, 9: 2,
                    10: 3, 11: 1, 12: 2, 13: 4, 14: 2, 15: 2, 16: 2, 17: 2, 18: 2,
                    19: 3, 20: 2, 21: 2, 22: 2, 23: 3, 24: 3, 25: 3, 26: 2, 27: 3,
                    28: 3, 29: 3, 30: 4, 31: 3, 32: 4, 33: 4}
    approx_karate_k_num = karate_k_num.copy()
    approx_karate_k_num[24] = 2
    approx_karate_k_num[25] = 2
    G = nx.karate_club_graph()
    k_comps = k_components(G)
    k_num = build_k_number_dict(k_comps)
    assert_in(k_num, (karate_k_num, approx_karate_k_num))


def test_example_1_detail_3_and_4():
    G = graph_example_1()
    result = k_components(G)
    # In this example graph there are 8 3-components, 4 with 15 nodes
    # and 4 with 5 nodes.
    assert_equal(len(result[3]), 8)
    assert_equal(len([c for c in result[3] if len(c) == 15]), 4)
    assert_equal(len([c for c in result[3] if len(c) == 5]), 4)
    # There are also 8 4-components all with 5 nodes.
    assert_equal(len(result[4]), 8)
    assert_true(all(len(c) == 5 for c in result[4]))
    # Finally check that the k-components detected have actually node
    # connectivity >= k.
    for k, components in result.items():
        if k < 3:
            continue
        for component in components:
            K = nx.node_connectivity(G.subgraph(component))
            assert_greater_equal(K, k)


@raises(nx.NetworkXNotImplemented)
def test_directed():
    G = nx.gnp_random_graph(10, 0.4, directed=True)
    kc = k_components(G)


def test_same():
    equal = {'A': 2, 'B': 2, 'C': 2}
    slightly_different = {'A': 2, 'B': 1, 'C': 2}
    different = {'A': 2, 'B': 8, 'C': 18}
    assert_true(_same(equal))
    assert_false(_same(slightly_different))
    assert_true(_same(slightly_different, tol=1))
    assert_false(_same(different))
    assert_false(_same(different, tol=4))


class TestAntiGraph:
    def setUp(self):
        self.Gnp = nx.gnp_random_graph(20, 0.8)
        self.Anp = _AntiGraph(nx.complement(self.Gnp))
        self.Gd = nx.davis_southern_women_graph()
        self.Ad = _AntiGraph(nx.complement(self.Gd))
        self.Gk = nx.karate_club_graph()
        self.Ak = _AntiGraph(nx.complement(self.Gk))
        self.GA = [(self.Gnp, self.Anp),
                   (self.Gd, self.Ad),
                   (self.Gk, self.Ak)]

    def test_size(self):
        for G, A in self.GA:
            n = G.order()
            s = len(list(G.edges())) + len(list(A.edges()))
            assert_true(s == (n * (n - 1)) / 2)

    def test_degree(self):
        for G, A in self.GA:
            assert_equal(G.degree(), A.degree())

    def test_core_number(self):
        for G, A in self.GA:
            assert_equal(nx.core_number(G), nx.core_number(A))

    def test_connected_components(self):
        for G, A in self.GA:
            gc = [set(c) for c in nx.connected_components(G)]
            ac = [set(c) for c in nx.connected_components(A)]
            for comp in ac:
                assert_true(comp in gc)

    def test_adj(self):
        for G, A in self.GA:
            for n, nbrs in G.adj.items():
                a_adj = sorted((n, sorted(ad)) for n, ad in A.adj.items())
                g_adj = sorted((n, sorted(ad)) for n, ad in G.adj.items())
                assert_equal(a_adj, g_adj)

    def test_adjacency(self):
        for G, A in self.GA:
            a_adj = list(A.adjacency())
            for n, nbrs in G.adjacency():
                assert_true((n, set(nbrs)) in a_adj)

    def test_neighbors(self):
        for G, A in self.GA:
            node = list(G.nodes())[0]
            assert_equal(set(G.neighbors(node)), set(A.neighbors(node)))

    def test_node_not_in_graph(self):
        for G, A in self.GA:
            node = 'non_existent_node'
            assert_raises(nx.NetworkXError, A.neighbors, node)
            assert_raises(nx.NetworkXError, G.neighbors, node)

    def test_degree(self):
        for G, A in self.GA:
            node = list(G.nodes())[0]
            nodes = list(G.nodes())[1:4]
            assert_equal(G.degree(node), A.degree(node))
            assert_equal(sum(d for n, d in G.degree()), sum(d for n, d in A.degree()))
            # AntiGraph is a ThinGraph, so all the weights are 1
            assert_equal(sum(d for n, d in A.degree()),
                         sum(d for n, d in A.degree(weight='weight')))
            assert_equal(sum(d for n, d in G.degree(nodes)),
                         sum(d for n, d in A.degree(nodes)))
