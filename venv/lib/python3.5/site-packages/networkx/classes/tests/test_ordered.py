from nose.tools import assert_equals
import networkx as nx


class SmokeTestOrdered(object):
    # Just test instantiation.
    def test_graph():
        G = nx.OrderedGraph()

    def test_digraph():
        G = nx.OrderedDiGraph()

    def test_multigraph():
        G = nx.OrderedMultiGraph()

    def test_multidigraph():
        G = nx.OrderedMultiDiGraph()


class TestOrderedFeatures(object):
    def setUp(self):
        self.G = nx.OrderedDiGraph()
        self.G.add_nodes_from([1, 2, 3])
        self.G.add_edges_from([(2, 3), (1, 3)])

    def test_subgraph_order(self):
        G = self.G
        G_sub = G.subgraph([1, 2, 3])
        assert_equals(list(G.nodes), list(G_sub.nodes))
        assert_equals(list(G.edges), list(G_sub.edges))
        assert_equals(list(G.pred[3]), list(G_sub.pred[3]))
        assert_equals([2, 1], list(G_sub.pred[3]))
        assert_equals([], list(G_sub.succ[3]))

        G_sub = nx.induced_subgraph(G, [1, 2, 3])
        assert_equals(list(G.nodes), list(G_sub.nodes))
        assert_equals(list(G.edges), list(G_sub.edges))
        assert_equals(list(G.pred[3]), list(G_sub.pred[3]))
        assert_equals([2, 1], list(G_sub.pred[3]))
        assert_equals([], list(G_sub.succ[3]))
