""" Tests for subgraphs attributes
"""
from copy import deepcopy
from nose.tools import assert_equal
import networkx as nx

# deprecated in 2.1 for removal in 2.2


class TestSubgraphAttributesDicts:

    def setUp(self):
        self.undirected = [
            nx.connected_component_subgraphs,
            nx.biconnected_component_subgraphs,
        ]
        self.directed = [
            nx.weakly_connected_component_subgraphs,
            nx.strongly_connected_component_subgraphs,
            nx.attracting_component_subgraphs,
        ]
        self.subgraph_funcs = self.undirected + self.directed

        self.D = nx.DiGraph()
        self.D.add_edge(1, 2, eattr='red')
        self.D.add_edge(2, 1, eattr='red')
        self.D.nodes[1]['nattr'] = 'blue'
        self.D.graph['gattr'] = 'green'

        self.G = nx.Graph()
        self.G.add_edge(1, 2, eattr='red')
        self.G.nodes[1]['nattr'] = 'blue'
        self.G.graph['gattr'] = 'green'

    def test_subgraphs_default_copy_behavior(self):
        # Test the default behavior of subgraph functions
        # For the moment (1.10) the default is to copy
        for subgraph_func in self.subgraph_funcs:
            G = deepcopy(self.G if subgraph_func in self.undirected else self.D)
            SG = list(subgraph_func(G))[0]
            assert_equal(SG[1][2]['eattr'], 'red')
            assert_equal(SG.nodes[1]['nattr'], 'blue')
            assert_equal(SG.graph['gattr'], 'green')
            SG[1][2]['eattr'] = 'foo'
            assert_equal(G[1][2]['eattr'], 'red')
            assert_equal(SG[1][2]['eattr'], 'foo')
            SG.nodes[1]['nattr'] = 'bar'
            assert_equal(G.nodes[1]['nattr'], 'blue')
            assert_equal(SG.nodes[1]['nattr'], 'bar')
            SG.graph['gattr'] = 'baz'
            assert_equal(G.graph['gattr'], 'green')
            assert_equal(SG.graph['gattr'], 'baz')

    def test_subgraphs_copy(self):
        for subgraph_func in self.subgraph_funcs:
            test_graph = self.G if subgraph_func in self.undirected else self.D
            G = deepcopy(test_graph)
            SG = list(subgraph_func(G, copy=True))[0]
            assert_equal(SG[1][2]['eattr'], 'red')
            assert_equal(SG.nodes[1]['nattr'], 'blue')
            assert_equal(SG.graph['gattr'], 'green')
            SG[1][2]['eattr'] = 'foo'
            assert_equal(G[1][2]['eattr'], 'red')
            assert_equal(SG[1][2]['eattr'], 'foo')
            SG.nodes[1]['nattr'] = 'bar'
            assert_equal(G.nodes[1]['nattr'], 'blue')
            assert_equal(SG.nodes[1]['nattr'], 'bar')
            SG.graph['gattr'] = 'baz'
            assert_equal(G.graph['gattr'], 'green')
            assert_equal(SG.graph['gattr'], 'baz')

    def test_subgraphs_no_copy(self):
        for subgraph_func in self.subgraph_funcs:
            G = deepcopy(self.G if subgraph_func in self.undirected else self.D)
            SG = list(subgraph_func(G, copy=False))[0]
            assert_equal(SG[1][2]['eattr'], 'red')
            assert_equal(SG.nodes[1]['nattr'], 'blue')
            assert_equal(SG.graph['gattr'], 'green')
            SG[1][2]['eattr'] = 'foo'
            assert_equal(G[1][2]['eattr'], 'foo')
            assert_equal(SG[1][2]['eattr'], 'foo')
            SG.nodes[1]['nattr'] = 'bar'
            assert_equal(G.nodes[1]['nattr'], 'bar')
            assert_equal(SG.nodes[1]['nattr'], 'bar')
            SG.graph['gattr'] = 'baz'
            assert_equal(G.graph['gattr'], 'baz')
            assert_equal(SG.graph['gattr'], 'baz')
