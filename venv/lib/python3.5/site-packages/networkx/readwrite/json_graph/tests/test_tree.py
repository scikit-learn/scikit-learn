import json
from nose.tools import assert_equal, assert_true, raises
import networkx as nx
from networkx.readwrite.json_graph import *


class TestTree:

    def test_graph(self):
        G = nx.DiGraph()
        G.add_nodes_from([1, 2, 3], color='red')
        G.add_edge(1, 2, foo=7)
        G.add_edge(1, 3, foo=10)
        G.add_edge(3, 4, foo=10)
        H = tree_graph(tree_data(G, 1))
        nx.is_isomorphic(G, H)

    def test_graph_attributes(self):
        G = nx.DiGraph()
        G.add_nodes_from([1, 2, 3], color='red')
        G.add_edge(1, 2, foo=7)
        G.add_edge(1, 3, foo=10)
        G.add_edge(3, 4, foo=10)
        H = tree_graph(tree_data(G, 1))
        assert_equal(H.nodes[1]['color'], 'red')

        d = json.dumps(tree_data(G, 1))
        H = tree_graph(json.loads(d))
        assert_equal(H.nodes[1]['color'], 'red')

    @raises(nx.NetworkXError)
    def test_exception(self):
        G = nx.MultiDiGraph()
        G.add_node(0)
        attrs = dict(id='node', children='node')
        tree_data(G, 0, attrs)
