import json
from nose.tools import assert_equal, assert_true, raises
import networkx as nx
from networkx.readwrite.json_graph import *


class TestCytoscape:

    def test_graph(self):
        G = nx.path_graph(4)
        H = cytoscape_graph(cytoscape_data(G))
        nx.is_isomorphic(G, H)

    def test_graph_attributes(self):
        G = nx.path_graph(4)
        G.add_node(1, color='red')
        G.add_edge(1, 2, width=7)
        G.graph['foo'] = 'bar'
        G.graph[1] = 'one'
        G.add_node(3, name="node", id="123")

        H = cytoscape_graph(cytoscape_data(G))
        assert_equal(H.graph['foo'], 'bar')
        assert_equal(H.nodes[1]['color'], 'red')
        assert_equal(H[1][2]['width'], 7)
        assert_equal(H.nodes[3]['name'], 'node')
        assert_equal(H.nodes[3]['id'], '123')

        d = json.dumps(cytoscape_data(G))
        H = cytoscape_graph(json.loads(d))
        assert_equal(H.graph['foo'], 'bar')
        assert_equal(H.graph[1], 'one')
        assert_equal(H.nodes[1]['color'], 'red')
        assert_equal(H[1][2]['width'], 7)
        assert_equal(H.nodes[3]['name'], 'node')
        assert_equal(H.nodes[3]['id'], '123')

    def test_digraph(self):
        G = nx.DiGraph()
        nx.add_path(G, [1, 2, 3])
        H = cytoscape_graph(cytoscape_data(G))
        assert_true(H.is_directed())
        nx.is_isomorphic(G, H)

    def test_multidigraph(self):
        G = nx.MultiDiGraph()
        nx.add_path(G, [1, 2, 3])
        H = cytoscape_graph(cytoscape_data(G))
        assert_true(H.is_directed())
        assert_true(H.is_multigraph())

    def test_multigraph(self):
        G = nx.MultiGraph()
        G.add_edge(1, 2, key='first')
        G.add_edge(1, 2, key='second', color='blue')
        H = cytoscape_graph(cytoscape_data(G))
        assert_true(nx.is_isomorphic(G, H))
        assert_equal(H[1][2]['second']['color'], 'blue')

    @raises(nx.NetworkXError)
    def test_exception(self):
        G = nx.MultiDiGraph()
        attrs = dict(name='node', ident='node')
        cytoscape_data(G, attrs)
