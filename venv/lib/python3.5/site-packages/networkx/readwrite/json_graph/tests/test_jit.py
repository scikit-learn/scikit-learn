import json
from nose.tools import assert_true, assert_false, assert_raises
import networkx as nx
from networkx.readwrite.json_graph import jit_data, jit_graph


class TestJIT(object):
    def test_jit(self):
        G = nx.Graph()
        G.add_node('Node1', node_data='foobar')
        G.add_node('Node3', node_data='bar')
        G.add_node('Node4')
        G.add_edge('Node1', 'Node2', weight=9, something='isSomething')
        G.add_edge('Node2', 'Node3', weight=4, something='isNotSomething')
        G.add_edge('Node1', 'Node2')
        d = jit_data(G)
        K = jit_graph(json.loads(d))
        assert_true(nx.is_isomorphic(G, K))

    def test_jit_2(self):
        G = nx.Graph()
        G.add_node(1, node_data=3)
        G.add_node(3, node_data=0)
        G.add_edge(1, 2, weight=9, something=0)
        G.add_edge(2, 3, weight=4, something=3)
        G.add_edge(1, 2)
        d = jit_data(G)
        K = jit_graph(json.loads(d))
        assert_true(nx.is_isomorphic(G, K))

    def test_jit_directed(self):
        G = nx.DiGraph()
        G.add_node(1, node_data=3)
        G.add_node(3, node_data=0)
        G.add_edge(1, 2, weight=9, something=0)
        G.add_edge(2, 3, weight=4, something=3)
        G.add_edge(1, 2)
        d = jit_data(G)
        K = jit_graph(json.loads(d), create_using=nx.DiGraph())
        assert_true(nx.is_isomorphic(G, K))

    def test_jit_multi_directed(self):
        G = nx.MultiDiGraph()
        G.add_node(1, node_data=3)
        G.add_node(3, node_data=0)
        G.add_edge(1, 2, weight=9, something=0)
        G.add_edge(2, 3, weight=4, something=3)
        G.add_edge(1, 2)
        assert_raises(nx.NetworkXNotImplemented, jit_data, G)

        H = nx.DiGraph(G)
        d = jit_data(H)
        K = jit_graph(json.loads(d), create_using=nx.MultiDiGraph())
        assert_true(nx.is_isomorphic(H, K))
        K.add_edge(1, 2)
        assert_false(nx.is_isomorphic(H, K))
        assert_true(nx.is_isomorphic(G, K))
