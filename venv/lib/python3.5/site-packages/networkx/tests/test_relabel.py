#!/usr/bin/env python
from nose.tools import *
from networkx import *
from networkx.convert import *
from networkx.algorithms.operators import *
from networkx.generators.classic import barbell_graph, cycle_graph
from networkx.testing import *


class TestRelabel():
    def test_convert_node_labels_to_integers(self):
        # test that empty graph converts fine for all options
        G = empty_graph()
        H = convert_node_labels_to_integers(G, 100)
        assert_equal(list(H.nodes()), [])
        assert_equal(list(H.edges()), [])

        for opt in ["default", "sorted", "increasing degree", "decreasing degree"]:
            G = empty_graph()
            H = convert_node_labels_to_integers(G, 100, ordering=opt)
            assert_equal(list(H.nodes()), [])
            assert_equal(list(H.edges()), [])

        G = empty_graph()
        G.add_edges_from([('A', 'B'), ('A', 'C'), ('B', 'C'), ('C', 'D')])
        H = convert_node_labels_to_integers(G)
        degH = (d for n, d in H.degree())
        degG = (d for n, d in G.degree())
        assert_equal(sorted(degH), sorted(degG))

        H = convert_node_labels_to_integers(G, 1000)
        degH = (d for n, d in H.degree())
        degG = (d for n, d in G.degree())
        assert_equal(sorted(degH), sorted(degG))
        assert_nodes_equal(H.nodes(), [1000, 1001, 1002, 1003])

        H = convert_node_labels_to_integers(G, ordering="increasing degree")
        degH = (d for n, d in H.degree())
        degG = (d for n, d in G.degree())
        assert_equal(sorted(degH), sorted(degG))
        assert_equal(degree(H, 0), 1)
        assert_equal(degree(H, 1), 2)
        assert_equal(degree(H, 2), 2)
        assert_equal(degree(H, 3), 3)

        H = convert_node_labels_to_integers(G, ordering="decreasing degree")
        degH = (d for n, d in H.degree())
        degG = (d for n, d in G.degree())
        assert_equal(sorted(degH), sorted(degG))
        assert_equal(degree(H, 0), 3)
        assert_equal(degree(H, 1), 2)
        assert_equal(degree(H, 2), 2)
        assert_equal(degree(H, 3), 1)

        H = convert_node_labels_to_integers(G, ordering="increasing degree",
                                            label_attribute='label')
        degH = (d for n, d in H.degree())
        degG = (d for n, d in G.degree())
        assert_equal(sorted(degH), sorted(degG))
        assert_equal(degree(H, 0), 1)
        assert_equal(degree(H, 1), 2)
        assert_equal(degree(H, 2), 2)
        assert_equal(degree(H, 3), 3)

        # check mapping
        assert_equal(H.nodes[3]['label'], 'C')
        assert_equal(H.nodes[0]['label'], 'D')
        assert_true(H.nodes[1]['label'] == 'A' or H.nodes[2]['label'] == 'A')
        assert_true(H.nodes[1]['label'] == 'B' or H.nodes[2]['label'] == 'B')

    def test_convert_to_integers2(self):
        G = empty_graph()
        G.add_edges_from([('C', 'D'), ('A', 'B'), ('A', 'C'), ('B', 'C')])
        H = convert_node_labels_to_integers(G, ordering="sorted")
        degH = (d for n, d in H.degree())
        degG = (d for n, d in G.degree())
        assert_equal(sorted(degH), sorted(degG))

        H = convert_node_labels_to_integers(G, ordering="sorted",
                                            label_attribute='label')
        assert_equal(H.nodes[0]['label'], 'A')
        assert_equal(H.nodes[1]['label'], 'B')
        assert_equal(H.nodes[2]['label'], 'C')
        assert_equal(H.nodes[3]['label'], 'D')

    @raises(nx.NetworkXError)
    def test_convert_to_integers_raise(self):
        G = nx.Graph()
        H = convert_node_labels_to_integers(G, ordering="increasing age")

    def test_relabel_nodes_copy(self):
        G = empty_graph()
        G.add_edges_from([('A', 'B'), ('A', 'C'), ('B', 'C'), ('C', 'D')])
        mapping = {'A': 'aardvark', 'B': 'bear', 'C': 'cat', 'D': 'dog'}
        H = relabel_nodes(G, mapping)
        assert_nodes_equal(H.nodes(), ['aardvark', 'bear', 'cat', 'dog'])

    def test_relabel_nodes_function(self):
        G = empty_graph()
        G.add_edges_from([('A', 'B'), ('A', 'C'), ('B', 'C'), ('C', 'D')])
        # function mapping no longer encouraged but works

        def mapping(n):
            return ord(n)
        H = relabel_nodes(G, mapping)
        assert_nodes_equal(H.nodes(), [65, 66, 67, 68])

    def test_relabel_nodes_graph(self):
        G = Graph([('A', 'B'), ('A', 'C'), ('B', 'C'), ('C', 'D')])
        mapping = {'A': 'aardvark', 'B': 'bear', 'C': 'cat', 'D': 'dog'}
        H = relabel_nodes(G, mapping)
        assert_nodes_equal(H.nodes(), ['aardvark', 'bear', 'cat', 'dog'])

    def test_relabel_nodes_orderedgraph(self):
        G = OrderedGraph()
        G.add_nodes_from([1, 2, 3])
        G.add_edges_from([(1, 3), (2, 3)])
        mapping = {1: 'a', 2: 'b', 3: 'c'}
        H = relabel_nodes(G, mapping)
        assert list(H.nodes) == ['a', 'b', 'c']

    def test_relabel_nodes_digraph(self):
        G = DiGraph([('A', 'B'), ('A', 'C'), ('B', 'C'), ('C', 'D')])
        mapping = {'A': 'aardvark', 'B': 'bear', 'C': 'cat', 'D': 'dog'}
        H = relabel_nodes(G, mapping, copy=False)
        assert_nodes_equal(H.nodes(), ['aardvark', 'bear', 'cat', 'dog'])

    def test_relabel_nodes_multigraph(self):
        G = MultiGraph([('a', 'b'), ('a', 'b')])
        mapping = {'a': 'aardvark', 'b': 'bear'}
        G = relabel_nodes(G, mapping, copy=False)
        assert_nodes_equal(G.nodes(), ['aardvark', 'bear'])
        assert_edges_equal(G.edges(), [('aardvark', 'bear'), ('aardvark', 'bear')])

    def test_relabel_nodes_multidigraph(self):
        G = MultiDiGraph([('a', 'b'), ('a', 'b')])
        mapping = {'a': 'aardvark', 'b': 'bear'}
        G = relabel_nodes(G, mapping, copy=False)
        assert_nodes_equal(G.nodes(), ['aardvark', 'bear'])
        assert_edges_equal(G.edges(), [('aardvark', 'bear'), ('aardvark', 'bear')])

    def test_relabel_isolated_nodes_to_same(self):
        G = Graph()
        G.add_nodes_from(range(4))
        mapping = {1: 1}
        H = relabel_nodes(G, mapping, copy=False)
        assert_nodes_equal(H.nodes(), list(range(4)))

    @raises(KeyError)
    def test_relabel_nodes_missing(self):
        G = Graph([('A', 'B'), ('A', 'C'), ('B', 'C'), ('C', 'D')])
        mapping = {0: 'aardvark'}
        G = relabel_nodes(G, mapping, copy=False)

    def test_relabel_copy_name(self):
        G = Graph()
        H = relabel_nodes(G, {}, copy=True)
        assert_equal(H.graph, G.graph)
        H = relabel_nodes(G, {}, copy=False)
        assert_equal(H.graph, G.graph)
        G.name = "first"
        H = relabel_nodes(G, {}, copy=True)
        assert_equal(H.graph, G.graph)
        H = relabel_nodes(G, {}, copy=False)
        assert_equal(H.graph, G.graph)

    def test_relabel_toposort(self):
        K5 = nx.complete_graph(4)
        G = nx.complete_graph(4)
        G = nx.relabel_nodes(G, dict([(i, i + 1) for i in range(4)]), copy=False)
        nx.is_isomorphic(K5, G)
        G = nx.complete_graph(4)
        G = nx.relabel_nodes(G, dict([(i, i - 1) for i in range(4)]), copy=False)
        nx.is_isomorphic(K5, G)

    def test_relabel_selfloop(self):
        G = nx.DiGraph([(1, 1), (1, 2), (2, 3)])
        G = nx.relabel_nodes(G, {1: 'One', 2: 'Two', 3: 'Three'}, copy=False)
        assert_nodes_equal(G.nodes(), ['One', 'Three', 'Two'])
        G = nx.MultiDiGraph([(1, 1), (1, 2), (2, 3)])
        G = nx.relabel_nodes(G, {1: 'One', 2: 'Two', 3: 'Three'}, copy=False)
        assert_nodes_equal(G.nodes(), ['One', 'Three', 'Two'])
        G = nx.MultiDiGraph([(1, 1)])
        G = nx.relabel_nodes(G, {1: 0}, copy=False)
        assert_nodes_equal(G.nodes(), [0])
