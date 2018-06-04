#!/usr/bin/env python
from nose.tools import *
import networkx as nx
from networkx import convert_node_labels_to_integers as cnlti
from networkx import NetworkXNotImplemented


class TestConnected:

    def setUp(self):
        G1 = cnlti(nx.grid_2d_graph(2, 2), first_label=0, ordering="sorted")
        G2 = cnlti(nx.lollipop_graph(3, 3), first_label=4, ordering="sorted")
        G3 = cnlti(nx.house_graph(), first_label=10, ordering="sorted")
        self.G = nx.union(G1, G2)
        self.G = nx.union(self.G, G3)
        self.DG = nx.DiGraph([(1, 2), (1, 3), (2, 3)])
        self.grid = cnlti(nx.grid_2d_graph(4, 4), first_label=1)

        self.gc = []
        G = nx.DiGraph()
        G.add_edges_from([(1, 2), (2, 3), (2, 8), (3, 4), (3, 7), (4, 5),
                          (5, 3), (5, 6), (7, 4), (7, 6), (8, 1), (8, 7)])
        C = [[3, 4, 5, 7], [1, 2, 8], [6]]
        self.gc.append((G, C))

        G = nx.DiGraph()
        G.add_edges_from([(1, 2), (1, 3), (1, 4), (4, 2), (3, 4), (2, 3)])
        C = [[2, 3, 4], [1]]
        self.gc.append((G, C))

        G = nx.DiGraph()
        G.add_edges_from([(1, 2), (2, 3), (3, 2), (2, 1)])
        C = [[1, 2, 3]]
        self.gc.append((G, C))

        # Eppstein's tests
        G = nx.DiGraph({0: [1], 1: [2, 3], 2: [4, 5], 3: [4, 5], 4: [6], 5: [], 6: []})
        C = [[0], [1], [2], [3], [4], [5], [6]]
        self.gc.append((G, C))

        G = nx.DiGraph({0: [1], 1: [2, 3, 4], 2: [0, 3], 3: [4], 4: [3]})
        C = [[0, 1, 2], [3, 4]]
        self.gc.append((G, C))

        G = nx.DiGraph()
        C = []
        self.gc.append((G, C))

    def test_connected_components(self):
        cc = nx.connected_components
        G = self.G
        C = {
            frozenset([0, 1, 2, 3]),
            frozenset([4, 5, 6, 7, 8, 9]),
            frozenset([10, 11, 12, 13, 14])
        }
        assert_equal({frozenset(g) for g in cc(G)}, C)

    def test_number_connected_components(self):
        ncc = nx.number_connected_components
        assert_equal(ncc(self.G), 3)

    def test_number_connected_components2(self):
        ncc = nx.number_connected_components
        assert_equal(ncc(self.grid), 1)

    def test_connected_components2(self):
        cc = nx.connected_components
        G = self.grid
        C = {frozenset([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16])}
        assert_equal({frozenset(g) for g in cc(G)}, C)

    def test_node_connected_components(self):
        ncc = nx.node_connected_component
        G = self.grid
        C = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16}
        assert_equal(ncc(G, 1), C)

    # deprecated
    def test_connected_component_subgraphs(self):
        wcc = nx.weakly_connected_component_subgraphs
        cc = nx.connected_component_subgraphs
        for G, C in self.gc:
            U = G.to_undirected()
            w = {frozenset(g) for g in wcc(G)}
            c = {frozenset(g) for g in cc(U)}
            assert_equal(w, c)

    def test_is_connected(self):
        assert_true(nx.is_connected(self.grid))
        G = nx.Graph()
        G.add_nodes_from([1, 2])
        assert_false(nx.is_connected(G))

    def test_connected_raise(self):
        assert_raises(NetworkXNotImplemented, nx.connected_components, self.DG)
        assert_raises(NetworkXNotImplemented, nx.number_connected_components, self.DG)
        assert_raises(NetworkXNotImplemented, nx.node_connected_component, self.DG, 1)
        assert_raises(NetworkXNotImplemented, nx.is_connected, self.DG)
        assert_raises(nx.NetworkXPointlessConcept, nx.is_connected, nx.Graph())
        # deprecated
        assert_raises(NetworkXNotImplemented, nx.connected_component_subgraphs, self.DG)
