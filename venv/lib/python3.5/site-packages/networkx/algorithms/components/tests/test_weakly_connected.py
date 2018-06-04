#!/usr/bin/env python
from nose.tools import *
import networkx as nx
from networkx import NetworkXNotImplemented


class TestWeaklyConnected:

    def setUp(self):
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

    def test_weakly_connected_components(self):
        for G, C in self.gc:
            U = G.to_undirected()
            w = {frozenset(g) for g in nx.weakly_connected_components(G)}
            c = {frozenset(g) for g in nx.connected_components(U)}
            assert_equal(w, c)

    def test_number_weakly_connected_components(self):
        for G, C in self.gc:
            U = G.to_undirected()
            w = nx.number_weakly_connected_components(G)
            c = nx.number_connected_components(U)
            assert_equal(w, c)

    # deprecated
    def test_weakly_connected_component_subgraphs(self):
        wcc = nx.weakly_connected_component_subgraphs
        cc = nx.connected_component_subgraphs
        for G, C in self.gc:
            U = G.to_undirected()
            w = {frozenset(g) for g in wcc(G)}
            c = {frozenset(g) for g in cc(U)}
            assert_equal(w, c)

    def test_is_weakly_connected(self):
        for G, C in self.gc:
            U = G.to_undirected()
            assert_equal(nx.is_weakly_connected(G), nx.is_connected(U))

    def test_null_graph(self):
        G = nx.DiGraph()
        assert_equal(list(nx.weakly_connected_components(G)), [])
        assert_equal(nx.number_weakly_connected_components(G), 0)
        assert_raises(nx.NetworkXPointlessConcept, nx.is_weakly_connected, G)

    def test_connected_raise(self):
        G = nx.Graph()
        assert_raises(NetworkXNotImplemented, nx.weakly_connected_components, G)
        assert_raises(NetworkXNotImplemented, nx.number_weakly_connected_components, G)
        assert_raises(NetworkXNotImplemented, nx.is_weakly_connected, G)
        # deprecated
        assert_raises(NetworkXNotImplemented, nx.weakly_connected_component_subgraphs, G)
