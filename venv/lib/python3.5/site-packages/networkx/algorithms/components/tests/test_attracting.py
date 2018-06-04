#!/usr/bin/env python
from nose.tools import *
import networkx as nx
from networkx import NetworkXNotImplemented


class TestAttractingComponents(object):
    def setUp(self):
        self.G1 = nx.DiGraph()
        self.G1.add_edges_from([(5, 11), (11, 2), (11, 9), (11, 10),
                                (7, 11), (7, 8), (8, 9), (3, 8), (3, 10)])
        self.G2 = nx.DiGraph()
        self.G2.add_edges_from([(0, 1), (0, 2), (1, 1), (1, 2), (2, 1)])

        self.G3 = nx.DiGraph()
        self.G3.add_edges_from([(0, 1), (1, 2), (2, 1), (0, 3), (3, 4), (4, 3)])

        self.G4 = nx.DiGraph()

    def test_attracting_components(self):
        ac = list(nx.attracting_components(self.G1))
        assert_true({2} in ac)
        assert_true({9} in ac)
        assert_true({10} in ac)

        ac = list(nx.attracting_components(self.G2))
        ac = [tuple(sorted(x)) for x in ac]
        assert_true(ac == [(1, 2)])

        ac = list(nx.attracting_components(self.G3))
        ac = [tuple(sorted(x)) for x in ac]
        assert_true((1, 2) in ac)
        assert_true((3, 4) in ac)
        assert_equal(len(ac), 2)

        ac = list(nx.attracting_components(self.G4))
        assert_equal(ac, [])

    def test_number_attacting_components(self):
        assert_equal(nx.number_attracting_components(self.G1), 3)
        assert_equal(nx.number_attracting_components(self.G2), 1)
        assert_equal(nx.number_attracting_components(self.G3), 2)
        assert_equal(nx.number_attracting_components(self.G4), 0)

    def test_is_attracting_component(self):
        assert_false(nx.is_attracting_component(self.G1))
        assert_false(nx.is_attracting_component(self.G2))
        assert_false(nx.is_attracting_component(self.G3))
        g2 = self.G3.subgraph([1, 2])
        assert_true(nx.is_attracting_component(g2))
        assert_false(nx.is_attracting_component(self.G4))

    def test_connected_raise(self):
        G = nx.Graph()
        assert_raises(NetworkXNotImplemented, nx.attracting_components, G)
        assert_raises(NetworkXNotImplemented, nx.number_attracting_components, G)
        assert_raises(NetworkXNotImplemented, nx.is_attracting_component, G)
        # deprecated
        assert_raises(NetworkXNotImplemented, nx.attracting_component_subgraphs, G)
