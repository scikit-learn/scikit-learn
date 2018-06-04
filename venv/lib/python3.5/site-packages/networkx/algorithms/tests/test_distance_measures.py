#!/usr/bin/env python
from nose.tools import *
import networkx


class TestDistance:

    def setUp(self):
        G = networkx.Graph()
        from networkx import convert_node_labels_to_integers as cnlti
        G = cnlti(networkx.grid_2d_graph(4, 4), first_label=1, ordering="sorted")
        self.G = G

    def test_eccentricity(self):
        assert_equal(networkx.eccentricity(self.G, 1), 6)
        e = networkx.eccentricity(self.G)
        assert_equal(e[1], 6)
        sp = dict(networkx.shortest_path_length(self.G))
        e = networkx.eccentricity(self.G, sp=sp)
        assert_equal(e[1], 6)
        e = networkx.eccentricity(self.G, v=1)
        assert_equal(e, 6)
        # This behavior changed in version 1.8 (ticket #739)
        e = networkx.eccentricity(self.G, v=[1, 1])
        assert_equal(e[1], 6)
        e = networkx.eccentricity(self.G, v=[1, 2])
        assert_equal(e[1], 6)
        # test against graph with one node
        G = networkx.path_graph(1)
        e = networkx.eccentricity(G)
        assert_equal(e[0], 0)
        e = networkx.eccentricity(G, v=0)
        assert_equal(e, 0)
        assert_raises(networkx.NetworkXError, networkx.eccentricity, G, 1)
        # test against empty graph
        G = networkx.empty_graph()
        e = networkx.eccentricity(G)
        assert_equal(e, {})

    def test_diameter(self):
        assert_equal(networkx.diameter(self.G), 6)

    def test_radius(self):
        assert_equal(networkx.radius(self.G), 4)

    def test_periphery(self):
        assert_equal(set(networkx.periphery(self.G)), set([1, 4, 13, 16]))

    def test_center(self):
        assert_equal(set(networkx.center(self.G)), set([6, 7, 10, 11]))

    def test_bound_diameter(self):
        assert_equal(networkx.diameter(self.G, usebounds=True), 6)

    def test_bound_radius(self):
        assert_equal(networkx.radius(self.G, usebounds=True), 4)

    def test_bound_periphery(self):
        assert_equal(set(networkx.periphery(self.G, usebounds=True)), set([1, 4, 13, 16]))

    def test_bound_center(self):
        assert_equal(set(networkx.center(self.G, usebounds=True)), set([6, 7, 10, 11]))

    def test_radius_exception(self):
        G = networkx.Graph()
        G.add_edge(1, 2)
        G.add_edge(3, 4)
        assert_raises(networkx.NetworkXError, networkx.diameter, G)

    @raises(networkx.NetworkXError)
    def test_eccentricity_infinite(self):
        G = networkx.Graph([(1, 2), (3, 4)])
        e = networkx.eccentricity(G)

    @raises(networkx.NetworkXError)
    def test_eccentricity_undirected_not_connected(self):
        G = networkx.Graph([(1, 2), (3, 4)])
        e = networkx.eccentricity(G, sp=1)

    @raises(networkx.NetworkXError)
    def test_eccentricity_directed_weakly_connected(self):
        DG = networkx.DiGraph([(1, 2), (1, 3)])
        networkx.eccentricity(DG)
