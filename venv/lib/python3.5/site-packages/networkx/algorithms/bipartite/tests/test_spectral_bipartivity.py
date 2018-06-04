# -*- coding: utf-8 -*-
from nose import SkipTest
from nose.tools import *
import networkx as nx
from networkx.algorithms.bipartite import spectral_bipartivity as sb

# Examples from Figure 1
# E. Estrada and J. A. Rodríguez-Velázquez, "Spectral measures of
# bipartivity in complex networks", PhysRev E 72, 046105 (2005)


class TestSpectralBipartivity(object):
    @classmethod
    def setupClass(cls):
        global scipy
        global assert_equal
        global assert_almost_equal
        try:
            import scipy.linalg
        except ImportError:
            raise SkipTest('SciPy not available.')

    def test_star_like(self):
        # star-like

        G = nx.star_graph(2)
        G.add_edge(1, 2)
        assert_almost_equal(sb(G), 0.843, places=3)

        G = nx.star_graph(3)
        G.add_edge(1, 2)
        assert_almost_equal(sb(G), 0.871, places=3)

        G = nx.star_graph(4)
        G.add_edge(1, 2)
        assert_almost_equal(sb(G), 0.890, places=3)

    def k23_like(self):
        # K2,3-like
        G = nx.complete_bipartite_graph(2, 3)
        G.add_edge(0, 1)
        assert_almost_equal(sb(G), 0.769, places=3)

        G = nx.complete_bipartite_graph(2, 3)
        G.add_edge(2, 4)
        assert_almost_equal(sb(G), 0.829, places=3)

        G = nx.complete_bipartite_graph(2, 3)
        G.add_edge(2, 4)
        G.add_edge(3, 4)
        assert_almost_equal(sb(G), 0.731, places=3)

        G = nx.complete_bipartite_graph(2, 3)
        G.add_edge(0, 1)
        G.add_edge(2, 4)
        assert_almost_equal(sb(G), 0.692, places=3)

        G = nx.complete_bipartite_graph(2, 3)
        G.add_edge(2, 4)
        G.add_edge(3, 4)
        G.add_edge(0, 1)
        assert_almost_equal(sb(G), 0.645, places=3)

        G = nx.complete_bipartite_graph(2, 3)
        G.add_edge(2, 4)
        G.add_edge(3, 4)
        G.add_edge(2, 3)
        assert_almost_equal(sb(G), 0.645, places=3)

        G = nx.complete_bipartite_graph(2, 3)
        G.add_edge(2, 4)
        G.add_edge(3, 4)
        G.add_edge(2, 3)
        G.add_edge(0, 1)
        assert_almost_equal(sb(G), 0.597, places=3)

    def test_single_nodes(self):

        # single nodes
        G = nx.complete_bipartite_graph(2, 3)
        G.add_edge(2, 4)
        sbn = sb(G, nodes=[1, 2])
        assert_almost_equal(sbn[1], 0.85, places=2)
        assert_almost_equal(sbn[2], 0.77, places=2)

        G = nx.complete_bipartite_graph(2, 3)
        G.add_edge(0, 1)
        sbn = sb(G, nodes=[1, 2])
        assert_almost_equal(sbn[1], 0.73, places=2)
        assert_almost_equal(sbn[2], 0.82, places=2)
