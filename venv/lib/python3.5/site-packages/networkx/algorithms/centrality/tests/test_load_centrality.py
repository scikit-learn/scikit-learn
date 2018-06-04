#!/usr/bin/env python
from nose.tools import *
import networkx as nx


class TestLoadCentrality:

    def setUp(self):

        G = nx.Graph()
        G.add_edge(0, 1, weight=3)
        G.add_edge(0, 2, weight=2)
        G.add_edge(0, 3, weight=6)
        G.add_edge(0, 4, weight=4)
        G.add_edge(1, 3, weight=5)
        G.add_edge(1, 5, weight=5)
        G.add_edge(2, 4, weight=1)
        G.add_edge(3, 4, weight=2)
        G.add_edge(3, 5, weight=1)
        G.add_edge(4, 5, weight=4)
        self.G = G
        self.exact_weighted = {0: 4.0, 1: 0.0, 2: 8.0, 3: 6.0, 4: 8.0, 5: 0.0}
        self.K = nx.krackhardt_kite_graph()
        self.P3 = nx.path_graph(3)
        self.P4 = nx.path_graph(4)
        self.K5 = nx.complete_graph(5)

        self.C4 = nx.cycle_graph(4)
        self.T = nx.balanced_tree(r=2, h=2)
        self.Gb = nx.Graph()
        self.Gb.add_edges_from([(0, 1), (0, 2), (1, 3), (2, 3),
                                (2, 4), (4, 5), (3, 5)])
        self.F = nx.florentine_families_graph()
        self.D = nx.cycle_graph(3, create_using=nx.DiGraph())
        self.D.add_edges_from([(3, 0), (4, 3)])

    def test_not_strongly_connected(self):
        b = nx.load_centrality(self.D)
        result = {0: 5. / 12,
                  1: 1. / 4,
                  2: 1. / 12,
                  3: 1. / 4,
                  4: 0.000}
        for n in sorted(self.D):
            assert_almost_equal(result[n], b[n], places=3)
            assert_almost_equal(result[n], nx.load_centrality(self.D, n), places=3)

    def test_weighted_load(self):
        b = nx.load_centrality(self.G, weight='weight', normalized=False)
        for n in sorted(self.G):
            assert_equal(b[n], self.exact_weighted[n])

    def test_k5_load(self):
        G = self.K5
        c = nx.load_centrality(G)
        d = {0: 0.000,
             1: 0.000,
             2: 0.000,
             3: 0.000,
             4: 0.000}
        for n in sorted(G):
            assert_almost_equal(c[n], d[n], places=3)

    def test_p3_load(self):
        G = self.P3
        c = nx.load_centrality(G)
        d = {0: 0.000,
             1: 1.000,
             2: 0.000}
        for n in sorted(G):
            assert_almost_equal(c[n], d[n], places=3)
        c = nx.load_centrality(G, v=1)
        assert_almost_equal(c, 1.0)
        c = nx.load_centrality(G, v=1, normalized=True)
        assert_almost_equal(c, 1.0)

    def test_p2_load(self):
        G = nx.path_graph(2)
        c = nx.load_centrality(G)
        d = {0: 0.000,
             1: 0.000}
        for n in sorted(G):
            assert_almost_equal(c[n], d[n], places=3)

    def test_krackhardt_load(self):
        G = self.K
        c = nx.load_centrality(G)
        d = {0: 0.023,
             1: 0.023,
             2: 0.000,
             3: 0.102,
             4: 0.000,
             5: 0.231,
             6: 0.231,
             7: 0.389,
             8: 0.222,
             9: 0.000}
        for n in sorted(G):
            assert_almost_equal(c[n], d[n], places=3)

    def test_florentine_families_load(self):
        G = self.F
        c = nx.load_centrality(G)
        d = {'Acciaiuoli':    0.000,
             'Albizzi':       0.211,
             'Barbadori':     0.093,
             'Bischeri':      0.104,
             'Castellani':    0.055,
             'Ginori':        0.000,
             'Guadagni':      0.251,
             'Lamberteschi':  0.000,
             'Medici':        0.522,
             'Pazzi':         0.000,
             'Peruzzi':       0.022,
             'Ridolfi':       0.117,
             'Salviati':      0.143,
             'Strozzi':       0.106,
             'Tornabuoni':    0.090}
        for n in sorted(G):
            assert_almost_equal(c[n], d[n], places=3)

    def test_unnormalized_k5_load(self):
        G = self.K5
        c = nx.load_centrality(G, normalized=False)
        d = {0: 0.000,
             1: 0.000,
             2: 0.000,
             3: 0.000,
             4: 0.000}
        for n in sorted(G):
            assert_almost_equal(c[n], d[n], places=3)

    def test_unnormalized_p3_load(self):
        G = self.P3
        c = nx.load_centrality(G, normalized=False)
        d = {0: 0.000,
             1: 2.000,
             2: 0.000}
        for n in sorted(G):
            assert_almost_equal(c[n], d[n], places=3)

    def test_unnormalized_krackhardt_load(self):
        G = self.K
        c = nx.load_centrality(G, normalized=False)
        d = {0: 1.667,
             1: 1.667,
             2: 0.000,
             3: 7.333,
             4: 0.000,
             5: 16.667,
             6: 16.667,
             7: 28.000,
             8: 16.000,
             9: 0.000}

        for n in sorted(G):
            assert_almost_equal(c[n], d[n], places=3)

    def test_unnormalized_florentine_families_load(self):
        G = self.F
        c = nx.load_centrality(G, normalized=False)

        d = {'Acciaiuoli':  0.000,
             'Albizzi':    38.333,
             'Barbadori':  17.000,
             'Bischeri':   19.000,
             'Castellani': 10.000,
             'Ginori':     0.000,
             'Guadagni':   45.667,
             'Lamberteschi': 0.000,
             'Medici':     95.000,
             'Pazzi':      0.000,
             'Peruzzi':    4.000,
             'Ridolfi':    21.333,
             'Salviati':   26.000,
             'Strozzi':    19.333,
             'Tornabuoni': 16.333}
        for n in sorted(G):
            assert_almost_equal(c[n], d[n], places=3)

    def test_load_betweenness_difference(self):
        # Difference Between Load and Betweenness
        # --------------------------------------- The smallest graph
        # that shows the difference between load and betweenness is
        # G=ladder_graph(3) (Graph B below)

        # Graph A and B are from Tao Zhou, Jian-Guo Liu, Bing-Hong
        # Wang: Comment on "Scientific collaboration
        # networks. II. Shortest paths, weighted networks, and
        # centrality". https://arxiv.org/pdf/physics/0511084

        # Notice that unlike here, their calculation adds to 1 to the
        # betweennes of every node i for every path from i to every
        # other node.  This is exactly what it should be, based on
        # Eqn. (1) in their paper: the eqn is B(v) = \sum_{s\neq t,
        # s\neq v}{\frac{\sigma_{st}(v)}{\sigma_{st}}}, therefore,
        # they allow v to be the target node.

        # We follow Brandes 2001, who follows Freeman 1977 that make
        # the sum for betweenness of v exclude paths where v is either
        # the source or target node.  To agree with their numbers, we
        # must additionally, remove edge (4,8) from the graph, see AC
        # example following (there is a mistake in the figure in their
        # paper - personal communication).

        # A = nx.Graph()
        # A.add_edges_from([(0,1), (1,2), (1,3), (2,4),
        #                  (3,5), (4,6), (4,7), (4,8),
        #                  (5,8), (6,9), (7,9), (8,9)])
        B = nx.Graph()  # ladder_graph(3)
        B.add_edges_from([(0, 1), (0, 2), (1, 3), (2, 3), (2, 4), (4, 5), (3, 5)])
        c = nx.load_centrality(B, normalized=False)
        d = {0: 1.750,
             1: 1.750,
             2: 6.500,
             3: 6.500,
             4: 1.750,
             5: 1.750}
        for n in sorted(B):
            assert_almost_equal(c[n], d[n], places=3)

    def test_c4_edge_load(self):
        G = self.C4
        c = nx.edge_load_centrality(G)
        d = {(0, 1): 6.000,
             (0, 3): 6.000,
             (1, 2): 6.000,
             (2, 3): 6.000}
        for n in G.edges():
            assert_almost_equal(c[n], d[n], places=3)

    def test_p4_edge_load(self):
        G = self.P4
        c = nx.edge_load_centrality(G)
        d = {(0, 1): 6.000,
             (1, 2): 8.000,
             (2, 3): 6.000}
        for n in G.edges():
            assert_almost_equal(c[n], d[n], places=3)

    def test_k5_edge_load(self):
        G = self.K5
        c = nx.edge_load_centrality(G)
        d = {(0, 1): 5.000,
             (0, 2): 5.000,
             (0, 3): 5.000,
             (0, 4): 5.000,
             (1, 2): 5.000,
             (1, 3): 5.000,
             (1, 4): 5.000,
             (2, 3): 5.000,
             (2, 4): 5.000,
             (3, 4): 5.000}
        for n in G.edges():
            assert_almost_equal(c[n], d[n], places=3)

    def test_tree_edge_load(self):
        G = self.T
        c = nx.edge_load_centrality(G)
        d = {(0, 1): 24.000,
             (0, 2): 24.000,
             (1, 3): 12.000,
             (1, 4): 12.000,
             (2, 5): 12.000,
             (2, 6): 12.000}
        for n in G.edges():
            assert_almost_equal(c[n], d[n], places=3)
