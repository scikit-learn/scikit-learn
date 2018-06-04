#!/usr/bin/env python
from nose.tools import *
import networkx as nx


def weighted_G():
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

    return G


class TestBetweennessCentrality(object):

    def test_K5(self):
        """Betweenness centrality: K5"""
        G = nx.complete_graph(5)
        b = nx.betweenness_centrality(G,
                                      weight=None,
                                      normalized=False)
        b_answer = {0: 0.0, 1: 0.0, 2: 0.0, 3: 0.0, 4: 0.0}
        for n in sorted(G):
            assert_almost_equal(b[n], b_answer[n])

    def test_K5_endpoints(self):
        """Betweenness centrality: K5 endpoints"""
        G = nx.complete_graph(5)
        b = nx.betweenness_centrality(G,
                                      weight=None,
                                      normalized=False,
                                      endpoints=True)
        b_answer = {0: 4.0, 1: 4.0, 2: 4.0, 3: 4.0, 4: 4.0}
        for n in sorted(G):
            assert_almost_equal(b[n], b_answer[n])

    def test_P3_normalized(self):
        """Betweenness centrality: P3 normalized"""
        G = nx.path_graph(3)
        b = nx.betweenness_centrality(G,
                                      weight=None,
                                      normalized=True)
        b_answer = {0: 0.0, 1: 1.0, 2: 0.0}
        for n in sorted(G):
            assert_almost_equal(b[n], b_answer[n])

    def test_P3(self):
        """Betweenness centrality: P3"""
        G = nx.path_graph(3)
        b_answer = {0: 0.0, 1: 1.0, 2: 0.0}
        b = nx.betweenness_centrality(G,
                                      weight=None,
                                      normalized=False)
        for n in sorted(G):
            assert_almost_equal(b[n], b_answer[n])

    def test_P3_endpoints(self):
        """Betweenness centrality: P3 endpoints"""
        G = nx.path_graph(3)
        b_answer = {0: 2.0, 1: 3.0, 2: 2.0}
        b = nx.betweenness_centrality(G,
                                      weight=None,
                                      normalized=False,
                                      endpoints=True)
        for n in sorted(G):
            assert_almost_equal(b[n], b_answer[n])

    def test_krackhardt_kite_graph(self):
        """Betweenness centrality: Krackhardt kite graph"""
        G = nx.krackhardt_kite_graph()
        b_answer = {0: 1.667, 1: 1.667, 2: 0.000, 3: 7.333, 4: 0.000,
                    5: 16.667, 6: 16.667, 7: 28.000, 8: 16.000, 9: 0.000}
        for b in b_answer:
            b_answer[b] /= 2.0
        b = nx.betweenness_centrality(G,
                                      weight=None,
                                      normalized=False)

        for n in sorted(G):
            assert_almost_equal(b[n], b_answer[n], places=3)

    def test_krackhardt_kite_graph_normalized(self):
        """Betweenness centrality: Krackhardt kite graph normalized"""
        G = nx.krackhardt_kite_graph()
        b_answer = {0: 0.023, 1: 0.023, 2: 0.000, 3: 0.102, 4: 0.000,
                    5: 0.231, 6: 0.231, 7: 0.389, 8: 0.222, 9: 0.000}
        b = nx.betweenness_centrality(G,
                                      weight=None,
                                      normalized=True)

        for n in sorted(G):
            assert_almost_equal(b[n], b_answer[n], places=3)

    def test_florentine_families_graph(self):
        """Betweenness centrality: Florentine families graph"""
        G = nx.florentine_families_graph()
        b_answer =\
            {'Acciaiuoli':    0.000,
             'Albizzi':       0.212,
             'Barbadori':     0.093,
             'Bischeri':      0.104,
             'Castellani':    0.055,
             'Ginori':        0.000,
             'Guadagni':      0.255,
             'Lamberteschi':  0.000,
             'Medici':        0.522,
             'Pazzi':         0.000,
             'Peruzzi':       0.022,
             'Ridolfi':       0.114,
             'Salviati':      0.143,
             'Strozzi':       0.103,
             'Tornabuoni':    0.092}

        b = nx.betweenness_centrality(G,
                                      weight=None,
                                      normalized=True)
        for n in sorted(G):
            assert_almost_equal(b[n], b_answer[n], places=3)

    def test_ladder_graph(self):
        """Betweenness centrality: Ladder graph"""
        G = nx.Graph()  # ladder_graph(3)
        G.add_edges_from([(0, 1), (0, 2), (1, 3), (2, 3),
                          (2, 4), (4, 5), (3, 5)])
        b_answer = {0: 1.667, 1: 1.667, 2: 6.667,
                    3: 6.667, 4: 1.667, 5: 1.667}
        for b in b_answer:
            b_answer[b] /= 2.0
        b = nx.betweenness_centrality(G,
                                      weight=None,
                                      normalized=False)
        for n in sorted(G):
            assert_almost_equal(b[n], b_answer[n], places=3)

    def test_disconnected_path(self):
        """Betweenness centrality: disconnected path"""
        G = nx.Graph()
        nx.add_path(G, [0, 1, 2])
        nx.add_path(G, [3, 4, 5, 6])
        b_answer = {0: 0, 1: 1, 2: 0, 3: 0, 4: 2, 5: 2, 6: 0}
        b = nx.betweenness_centrality(G,
                                      weight=None,
                                      normalized=False)
        for n in sorted(G):
            assert_almost_equal(b[n], b_answer[n])

    def test_disconnected_path_endpoints(self):
        """Betweenness centrality: disconnected path endpoints"""
        G = nx.Graph()
        nx.add_path(G, [0, 1, 2])
        nx.add_path(G, [3, 4, 5, 6])
        b_answer = {0: 2, 1: 3, 2: 2, 3: 3, 4: 5, 5: 5, 6: 3}
        b = nx.betweenness_centrality(G,
                                      weight=None,
                                      normalized=False,
                                      endpoints=True)
        for n in sorted(G):
            assert_almost_equal(b[n], b_answer[n])

    def test_directed_path(self):
        """Betweenness centrality: directed path"""
        G = nx.DiGraph()
        nx.add_path(G, [0, 1, 2])
        b = nx.betweenness_centrality(G,
                                      weight=None,
                                      normalized=False)
        b_answer = {0: 0.0, 1: 1.0, 2: 0.0}
        for n in sorted(G):
            assert_almost_equal(b[n], b_answer[n])

    def test_directed_path_normalized(self):
        """Betweenness centrality: directed path normalized"""
        G = nx.DiGraph()
        nx.add_path(G, [0, 1, 2])
        b = nx.betweenness_centrality(G,
                                      weight=None,
                                      normalized=True)
        b_answer = {0: 0.0, 1: 0.5, 2: 0.0}
        for n in sorted(G):
            assert_almost_equal(b[n], b_answer[n])


class TestWeightedBetweennessCentrality(object):

    def test_K5(self):
        """Weighted betweenness centrality: K5"""
        G = nx.complete_graph(5)
        b = nx.betweenness_centrality(G,
                                      weight='weight',
                                      normalized=False)
        b_answer = {0: 0.0, 1: 0.0, 2: 0.0, 3: 0.0, 4: 0.0}
        for n in sorted(G):
            assert_almost_equal(b[n], b_answer[n])

    def test_P3_normalized(self):
        """Weighted betweenness centrality: P3 normalized"""
        G = nx.path_graph(3)
        b = nx.betweenness_centrality(G,
                                      weight='weight',
                                      normalized=True)
        b_answer = {0: 0.0, 1: 1.0, 2: 0.0}
        for n in sorted(G):
            assert_almost_equal(b[n], b_answer[n])

    def test_P3(self):
        """Weighted betweenness centrality: P3"""
        G = nx.path_graph(3)
        b_answer = {0: 0.0, 1: 1.0, 2: 0.0}
        b = nx.betweenness_centrality(G,
                                      weight='weight',
                                      normalized=False)
        for n in sorted(G):
            assert_almost_equal(b[n], b_answer[n])

    def test_krackhardt_kite_graph(self):
        """Weighted betweenness centrality: Krackhardt kite graph"""
        G = nx.krackhardt_kite_graph()
        b_answer = {0: 1.667, 1: 1.667, 2: 0.000, 3: 7.333, 4: 0.000,
                    5: 16.667, 6: 16.667, 7: 28.000, 8: 16.000, 9: 0.000}
        for b in b_answer:
            b_answer[b] /= 2.0

        b = nx.betweenness_centrality(G,
                                      weight='weight',
                                      normalized=False)

        for n in sorted(G):
            assert_almost_equal(b[n], b_answer[n], places=3)

    def test_krackhardt_kite_graph_normalized(self):
        """Weighted betweenness centrality: 
        Krackhardt kite graph normalized
        """
        G = nx.krackhardt_kite_graph()
        b_answer = {0: 0.023, 1: 0.023, 2: 0.000, 3: 0.102, 4: 0.000,
                    5: 0.231, 6: 0.231, 7: 0.389, 8: 0.222, 9: 0.000}
        b = nx.betweenness_centrality(G,
                                      weight='weight',
                                      normalized=True)

        for n in sorted(G):
            assert_almost_equal(b[n], b_answer[n], places=3)

    def test_florentine_families_graph(self):
        """Weighted betweenness centrality: 
        Florentine families graph"""
        G = nx.florentine_families_graph()
        b_answer =\
            {'Acciaiuoli':    0.000,
             'Albizzi':       0.212,
             'Barbadori':     0.093,
             'Bischeri':      0.104,
             'Castellani':    0.055,
             'Ginori':        0.000,
             'Guadagni':      0.255,
             'Lamberteschi':  0.000,
             'Medici':        0.522,
             'Pazzi':         0.000,
             'Peruzzi':       0.022,
             'Ridolfi':       0.114,
             'Salviati':      0.143,
             'Strozzi':       0.103,
             'Tornabuoni':    0.092}

        b = nx.betweenness_centrality(G,
                                      weight='weight',
                                      normalized=True)
        for n in sorted(G):
            assert_almost_equal(b[n], b_answer[n], places=3)

    def test_ladder_graph(self):
        """Weighted betweenness centrality: Ladder graph"""
        G = nx.Graph()  # ladder_graph(3)
        G.add_edges_from([(0, 1), (0, 2), (1, 3), (2, 3),
                          (2, 4), (4, 5), (3, 5)])
        b_answer = {0: 1.667, 1: 1.667, 2: 6.667,
                    3: 6.667, 4: 1.667, 5: 1.667}
        for b in b_answer:
            b_answer[b] /= 2.0
        b = nx.betweenness_centrality(G,
                                      weight='weight',
                                      normalized=False)
        for n in sorted(G):
            assert_almost_equal(b[n], b_answer[n], places=3)

    def test_G(self):
        """Weighted betweenness centrality: G"""
        G = weighted_G()
        b_answer = {0: 2.0, 1: 0.0, 2: 4.0, 3: 3.0, 4: 4.0, 5: 0.0}
        b = nx.betweenness_centrality(G,
                                      weight='weight',
                                      normalized=False)
        for n in sorted(G):
            assert_almost_equal(b[n], b_answer[n])

    def test_G2(self):
        """Weighted betweenness centrality: G2"""
        G = nx.DiGraph()
        G.add_weighted_edges_from([('s', 'u', 10), ('s', 'x', 5),
                                   ('u', 'v', 1), ('u', 'x', 2),
                                   ('v', 'y', 1), ('x', 'u', 3),
                                   ('x', 'v', 5), ('x', 'y', 2),
                                   ('y', 's', 7), ('y', 'v', 6)])

        b_answer = {'y': 5.0, 'x': 5.0, 's': 4.0, 'u': 2.0, 'v': 2.0}

        b = nx.betweenness_centrality(G,
                                      weight='weight',
                                      normalized=False)
        for n in sorted(G):
            assert_almost_equal(b[n], b_answer[n])


class TestEdgeBetweennessCentrality(object):

    def test_K5(self):
        """Edge betweenness centrality: K5"""
        G = nx.complete_graph(5)
        b = nx.edge_betweenness_centrality(G, weight=None, normalized=False)
        b_answer = dict.fromkeys(G.edges(), 1)
        for n in sorted(G.edges()):
            assert_almost_equal(b[n], b_answer[n])

    def test_normalized_K5(self):
        """Edge betweenness centrality: K5"""
        G = nx.complete_graph(5)
        b = nx.edge_betweenness_centrality(G, weight=None, normalized=True)
        b_answer = dict.fromkeys(G.edges(), 1 / 10.0)
        for n in sorted(G.edges()):
            assert_almost_equal(b[n], b_answer[n])

    def test_C4(self):
        """Edge betweenness centrality: C4"""
        G = nx.cycle_graph(4)
        b = nx.edge_betweenness_centrality(G, weight=None, normalized=True)
        b_answer = {(0, 1): 2, (0, 3): 2, (1, 2): 2, (2, 3): 2}
        for n in sorted(G.edges()):
            assert_almost_equal(b[n], b_answer[n] / 6.0)

    def test_P4(self):
        """Edge betweenness centrality: P4"""
        G = nx.path_graph(4)
        b = nx.edge_betweenness_centrality(G, weight=None, normalized=False)
        b_answer = {(0, 1): 3, (1, 2): 4, (2, 3): 3}
        for n in sorted(G.edges()):
            assert_almost_equal(b[n], b_answer[n])

    def test_normalized_P4(self):
        """Edge betweenness centrality: P4"""
        G = nx.path_graph(4)
        b = nx.edge_betweenness_centrality(G, weight=None, normalized=True)
        b_answer = {(0, 1): 3, (1, 2): 4, (2, 3): 3}
        for n in sorted(G.edges()):
            assert_almost_equal(b[n], b_answer[n] / 6.0)

    def test_balanced_tree(self):
        """Edge betweenness centrality: balanced tree"""
        G = nx.balanced_tree(r=2, h=2)
        b = nx.edge_betweenness_centrality(G, weight=None, normalized=False)
        b_answer = {(0, 1): 12, (0, 2): 12,
                    (1, 3): 6, (1, 4): 6, (2, 5): 6, (2, 6): 6}
        for n in sorted(G.edges()):
            assert_almost_equal(b[n], b_answer[n])


class TestWeightedEdgeBetweennessCentrality(object):

    def test_K5(self):
        """Edge betweenness centrality: K5"""
        G = nx.complete_graph(5)
        b = nx.edge_betweenness_centrality(G, weight='weight', normalized=False)
        b_answer = dict.fromkeys(G.edges(), 1)
        for n in sorted(G.edges()):
            assert_almost_equal(b[n], b_answer[n])

    def test_C4(self):
        """Edge betweenness centrality: C4"""
        G = nx.cycle_graph(4)
        b = nx.edge_betweenness_centrality(G, weight='weight', normalized=False)
        b_answer = {(0, 1): 2, (0, 3): 2, (1, 2): 2, (2, 3): 2}
        for n in sorted(G.edges()):
            assert_almost_equal(b[n], b_answer[n])

    def test_P4(self):
        """Edge betweenness centrality: P4"""
        G = nx.path_graph(4)
        b = nx.edge_betweenness_centrality(G, weight='weight', normalized=False)
        b_answer = {(0, 1): 3, (1, 2): 4, (2, 3): 3}
        for n in sorted(G.edges()):
            assert_almost_equal(b[n], b_answer[n])

    def test_balanced_tree(self):
        """Edge betweenness centrality: balanced tree"""
        G = nx.balanced_tree(r=2, h=2)
        b = nx.edge_betweenness_centrality(G, weight='weight', normalized=False)
        b_answer = {(0, 1): 12, (0, 2): 12,
                    (1, 3): 6, (1, 4): 6, (2, 5): 6, (2, 6): 6}
        for n in sorted(G.edges()):
            assert_almost_equal(b[n], b_answer[n])

    def test_weighted_graph(self):
        eList = [(0, 1, 5), (0, 2, 4), (0, 3, 3),
                 (0, 4, 2), (1, 2, 4), (1, 3, 1),
                 (1, 4, 3), (2, 4, 5), (3, 4, 4)]
        G = nx.Graph()
        G.add_weighted_edges_from(eList)
        b = nx.edge_betweenness_centrality(G, weight='weight', normalized=False)
        b_answer = {(0, 1): 0.0,
                    (0, 2): 1.0,
                    (0, 3): 2.0,
                    (0, 4): 1.0,
                    (1, 2): 2.0,
                    (1, 3): 3.5,
                    (1, 4): 1.5,
                    (2, 4): 1.0,
                    (3, 4): 0.5}

        for n in sorted(G.edges()):
            assert_almost_equal(b[n], b_answer[n])

    def test_normalized_weighted_graph(self):
        eList = [(0, 1, 5), (0, 2, 4), (0, 3, 3),
                 (0, 4, 2), (1, 2, 4), (1, 3, 1),
                 (1, 4, 3), (2, 4, 5), (3, 4, 4)]
        G = nx.Graph()
        G.add_weighted_edges_from(eList)
        b = nx.edge_betweenness_centrality(G, weight='weight', normalized=True)
        b_answer = {(0, 1): 0.0,
                    (0, 2): 1.0,
                    (0, 3): 2.0,
                    (0, 4): 1.0,
                    (1, 2): 2.0,
                    (1, 3): 3.5,
                    (1, 4): 1.5,
                    (2, 4): 1.0,
                    (3, 4): 0.5}

        norm = len(G) * (len(G) - 1) / 2.0
        for n in sorted(G.edges()):
            assert_almost_equal(b[n], b_answer[n] / norm)
