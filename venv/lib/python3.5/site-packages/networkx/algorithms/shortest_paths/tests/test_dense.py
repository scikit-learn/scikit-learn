#!/usr/bin/env python
from nose.tools import *
import networkx as nx


class TestFloyd:
    def setUp(self):
        pass

    def test_floyd_warshall_predecessor_and_distance(self):
        XG = nx.DiGraph()
        XG.add_weighted_edges_from([('s', 'u', 10), ('s', 'x', 5),
                                    ('u', 'v', 1), ('u', 'x', 2),
                                    ('v', 'y', 1), ('x', 'u', 3),
                                    ('x', 'v', 5), ('x', 'y', 2),
                                    ('y', 's', 7), ('y', 'v', 6)])
        path, dist = nx.floyd_warshall_predecessor_and_distance(XG)
        assert_equal(dist['s']['v'], 9)
        assert_equal(path['s']['v'], 'u')
        assert_equal(dist,
                     {'y': {'y': 0, 'x': 12, 's': 7, 'u': 15, 'v': 6},
                      'x': {'y': 2, 'x': 0, 's': 9, 'u': 3, 'v': 4},
                      's': {'y': 7, 'x': 5, 's': 0, 'u': 8, 'v': 9},
                      'u': {'y': 2, 'x': 2, 's': 9, 'u': 0, 'v': 1},
                      'v': {'y': 1, 'x': 13, 's': 8, 'u': 16, 'v': 0}})

        GG = XG.to_undirected()
        # make sure we get lower weight
        # to_undirected might choose either edge with weight 2 or weight 3
        GG['u']['x']['weight'] = 2
        path, dist = nx.floyd_warshall_predecessor_and_distance(GG)
        assert_equal(dist['s']['v'], 8)
        # skip this test, could be alternate path s-u-v
#        assert_equal(path['s']['v'],'y')

        G = nx.DiGraph()  # no weights
        G.add_edges_from([('s', 'u'), ('s', 'x'),
                          ('u', 'v'), ('u', 'x'),
                          ('v', 'y'), ('x', 'u'),
                          ('x', 'v'), ('x', 'y'),
                          ('y', 's'), ('y', 'v')])
        path, dist = nx.floyd_warshall_predecessor_and_distance(G)
        assert_equal(dist['s']['v'], 2)
        # skip this test, could be alternate path s-u-v
 #       assert_equal(path['s']['v'],'x')

        # alternate interface
        dist = nx.floyd_warshall(G)
        assert_equal(dist['s']['v'], 2)

    def test_cycle(self):
        path, dist = nx.floyd_warshall_predecessor_and_distance(nx.cycle_graph(7))
        assert_equal(dist[0][3], 3)
        assert_equal(path[0][3], 2)
        assert_equal(dist[0][4], 3)

    def test_weighted(self):
        XG3 = nx.Graph()
        XG3.add_weighted_edges_from([[0, 1, 2], [1, 2, 12], [2, 3, 1],
                                     [3, 4, 5], [4, 5, 1], [5, 0, 10]])
        path, dist = nx.floyd_warshall_predecessor_and_distance(XG3)
        assert_equal(dist[0][3], 15)
        assert_equal(path[0][3], 2)

    def test_weighted2(self):
        XG4 = nx.Graph()
        XG4.add_weighted_edges_from([[0, 1, 2], [1, 2, 2], [2, 3, 1],
                                     [3, 4, 1], [4, 5, 1], [5, 6, 1],
                                     [6, 7, 1], [7, 0, 1]])
        path, dist = nx.floyd_warshall_predecessor_and_distance(XG4)
        assert_equal(dist[0][2], 4)
        assert_equal(path[0][2], 1)

    def test_weight_parameter(self):
        XG4 = nx.Graph()
        XG4.add_edges_from([(0, 1, {'heavy': 2}), (1, 2, {'heavy': 2}),
                            (2, 3, {'heavy': 1}), (3, 4, {'heavy': 1}),
                            (4, 5, {'heavy': 1}), (5, 6, {'heavy': 1}),
                            (6, 7, {'heavy': 1}), (7, 0, {'heavy': 1})])
        path, dist = nx.floyd_warshall_predecessor_and_distance(XG4,
                                                                weight='heavy')
        assert_equal(dist[0][2], 4)
        assert_equal(path[0][2], 1)

    def test_zero_distance(self):
        XG = nx.DiGraph()
        XG.add_weighted_edges_from([('s', 'u', 10), ('s', 'x', 5),
                                    ('u', 'v', 1), ('u', 'x', 2),
                                    ('v', 'y', 1), ('x', 'u', 3),
                                    ('x', 'v', 5), ('x', 'y', 2),
                                    ('y', 's', 7), ('y', 'v', 6)])
        path, dist = nx.floyd_warshall_predecessor_and_distance(XG)

        for u in XG:
            assert_equal(dist[u][u], 0)

        GG = XG.to_undirected()
        # make sure we get lower weight
        # to_undirected might choose either edge with weight 2 or weight 3
        GG['u']['x']['weight'] = 2
        path, dist = nx.floyd_warshall_predecessor_and_distance(GG)

        for u in GG:
            dist[u][u] = 0

    def test_zero_weight(self):
        G = nx.DiGraph()
        edges = [(1, 2, -2), (2, 3, -4), (1, 5, 1), (5, 4, 0), (4, 3, -5), (2, 5, -7)]
        G.add_weighted_edges_from(edges)
        dist = nx.floyd_warshall(G)
        assert_equal(dist[1][3], -14)

        G = nx.MultiDiGraph()
        edges.append((2, 5, -7))
        G.add_weighted_edges_from(edges)
        dist = nx.floyd_warshall(G)
        assert_equal(dist[1][3], -14)
