from nose.tools import assert_equal
from nose.tools import assert_raises
from nose.tools import raises

from math import sqrt
from random import random, choice

import networkx as nx
from networkx.utils import pairwise


def dist(a, b):
    """Returns the Euclidean distance between points `a` and `b`."""
    return sqrt(sum((x1 - x2) ** 2 for x1, x2 in zip(a, b)))


class TestAStar:

    def setUp(self):
        edges = [('s', 'u', 10), ('s', 'x', 5), ('u', 'v', 1), ('u', 'x', 2),
                 ('v', 'y', 1), ('x', 'u', 3), ('x', 'v', 5), ('x', 'y', 2),
                 ('y', 's', 7), ('y', 'v', 6)]
        self.XG = nx.DiGraph()
        self.XG.add_weighted_edges_from(edges)

    def test_random_graph(self):
        """Tests that the A* shortest path agrees with Dijkstra's
        shortest path for a random graph.

        """

        G = nx.Graph()

        points = [(random(), random()) for _ in range(100)]

        # Build a path from points[0] to points[-1] to be sure it exists
        for p1, p2 in pairwise(points):
            G.add_edge(p1, p2, weight=dist(p1, p2))

        # Add other random edges
        for _ in range(100):
            p1, p2 = choice(points), choice(points)
            G.add_edge(p1, p2, weight=dist(p1, p2))

        path = nx.astar_path(G, points[0], points[-1], dist)
        assert_equal(path, nx.dijkstra_path(G, points[0], points[-1]))

    def test_astar_directed(self):
        assert_equal(nx.astar_path(self.XG, 's', 'v'), ['s', 'x', 'u', 'v'])
        assert_equal(nx.astar_path_length(self.XG, 's', 'v'), 9)

    def test_astar_multigraph(self):
        G = nx.MultiDiGraph(self.XG)
        assert_raises(nx.NetworkXNotImplemented, nx.astar_path, G, 's', 'v')
        assert_raises(nx.NetworkXNotImplemented, nx.astar_path_length,
                      G, 's', 'v')

    def test_astar_undirected(self):
        GG = self.XG.to_undirected()
        # make sure we get lower weight
        # to_undirected might choose either edge with weight 2 or weight 3
        GG['u']['x']['weight'] = 2
        GG['y']['v']['weight'] = 2
        assert_equal(nx.astar_path(GG, 's', 'v'), ['s', 'x', 'u', 'v'])
        assert_equal(nx.astar_path_length(GG, 's', 'v'), 8)

    def test_astar_directed2(self):
        XG2 = nx.DiGraph()
        edges = [(1, 4, 1), (4, 5, 1), (5, 6, 1), (6, 3, 1), (1, 3, 50),
                 (1, 2, 100), (2, 3, 100)]
        XG2.add_weighted_edges_from(edges)
        assert_equal(nx.astar_path(XG2, 1, 3), [1, 4, 5, 6, 3])

    def test_astar_undirected2(self):
        XG3 = nx.Graph()
        edges = [(0, 1, 2), (1, 2, 12), (2, 3, 1), (3, 4, 5), (4, 5, 1),
                 (5, 0, 10)]
        XG3.add_weighted_edges_from(edges)
        assert_equal(nx.astar_path(XG3, 0, 3), [0, 1, 2, 3])
        assert_equal(nx.astar_path_length(XG3, 0, 3), 15)

    def test_astar_undirected3(self):
        XG4 = nx.Graph()
        edges = [(0, 1, 2), (1, 2, 2), (2, 3, 1), (3, 4, 1), (4, 5, 1),
                 (5, 6, 1), (6, 7, 1), (7, 0, 1)]
        XG4.add_weighted_edges_from(edges)
        assert_equal(nx.astar_path(XG4, 0, 2), [0, 1, 2])
        assert_equal(nx.astar_path_length(XG4, 0, 2), 4)

# >>> MXG4=NX.MultiGraph(XG4)
# >>> MXG4.add_edge(0,1,3)
# >>> NX.dijkstra_path(MXG4,0,2)
# [0, 1, 2]

    def test_astar_w1(self):
        G = nx.DiGraph()
        G.add_edges_from([('s', 'u'), ('s', 'x'), ('u', 'v'), ('u', 'x'),
                          ('v', 'y'), ('x', 'u'), ('x', 'w'), ('w', 'v'),
                          ('x', 'y'), ('y', 's'), ('y', 'v')])
        assert_equal(nx.astar_path(G, 's', 'v'), ['s', 'u', 'v'])
        assert_equal(nx.astar_path_length(G, 's', 'v'), 2)

    @raises(nx.NodeNotFound)
    def test_astar_nopath(self):
        nx.astar_path(self.XG, 's', 'moon')

    def test_cycle(self):
        C = nx.cycle_graph(7)
        assert_equal(nx.astar_path(C, 0, 3), [0, 1, 2, 3])
        assert_equal(nx.dijkstra_path(C, 0, 4), [0, 6, 5, 4])

    def test_unorderable_nodes(self):
        """Tests that A* accomodates nodes that are not orderable.

        For more information, see issue #554.

        """
        # TODO In Python 3, instances of the `object` class are
        # unorderable by default, so we wouldn't need to define our own
        # class here, we could just instantiate an instance of the
        # `object` class. However, we still support Python 2; when
        # support for Python 2 is dropped, this test can be simplified
        # by replacing `Unorderable()` by `object()`.
        class Unorderable(object):

            def __le__(self):
                raise NotImplemented

            def __ge__(self):
                raise NotImplemented

        # Create the cycle graph on four nodes, with nodes represented
        # as (unorderable) Python objects.
        nodes = [Unorderable() for n in range(4)]
        G = nx.Graph()
        G.add_edges_from(pairwise(nodes, cyclic=True))
        path = nx.astar_path(G, nodes[0], nodes[2])
        assert_equal(len(path), 3)
