from __future__ import division

from nose.tools import assert_almost_equal
from nose.tools import assert_equal
from nose.tools import assert_false
from nose.tools import assert_true
from nose.tools import assert_raises
from nose.tools import ok_
from nose.tools import raises

import networkx as nx


def validate_grid_path(r, c, s, t, p):
    ok_(isinstance(p, list))
    assert_equal(p[0], s)
    assert_equal(p[-1], t)
    s = ((s - 1) // c, (s - 1) % c)
    t = ((t - 1) // c, (t - 1) % c)
    assert_equal(len(p), abs(t[0] - s[0]) + abs(t[1] - s[1]) + 1)
    p = [((u - 1) // c, (u - 1) % c) for u in p]
    for u in p:
        ok_(0 <= u[0] < r)
        ok_(0 <= u[1] < c)
    for u, v in zip(p[:-1], p[1:]):
        ok_((abs(v[0] - u[0]), abs(v[1] - u[1])) in [(0, 1), (1, 0)])


class TestGenericPath:

    def setUp(self):
        from networkx import convert_node_labels_to_integers as cnlti
        self.grid = cnlti(nx.grid_2d_graph(4, 4), first_label=1, ordering="sorted")
        self.cycle = nx.cycle_graph(7)
        self.directed_cycle = nx.cycle_graph(7, create_using=nx.DiGraph())

    def test_shortest_path(self):
        assert_equal(nx.shortest_path(self.cycle, 0, 3), [0, 1, 2, 3])
        assert_equal(nx.shortest_path(self.cycle, 0, 4), [0, 6, 5, 4])
        validate_grid_path(4, 4, 1, 12, nx.shortest_path(self.grid, 1, 12))
        assert_equal(nx.shortest_path(self.directed_cycle, 0, 3), [0, 1, 2, 3])
        # now with weights
        assert_equal(nx.shortest_path(self.cycle, 0, 3, weight='weight'), [0, 1, 2, 3])
        assert_equal(nx.shortest_path(self.cycle, 0, 4, weight='weight'), [0, 6, 5, 4])
        validate_grid_path(4, 4, 1, 12, nx.shortest_path(self.grid, 1, 12, weight='weight'))
        assert_equal(nx.shortest_path(self.directed_cycle, 0, 3, weight='weight'),
                     [0, 1, 2, 3])

    def test_shortest_path_target(self):
        sp = nx.shortest_path(nx.path_graph(3), target=1)
        assert_equal(sp, {0: [0, 1], 1: [1], 2: [2, 1]})

    def test_shortest_path_length(self):
        assert_equal(nx.shortest_path_length(self.cycle, 0, 3), 3)
        assert_equal(nx.shortest_path_length(self.grid, 1, 12), 5)
        assert_equal(nx.shortest_path_length(self.directed_cycle, 0, 4), 4)
        # now with weights
        assert_equal(nx.shortest_path_length(self.cycle, 0, 3, weight='weight'), 3)
        assert_equal(nx.shortest_path_length(self.grid, 1, 12, weight='weight'), 5)
        assert_equal(nx.shortest_path_length(self.directed_cycle, 0, 4, weight='weight'), 4)

    def test_shortest_path_length_target(self):
        sp = dict(nx.shortest_path_length(nx.path_graph(3), target=1))
        assert_equal(sp[0], 1)
        assert_equal(sp[1], 0)
        assert_equal(sp[2], 1)

    def test_single_source_shortest_path(self):
        p = nx.shortest_path(self.cycle, 0)
        assert_equal(p[3], [0, 1, 2, 3])
        assert_equal(p, nx.single_source_shortest_path(self.cycle, 0))
        p = nx.shortest_path(self.grid, 1)
        validate_grid_path(4, 4, 1, 12, p[12])
        # now with weights
        p = nx.shortest_path(self.cycle, 0, weight='weight')
        assert_equal(p[3], [0, 1, 2, 3])
        assert_equal(p, nx.single_source_dijkstra_path(self.cycle, 0))
        p = nx.shortest_path(self.grid, 1, weight='weight')
        validate_grid_path(4, 4, 1, 12, p[12])

    def test_single_source_shortest_path_length(self):
        l = dict(nx.shortest_path_length(self.cycle, 0))
        assert_equal(l, {0: 0, 1: 1, 2: 2, 3: 3, 4: 3, 5: 2, 6: 1})
        assert_equal(l, dict(nx.single_source_shortest_path_length(self.cycle, 0)))
        l = dict(nx.shortest_path_length(self.grid, 1))
        assert_equal(l[16], 6)
        # now with weights
        l = dict(nx.shortest_path_length(self.cycle, 0, weight='weight'))
        assert_equal(l, {0: 0, 1: 1, 2: 2, 3: 3, 4: 3, 5: 2, 6: 1})
        assert_equal(l, dict(nx.single_source_dijkstra_path_length(
            self.cycle, 0)))
        l = dict(nx.shortest_path_length(self.grid, 1, weight='weight'))
        assert_equal(l[16], 6)

    def test_all_pairs_shortest_path(self):
        p = nx.shortest_path(self.cycle)
        assert_equal(p[0][3], [0, 1, 2, 3])
        assert_equal(p, dict(nx.all_pairs_shortest_path(self.cycle)))
        p = nx.shortest_path(self.grid)
        validate_grid_path(4, 4, 1, 12, p[1][12])
        # now with weights
        p = nx.shortest_path(self.cycle, weight='weight')
        assert_equal(p[0][3], [0, 1, 2, 3])
        assert_equal(p, dict(nx.all_pairs_dijkstra_path(self.cycle)))
        p = nx.shortest_path(self.grid, weight='weight')
        validate_grid_path(4, 4, 1, 12, p[1][12])

    def test_all_pairs_shortest_path_length(self):
        l = dict(nx.shortest_path_length(self.cycle))
        assert_equal(l[0], {0: 0, 1: 1, 2: 2, 3: 3, 4: 3, 5: 2, 6: 1})
        assert_equal(l, dict(nx.all_pairs_shortest_path_length(self.cycle)))
        l = dict(nx.shortest_path_length(self.grid))
        assert_equal(l[1][16], 6)
        # now with weights
        l = dict(nx.shortest_path_length(self.cycle, weight='weight'))
        assert_equal(l[0], {0: 0, 1: 1, 2: 2, 3: 3, 4: 3, 5: 2, 6: 1})
        assert_equal(l, dict(nx.all_pairs_dijkstra_path_length(self.cycle)))
        l = dict(nx.shortest_path_length(self.grid, weight='weight'))
        assert_equal(l[1][16], 6)

    def test_has_path(self):
        G = nx.Graph()
        nx.add_path(G, range(3))
        nx.add_path(G, range(3, 5))
        assert_true(nx.has_path(G, 0, 2))
        assert_false(nx.has_path(G, 0, 4))

    def test_all_shortest_paths(self):
        G = nx.Graph()
        nx.add_path(G, [0, 1, 2, 3])
        nx.add_path(G, [0, 10, 20, 3])
        assert_equal([[0, 1, 2, 3], [0, 10, 20, 3]],
                     sorted(nx.all_shortest_paths(G, 0, 3)))

    @raises(nx.NetworkXNoPath)
    def test_all_shortest_paths_raise(self):
        G = nx.path_graph(4)
        G.add_node(4)
        paths = list(nx.all_shortest_paths(G, 0, 4))


class TestAverageShortestPathLength(object):

    def test_cycle_graph(self):
        l = nx.average_shortest_path_length(nx.cycle_graph(7))
        assert_almost_equal(l, 2)

    def test_path_graph(self):
        l = nx.average_shortest_path_length(nx.path_graph(5))
        assert_almost_equal(l, 2)

    def test_weighted(self):
        G = nx.Graph()
        nx.add_cycle(G, range(7), weight=2)
        l = nx.average_shortest_path_length(G, weight='weight')
        assert_almost_equal(l, 4)
        G = nx.Graph()
        nx.add_path(G, range(5), weight=2)
        l = nx.average_shortest_path_length(G, weight='weight')
        assert_almost_equal(l, 4)

    def test_disconnected(self):
        g = nx.Graph()
        g.add_nodes_from(range(3))
        g.add_edge(0, 1)
        assert_raises(nx.NetworkXError, nx.average_shortest_path_length, g)
        g = g.to_directed()
        assert_raises(nx.NetworkXError, nx.average_shortest_path_length, g)

    def test_trivial_graph(self):
        """Tests that the trivial graph has average path length zero,
        since there is exactly one path of length zero in the trivial
        graph.

        For more information, see issue #1960.

        """
        G = nx.trivial_graph()
        assert_equal(nx.average_shortest_path_length(G), 0)

    @raises(nx.NetworkXPointlessConcept)
    def test_null_graph(self):
        nx.average_shortest_path_length(nx.null_graph())
