#!/usr/bin/env python
from nose.tools import *
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


class TestUnweightedPath:

    def setUp(self):
        from networkx import convert_node_labels_to_integers as cnlti
        self.grid = cnlti(nx.grid_2d_graph(4, 4), first_label=1, ordering="sorted")
        self.cycle = nx.cycle_graph(7)
        self.directed_cycle = nx.cycle_graph(7, create_using=nx.DiGraph())

    def test_bidirectional_shortest_path(self):
        assert_equal(nx.bidirectional_shortest_path(self.cycle, 0, 3),
                     [0, 1, 2, 3])
        assert_equal(nx.bidirectional_shortest_path(self.cycle, 0, 4),
                     [0, 6, 5, 4])
        validate_grid_path(4, 4, 1, 12, nx.bidirectional_shortest_path(self.grid, 1, 12))
        assert_equal(nx.bidirectional_shortest_path(self.directed_cycle, 0, 3),
                     [0, 1, 2, 3])

    def test_shortest_path_length(self):
        assert_equal(nx.shortest_path_length(self.cycle, 0, 3), 3)
        assert_equal(nx.shortest_path_length(self.grid, 1, 12), 5)
        assert_equal(nx.shortest_path_length(self.directed_cycle, 0, 4), 4)
        # now with weights
        assert_equal(nx.shortest_path_length(self.cycle, 0, 3, weight=True), 3)
        assert_equal(nx.shortest_path_length(self.grid, 1, 12, weight=True), 5)
        assert_equal(nx.shortest_path_length(self.directed_cycle, 0, 4, weight=True), 4)

    def test_single_source_shortest_path(self):
        p = nx.single_source_shortest_path(self.directed_cycle, 3)
        assert_equal(p[0], [3, 4, 5, 6, 0])
        p = nx.single_source_shortest_path(self.cycle, 0)
        assert_equal(p[3], [0, 1, 2, 3])
        p = nx.single_source_shortest_path(self.cycle, 0, cutoff=0)
        assert_equal(p, {0: [0]})

    def test_single_source_shortest_path_length(self):
        pl = nx.single_source_shortest_path_length
        lengths = {0: 0, 1: 1, 2: 2, 3: 3, 4: 3, 5: 2, 6: 1}
        assert_equal(dict(pl(self.cycle, 0)), lengths)
        lengths = {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6}
        assert_equal(dict(pl(self.directed_cycle, 0)), lengths)

    def test_single_target_shortest_path(self):
        p = nx.single_target_shortest_path(self.directed_cycle, 0)
        assert_equal(p[3], [3, 4, 5, 6, 0])
        p = nx.single_target_shortest_path(self.cycle, 0)
        assert_equal(p[3], [3, 2, 1, 0])
        p = nx.single_target_shortest_path(self.cycle, 0, cutoff=0)
        assert_equal(p, {0: [0]})

    def test_single_target_shortest_path_length(self):
        pl = nx.single_target_shortest_path_length
        lengths = {0: 0, 1: 1, 2: 2, 3: 3, 4: 3, 5: 2, 6: 1}
        assert_equal(dict(pl(self.cycle, 0)), lengths)
        lengths = {0: 0, 1: 6, 2: 5, 3: 4, 4: 3, 5: 2, 6: 1}
        assert_equal(dict(pl(self.directed_cycle, 0)), lengths)

    def test_all_pairs_shortest_path(self):
        p = dict(nx.all_pairs_shortest_path(self.cycle))
        assert_equal(p[0][3], [0, 1, 2, 3])
        p = dict(nx.all_pairs_shortest_path(self.grid))
        validate_grid_path(4, 4, 1, 12, p[1][12])

    def test_all_pairs_shortest_path_length(self):
        l = dict(nx.all_pairs_shortest_path_length(self.cycle))
        assert_equal(l[0], {0: 0, 1: 1, 2: 2, 3: 3, 4: 3, 5: 2, 6: 1})
        l = dict(nx.all_pairs_shortest_path_length(self.grid))
        assert_equal(l[1][16], 6)

    def test_predecessor_path(self):
        G = nx.path_graph(4)
        assert_equal(nx.predecessor(G, 0), {0: [], 1: [0], 2: [1], 3: [2]})
        assert_equal(nx.predecessor(G, 0, 3), [2])

    def test_predecessor_cycle(self):
        G = nx.cycle_graph(4)
        pred = nx.predecessor(G, 0)
        assert_equal(pred[0], [])
        assert_equal(pred[1], [0])
        assert_true(pred[2] in [[1, 3], [3, 1]])
        assert_equal(pred[3], [0])

    def test_predecessor_cutoff(self):
        G = nx.path_graph(4)
        p = nx.predecessor(G, 0, 3)
        assert_false(4 in p)

    def test_predecessor_target(self):
        G = nx.path_graph(4)
        p = nx.predecessor(G, 0, 3)
        assert_equal(p, [2])
        p = nx.predecessor(G, 0, 3, cutoff=2)
        assert_equal(p, [])
        p, s = nx.predecessor(G, 0, 3, return_seen=True)
        assert_equal(p, [2])
        assert_equal(s, 3)
        p, s = nx.predecessor(G, 0, 3, cutoff=2, return_seen=True)
        assert_equal(p, [])
        assert_equal(s, -1)
