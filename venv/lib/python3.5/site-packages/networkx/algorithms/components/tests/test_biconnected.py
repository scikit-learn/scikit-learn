#!/usr/bin/env python
from nose.tools import *
import networkx as nx
from networkx.algorithms.components import biconnected
from networkx import NetworkXNotImplemented


def assert_components_edges_equal(x, y):
    sx = {frozenset([frozenset(e) for e in c]) for c in x}
    sy = {frozenset([frozenset(e) for e in c]) for c in y}
    assert_equal(sx, sy)


def assert_components_equal(x, y):
    sx = {frozenset(c) for c in x}
    sy = {frozenset(c) for c in y}
    assert_equal(sx, sy)


def test_barbell():
    G = nx.barbell_graph(8, 4)
    nx.add_path(G, [7, 20, 21, 22])
    nx.add_cycle(G, [22, 23, 24, 25])
    pts = set(nx.articulation_points(G))
    assert_equal(pts, {7, 8, 9, 10, 11, 12, 20, 21, 22})

    answer = [
        {12, 13, 14, 15, 16, 17, 18, 19},
        {0, 1, 2, 3, 4, 5, 6, 7},
        {22, 23, 24, 25},
        {11, 12},
        {10, 11},
        {9, 10},
        {8, 9},
        {7, 8},
        {21, 22},
        {20, 21},
        {7, 20},
    ]
    assert_components_equal(list(nx.biconnected_components(G)), answer)

    G.add_edge(2, 17)
    pts = set(nx.articulation_points(G))
    assert_equal(pts, {7, 20, 21, 22})


def test_articulation_points_repetitions():
    G = nx.Graph()
    G.add_edges_from([(0, 1), (1, 2), (1, 3)])
    assert_equal(list(nx.articulation_points(G)), [1])


def test_articulation_points_cycle():
    G = nx.cycle_graph(3)
    nx.add_cycle(G, [1, 3, 4])
    pts = set(nx.articulation_points(G))
    assert_equal(pts, {1})


def test_is_biconnected():
    G = nx.cycle_graph(3)
    assert_true(nx.is_biconnected(G))
    nx.add_cycle(G, [1, 3, 4])
    assert_false(nx.is_biconnected(G))


def test_empty_is_biconnected():
    G = nx.empty_graph(5)
    assert_false(nx.is_biconnected(G))
    G.add_edge(0, 1)
    assert_false(nx.is_biconnected(G))


def test_biconnected_components_cycle():
    G = nx.cycle_graph(3)
    nx.add_cycle(G, [1, 3, 4])
    answer = [{0, 1, 2}, {1, 3, 4}]
    assert_components_equal(list(nx.biconnected_components(G)), answer)

# deprecated


def test_biconnected_component_subgraphs_cycle():
    G = nx.cycle_graph(3)
    nx.add_cycle(G, [1, 3, 4, 5])
    Gc = set(nx.biconnected_component_subgraphs(G))
    assert_equal(len(Gc), 2)
    g1, g2 = Gc
    if 0 in g1:
        assert_true(nx.is_isomorphic(g1, nx.Graph([(0, 1), (0, 2), (1, 2)])))
        assert_true(nx.is_isomorphic(g2, nx.Graph([(1, 3), (1, 5), (3, 4), (4, 5)])))
    else:
        assert_true(nx.is_isomorphic(g1, nx.Graph([(1, 3), (1, 5), (3, 4), (4, 5)])))
        assert_true(nx.is_isomorphic(g2, nx.Graph([(0, 1), (0, 2), (1, 2)])))


def test_biconnected_components1():
    # graph example from
    # http://www.ibluemojo.com/school/articul_algorithm.html
    edges = [
        (0, 1), (0, 5), (0, 6), (0, 14), (1, 5), (1, 6), (1, 14), (2, 4),
        (2, 10), (3, 4), (3, 15), (4, 6), (4, 7), (4, 10), (5, 14), (6, 14),
        (7, 9), (8, 9), (8, 12), (8, 13), (10, 15), (11, 12), (11, 13), (12, 13)
    ]
    G = nx.Graph(edges)
    pts = set(nx.articulation_points(G))
    assert_equal(pts, {4, 6, 7, 8, 9})
    comps = list(nx.biconnected_component_edges(G))
    answer = [
        [(3, 4), (15, 3), (10, 15), (10, 4), (2, 10), (4, 2)],
        [(13, 12), (13, 8), (11, 13), (12, 11), (8, 12)],
        [(9, 8)],
        [(7, 9)],
        [(4, 7)],
        [(6, 4)],
        [(14, 0), (5, 1), (5, 0), (14, 5), (14, 1), (6, 14), (6, 0), (1, 6), (0, 1)],
    ]
    assert_components_edges_equal(comps, answer)


def test_biconnected_components2():
    G = nx.Graph()
    nx.add_cycle(G, 'ABC')
    nx.add_cycle(G, 'CDE')
    nx.add_cycle(G, 'FIJHG')
    nx.add_cycle(G, 'GIJ')
    G.add_edge('E', 'G')
    comps = list(nx.biconnected_component_edges(G))
    answer = [
        [tuple('GF'), tuple('FI'), tuple('IG'), tuple('IJ'),
         tuple('JG'), tuple('JH'), tuple('HG')],
        [tuple('EG')],
        [tuple('CD'), tuple('DE'), tuple('CE')],
        [tuple('AB'), tuple('BC'), tuple('AC')]
    ]
    assert_components_edges_equal(comps, answer)


def test_biconnected_davis():
    D = nx.davis_southern_women_graph()
    bcc = list(nx.biconnected_components(D))[0]
    assert_true(set(D) == bcc)  # All nodes in a giant bicomponent
    # So no articulation points
    assert_equal(len(list(nx.articulation_points(D))), 0)


def test_biconnected_karate():
    K = nx.karate_club_graph()
    answer = [{0, 1, 2, 3, 7, 8, 9, 12, 13, 14, 15, 17, 18, 19,
               20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33},
              {0, 4, 5, 6, 10, 16},
              {0, 11}]
    bcc = list(nx.biconnected_components(K))
    assert_components_equal(bcc, answer)
    assert_equal(set(nx.articulation_points(K)), {0})


def test_biconnected_eppstein():
    # tests from http://www.ics.uci.edu/~eppstein/PADS/Biconnectivity.py
    G1 = nx.Graph({
        0: [1, 2, 5],
        1: [0, 5],
        2: [0, 3, 4],
        3: [2, 4, 5, 6],
        4: [2, 3, 5, 6],
        5: [0, 1, 3, 4],
        6: [3, 4],
    })
    G2 = nx.Graph({
        0: [2, 5],
        1: [3, 8],
        2: [0, 3, 5],
        3: [1, 2, 6, 8],
        4: [7],
        5: [0, 2],
        6: [3, 8],
        7: [4],
        8: [1, 3, 6],
    })
    assert_true(nx.is_biconnected(G1))
    assert_false(nx.is_biconnected(G2))
    answer_G2 = [{1, 3, 6, 8}, {0, 2, 5}, {2, 3}, {4, 7}]
    bcc = list(nx.biconnected_components(G2))
    assert_components_equal(bcc, answer_G2)


def test_null_graph():
    G = nx.Graph()
    assert_false(nx.is_biconnected(G))
    assert_equal(list(nx.biconnected_components(G)), [])
    assert_equal(list(nx.biconnected_component_edges(G)), [])
    assert_equal(list(nx.articulation_points(G)), [])


def test_connected_raise():
    DG = nx.DiGraph()
    assert_raises(NetworkXNotImplemented, nx.biconnected_components, DG)
    assert_raises(NetworkXNotImplemented, nx.biconnected_component_edges, DG)
    assert_raises(NetworkXNotImplemented, nx.articulation_points, DG)
    assert_raises(NetworkXNotImplemented, nx.is_biconnected, DG)
    # deprecated
    assert_raises(NetworkXNotImplemented, nx.biconnected_component_subgraphs, DG)
