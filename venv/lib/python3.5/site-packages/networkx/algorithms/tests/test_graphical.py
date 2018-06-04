#!/usr/bin/env python
from nose.tools import *
from nose import SkipTest
import networkx as nx


def test_valid_degree_sequence1():
    n = 100
    p = .3
    for i in range(10):
        G = nx.erdos_renyi_graph(n, p)
        deg = (d for n, d in G.degree())
        assert_true(nx.is_graphical(deg, method='eg'))
        assert_true(nx.is_graphical(deg, method='hh'))


def test_valid_degree_sequence2():
    n = 100
    for i in range(10):
        G = nx.barabasi_albert_graph(n, 1)
        deg = (d for n, d in G.degree())
        assert_true(nx.is_graphical(deg, method='eg'))
        assert_true(nx.is_graphical(deg, method='hh'))


@raises(nx.NetworkXException)
def test_string_input():
    a = nx.is_graphical([], 'foo')


def test_negative_input():
    assert_false(nx.is_graphical([-1], 'hh'))
    assert_false(nx.is_graphical([-1], 'eg'))
    assert_false(nx.is_graphical([72.5], 'eg'))


class TestAtlas(object):
    @classmethod
    def setupClass(cls):
        global atlas
        import platform
        if platform.python_implementation() == 'Jython':
            raise SkipTest('graph atlas not available under Jython.')
        import networkx.generators.atlas as atlas

    def setUp(self):
        self.GAG = atlas.graph_atlas_g()

    def test_atlas(self):
        for graph in self.GAG:
            deg = (d for n, d in graph.degree())
            assert_true(nx.is_graphical(deg, method='eg'))
            assert_true(nx.is_graphical(deg, method='hh'))


def test_small_graph_true():
    z = [5, 3, 3, 3, 3, 2, 2, 2, 1, 1, 1]
    assert_true(nx.is_graphical(z, method='hh'))
    assert_true(nx.is_graphical(z, method='eg'))
    z = [10, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2]
    assert_true(nx.is_graphical(z, method='hh'))
    assert_true(nx.is_graphical(z, method='eg'))
    z = [1, 1, 1, 1, 1, 2, 2, 2, 3, 4]
    assert_true(nx.is_graphical(z, method='hh'))
    assert_true(nx.is_graphical(z, method='eg'))


def test_small_graph_false():
    z = [1000, 3, 3, 3, 3, 2, 2, 2, 1, 1, 1]
    assert_false(nx.is_graphical(z, method='hh'))
    assert_false(nx.is_graphical(z, method='eg'))
    z = [6, 5, 4, 4, 2, 1, 1, 1]
    assert_false(nx.is_graphical(z, method='hh'))
    assert_false(nx.is_graphical(z, method='eg'))
    z = [1, 1, 1, 1, 1, 1, 2, 2, 2, 3, 4]
    assert_false(nx.is_graphical(z, method='hh'))
    assert_false(nx.is_graphical(z, method='eg'))


def test_directed_degree_sequence():
    # Test a range of valid directed degree sequences
    n, r = 100, 10
    p = 1.0 / r
    for i in range(r):
        G = nx.erdos_renyi_graph(n, p * (i + 1), None, True)
        din = (d for n, d in G.in_degree())
        dout = (d for n, d in G.out_degree())
        assert_true(nx.is_digraphical(din, dout))


def test_small_directed_sequences():
    dout = [5, 3, 3, 3, 3, 2, 2, 2, 1, 1, 1]
    din = [3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 1]
    assert_true(nx.is_digraphical(din, dout))
    # Test nongraphical directed sequence
    dout = [1000, 3, 3, 3, 3, 2, 2, 2, 1, 1, 1]
    din = [103, 102, 102, 102, 102, 102, 102, 102, 102, 102]
    assert_false(nx.is_digraphical(din, dout))
    # Test digraphical small sequence
    dout = [1, 1, 1, 1, 1, 2, 2, 2, 3, 4]
    din = [2, 2, 2, 2, 2, 2, 2, 2, 1, 1]
    assert_true(nx.is_digraphical(din, dout))
    # Test nonmatching sum
    din = [2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1]
    assert_false(nx.is_digraphical(din, dout))
    # Test for negative integer in sequence
    din = [2, 2, 2, -2, 2, 2, 2, 2, 1, 1, 4]
    assert_false(nx.is_digraphical(din, dout))


def test_multi_sequence():
    # Test nongraphical multi sequence
    seq = [1000, 3, 3, 3, 3, 2, 2, 2, 1, 1]
    assert_false(nx.is_multigraphical(seq))
    # Test small graphical multi sequence
    seq = [6, 5, 4, 4, 2, 1, 1, 1]
    assert_true(nx.is_multigraphical(seq))
    # Test for negative integer in sequence
    seq = [6, 5, 4, -4, 2, 1, 1, 1]
    assert_false(nx.is_multigraphical(seq))
    # Test for sequence with odd sum
    seq = [1, 1, 1, 1, 1, 1, 2, 2, 2, 3, 4]
    assert_false(nx.is_multigraphical(seq))


def test_pseudo_sequence():
    # Test small valid pseudo sequence
    seq = [1000, 3, 3, 3, 3, 2, 2, 2, 1, 1]
    assert_true(nx.is_pseudographical(seq))
    # Test for sequence with odd sum
    seq = [1000, 3, 3, 3, 3, 2, 2, 2, 1, 1, 1]
    assert_false(nx.is_pseudographical(seq))
    # Test for negative integer in sequence
    seq = [1000, 3, 3, 3, 3, 2, 2, -2, 1, 1]
    assert_false(nx.is_pseudographical(seq))
