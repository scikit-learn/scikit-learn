#!/usr/bin/env python
from nose.tools import *
import networkx as nx


def test_hierarchy_exception():
    G = nx.cycle_graph(5)
    assert_raises(nx.NetworkXError, nx.flow_hierarchy, G)


def test_hierarchy_cycle():
    G = nx.cycle_graph(5, create_using=nx.DiGraph())
    assert_equal(nx.flow_hierarchy(G), 0.0)


def test_hierarchy_tree():
    G = nx.full_rary_tree(2, 16, create_using=nx.DiGraph())
    assert_equal(nx.flow_hierarchy(G), 1.0)


def test_hierarchy_1():
    G = nx.DiGraph()
    G.add_edges_from([(0, 1), (1, 2), (2, 3), (3, 1), (3, 4), (0, 4)])
    assert_equal(nx.flow_hierarchy(G), 0.5)


def test_hierarchy_weight():
    G = nx.DiGraph()
    G.add_edges_from([(0, 1, {'weight': .3}),
                      (1, 2, {'weight': .1}),
                      (2, 3, {'weight': .1}),
                      (3, 1, {'weight': .1}),
                      (3, 4, {'weight': .3}),
                      (0, 4, {'weight': .3})])
    assert_equal(nx.flow_hierarchy(G, weight='weight'), .75)
