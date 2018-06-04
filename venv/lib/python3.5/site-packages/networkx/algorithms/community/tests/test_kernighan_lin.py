# -*- encoding: utf-8 -*-
# test_kernighan_lin.py - unit tests for Kernighan–Lin bipartition algorithm
#
# Copyright 2011 Ben Edwards <bedwards@cs.unm.edu>.
# Copyright 2011 Aric Hagberg <hagberg@lanl.gov>.
# Copyright 2015 NetworkX developers.
#
# This file is part of NetworkX.
#
# NetworkX is distributed under a BSD license; see LICENSE.txt for more
# information.
"""Unit tests for the :mod:`networkx.algorithms.community.kernighan_lin`
module.

"""
from nose.tools import assert_equal
from nose.tools import raises

import networkx as nx
from networkx.algorithms.community import kernighan_lin_bisection


def assert_partition_equal(x, y):
    assert_equal(set(map(frozenset, x)), set(map(frozenset, y)))


def test_partition():
    G = nx.barbell_graph(3, 0)
    C = kernighan_lin_bisection(G)
    assert_partition_equal(C, [{0, 1, 2}, {3, 4, 5}])


@raises(nx.NetworkXError)
def test_non_disjoint_partition():
    G = nx.barbell_graph(3, 0)
    partition = ({0, 1, 2}, {2, 3, 4, 5})
    kernighan_lin_bisection(G, partition)


@raises(nx.NetworkXError)
def test_too_many_blocks():
    G = nx.barbell_graph(3, 0)
    partition = ({0, 1}, {2}, {3, 4, 5})
    kernighan_lin_bisection(G, partition)


def test_multigraph():
    G = nx.cycle_graph(4)
    M = nx.MultiGraph(G.edges())
    M.add_edges_from(G.edges())
    M.remove_edge(1, 2)
    A, B = kernighan_lin_bisection(M)
    assert_partition_equal([A, B], [{0, 1}, {2, 3}])
