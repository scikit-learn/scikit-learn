# test_generators.py - unit tests for the community.generators module
#
# Copyright 2011 Ben Edwards <bedwards@cs.unm.edu>.
# Copyright 2011 Aric Hagberg <hagberg@lanl.gov>
# Copyright 2015 NetworkX developers.
#
# This file is part of NetworkX.
#
# NetworkX is distributed under a BSD license; see LICENSE.txt for more
# information.
"""Unit tests for the :mod:`networkx.algorithms.community.generators`
module.

"""
from nose.tools import assert_equal
from nose.tools import assert_true
from nose.tools import raises

import networkx as nx
from networkx.algorithms.community import LFR_benchmark_graph
from networkx.algorithms.community.community_utils import is_partition


def test_generator():
    n = 250
    tau1 = 3
    tau2 = 1.5
    mu = 0.1
    G = LFR_benchmark_graph(n, tau1, tau2, mu, average_degree=5,
                            min_community=20, seed=10)
    assert_equal(len(G), 250)
    C = {frozenset(G.nodes[v]['community']) for v in G}
    assert_true(is_partition(G.nodes(), C))
    # assert_equal([len(c) for c in C], [53, 12, 10, 15, 10])
    # assert_equal(len(G.edges()), 157)


@raises(nx.NetworkXError)
def test_invalid_tau1():
    n = 100
    tau1 = 2
    tau2 = 1
    mu = 0.1
    LFR_benchmark_graph(n, tau1, tau2, mu, min_degree=2)


@raises(nx.NetworkXError)
def test_invalid_tau2():
    n = 100
    tau1 = 1
    tau2 = 2
    mu = 0.1
    LFR_benchmark_graph(n, tau1, tau2, mu, min_degree=2)


@raises(nx.NetworkXError)
def test_mu_too_large():
    n = 100
    tau1 = 2
    tau2 = 2
    mu = 1.1
    LFR_benchmark_graph(n, tau1, tau2, mu, min_degree=2)


@raises(nx.NetworkXError)
def test_mu_too_small():
    n = 100
    tau1 = 2
    tau2 = 2
    mu = -1
    LFR_benchmark_graph(n, tau1, tau2, mu, min_degree=2)


@raises(nx.NetworkXError)
def test_both_degrees_none():
    n = 100
    tau1 = 2
    tau2 = 2
    mu = -1
    LFR_benchmark_graph(n, tau1, tau2, mu)


@raises(nx.NetworkXError)
def test_neither_degrees_none():
    n = 100
    tau1 = 2
    tau2 = 2
    mu = -1
    LFR_benchmark_graph(n, tau1, tau2, mu, min_degree=2, average_degree=5)
