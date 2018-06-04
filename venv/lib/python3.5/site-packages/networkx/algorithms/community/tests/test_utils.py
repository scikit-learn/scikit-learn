# test_utils.py - unit tests for the community utils module
#
# Copyright 2016 NetworkX developers.
#
# This file is part of NetworkX.
#
# NetworkX is distributed under a BSD license; see LICENSE.txt for more
# information.
"""Unit tests for the :mod:`networkx.algorithms.community.utils` module.

"""
from nose.tools import assert_false
from nose.tools import assert_true

import networkx as nx
from networkx.algorithms.community import is_partition


def test_is_partition():
    G = nx.empty_graph(3)
    assert_true(is_partition(G, [{0, 1}, {2}]))


def test_not_covering():
    G = nx.empty_graph(3)
    assert_false(is_partition(G, [{0}, {1}]))


def test_not_disjoint():
    G = nx.empty_graph(3)
    assert_false(is_partition(G, [{0, 1}, {1, 2}]))
