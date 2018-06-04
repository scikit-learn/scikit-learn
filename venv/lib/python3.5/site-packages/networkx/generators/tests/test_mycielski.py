# test_mycielski.py - unit tests for the mycielski module
#
# Copyright 2010, 2011, 2012, 2013, 2014, 2015 NetworkX developers.
#
# This file is part of NetworkX.
#
# NetworkX is distributed under a BSD license; see LICENSE.txt for more
# information.

"""Unit tests for the :mod:`networkx.generators.mycielski` module."""

from nose.tools import assert_true, assert_equal, raises
import networkx as nx
from networkx import *


class TestMycielski(object):

    def test_construction(self):
        G = nx.path_graph(2)
        M = mycielskian(G)
        assert_true(is_isomorphic(M, cycle_graph(5)))

    def test_size(self):
        G = nx.path_graph(2)
        M = mycielskian(G, 2)
        assert_equal(len(M), 11)
        assert_equal(M.size(), 20)

    def test_mycielski_graph_generator(self):
        G = mycielski_graph(1)
        assert_true(is_isomorphic(G, nx.empty_graph(1)))
        G = mycielski_graph(2)
        assert_true(is_isomorphic(G, nx.path_graph(2)))
        G = mycielski_graph(3)
        assert_true(is_isomorphic(G, cycle_graph(5)))
        G = mycielski_graph(4)
        assert_true(is_isomorphic(G, mycielskian(cycle_graph(5))))
