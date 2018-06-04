# -*- encoding: utf-8 -*-
# test_duplication.py - unit tests for the generators.duplication module
#
# Copyright 2010-2018 NetworkX developers.
#
# This file is part of NetworkX.
#
# NetworkX is distributed under a BSD license; see LICENSE.txt for more
# information.
"""Unit tests for the :mod:`networkx.generators.duplication` module.

"""
from nose.tools import assert_equal
from nose.tools import assert_raises
from nose.tools import raises

from networkx.exception import NetworkXError
from networkx.generators.duplication import duplication_divergence_graph
from networkx.generators.duplication import partial_duplication_graph


class TestDuplicationDivergenceGraph(object):
    """Unit tests for the
    :func:`networkx.generators.duplication.duplication_divergence_graph`
    function.

    """

    def test_final_size(self):
        G = duplication_divergence_graph(3, 1)
        assert_equal(len(G), 3)

    @raises(NetworkXError)
    def test_probability_too_large(self):
        duplication_divergence_graph(3, 2)

    @raises(NetworkXError)
    def test_probability_too_small(self):
        duplication_divergence_graph(3, -1)


class TestPartialDuplicationGraph(object):
    """Unit tests for the
    :func:`networkx.generators.duplication.partial_duplication_graph`
    function.

    """

    def test_final_size(self):
        N = 10
        n = 5
        p = 0.5
        q = 0.5
        G = partial_duplication_graph(N, n, p, q)
        assert_equal(len(G), N)

    def test_initial_clique_size(self):
        N = 10
        n = 10
        p = 0.5
        q = 0.5
        G = partial_duplication_graph(N, n, p, q)
        assert_equal(len(G), n)

    @raises(NetworkXError)
    def test_invalid_initial_size(self):
        N = 5
        n = 10
        p = 0.5
        q = 0.5
        G = partial_duplication_graph(N, n, p, q)
        assert_equal(len(G), n)

    def test_invalid_probabilities(self):
        N = 1
        n = 1
        for p, q in [(0.5, 2), (0.5, -1), (2, 0.5), (-1, 0.5)]:
            args = (N, n, p, q)
            assert_raises(NetworkXError, partial_duplication_graph, *args)
