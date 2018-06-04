# -*- encoding: utf-8 -*-
# test_random_graphs.py - unit tests for random graph generators
#
# Copyright 2010-2018 NetworkX developers.
#
# This file is part of NetworkX.
#
# NetworkX is distributed under a BSD license; see LICENSE.txt for more
# information.
"""Unit tests for the :mod:`networkx.generators.random_graphs` module.

"""
from nose.tools import assert_almost_equal
from nose.tools import assert_greater
from nose.tools import assert_less
from nose.tools import assert_equal
from nose.tools import assert_raises
from nose.tools import assert_true

from networkx.exception import NetworkXError
from networkx.generators.random_graphs import barabasi_albert_graph
from networkx.generators.random_graphs import extended_barabasi_albert_graph
from networkx.generators.random_graphs import binomial_graph
from networkx.generators.random_graphs import connected_watts_strogatz_graph
from networkx.generators.random_graphs import dense_gnm_random_graph
from networkx.generators.random_graphs import erdos_renyi_graph
from networkx.generators.random_graphs import fast_gnp_random_graph
from networkx.generators.random_graphs import gnm_random_graph
from networkx.generators.random_graphs import gnp_random_graph
from networkx.generators.random_graphs import newman_watts_strogatz_graph
from networkx.generators.random_graphs import powerlaw_cluster_graph
from networkx.generators.random_graphs import random_kernel_graph
from networkx.generators.random_graphs import random_lobster
from networkx.generators.random_graphs import random_regular_graph
from networkx.generators.random_graphs import random_shell_graph
from networkx.generators.random_graphs import watts_strogatz_graph


class TestGeneratorsRandom(object):

    def smoke_test_random_graph(self):
        seed = 42
        G = gnp_random_graph(100, 0.25, seed)
        G = binomial_graph(100, 0.25, seed)
        G = erdos_renyi_graph(100, 0.25, seed)
        G = fast_gnp_random_graph(100, 0.25, seed)
        G = gnm_random_graph(100, 20, seed)
        G = dense_gnm_random_graph(100, 20, seed)

        G = watts_strogatz_graph(10, 2, 0.25, seed)
        assert_equal(len(G), 10)
        assert_equal(G.number_of_edges(), 10)

        G = connected_watts_strogatz_graph(10, 2, 0.1, seed)
        assert_equal(len(G), 10)
        assert_equal(G.number_of_edges(), 10)

        G = watts_strogatz_graph(10, 4, 0.25, seed)
        assert_equal(len(G), 10)
        assert_equal(G.number_of_edges(), 20)

        G = newman_watts_strogatz_graph(10, 2, 0.0, seed)
        assert_equal(len(G), 10)
        assert_equal(G.number_of_edges(), 10)

        G = newman_watts_strogatz_graph(10, 4, 0.25, seed)
        assert_equal(len(G), 10)
        assert_true(G.number_of_edges() >= 20)

        G = barabasi_albert_graph(100, 1, seed)
        G = barabasi_albert_graph(100, 3, seed)
        assert_equal(G.number_of_edges(), (97 * 3))

        G = extended_barabasi_albert_graph(100, 1, 0, 0, seed)
        assert_equal(G.number_of_edges(), 99)
        G = extended_barabasi_albert_graph(100, 3, 0, 0, seed)
        assert_equal(G.number_of_edges(), 97 * 3)
        G = extended_barabasi_albert_graph(100, 1, 0, 0.5, seed)
        assert_equal(G.number_of_edges(), 99)
        G = extended_barabasi_albert_graph(100, 2, 0.5, 0, seed)
        assert_greater(G.number_of_edges(), 100 * 3)
        assert_less(G.number_of_edges(), 100 * 4)

        G = extended_barabasi_albert_graph(100, 2, 0.3, 0.3, seed)
        assert_greater(G.number_of_edges(), 100 * 2)
        assert_less(G.number_of_edges(), 100 * 4)

        G = powerlaw_cluster_graph(100, 1, 1.0, seed)
        G = powerlaw_cluster_graph(100, 3, 0.0, seed)
        assert_equal(G.number_of_edges(), (97 * 3))

        G = random_regular_graph(10, 20, seed)

        assert_raises(NetworkXError, random_regular_graph, 3, 21)

        constructor = [(10, 20, 0.8), (20, 40, 0.8)]
        G = random_shell_graph(constructor, seed)

        G = random_lobster(10, 0.1, 0.5, seed)

    def test_extended_barabasi_albert(self, m=2):
        """
        Tests that the extended BA random graph generated behaves consistenly.

        Tests the exceptions are raised as expected.

        The graphs generation are repeated several times to prevent lucky-shots

        """
        seed = 42
        repeats = 2
        BA_model = barabasi_albert_graph(100, m, seed)
        BA_model_edges = BA_model.number_of_edges()

        while repeats:
            repeats -= 1

            # This behaves just like BA, the number of edges must be the same
            G1 = extended_barabasi_albert_graph(100, m, 0, 0, seed)
            assert_equal(G1.size(), BA_model_edges)

            # More than twice more edges should have been added
            G1 = extended_barabasi_albert_graph(100, m, 0.8, 0, seed)
            assert_greater(G1.size(), BA_model_edges * 2)

            # Only edge rewiring, so the number of edges less than original
            G2 = extended_barabasi_albert_graph(100, m, 0, 0.8, seed)
            assert_equal(G2.size(), BA_model_edges)

            # Mixed scenario: less edges than G1 and more edges than G2
            G3 = extended_barabasi_albert_graph(100, m, 0.3, 0.3, seed)
            assert_greater(G3.size(), G2.size())
            assert_less(G3.size(), G1.size())

        # Testing exceptions
        ebag = extended_barabasi_albert_graph
        assert_raises(NetworkXError, ebag, m, m, 0, 0)
        assert_raises(NetworkXError, ebag, 1, 0.5, 0, 0)
        assert_raises(NetworkXError, ebag, 100, 2, 0.5, 0.5)

    def test_random_zero_regular_graph(self):
        """Tests that a 0-regular graph has the correct number of nodes and
        edges.

        """
        seed = 42
        G = random_regular_graph(0, 10)
        assert_equal(len(G), 10)
        assert_equal(sum(1 for _ in G.edges()), 0)

    def test_gnp(self):
        for generator in [gnp_random_graph, binomial_graph, erdos_renyi_graph,
                          fast_gnp_random_graph]:
            G = generator(10, -1.1)
            assert_equal(len(G), 10)
            assert_equal(sum(1 for _ in G.edges()), 0)

            G = generator(10, 0.1)
            assert_equal(len(G), 10)

            G = generator(10, 0.1, seed=42)
            assert_equal(len(G), 10)

            G = generator(10, 1.1)
            assert_equal(len(G), 10)
            assert_equal(sum(1 for _ in G.edges()), 45)

            G = generator(10, -1.1, directed=True)
            assert_true(G.is_directed())
            assert_equal(len(G), 10)
            assert_equal(sum(1 for _ in G.edges()), 0)

            G = generator(10, 0.1, directed=True)
            assert_true(G.is_directed())
            assert_equal(len(G), 10)

            G = generator(10, 1.1, directed=True)
            assert_true(G.is_directed())
            assert_equal(len(G), 10)
            assert_equal(sum(1 for _ in G.edges()), 90)

            # assert that random graphs generate all edges for p close to 1
            edges = 0
            runs = 100
            for i in range(runs):
                edges += sum(1 for _ in generator(10, 0.99999, directed=True).edges())
            assert_almost_equal(edges / float(runs), 90, delta=runs * 2.0 / 100)

    def test_gnm(self):
        G = gnm_random_graph(10, 3)
        assert_equal(len(G), 10)
        assert_equal(sum(1 for _ in G.edges()), 3)

        G = gnm_random_graph(10, 3, seed=42)
        assert_equal(len(G), 10)
        assert_equal(sum(1 for _ in G.edges()), 3)

        G = gnm_random_graph(10, 100)
        assert_equal(len(G), 10)
        assert_equal(sum(1 for _ in G.edges()), 45)

        G = gnm_random_graph(10, 100, directed=True)
        assert_equal(len(G), 10)
        assert_equal(sum(1 for _ in G.edges()), 90)

        G = gnm_random_graph(10, -1.1)
        assert_equal(len(G), 10)
        assert_equal(sum(1 for _ in G.edges()), 0)

    def test_watts_strogatz_big_k(self):
        assert_raises(NetworkXError, watts_strogatz_graph, 10, 10, 0.25)
        assert_raises(NetworkXError, newman_watts_strogatz_graph, 10, 10, 0.25)
        # could create an infinite loop, now doesn't
        # infinite loop used to occur when a node has degree n-1 and needs to rewire
        watts_strogatz_graph(10, 9, 0.25, seed=0)
        newman_watts_strogatz_graph(10, 9, 0.5, seed=0)

    def test_random_kernel_graph(self):
        def integral(u, w, z):
            return c * (z - w)

        def root(u, w, r):
            return r / c + w
        c = 1
        graph = random_kernel_graph(1000, integral, root)
        assert_equal(len(graph), 1000)
