from itertools import groupby

from nose.tools import assert_equal
from nose.tools import assert_less_equal
from nose.tools import raises

import networkx as nx
from networkx.testing import *
from networkx import graph_atlas
from networkx import graph_atlas_g
from networkx.generators.atlas import NUM_GRAPHS
from networkx.utils import pairwise


class TestAtlasGraph(object):
    """Unit tests for the :func:`~networkx.graph_atlas` function."""

    @raises(ValueError)
    def test_index_too_small(self):
        graph_atlas(-1)

    @raises(ValueError)
    def test_index_too_large(self):
        graph_atlas(NUM_GRAPHS)

    def test_graph(self):
        G = graph_atlas(6)
        assert_nodes_equal(G.nodes(), range(3))
        assert_edges_equal(G.edges(), [(0, 1), (0, 2)])


class TestAtlasGraphG(object):
    """Unit tests for the :func:`~networkx.graph_atlas_g` function."""

    def setUp(self):
        self.GAG = graph_atlas_g()

    def test_sizes(self):
        G = self.GAG[0]
        assert_equal(G.number_of_nodes(), 0)
        assert_equal(G.number_of_edges(), 0)

        G = self.GAG[7]
        assert_equal(G.number_of_nodes(), 3)
        assert_equal(G.number_of_edges(), 3)

    def test_names(self):
        for i, G in enumerate(self.GAG):
            assert_equal(int(G.name[1:]), i)

    def test_nondecreasing_nodes(self):
        # check for nondecreasing number of nodes
        for n1, n2 in pairwise(map(len, self.GAG)):
            assert_less_equal(n2, n1 + 1)

    def test_nondecreasing_edges(self):
        # check for nondecreasing number of edges (for fixed number of
        # nodes)
        for n, group in groupby(self.GAG, key=nx.number_of_nodes):
            for m1, m2 in pairwise(map(nx.number_of_edges, group)):
                assert_less_equal(m2, m1 + 1)

    def test_nondecreasing_degree_sequence(self):
        # Check for lexicographically nondecreasing degree sequences
        # (for fixed number of nodes and edges).
        #
        # There are three exceptions to this rule in the order given in
        # the "Atlas of Graphs" book, so we need to manually exclude
        # those.
        exceptions = [('G55', 'G56'), ('G1007', 'G1008'), ('G1012', 'G1013')]
        for n, group in groupby(self.GAG, key=nx.number_of_nodes):
            for m, group in groupby(group, key=nx.number_of_edges):
                for G1, G2 in pairwise(group):
                    if (G1.name, G2.name) in exceptions:
                        continue
                    d1 = sorted(d for v, d in G1.degree())
                    d2 = sorted(d for v, d in G2.degree())
                    assert_less_equal(d1, d2)
