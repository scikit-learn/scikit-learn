#!/usr/bin/env python

from nose.tools import *
from networkx import *
from networkx.algorithms.bipartite.generators import *

"""Generators - Bipartite
----------------------
"""


class TestGeneratorsBipartite():
    def test_complete_bipartite_graph(self):
        G = complete_bipartite_graph(0, 0)
        assert_true(is_isomorphic(G, null_graph()))

        for i in [1, 5]:
            G = complete_bipartite_graph(i, 0)
            assert_true(is_isomorphic(G, empty_graph(i)))
            G = complete_bipartite_graph(0, i)
            assert_true(is_isomorphic(G, empty_graph(i)))

        G = complete_bipartite_graph(2, 2)
        assert_true(is_isomorphic(G, cycle_graph(4)))

        G = complete_bipartite_graph(1, 5)
        assert_true(is_isomorphic(G, star_graph(5)))

        G = complete_bipartite_graph(5, 1)
        assert_true(is_isomorphic(G, star_graph(5)))

        # complete_bipartite_graph(m1,m2) is a connected graph with
        # m1+m2 nodes and  m1*m2 edges
        for m1, m2 in [(5, 11), (7, 3)]:
            G = complete_bipartite_graph(m1, m2)
            assert_equal(number_of_nodes(G), m1 + m2)
            assert_equal(number_of_edges(G), m1 * m2)

        assert_raises(networkx.exception.NetworkXError,
                      complete_bipartite_graph, 7, 3, create_using=DiGraph())

        mG = complete_bipartite_graph(7, 3, create_using=MultiGraph())
        assert_equal(sorted(mG.edges()), sorted(G.edges()))

        # specify nodes rather than number of nodes
        G = complete_bipartite_graph([1, 2], ['a', 'b'])
        has_edges = G.has_edge(1, 'a') & G.has_edge(1, 'b') &\
            G.has_edge(2, 'a') & G.has_edge(2, 'b')
        assert_true(has_edges)
        assert_equal(G.size(), 4)

    def test_configuration_model(self):
        aseq = [3, 3, 3, 3]
        bseq = [2, 2, 2, 2, 2]
        assert_raises(networkx.exception.NetworkXError,
                      configuration_model, aseq, bseq)

        aseq = [3, 3, 3, 3]
        bseq = [2, 2, 2, 2, 2, 2]
        G = configuration_model(aseq, bseq)
        assert_equal(sorted(d for n, d in G.degree()),
                     [2, 2, 2, 2, 2, 2, 3, 3, 3, 3])

        aseq = [2, 2, 2, 2, 2, 2]
        bseq = [3, 3, 3, 3]
        G = configuration_model(aseq, bseq)
        assert_equal(sorted(d for n, d in G.degree()),
                     [2, 2, 2, 2, 2, 2, 3, 3, 3, 3])

        aseq = [2, 2, 2, 1, 1, 1]
        bseq = [3, 3, 3]
        G = configuration_model(aseq, bseq)
        assert_equal(sorted(d for n, d in G.degree()),
                     [1, 1, 1, 2, 2, 2, 3, 3, 3])

        GU = project(Graph(G), range(len(aseq)))
        assert_equal(GU.number_of_nodes(), 6)

        GD = project(Graph(G), range(len(aseq), len(aseq) + len(bseq)))
        assert_equal(GD.number_of_nodes(), 3)

        assert_raises(networkx.exception.NetworkXError,
                      configuration_model, aseq, bseq,
                      create_using=DiGraph())

    def test_havel_hakimi_graph(self):
        aseq = [3, 3, 3, 3]
        bseq = [2, 2, 2, 2, 2]
        assert_raises(networkx.exception.NetworkXError,
                      havel_hakimi_graph, aseq, bseq)

        bseq = [2, 2, 2, 2, 2, 2]
        G = havel_hakimi_graph(aseq, bseq)
        assert_equal(sorted(d for n, d in G.degree()),
                     [2, 2, 2, 2, 2, 2, 3, 3, 3, 3])

        aseq = [2, 2, 2, 2, 2, 2]
        bseq = [3, 3, 3, 3]
        G = havel_hakimi_graph(aseq, bseq)
        assert_equal(sorted(d for n, d in G.degree()),
                     [2, 2, 2, 2, 2, 2, 3, 3, 3, 3])

        GU = project(Graph(G), range(len(aseq)))
        assert_equal(GU.number_of_nodes(), 6)

        GD = project(Graph(G), range(len(aseq), len(aseq) + len(bseq)))
        assert_equal(GD.number_of_nodes(), 4)
        assert_raises(networkx.exception.NetworkXError,
                      havel_hakimi_graph, aseq, bseq,
                      create_using=DiGraph())

    def test_reverse_havel_hakimi_graph(self):
        aseq = [3, 3, 3, 3]
        bseq = [2, 2, 2, 2, 2]
        assert_raises(networkx.exception.NetworkXError,
                      reverse_havel_hakimi_graph, aseq, bseq)

        bseq = [2, 2, 2, 2, 2, 2]
        G = reverse_havel_hakimi_graph(aseq, bseq)
        assert_equal(sorted(d for n, d in G.degree()),
                     [2, 2, 2, 2, 2, 2, 3, 3, 3, 3])

        aseq = [2, 2, 2, 2, 2, 2]
        bseq = [3, 3, 3, 3]
        G = reverse_havel_hakimi_graph(aseq, bseq)
        assert_equal(sorted(d for n, d in G.degree()),
                     [2, 2, 2, 2, 2, 2, 3, 3, 3, 3])

        aseq = [2, 2, 2, 1, 1, 1]
        bseq = [3, 3, 3]
        G = reverse_havel_hakimi_graph(aseq, bseq)
        assert_equal(sorted(d for n, d in G.degree()),
                     [1, 1, 1, 2, 2, 2, 3, 3, 3])

        GU = project(Graph(G), range(len(aseq)))
        assert_equal(GU.number_of_nodes(), 6)

        GD = project(Graph(G), range(len(aseq), len(aseq) + len(bseq)))
        assert_equal(GD.number_of_nodes(), 3)
        assert_raises(networkx.exception.NetworkXError,
                      reverse_havel_hakimi_graph, aseq, bseq,
                      create_using=DiGraph())

    def test_alternating_havel_hakimi_graph(self):
        aseq = [3, 3, 3, 3]
        bseq = [2, 2, 2, 2, 2]
        assert_raises(networkx.exception.NetworkXError,
                      alternating_havel_hakimi_graph, aseq, bseq)

        bseq = [2, 2, 2, 2, 2, 2]
        G = alternating_havel_hakimi_graph(aseq, bseq)
        assert_equal(sorted(d for n, d in G.degree()),
                     [2, 2, 2, 2, 2, 2, 3, 3, 3, 3])

        aseq = [2, 2, 2, 2, 2, 2]
        bseq = [3, 3, 3, 3]
        G = alternating_havel_hakimi_graph(aseq, bseq)
        assert_equal(sorted(d for n, d in G.degree()),
                     [2, 2, 2, 2, 2, 2, 3, 3, 3, 3])

        aseq = [2, 2, 2, 1, 1, 1]
        bseq = [3, 3, 3]
        G = alternating_havel_hakimi_graph(aseq, bseq)
        assert_equal(sorted(d for n, d in G.degree()),
                     [1, 1, 1, 2, 2, 2, 3, 3, 3])

        GU = project(Graph(G), range(len(aseq)))
        assert_equal(GU.number_of_nodes(), 6)

        GD = project(Graph(G), range(len(aseq), len(aseq) + len(bseq)))
        assert_equal(GD.number_of_nodes(), 3)

        assert_raises(networkx.exception.NetworkXError,
                      alternating_havel_hakimi_graph, aseq, bseq,
                      create_using=DiGraph())

    def test_preferential_attachment(self):
        aseq = [3, 2, 1, 1]
        G = preferential_attachment_graph(aseq, 0.5)
        assert_raises(networkx.exception.NetworkXError,
                      preferential_attachment_graph, aseq, 0.5,
                      create_using=DiGraph())

    def test_random_graph(self):
        n = 10
        m = 20
        G = random_graph(n, m, 0.9)
        assert_equal(len(G), 30)
        assert_true(is_bipartite(G))
        X, Y = nx.algorithms.bipartite.sets(G)
        assert_equal(set(range(n)), X)
        assert_equal(set(range(n, n + m)), Y)

    def test_random_graph(self):
        n = 10
        m = 20
        G = random_graph(n, m, 0.9, directed=True)
        assert_equal(len(G), 30)
        assert_true(is_bipartite(G))
        X, Y = nx.algorithms.bipartite.sets(G)
        assert_equal(set(range(n)), X)
        assert_equal(set(range(n, n + m)), Y)

    def test_gnmk_random_graph(self):
        n = 10
        m = 20
        edges = 200
        G = gnmk_random_graph(n, m, edges)
        assert_equal(len(G), 30)
        assert_true(is_bipartite(G))
        X, Y = nx.algorithms.bipartite.sets(G)
        print(X)
        assert_equal(set(range(n)), X)
        assert_equal(set(range(n, n + m)), Y)
        assert_equal(edges, len(list(G.edges())))
