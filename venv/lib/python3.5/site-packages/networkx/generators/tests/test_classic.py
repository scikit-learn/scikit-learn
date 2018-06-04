#!/usr/bin/env python
"""
====================
Generators - Classic
====================

Unit tests for various classic graph generators in generators/classic.py
"""
import itertools

from nose.tools import *
import networkx as nx
from networkx import *
from networkx.algorithms.isomorphism.isomorph import graph_could_be_isomorphic
from networkx.testing import assert_edges_equal
from networkx.testing import assert_nodes_equal

is_isomorphic = graph_could_be_isomorphic


class TestGeneratorClassic():
    def test_balanced_tree(self):
        # balanced_tree(r,h) is a tree with (r**(h+1)-1)/(r-1) edges
        for r, h in [(2, 2), (3, 3), (6, 2)]:
            t = balanced_tree(r, h)
            order = t.order()
            assert_true(order == (r**(h + 1) - 1) / (r - 1))
            assert_true(is_connected(t))
            assert_true(t.size() == order - 1)
            dh = degree_histogram(t)
            assert_equal(dh[0], 0)  # no nodes of 0
            assert_equal(dh[1], r**h)  # nodes of degree 1 are leaves
            assert_equal(dh[r], 1)  # root is degree r
            assert_equal(dh[r + 1], order - r**h - 1)  # everyone else is degree r+1
            assert_equal(len(dh), r + 2)

    def test_balanced_tree_star(self):
        # balanced_tree(r,1) is the r-star
        t = balanced_tree(r=2, h=1)
        assert_true(is_isomorphic(t, star_graph(2)))
        t = balanced_tree(r=5, h=1)
        assert_true(is_isomorphic(t, star_graph(5)))
        t = balanced_tree(r=10, h=1)
        assert_true(is_isomorphic(t, star_graph(10)))

    def test_balanced_tree_path(self):
        """Tests that the balanced tree with branching factor one is the
        path graph.

        """
        # A tree of height four has five levels.
        T = balanced_tree(1, 4)
        P = path_graph(5)
        assert_true(is_isomorphic(T, P))

    def test_full_rary_tree(self):
        r = 2
        n = 9
        t = full_rary_tree(r, n)
        assert_equal(t.order(), n)
        assert_true(is_connected(t))
        dh = degree_histogram(t)
        assert_equal(dh[0], 0)  # no nodes of 0
        assert_equal(dh[1], 5)  # nodes of degree 1 are leaves
        assert_equal(dh[r], 1)  # root is degree r
        assert_equal(dh[r + 1], 9 - 5 - 1)  # everyone else is degree r+1
        assert_equal(len(dh), r + 2)

    def test_full_rary_tree_balanced(self):
        t = full_rary_tree(2, 15)
        th = balanced_tree(2, 3)
        assert_true(is_isomorphic(t, th))

    def test_full_rary_tree_path(self):
        t = full_rary_tree(1, 10)
        assert_true(is_isomorphic(t, path_graph(10)))

    def test_full_rary_tree_empty(self):
        t = full_rary_tree(0, 10)
        assert_true(is_isomorphic(t, empty_graph(10)))
        t = full_rary_tree(3, 0)
        assert_true(is_isomorphic(t, empty_graph(0)))

    def test_full_rary_tree_3_20(self):
        t = full_rary_tree(3, 20)
        assert_equal(t.order(), 20)

    def test_barbell_graph(self):
        # number of nodes = 2*m1 + m2 (2 m1-complete graphs + m2-path + 2 edges)
        # number of edges = 2*(number_of_edges(m1-complete graph) + m2 + 1
        m1 = 3
        m2 = 5
        b = barbell_graph(m1, m2)
        assert_true(number_of_nodes(b) == 2 * m1 + m2)
        assert_true(number_of_edges(b) == m1 * (m1 - 1) + m2 + 1)

        m1 = 4
        m2 = 10
        b = barbell_graph(m1, m2)
        assert_true(number_of_nodes(b) == 2 * m1 + m2)
        assert_true(number_of_edges(b) == m1 * (m1 - 1) + m2 + 1)

        m1 = 3
        m2 = 20
        b = barbell_graph(m1, m2)
        assert_true(number_of_nodes(b) == 2 * m1 + m2)
        assert_true(number_of_edges(b) == m1 * (m1 - 1) + m2 + 1)

        # Raise NetworkXError if m1<2
        m1 = 1
        m2 = 20
        assert_raises(networkx.exception.NetworkXError, barbell_graph, m1, m2)

        # Raise NetworkXError if m2<0
        m1 = 5
        m2 = -2
        assert_raises(networkx.exception.NetworkXError, barbell_graph, m1, m2)

        # barbell_graph(2,m) = path_graph(m+4)
        m1 = 2
        m2 = 5
        b = barbell_graph(m1, m2)
        assert_true(is_isomorphic(b, path_graph(m2 + 4)))

        m1 = 2
        m2 = 10
        b = barbell_graph(m1, m2)
        assert_true(is_isomorphic(b, path_graph(m2 + 4)))

        m1 = 2
        m2 = 20
        b = barbell_graph(m1, m2)
        assert_true(is_isomorphic(b, path_graph(m2 + 4)))

        assert_raises(networkx.exception.NetworkXError, barbell_graph, m1, m2,
                      create_using=DiGraph())

        mb = barbell_graph(m1, m2, create_using=MultiGraph())
        assert_edges_equal(mb.edges(), b.edges())

    def test_complete_graph(self):
        # complete_graph(m) is a connected graph with
        # m nodes and  m*(m+1)/2 edges
        for m in [0, 1, 3, 5]:
            g = complete_graph(m)
            assert_true(number_of_nodes(g) == m)
            assert_true(number_of_edges(g) == m * (m - 1) // 2)

        mg = complete_graph(m, create_using=MultiGraph())
        assert_edges_equal(mg.edges(), g.edges())

        g = complete_graph("abc")
        assert_nodes_equal(g.nodes(), ['a', 'b', 'c'])
        assert_equal(g.size(), 3)

    def test_complete_digraph(self):
        # complete_graph(m) is a connected graph with
        # m nodes and  m*(m+1)/2 edges
        for m in [0, 1, 3, 5]:
            g = complete_graph(m, create_using=nx.DiGraph())
            assert_true(number_of_nodes(g) == m)
            assert_true(number_of_edges(g) == m * (m - 1))

        g = complete_graph("abc", create_using=nx.DiGraph())
        assert_equal(len(g), 3)
        assert_equal(g.size(), 6)
        assert_true(g.is_directed())

    def test_circular_ladder_graph(self):
        G = circular_ladder_graph(5)
        assert_raises(networkx.exception.NetworkXError, circular_ladder_graph,
                      5, create_using=DiGraph())
        mG = circular_ladder_graph(5, create_using=MultiGraph())
        assert_edges_equal(mG.edges(), G.edges())

    def test_circulant_graph(self):
        # Ci_n(1) is the cycle graph for all n
        Ci6_1 = circulant_graph(6, [1])
        C6 = cycle_graph(6)
        assert_edges_equal(Ci6_1.edges(), C6.edges())

        # Ci_n(1, 2, ..., n div 2) is the complete graph for all n
        Ci7 = circulant_graph(7, [1, 2, 3])
        K7 = complete_graph(7)
        assert_edges_equal(Ci7.edges(), K7.edges())

        # Ci_6(1, 3) is K_3,3 i.e. the utility graph
        Ci6_1_3 = circulant_graph(6, [1, 3])
        K3_3 = complete_bipartite_graph(3, 3)
        assert_true(is_isomorphic(Ci6_1_3, K3_3))

    def test_cycle_graph(self):
        G = cycle_graph(4)
        assert_edges_equal(G.edges(), [(0, 1), (0, 3), (1, 2), (2, 3)])
        mG = cycle_graph(4, create_using=MultiGraph())
        assert_edges_equal(mG.edges(), [(0, 1), (0, 3), (1, 2), (2, 3)])
        G = cycle_graph(4, create_using=DiGraph())
        assert_false(G.has_edge(2, 1))
        assert_true(G.has_edge(1, 2))
        assert_true(G.is_directed())

        G = cycle_graph("abc")
        assert_equal(len(G), 3)
        assert_equal(G.size(), 3)
        g = cycle_graph("abc", nx.DiGraph())
        assert_equal(len(g), 3)
        assert_equal(g.size(), 3)
        assert_true(g.is_directed())

    def test_dorogovtsev_goltsev_mendes_graph(self):
        G = dorogovtsev_goltsev_mendes_graph(0)
        assert_edges_equal(G.edges(), [(0, 1)])
        assert_nodes_equal(list(G), [0, 1])
        G = dorogovtsev_goltsev_mendes_graph(1)
        assert_edges_equal(G.edges(), [(0, 1), (0, 2), (1, 2)])
        assert_equal(average_clustering(G), 1.0)
        assert_equal(sorted(triangles(G).values()), [1, 1, 1])
        G = dorogovtsev_goltsev_mendes_graph(10)
        assert_equal(number_of_nodes(G), 29526)
        assert_equal(number_of_edges(G), 59049)
        assert_equal(G.degree(0), 1024)
        assert_equal(G.degree(1), 1024)
        assert_equal(G.degree(2), 1024)

        assert_raises(networkx.exception.NetworkXError,
                      dorogovtsev_goltsev_mendes_graph, 7,
                      create_using=DiGraph())
        assert_raises(networkx.exception.NetworkXError,
                      dorogovtsev_goltsev_mendes_graph, 7,
                      create_using=MultiGraph())

    def test_empty_graph(self):
        G = empty_graph()
        assert_equal(number_of_nodes(G), 0)
        G = empty_graph(42)
        assert_equal(number_of_nodes(G), 42)
        assert_equal(number_of_edges(G), 0)

        G = empty_graph("abc")
        assert_equal(len(G), 3)
        assert_equal(G.size(), 0)

        # create empty digraph
        G = empty_graph(42, create_using=DiGraph(name="duh"))
        assert_equal(number_of_nodes(G), 42)
        assert_equal(number_of_edges(G), 0)
        assert_true(isinstance(G, DiGraph))

        # create empty multigraph
        G = empty_graph(42, create_using=MultiGraph(name="duh"))
        assert_equal(number_of_nodes(G), 42)
        assert_equal(number_of_edges(G), 0)
        assert_true(isinstance(G, MultiGraph))

        # create empty graph from another
        pete = petersen_graph()
        G = empty_graph(42, create_using=pete)
        assert_equal(number_of_nodes(G), 42)
        assert_equal(number_of_edges(G), 0)
        assert_true(isinstance(G, Graph))

    def test_ladder_graph(self):
        for i, G in [(0, empty_graph(0)), (1, path_graph(2)),
                     (2, hypercube_graph(2)), (10, grid_graph([2, 10]))]:
            assert_true(is_isomorphic(ladder_graph(i), G))

        assert_raises(networkx.exception.NetworkXError,
                      ladder_graph, 2, create_using=DiGraph())

        g = ladder_graph(2)
        mg = ladder_graph(2, create_using=MultiGraph())
        assert_edges_equal(mg.edges(), g.edges())

    def test_lollipop_graph(self):
        # number of nodes = m1 + m2
        # number of edges = number_of_edges(complete_graph(m1)) + m2
        for m1, m2 in [(3, 5), (4, 10), (3, 20)]:
            b = lollipop_graph(m1, m2)
            assert_equal(number_of_nodes(b), m1 + m2)
            assert_equal(number_of_edges(b), m1 * (m1 - 1) / 2 + m2)

        # Raise NetworkXError if m<2
        assert_raises(networkx.exception.NetworkXError,
                      lollipop_graph, 1, 20)

        # Raise NetworkXError if n<0
        assert_raises(networkx.exception.NetworkXError,
                      lollipop_graph, 5, -2)

        # lollipop_graph(2,m) = path_graph(m+2)
        for m1, m2 in [(2, 5), (2, 10), (2, 20)]:
            b = lollipop_graph(m1, m2)
            assert_true(is_isomorphic(b, path_graph(m2 + 2)))

        assert_raises(networkx.exception.NetworkXError,
                      lollipop_graph, m1, m2, create_using=DiGraph())

        mb = lollipop_graph(m1, m2, create_using=MultiGraph())
        assert_edges_equal(mb.edges(), b.edges())

        g = lollipop_graph([1, 2, 3, 4], "abc")
        assert_equal(len(g), 7)
        assert_equal(g.size(), 9)

    def test_null_graph(self):
        assert_equal(number_of_nodes(null_graph()), 0)

    def test_path_graph(self):
        p = path_graph(0)
        assert_true(is_isomorphic(p, null_graph()))

        p = path_graph(1)
        assert_true(is_isomorphic(p, empty_graph(1)))

        p = path_graph(10)
        assert_true(is_connected(p))
        assert_equal(sorted(d for n, d in p.degree()),
                     [1, 1, 2, 2, 2, 2, 2, 2, 2, 2])
        assert_equal(p.order() - 1, p.size())

        dp = path_graph(3, create_using=DiGraph())
        assert_true(dp.has_edge(0, 1))
        assert_false(dp.has_edge(1, 0))

        mp = path_graph(10, create_using=MultiGraph())
        assert_edges_equal(mp.edges(), p.edges())

        G = path_graph("abc")
        assert_equal(len(G), 3)
        assert_equal(G.size(), 2)
        g = path_graph("abc", nx.DiGraph())
        assert_equal(len(g), 3)
        assert_equal(g.size(), 2)
        assert_true(g.is_directed())

    def test_star_graph(self):
        assert_true(is_isomorphic(star_graph(0), empty_graph(1)))
        assert_true(is_isomorphic(star_graph(1), path_graph(2)))
        assert_true(is_isomorphic(star_graph(2), path_graph(3)))
        assert_true(is_isomorphic(star_graph(5), nx.complete_bipartite_graph(1, 5)))

        s = star_graph(10)
        assert_equal(sorted(d for n, d in s.degree()),
                     [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 10])

        assert_raises(networkx.exception.NetworkXError,
                      star_graph, 10, create_using=DiGraph())

        ms = star_graph(10, create_using=MultiGraph())
        assert_edges_equal(ms.edges(), s.edges())

        G = star_graph("abcdefg")
        assert_equal(len(G), 7)
        assert_equal(G.size(), 6)

    def test_trivial_graph(self):
        assert_equal(number_of_nodes(trivial_graph()), 1)

    def test_turan_graph(self):
        assert_equal(number_of_edges(turan_graph(13, 4)), 63)
        assert_true(is_isomorphic(turan_graph(13, 4), complete_multipartite_graph(3, 4, 3, 3)))

    def test_wheel_graph(self):
        for n, G in [(0, null_graph()), (1, empty_graph(1)),
                     (2, path_graph(2)), (3, complete_graph(3)),
                     (4, complete_graph(4))]:
            g = wheel_graph(n)
            assert_true(is_isomorphic(g, G))

        g = wheel_graph(10)
        assert_equal(sorted(d for n, d in g.degree()),
                     [3, 3, 3, 3, 3, 3, 3, 3, 3, 9])

        assert_raises(networkx.exception.NetworkXError,
                      wheel_graph, 10, create_using=DiGraph())

        mg = wheel_graph(10, create_using=MultiGraph())
        assert_edges_equal(mg.edges(), g.edges())

        G = wheel_graph("abc")
        assert_equal(len(G), 3)
        assert_equal(G.size(), 3)

    def test_complete_0_partite_graph(self):
        """Tests that the complete 0-partite graph is the null graph."""
        G = nx.complete_multipartite_graph()
        H = nx.null_graph()
        assert_nodes_equal(G, H)
        assert_edges_equal(G.edges(), H.edges())

    def test_complete_1_partite_graph(self):
        """Tests that the complete 1-partite graph is the empty graph."""
        G = nx.complete_multipartite_graph(3)
        H = nx.empty_graph(3)
        assert_nodes_equal(G, H)
        assert_edges_equal(G.edges(), H.edges())

    def test_complete_2_partite_graph(self):
        """Tests that the complete 2-partite graph is the complete bipartite
        graph.

        """
        G = nx.complete_multipartite_graph(2, 3)
        H = nx.complete_bipartite_graph(2, 3)
        assert_nodes_equal(G, H)
        assert_edges_equal(G.edges(), H.edges())

    def test_complete_multipartite_graph(self):
        """Tests for generating the complete multipartite graph."""
        G = nx.complete_multipartite_graph(2, 3, 4)
        blocks = [(0, 1), (2, 3, 4), (5, 6, 7, 8)]
        # Within each block, no two vertices should be adjacent.
        for block in blocks:
            for u, v in itertools.combinations_with_replacement(block, 2):
                assert_true(v not in G[u])
                assert_equal(G.nodes[u], G.nodes[v])
        # Across blocks, all vertices should be adjacent.
        for (block1, block2) in itertools.combinations(blocks, 2):
            for u, v in itertools.product(block1, block2):
                assert_true(v in G[u])
                assert_not_equal(G.nodes[u], G.nodes[v])
