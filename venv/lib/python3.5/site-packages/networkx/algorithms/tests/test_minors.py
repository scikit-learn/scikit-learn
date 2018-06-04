# test_minors.py - unit tests for the minors module
#
# Copyright 2015 Jeffrey Finkelstein <jeffrey.finkelstein@gmail.com>.
#
# This file is part of NetworkX.
#
# NetworkX is distributed under a BSD license; see LICENSE.txt for more
# information.
"""Unit tests for the :mod:`networkx.algorithms.minors` module."""
from nose.tools import assert_equal
from nose.tools import assert_true
from nose.tools import raises

import networkx as nx
from networkx.testing.utils import *
from networkx.utils import arbitrary_element


class TestQuotient(object):
    """Unit tests for computing quotient graphs."""

    def test_quotient_graph_complete_multipartite(self):
        """Tests that the quotient graph of the complete *n*-partite graph
        under the "same neighbors" node relation is the complete graph on *n*
        nodes.

        """
        G = nx.complete_multipartite_graph(2, 3, 4)
        # Two nodes are equivalent if they are not adjacent but have the same
        # neighbor set.

        def same_neighbors(u, v):
            return (u not in G[v] and v not in G[u] and G[u] == G[v])

        expected = nx.complete_graph(3)
        actual = nx.quotient_graph(G, same_neighbors)
        # It won't take too long to run a graph isomorphism algorithm on such
        # small graphs.
        assert_true(nx.is_isomorphic(expected, actual))

    def test_quotient_graph_complete_bipartite(self):
        """Tests that the quotient graph of the complete bipartite graph under
        the "same neighbors" node relation is `K_2`.

        """
        G = nx.complete_bipartite_graph(2, 3)
        # Two nodes are equivalent if they are not adjacent but have the same
        # neighbor set.

        def same_neighbors(u, v):
            return (u not in G[v] and v not in G[u] and G[u] == G[v])

        expected = nx.complete_graph(2)
        actual = nx.quotient_graph(G, same_neighbors)
        # It won't take too long to run a graph isomorphism algorithm on such
        # small graphs.
        assert_true(nx.is_isomorphic(expected, actual))

    def test_quotient_graph_edge_relation(self):
        """Tests for specifying an alternate edge relation for the quotient
        graph.

        """
        G = nx.path_graph(5)

        def identity(u, v):
            return u == v

        def same_parity(b, c):
            return (arbitrary_element(b) % 2 == arbitrary_element(c) % 2)

        actual = nx.quotient_graph(G, identity, same_parity)
        expected = nx.Graph()
        expected.add_edges_from([(0, 2), (0, 4), (2, 4)])
        expected.add_edge(1, 3)
        assert_true(nx.is_isomorphic(actual, expected))

    def test_condensation_as_quotient(self):
        """This tests that the condensation of a graph can be viewed as the
        quotient graph under the "in the same connected component" equivalence
        relation.

        """
        # This example graph comes from the file `test_strongly_connected.py`.
        G = nx.DiGraph()
        G.add_edges_from([(1, 2), (2, 3), (2, 11), (2, 12), (3, 4), (4, 3),
                          (4, 5), (5, 6), (6, 5), (6, 7), (7, 8), (7, 9),
                          (7, 10), (8, 9), (9, 7), (10, 6), (11, 2), (11, 4),
                          (11, 6), (12, 6), (12, 11)])
        scc = list(nx.strongly_connected_components(G))
        C = nx.condensation(G, scc)
        component_of = C.graph['mapping']
        # Two nodes are equivalent if they are in the same connected component.

        def same_component(u, v):
            return component_of[u] == component_of[v]

        Q = nx.quotient_graph(G, same_component)
        assert_true(nx.is_isomorphic(C, Q))

    def test_path(self):
        G = nx.path_graph(6)
        partition = [{0, 1}, {2, 3}, {4, 5}]
        M = nx.quotient_graph(G, partition, relabel=True)
        assert_nodes_equal(M, [0, 1, 2])
        assert_edges_equal(M.edges(), [(0, 1), (1, 2)])
        for n in M:
            assert_equal(M.nodes[n]['nedges'], 1)
            assert_equal(M.nodes[n]['nnodes'], 2)
            assert_equal(M.nodes[n]['density'], 1)

    def test_multigraph_path(self):
        G = nx.MultiGraph(nx.path_graph(6))
        partition = [{0, 1}, {2, 3}, {4, 5}]
        M = nx.quotient_graph(G, partition, relabel=True)
        assert_nodes_equal(M, [0, 1, 2])
        assert_edges_equal(M.edges(), [(0, 1), (1, 2)])
        for n in M:
            assert_equal(M.nodes[n]['nedges'], 1)
            assert_equal(M.nodes[n]['nnodes'], 2)
            assert_equal(M.nodes[n]['density'], 1)

    def test_directed_path(self):
        G = nx.DiGraph()
        nx.add_path(G, range(6))
        partition = [{0, 1}, {2, 3}, {4, 5}]
        M = nx.quotient_graph(G, partition, relabel=True)
        assert_nodes_equal(M, [0, 1, 2])
        assert_edges_equal(M.edges(), [(0, 1), (1, 2)])
        for n in M:
            assert_equal(M.nodes[n]['nedges'], 1)
            assert_equal(M.nodes[n]['nnodes'], 2)
            assert_equal(M.nodes[n]['density'], 0.5)

    def test_directed_multigraph_path(self):
        G = nx.MultiDiGraph()
        nx.add_path(G, range(6))
        partition = [{0, 1}, {2, 3}, {4, 5}]
        M = nx.quotient_graph(G, partition, relabel=True)
        assert_nodes_equal(M, [0, 1, 2])
        assert_edges_equal(M.edges(), [(0, 1), (1, 2)])
        for n in M:
            assert_equal(M.nodes[n]['nedges'], 1)
            assert_equal(M.nodes[n]['nnodes'], 2)
            assert_equal(M.nodes[n]['density'], 0.5)

    @raises(nx.NetworkXException)
    def test_overlapping_blocks(self):
        G = nx.path_graph(6)
        partition = [{0, 1, 2}, {2, 3}, {4, 5}]
        nx.quotient_graph(G, partition)

    def test_weighted_path(self):
        G = nx.path_graph(6)
        for i in range(5):
            G[i][i + 1]['weight'] = i + 1
        partition = [{0, 1}, {2, 3}, {4, 5}]
        M = nx.quotient_graph(G, partition, relabel=True)
        assert_nodes_equal(M, [0, 1, 2])
        assert_edges_equal(M.edges(), [(0, 1), (1, 2)])
        assert_equal(M[0][1]['weight'], 2)
        assert_equal(M[1][2]['weight'], 4)
        for n in M:
            assert_equal(M.nodes[n]['nedges'], 1)
            assert_equal(M.nodes[n]['nnodes'], 2)
            assert_equal(M.nodes[n]['density'], 1)

    def test_barbell(self):
        G = nx.barbell_graph(3, 0)
        partition = [{0, 1, 2}, {3, 4, 5}]
        M = nx.quotient_graph(G, partition, relabel=True)
        assert_nodes_equal(M, [0, 1])
        assert_edges_equal(M.edges(), [(0, 1)])
        for n in M:
            assert_equal(M.nodes[n]['nedges'], 3)
            assert_equal(M.nodes[n]['nnodes'], 3)
            assert_equal(M.nodes[n]['density'], 1)

    def test_barbell_plus(self):
        G = nx.barbell_graph(3, 0)
        # Add an extra edge joining the bells.
        G.add_edge(0, 5)
        partition = [{0, 1, 2}, {3, 4, 5}]
        M = nx.quotient_graph(G, partition, relabel=True)
        assert_nodes_equal(M, [0, 1])
        assert_edges_equal(M.edges(), [(0, 1)])
        assert_equal(M[0][1]['weight'], 2)
        for n in M:
            assert_equal(M.nodes[n]['nedges'], 3)
            assert_equal(M.nodes[n]['nnodes'], 3)
            assert_equal(M.nodes[n]['density'], 1)

    def test_blockmodel(self):
        G = nx.path_graph(6)
        partition = [[0, 1], [2, 3], [4, 5]]
        M = nx.quotient_graph(G, partition, relabel=True)
        assert_nodes_equal(M.nodes(), [0, 1, 2])
        assert_edges_equal(M.edges(), [(0, 1), (1, 2)])
        for n in M.nodes():
            assert_equal(M.nodes[n]['nedges'], 1)
            assert_equal(M.nodes[n]['nnodes'], 2)
            assert_equal(M.nodes[n]['density'], 1.0)

    def test_multigraph_blockmodel(self):
        G = nx.MultiGraph(nx.path_graph(6))
        partition = [[0, 1], [2, 3], [4, 5]]
        M = nx.quotient_graph(G, partition,
                              create_using=nx.MultiGraph(), relabel=True)
        assert_nodes_equal(M.nodes(), [0, 1, 2])
        assert_edges_equal(M.edges(), [(0, 1), (1, 2)])
        for n in M.nodes():
            assert_equal(M.nodes[n]['nedges'], 1)
            assert_equal(M.nodes[n]['nnodes'], 2)
            assert_equal(M.nodes[n]['density'], 1.0)

    def test_quotient_graph_incomplete_partition(self):
        G = nx.path_graph(6)
        partition = []
        H = nx.quotient_graph(G, partition, relabel=True)
        assert_nodes_equal(H.nodes(), [])
        assert_edges_equal(H.edges(), [])

        partition = [[0, 1], [2, 3], [5]]
        H = nx.quotient_graph(G, partition, relabel=True)
        assert_nodes_equal(H.nodes(), [0, 1, 2])
        assert_edges_equal(H.edges(), [(0, 1)])


class TestContraction(object):
    """Unit tests for node and edge contraction functions."""

    def test_undirected_node_contraction(self):
        """Tests for node contraction in an undirected graph."""
        G = nx.cycle_graph(4)
        actual = nx.contracted_nodes(G, 0, 1)
        expected = nx.complete_graph(3)
        expected.add_edge(0, 0)
        assert_true(nx.is_isomorphic(actual, expected))

    def test_directed_node_contraction(self):
        """Tests for node contraction in a directed graph."""
        G = nx.DiGraph(nx.cycle_graph(4))
        actual = nx.contracted_nodes(G, 0, 1)
        expected = nx.DiGraph(nx.complete_graph(3))
        expected.add_edge(0, 0)
        expected.add_edge(0, 0)
        assert_true(nx.is_isomorphic(actual, expected))

    def test_create_multigraph(self):
        """Tests that using a MultiGraph creates multiple edges."""
        G = nx.path_graph(3, create_using=nx.MultiGraph())
        G.add_edge(0, 1)
        G.add_edge(0, 0)
        G.add_edge(0, 2)
        actual = nx.contracted_nodes(G, 0, 2)
        expected = nx.MultiGraph()
        expected.add_edge(0, 1)
        expected.add_edge(0, 1)
        expected.add_edge(0, 1)
        expected.add_edge(0, 0)
        expected.add_edge(0, 0)
        assert_edges_equal(actual.edges, expected.edges)

    def test_multigraph_keys(self):
        """Tests that multiedge keys are reset in new graph."""
        G = nx.path_graph(3, create_using=nx.MultiGraph())
        G.add_edge(0, 1, 5)
        G.add_edge(0, 0, 0)
        G.add_edge(0, 2, 5)
        actual = nx.contracted_nodes(G, 0, 2)
        expected = nx.MultiGraph()
        expected.add_edge(0, 1, 0)
        expected.add_edge(0, 1, 5)
        expected.add_edge(0, 1, 2)  # keyed as 2 b/c 2 edges already in G
        expected.add_edge(0, 0, 0)
        expected.add_edge(0, 0, 1)  # this comes from (0, 2, 5)
        assert_edges_equal(actual.edges, expected.edges)

    def test_node_attributes(self):
        """Tests that node contraction preserves node attributes."""
        G = nx.cycle_graph(4)
        # Add some data to the two nodes being contracted.
        G.nodes[0]['foo'] = 'bar'
        G.nodes[1]['baz'] = 'xyzzy'
        actual = nx.contracted_nodes(G, 0, 1)
        # We expect that contracting the nodes 0 and 1 in C_4 yields K_3, but
        # with nodes labeled 0, 2, and 3, and with a self-loop on 0.
        expected = nx.complete_graph(3)
        expected = nx.relabel_nodes(expected, {1: 2, 2: 3})
        expected.add_edge(0, 0)
        cdict = {1: {'baz': 'xyzzy'}}
        expected.nodes[0].update(dict(foo='bar', contraction=cdict))
        assert_true(nx.is_isomorphic(actual, expected))
        assert_equal(actual.nodes, expected.nodes)

    def test_without_self_loops(self):
        """Tests for node contraction without preserving self-loops."""
        G = nx.cycle_graph(4)
        actual = nx.contracted_nodes(G, 0, 1, self_loops=False)
        expected = nx.complete_graph(3)
        assert_true(nx.is_isomorphic(actual, expected))

    def test_contract_selfloop_graph(self):
        """Tests for node contraction when nodes have selfloops."""
        G = nx.cycle_graph(4)
        G.add_edge(0, 0)
        actual = nx.contracted_nodes(G, 0, 1)
        expected = nx.complete_graph([0, 2, 3])
        expected.add_edge(0, 0)
        expected.add_edge(0, 0)
        assert_edges_equal(actual.edges, expected.edges)
        actual = nx.contracted_nodes(G, 1, 0)
        expected = nx.complete_graph([1, 2, 3])
        expected.add_edge(1, 1)
        expected.add_edge(1, 1)
        assert_edges_equal(actual.edges, expected.edges)

    def test_undirected_edge_contraction(self):
        """Tests for edge contraction in an undirected graph."""
        G = nx.cycle_graph(4)
        actual = nx.contracted_edge(G, (0, 1))
        expected = nx.complete_graph(3)
        expected.add_edge(0, 0)
        assert_true(nx.is_isomorphic(actual, expected))

    @raises(ValueError)
    def test_nonexistent_edge(self):
        """Tests that attempting to contract a non-existent edge raises an
        exception.

        """
        G = nx.cycle_graph(4)
        nx.contracted_edge(G, (0, 2))
