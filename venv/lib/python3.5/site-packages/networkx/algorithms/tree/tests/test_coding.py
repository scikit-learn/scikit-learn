# -*- encoding: utf-8 -*-
# test_coding.py - unit tests for the coding module
#
# Copyright 2015-2018 NetworkX developers.
#
# This file is part of NetworkX.
#
# NetworkX is distributed under a BSD license; see LICENSE.txt for more
# information.
"""Unit tests for the :mod:`~networkx.algorithms.tree.coding` module."""
from itertools import product

from nose.tools import assert_equal
from nose.tools import assert_true
from nose.tools import raises

import networkx as nx
from networkx.testing import assert_nodes_equal
from networkx.testing import assert_edges_equal


class TestPruferSequence(object):
    """Unit tests for the Prüfer sequence encoding and decoding
    functions.

    """

    @raises(nx.NotATree)
    def test_nontree(self):
        G = nx.cycle_graph(3)
        nx.to_prufer_sequence(G)

    @raises(nx.NetworkXPointlessConcept)
    def test_null_graph(self):
        nx.to_prufer_sequence(nx.null_graph())

    @raises(nx.NetworkXPointlessConcept)
    def test_trivial_graph(self):
        nx.to_prufer_sequence(nx.trivial_graph())

    @raises(KeyError)
    def test_bad_integer_labels(self):
        T = nx.Graph(nx.utils.pairwise('abc'))
        nx.to_prufer_sequence(T)

    def test_encoding(self):
        """Tests for encoding a tree as a Prüfer sequence using the
        iterative strategy.

        """
        # Example from Wikipedia.
        tree = nx.Graph([(0, 3), (1, 3), (2, 3), (3, 4), (4, 5)])
        sequence = nx.to_prufer_sequence(tree)
        assert_equal(sequence, [3, 3, 3, 4])

    def test_decoding(self):
        """Tests for decoding a tree from a Prüfer sequence."""
        # Example from Wikipedia.
        sequence = [3, 3, 3, 4]
        tree = nx.from_prufer_sequence(sequence)
        assert_nodes_equal(list(tree), list(range(6)))
        edges = [(0, 3), (1, 3), (2, 3), (3, 4), (4, 5)]
        assert_edges_equal(list(tree.edges()), edges)

    def test_decoding2(self):
        # Example from "An Optimal Algorithm for Prufer Codes".
        sequence = [2, 4, 0, 1, 3, 3]
        tree = nx.from_prufer_sequence(sequence)
        assert_nodes_equal(list(tree), list(range(8)))
        edges = [(0, 1), (0, 4), (1, 3), (2, 4), (2, 5), (3, 6), (3, 7)]
        assert_edges_equal(list(tree.edges()), edges)

    def test_inverse(self):
        """Tests that the encoding and decoding functions are inverses.

        """
        for T in nx.nonisomorphic_trees(4):
            T2 = nx.from_prufer_sequence(nx.to_prufer_sequence(T))
            assert_nodes_equal(list(T), list(T2))
            assert_edges_equal(list(T.edges()), list(T2.edges()))

        for seq in product(range(4), repeat=2):
            seq2 = nx.to_prufer_sequence(nx.from_prufer_sequence(seq))
            assert_equal(list(seq), seq2)


class TestNestedTuple(object):
    """Unit tests for the nested tuple encoding and decoding functions.

    """

    @raises(nx.NotATree)
    def test_nontree(self):
        G = nx.cycle_graph(3)
        nx.to_nested_tuple(G, 0)

    @raises(nx.NodeNotFound)
    def test_unknown_root(self):
        G = nx.path_graph(2)
        nx.to_nested_tuple(G, 'bogus')

    def test_encoding(self):
        T = nx.full_rary_tree(2, 2 ** 3 - 1)
        expected = (((), ()), ((), ()))
        actual = nx.to_nested_tuple(T, 0)
        assert_nodes_equal(expected, actual)

    def test_canonical_form(self):
        T = nx.Graph()
        T.add_edges_from([(0, 1), (0, 2), (0, 3)])
        T.add_edges_from([(1, 4), (1, 5)])
        T.add_edges_from([(3, 6), (3, 7)])
        root = 0
        actual = nx.to_nested_tuple(T, root, canonical_form=True)
        expected = ((), ((), ()), ((), ()))
        assert_equal(actual, expected)

    def test_decoding(self):
        balanced = (((), ()), ((), ()))
        expected = nx.full_rary_tree(2, 2 ** 3 - 1)
        actual = nx.from_nested_tuple(balanced)
        assert_true(nx.is_isomorphic(expected, actual))

    def test_sensible_relabeling(self):
        balanced = (((), ()), ((), ()))
        T = nx.from_nested_tuple(balanced, sensible_relabeling=True)
        edges = [(0, 1), (0, 2), (1, 3), (1, 4), (2, 5), (2, 6)]
        assert_nodes_equal(list(T), list(range(2 ** 3 - 1)))
        assert_edges_equal(list(T.edges()), edges)
