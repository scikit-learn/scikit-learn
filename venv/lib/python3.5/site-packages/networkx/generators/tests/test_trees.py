from nose.tools import assert_equal, assert_true

import networkx as nx
from networkx.generators.trees import NIL
from networkx.utils import arbitrary_element


class TestPrefixTree(object):
    """Unit tests for the prefix tree generator function."""

    def test_basic(self):
        # This example is from the Wikipedia article "Trie"
        # <https://en.wikipedia.org/wiki/Trie>.
        strings = ['a', 'to', 'tea', 'ted', 'ten', 'i', 'in', 'inn']
        T, root = nx.prefix_tree(strings)

        def source_label(v): return T.node[v]['source']

        # First, we check that the tree has the expected
        # structure. Recall that each node that corresponds to one of
        # the input strings has an edge to the NIL node.
        #
        # Consider the three children at level 1 in the trie.
        a, i, t = sorted(T[root], key=source_label)
        # Check the 'a' branch.
        assert_equal(len(T[a]), 1)
        nil = arbitrary_element(T[a])
        assert_equal(len(T[nil]), 0)
        # Check the 'i' branch.
        assert_equal(len(T[i]), 2)
        nil, in_ = sorted(T[i], key=source_label)
        assert_equal(len(T[nil]), 0)
        assert_equal(len(T[in_]), 2)
        nil, inn = sorted(T[in_], key=source_label)
        assert_equal(len(T[nil]), 0)
        assert_equal(len(T[inn]), 1)
        nil = arbitrary_element(T[inn])
        assert_equal(len(T[nil]), 0)
        # Check the 't' branch.
        te, to = sorted(T[t], key=source_label)
        assert_equal(len(T[to]), 1)
        nil = arbitrary_element(T[to])
        assert_equal(len(T[nil]), 0)
        tea, ted, ten = sorted(T[te], key=source_label)
        assert_equal(len(T[tea]), 1)
        assert_equal(len(T[ted]), 1)
        assert_equal(len(T[ten]), 1)
        nil = arbitrary_element(T[tea])
        assert_equal(len(T[nil]), 0)
        nil = arbitrary_element(T[ted])
        assert_equal(len(T[nil]), 0)
        nil = arbitrary_element(T[ten])
        assert_equal(len(T[nil]), 0)

        # Next, we check that the "sources" of each of the nodes is the
        # rightmost letter in the string corresponding to the path to
        # that node.
        assert_equal(source_label(root), None)
        assert_equal(source_label(a), 'a')
        assert_equal(source_label(i), 'i')
        assert_equal(source_label(t), 't')
        assert_equal(source_label(in_), 'n')
        assert_equal(source_label(inn), 'n')
        assert_equal(source_label(to), 'o')
        assert_equal(source_label(te), 'e')
        assert_equal(source_label(tea), 'a')
        assert_equal(source_label(ted), 'd')
        assert_equal(source_label(ten), 'n')
        assert_equal(source_label(NIL), NIL)


def test_random_tree():
    """Tests that a random tree is in fact a tree."""
    T = nx.random_tree(10, seed=1234)
    assert_true(nx.is_tree(T))
