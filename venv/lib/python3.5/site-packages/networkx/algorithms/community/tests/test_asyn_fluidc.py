from nose.tools import assert_equal
from networkx import Graph
from networkx.algorithms.community.asyn_fluidc import *
import random


def test_single_node():
    test = Graph()

    test.add_node('a')

    # ground truth
    ground_truth = set([frozenset(['a'])])

    communities = asyn_fluidc(test, 1)
    result = {frozenset(c) for c in communities}
    assert_equal(result, ground_truth)


def test_two_nodes():
    test = Graph()

    test.add_edge('a', 'b')

    # ground truth
    ground_truth = set([frozenset(['a']), frozenset(['b'])])

    communities = asyn_fluidc(test, 2)
    result = {frozenset(c) for c in communities}
    assert_equal(result, ground_truth)


def test_two_clique_communities():
    random.seed(7)
    test = Graph()

    # c1
    test.add_edge('a', 'b')
    test.add_edge('a', 'c')
    test.add_edge('b', 'c')

    # connection
    test.add_edge('c', 'd')

    # c2
    test.add_edge('d', 'e')
    test.add_edge('d', 'f')
    test.add_edge('f', 'e')

    # ground truth
    ground_truth = set([frozenset(['a', 'c', 'b']),
                        frozenset(['e', 'd', 'f'])])

    communities = asyn_fluidc(test, 2)
    result = {frozenset(c) for c in communities}
    assert_equal(result, ground_truth)


def five_clique_ring():
    """Not auto-tested (not named test_...) due to cross-version seed issues"""
    random.seed(9)
    test = Graph()

    # c1
    test.add_edge('1a', '1b')
    test.add_edge('1a', '1c')
    test.add_edge('1a', '1d')
    test.add_edge('1b', '1c')
    test.add_edge('1b', '1d')
    test.add_edge('1c', '1d')

    # c2
    test.add_edge('2a', '2b')
    test.add_edge('2a', '2c')
    test.add_edge('2a', '2d')
    test.add_edge('2b', '2c')
    test.add_edge('2b', '2d')
    test.add_edge('2c', '2d')

    # c3
    test.add_edge('3a', '3b')
    test.add_edge('3a', '3c')
    test.add_edge('3a', '3d')
    test.add_edge('3b', '3c')
    test.add_edge('3b', '3d')
    test.add_edge('3c', '3d')

    # c4
    test.add_edge('4a', '4b')
    test.add_edge('4a', '4c')
    test.add_edge('4a', '4d')
    test.add_edge('4b', '4c')
    test.add_edge('4b', '4d')
    test.add_edge('4c', '4d')

    # c5
    test.add_edge('5a', '5b')
    test.add_edge('5a', '5c')
    test.add_edge('5a', '5d')
    test.add_edge('5b', '5c')
    test.add_edge('5b', '5d')
    test.add_edge('5c', '5d')

    # connections
    test.add_edge('1a', '2c')
    test.add_edge('2a', '3c')
    test.add_edge('3a', '4c')
    test.add_edge('4a', '5c')
    test.add_edge('5a', '1c')

    # ground truth
    ground_truth = set([frozenset(['1a', '1b', '1c', '1d']),
                        frozenset(['2a', '2b', '2c', '2d']),
                        frozenset(['3a', '3b', '3c', '3d']),
                        frozenset(['4a', '4b', '4c', '4d']),
                        frozenset(['5a', '5b', '5c', '5d'])])

    communities = asyn_fluidc(test, 5)
    result = {frozenset(c) for c in communities}
    assert_equal(result, ground_truth)
