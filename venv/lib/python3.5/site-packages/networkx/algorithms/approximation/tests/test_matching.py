from nose.tools import *
import networkx as nx
import networkx.algorithms.approximation as a


def test_min_maximal_matching():
    # smoke test
    G = nx.Graph()
    assert_equal(len(a.min_maximal_matching(G)), 0)
