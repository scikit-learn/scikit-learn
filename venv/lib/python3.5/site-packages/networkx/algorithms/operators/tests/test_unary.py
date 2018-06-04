from nose.tools import *
import networkx as nx
from networkx import *


def test_complement():
    null = null_graph()
    empty1 = empty_graph(1)
    empty10 = empty_graph(10)
    K3 = complete_graph(3)
    K5 = complete_graph(5)
    K10 = complete_graph(10)
    P2 = path_graph(2)
    P3 = path_graph(3)
    P5 = path_graph(5)
    P10 = path_graph(10)
    # complement of the complete graph is empty

    G = complement(K3)
    assert_true(is_isomorphic(G, empty_graph(3)))
    G = complement(K5)
    assert_true(is_isomorphic(G, empty_graph(5)))
    # for any G, G=complement(complement(G))
    P3cc = complement(complement(P3))
    assert_true(is_isomorphic(P3, P3cc))
    nullcc = complement(complement(null))
    assert_true(is_isomorphic(null, nullcc))
    b = bull_graph()
    bcc = complement(complement(b))
    assert_true(is_isomorphic(b, bcc))


def test_complement_2():
    G1 = nx.DiGraph()
    G1.add_edge('A', 'B')
    G1.add_edge('A', 'C')
    G1.add_edge('A', 'D')
    G1C = complement(G1)
    assert_equal(sorted(G1C.edges()),
                 [('B', 'A'), ('B', 'C'),
                  ('B', 'D'), ('C', 'A'), ('C', 'B'),
                  ('C', 'D'), ('D', 'A'), ('D', 'B'), ('D', 'C')])


def test_reverse1():
    # Other tests for reverse are done by the DiGraph and MultiDigraph.
    G1 = nx.Graph()
    assert_raises(nx.NetworkXError, nx.reverse, G1)
