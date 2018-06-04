from itertools import chain
import networkx as nx
from nose.tools import *


class TestIsSemiconnected(object):

    def test_undirected(self):
        assert_raises(nx.NetworkXNotImplemented, nx.is_semiconnected,
                      nx.Graph())
        assert_raises(nx.NetworkXNotImplemented, nx.is_semiconnected,
                      nx.MultiGraph())

    def test_empty(self):
        assert_raises(nx.NetworkXPointlessConcept, nx.is_semiconnected,
                      nx.DiGraph())
        assert_raises(nx.NetworkXPointlessConcept, nx.is_semiconnected,
                      nx.MultiDiGraph())

    def test_single_node_graph(self):
        G = nx.DiGraph()
        G.add_node(0)
        ok_(nx.is_semiconnected(G))

    def test_path(self):
        G = nx.path_graph(100, create_using=nx.DiGraph())
        ok_(nx.is_semiconnected(G))
        G.add_edge(100, 99)
        ok_(not nx.is_semiconnected(G))

    def test_cycle(self):
        G = nx.cycle_graph(100, create_using=nx.DiGraph())
        ok_(nx.is_semiconnected(G))
        G = nx.path_graph(100, create_using=nx.DiGraph())
        G.add_edge(0, 99)
        ok_(nx.is_semiconnected(G))

    def test_tree(self):
        G = nx.DiGraph()
        G.add_edges_from(chain.from_iterable([(i, 2 * i + 1), (i, 2 * i + 2)]
                                             for i in range(100)))
        ok_(not nx.is_semiconnected(G))

    def test_dumbbell(self):
        G = nx.cycle_graph(100, create_using=nx.DiGraph())
        G.add_edges_from((i + 100, (i + 1) % 100 + 100) for i in range(100))
        ok_(not nx.is_semiconnected(G))  # G is disconnected.
        G.add_edge(100, 99)
        ok_(nx.is_semiconnected(G))

    def test_alternating_path(self):
        G = nx.DiGraph(chain.from_iterable([(i, i - 1), (i, i + 1)]
                                           for i in range(0, 100, 2)))
        ok_(not nx.is_semiconnected(G))
