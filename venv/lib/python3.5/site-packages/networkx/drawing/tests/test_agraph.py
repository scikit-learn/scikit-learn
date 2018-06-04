"""Unit tests for PyGraphviz inteface."""
import os
import tempfile
from nose import SkipTest
from nose.tools import assert_true, assert_equal, assert_raises
from networkx.testing import assert_edges_equal, assert_nodes_equal

import networkx as nx


class TestAGraph(object):
    @classmethod
    def setupClass(cls):
        global pygraphviz
        try:
            import pygraphviz
        except ImportError:
            raise SkipTest('PyGraphviz not available.')

    def build_graph(self, G):
        edges = [('A', 'B'), ('A', 'C'), ('A', 'C'), ('B', 'C'), ('A', 'D')]
        G.add_edges_from(edges)
        G.add_node('E')
        G.graph['metal'] = 'bronze'
        return G

    def assert_equal(self, G1, G2):
        assert_nodes_equal(G1.nodes(), G2.nodes())
        assert_edges_equal(G1.edges(), G2.edges())
        assert_equal(G1.graph['metal'], G2.graph['metal'])

    def agraph_checks(self, G):
        G = self.build_graph(G)
        A = nx.nx_agraph.to_agraph(G)
        H = nx.nx_agraph.from_agraph(A)
        self.assert_equal(G, H)

        fname = tempfile.mktemp()
        nx.drawing.nx_agraph.write_dot(H, fname)
        Hin = nx.nx_agraph.read_dot(fname)
        os.unlink(fname)
        self.assert_equal(H, Hin)

        (fd, fname) = tempfile.mkstemp()
        with open(fname, 'w') as fh:
            nx.drawing.nx_agraph.write_dot(H, fh)

        with open(fname, 'r') as fh:
            Hin = nx.nx_agraph.read_dot(fh)
        os.unlink(fname)
        self.assert_equal(H, Hin)

    def test_from_agraph_name(self):
        G = nx.Graph(name='test')
        A = nx.nx_agraph.to_agraph(G)
        H = nx.nx_agraph.from_agraph(A)
        assert_equal(G.name, 'test')

    def test_undirected(self):
        self.agraph_checks(nx.Graph())

    def test_directed(self):
        self.agraph_checks(nx.DiGraph())

    def test_multi_undirected(self):
        self.agraph_checks(nx.MultiGraph())

    def test_multi_directed(self):
        self.agraph_checks(nx.MultiDiGraph())

    def test_view_pygraphviz(self):
        G = nx.Graph()  # "An empty graph cannot be drawn."
        assert_raises(nx.NetworkXException, nx.nx_agraph.view_pygraphviz, G)
        G = nx.barbell_graph(4, 6)
        nx.nx_agraph.view_pygraphviz(G)

    def test_view_pygraphviz_edgelable(self):
        G = nx.Graph()
        G.add_edge(1, 2, weight=7)
        G.add_edge(2, 3, weight=8)
        nx.nx_agraph.view_pygraphviz(G, edgelabel='weight')

    def test_from_agraph_name(self):
        G = nx.Graph(name='test')
        A = nx.nx_agraph.to_agraph(G)
        H = nx.nx_agraph.from_agraph(A)
        assert_equal(G.name, 'test')

    def test_graph_with_reserved_keywords(self):
        # test attribute/keyword clash case for #1582
        # node: n
        # edges: u,v
        G = nx.Graph()
        G = self.build_graph(G)
        G.node['E']['n'] = 'keyword'
        G.edges[('A', 'B')]['u'] = 'keyword'
        G.edges[('A', 'B')]['v'] = 'keyword'
        A = nx.nx_agraph.to_agraph(G)
