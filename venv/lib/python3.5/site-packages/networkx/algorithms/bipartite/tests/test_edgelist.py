"""
    Unit tests for bipartite edgelists.
"""
from nose.tools import assert_equal, assert_raises, assert_not_equal, raises
import io
import tempfile
import os

import networkx as nx
from networkx.testing import (assert_edges_equal, assert_nodes_equal,
                              assert_graphs_equal)
from networkx.algorithms import bipartite


class TestEdgelist:

    def setUp(self):
        self.G = nx.Graph(name="test")
        e = [('a', 'b'), ('b', 'c'), ('c', 'd'), ('d', 'e'), ('e', 'f'), ('a', 'f')]
        self.G.add_edges_from(e)
        self.G.add_nodes_from(['a', 'c', 'e'], bipartite=0)
        self.G.add_nodes_from(['b', 'd', 'f'], bipartite=1)
        self.G.add_node('g', bipartite=0)
        self.DG = nx.DiGraph(self.G)
        self.MG = nx.MultiGraph()
        self.MG.add_edges_from([(1, 2), (1, 2), (1, 2)])
        self.MG.add_node(1, bipartite=0)
        self.MG.add_node(2, bipartite=1)

    def test_read_edgelist_1(self):
        s = b"""\
# comment line
1 2
# comment line
2 3
"""
        bytesIO = io.BytesIO(s)
        G = bipartite.read_edgelist(bytesIO, nodetype=int)
        assert_edges_equal(G.edges(), [(1, 2), (2, 3)])

    def test_read_edgelist_3(self):
        s = b"""\
# comment line
1 2 {'weight':2.0}
# comment line
2 3 {'weight':3.0}
"""
        bytesIO = io.BytesIO(s)
        G = bipartite.read_edgelist(bytesIO, nodetype=int, data=False)
        assert_edges_equal(G.edges(), [(1, 2), (2, 3)])

        bytesIO = io.BytesIO(s)
        G = bipartite.read_edgelist(bytesIO, nodetype=int, data=True)
        assert_edges_equal(G.edges(data=True),
                           [(1, 2, {'weight': 2.0}), (2, 3, {'weight': 3.0})])

    def test_write_edgelist_1(self):
        fh = io.BytesIO()
        G = nx.Graph()
        G.add_edges_from([(1, 2), (2, 3)])
        G.add_node(1, bipartite=0)
        G.add_node(2, bipartite=1)
        G.add_node(3, bipartite=0)
        bipartite.write_edgelist(G, fh, data=False)
        fh.seek(0)
        assert_equal(fh.read(), b"1 2\n3 2\n")

    def test_write_edgelist_2(self):
        fh = io.BytesIO()
        G = nx.Graph()
        G.add_edges_from([(1, 2), (2, 3)])
        G.add_node(1, bipartite=0)
        G.add_node(2, bipartite=1)
        G.add_node(3, bipartite=0)
        bipartite.write_edgelist(G, fh, data=True)
        fh.seek(0)
        assert_equal(fh.read(), b"1 2 {}\n3 2 {}\n")

    def test_write_edgelist_3(self):
        fh = io.BytesIO()
        G = nx.Graph()
        G.add_edge(1, 2, weight=2.0)
        G.add_edge(2, 3, weight=3.0)
        G.add_node(1, bipartite=0)
        G.add_node(2, bipartite=1)
        G.add_node(3, bipartite=0)
        bipartite.write_edgelist(G, fh, data=True)
        fh.seek(0)
        assert_equal(fh.read(), b"1 2 {'weight': 2.0}\n3 2 {'weight': 3.0}\n")

    def test_write_edgelist_4(self):
        fh = io.BytesIO()
        G = nx.Graph()
        G.add_edge(1, 2, weight=2.0)
        G.add_edge(2, 3, weight=3.0)
        G.add_node(1, bipartite=0)
        G.add_node(2, bipartite=1)
        G.add_node(3, bipartite=0)
        bipartite.write_edgelist(G, fh, data=[('weight')])
        fh.seek(0)
        assert_equal(fh.read(), b"1 2 2.0\n3 2 3.0\n")

    def test_unicode(self):
        G = nx.Graph()
        try:  # Python 3.x
            name1 = chr(2344) + chr(123) + chr(6543)
            name2 = chr(5543) + chr(1543) + chr(324)
        except ValueError:  # Python 2.6+
            name1 = unichr(2344) + unichr(123) + unichr(6543)
            name2 = unichr(5543) + unichr(1543) + unichr(324)
        G.add_edge(name1, 'Radiohead', **{name2: 3})
        G.add_node(name1, bipartite=0)
        G.add_node('Radiohead', bipartite=1)
        fd, fname = tempfile.mkstemp()
        bipartite.write_edgelist(G, fname)
        H = bipartite.read_edgelist(fname)
        assert_graphs_equal(G, H)
        os.close(fd)
        os.unlink(fname)

    def test_latin1_issue(self):
        G = nx.Graph()
        try:  # Python 3.x
            name1 = chr(2344) + chr(123) + chr(6543)
            name2 = chr(5543) + chr(1543) + chr(324)
        except ValueError:  # Python 2.6+
            name1 = unichr(2344) + unichr(123) + unichr(6543)
            name2 = unichr(5543) + unichr(1543) + unichr(324)
        G.add_edge(name1, 'Radiohead', **{name2: 3})
        G.add_node(name1, bipartite=0)
        G.add_node('Radiohead', bipartite=1)
        fd, fname = tempfile.mkstemp()
        assert_raises(UnicodeEncodeError,
                      bipartite.write_edgelist,
                      G, fname, encoding='latin-1')
        os.close(fd)
        os.unlink(fname)

    def test_latin1(self):
        G = nx.Graph()
        try:  # Python 3.x
            blurb = chr(1245)  # just to trigger the exception
            name1 = 'Bj' + chr(246) + 'rk'
            name2 = chr(220) + 'ber'
        except ValueError:  # Python 2.6+
            name1 = 'Bj' + unichr(246) + 'rk'
            name2 = unichr(220) + 'ber'
        G.add_edge(name1, 'Radiohead', **{name2: 3})
        G.add_node(name1, bipartite=0)
        G.add_node('Radiohead', bipartite=1)
        fd, fname = tempfile.mkstemp()
        bipartite.write_edgelist(G, fname, encoding='latin-1')
        H = bipartite.read_edgelist(fname, encoding='latin-1')
        assert_graphs_equal(G, H)
        os.close(fd)
        os.unlink(fname)

    def test_edgelist_graph(self):
        G = self.G
        (fd, fname) = tempfile.mkstemp()
        bipartite.write_edgelist(G, fname)
        H = bipartite.read_edgelist(fname)
        H2 = bipartite.read_edgelist(fname)
        assert_not_equal(H, H2)  # they should be different graphs
        G.remove_node('g')  # isolated nodes are not written in edgelist
        assert_nodes_equal(list(H), list(G))
        assert_edges_equal(list(H.edges()), list(G.edges()))
        os.close(fd)
        os.unlink(fname)

    def test_edgelist_integers(self):
        G = nx.convert_node_labels_to_integers(self.G)
        (fd, fname) = tempfile.mkstemp()
        bipartite.write_edgelist(G, fname)
        H = bipartite.read_edgelist(fname, nodetype=int)
        # isolated nodes are not written in edgelist
        G.remove_nodes_from(list(nx.isolates(G)))
        assert_nodes_equal(list(H), list(G))
        assert_edges_equal(list(H.edges()), list(G.edges()))
        os.close(fd)
        os.unlink(fname)

    def test_edgelist_multigraph(self):
        G = self.MG
        (fd, fname) = tempfile.mkstemp()
        bipartite.write_edgelist(G, fname)
        H = bipartite.read_edgelist(fname, nodetype=int, create_using=nx.MultiGraph())
        H2 = bipartite.read_edgelist(fname, nodetype=int, create_using=nx.MultiGraph())
        assert_not_equal(H, H2)  # they should be different graphs
        assert_nodes_equal(list(H), list(G))
        assert_edges_equal(list(H.edges()), list(G.edges()))
        os.close(fd)
        os.unlink(fname)

    @raises(nx.NetworkXNotImplemented)
    def test_empty_digraph(self):
        bytesIO = io.BytesIO()
        bipartite.write_edgelist(nx.DiGraph(), bytesIO)

    @raises(AttributeError)
    def test_raise_attribute(self):
        G = nx.path_graph(4)
        bytesIO = io.BytesIO()
        bipartite.write_edgelist(G, bytesIO)
