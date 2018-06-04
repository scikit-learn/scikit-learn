from __future__ import division

from io import BytesIO
import tempfile
from unittest import TestCase

from nose.tools import assert_equal

import networkx as nx
import networkx.readwrite.graph6 as g6
from networkx.testing.utils import assert_edges_equal
from networkx.testing.utils import assert_graphs_equal
from networkx.testing.utils import assert_nodes_equal


class TestGraph6Utils(object):

    def test_n_data_n_conversion(self):
        for i in [0, 1, 42, 62, 63, 64, 258047, 258048, 7744773, 68719476735]:
            assert_equal(g6.data_to_n(g6.n_to_data(i))[0], i)
            assert_equal(g6.data_to_n(g6.n_to_data(i))[1], [])
            assert_equal(g6.data_to_n(g6.n_to_data(i) + [42, 43])[1],
                         [42, 43])


class TestFromGraph6Bytes(TestCase):

    def test_from_graph6_bytes(self):
        data = b'DF{'
        G = nx.from_graph6_bytes(data)
        assert_nodes_equal(G.nodes(), [0, 1, 2, 3, 4])
        assert_edges_equal(G.edges(),
                           [(0, 3), (0, 4), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)])

    def test_read_equals_from_bytes(self):
        data = b'DF{'
        G = nx.from_graph6_bytes(data)
        fh = BytesIO(data)
        Gin = nx.read_graph6(fh)
        assert_nodes_equal(G.nodes(), Gin.nodes())
        assert_edges_equal(G.edges(), Gin.edges())


class TestReadGraph6(TestCase):

    def test_read_many_graph6(self):
        """Test for reading many graphs from a file into a list."""
        data = b'DF{\nD`{\nDqK\nD~{\n'
        fh = BytesIO(data)
        glist = nx.read_graph6(fh)
        assert_equal(len(glist), 4)
        for G in glist:
            assert_equal(sorted(G), list(range(5)))


class TestWriteGraph6(TestCase):
    """Unit tests for writing a graph to a file in graph6 format."""

    def test_null_graph(self):
        result = BytesIO()
        nx.write_graph6(nx.null_graph(), result)
        self.assertEqual(result.getvalue(), b'>>graph6<<?\n')

    def test_trivial_graph(self):
        result = BytesIO()
        nx.write_graph6(nx.trivial_graph(), result)
        self.assertEqual(result.getvalue(), b'>>graph6<<@\n')

    def test_complete_graph(self):
        result = BytesIO()
        nx.write_graph6(nx.complete_graph(4), result)
        self.assertEqual(result.getvalue(), b'>>graph6<<C~\n')

    def test_large_complete_graph(self):
        result = BytesIO()
        nx.write_graph6(nx.complete_graph(67), result, header=False)
        self.assertEqual(result.getvalue(), b'~?@B' + b'~' * 368 + b'w\n')

    def test_no_header(self):
        result = BytesIO()
        nx.write_graph6(nx.complete_graph(4), result, header=False)
        self.assertEqual(result.getvalue(), b'C~\n')

    def test_complete_bipartite_graph(self):
        result = BytesIO()
        G = nx.complete_bipartite_graph(6, 9)
        nx.write_graph6(G, result, header=False)
        # The expected encoding here was verified by Sage.
        self.assertEqual(result.getvalue(), b'N??F~z{~Fw^_~?~?^_?\n')

    def no_directed_graphs(self):
        with self.assertRaises(nx.NetworkXNotImplemented):
            nx.write_graph6(nx.DiGraph(), BytesIO())

    def test_length(self):
        for i in list(range(13)) + [31, 47, 62, 63, 64, 72]:
            g = nx.random_graphs.gnm_random_graph(i, i * i // 4, seed=i)
            gstr = BytesIO()
            nx.write_graph6(g, gstr, header=False)
            # Strip the trailing newline.
            gstr = gstr.getvalue().rstrip()
            assert_equal(len(gstr),
                         ((i - 1) * i // 2 + 5) // 6 + (1 if i < 63 else 4))

    def test_roundtrip(self):
        for i in list(range(13)) + [31, 47, 62, 63, 64, 72]:
            G = nx.random_graphs.gnm_random_graph(i, i * i // 4, seed=i)
            f = BytesIO()
            nx.write_graph6(G, f)
            f.seek(0)
            H = nx.read_graph6(f)
            assert_nodes_equal(G.nodes(), H.nodes())
            assert_edges_equal(G.edges(), H.edges())

    def test_write_path(self):
        with tempfile.NamedTemporaryFile() as f:
            g6.write_graph6_file(nx.null_graph(), f)
            f.seek(0)
            self.assertEqual(f.read(), b'>>graph6<<?\n')

    def test_relabeling(self):
        G = nx.Graph([(0, 1)])
        assert_equal(g6.to_graph6_bytes(G), b'>>graph6<<A_\n')
        G = nx.Graph([(1, 2)])
        assert_equal(g6.to_graph6_bytes(G), b'>>graph6<<A_\n')
        G = nx.Graph([(1, 42)])
        assert_equal(g6.to_graph6_bytes(G), b'>>graph6<<A_\n')
