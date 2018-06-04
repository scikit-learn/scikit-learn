"""Unit tests for pydot drawing functions."""
try:
    try:
        from cStringIO import StringIO
    except ImportError:
        from StringIO import StringIO
except ImportError:
    from io import StringIO
import sys
import tempfile
from nose.tools import assert_equal, assert_is_instance, assert_true
import networkx as nx
from networkx.testing import assert_graphs_equal


class TestPydot(object):
    @classmethod
    def setupClass(cls):
        '''
        Fixture defining the `pydot` global to be the `pydot` module if both
        importable and of sufficient version _or_ skipping this test.
        '''
        global pydot
        pydot = nx.nx_pydot.setup_module(sys.modules[__name__])
        assert pydot is not None

    def pydot_checks(self, G, prog):
        '''
        Validate :mod:`pydot`-based usage of the passed NetworkX graph with the
        passed basename of an external GraphViz command (e.g., `dot`, `neato`).
        '''

        # Set the name of this graph to... "G". Failing to do so will
        # subsequently trip an assertion expecting this name.
        G.graph['name'] = 'G'

        # Add arbitrary nodes and edges to the passed empty graph.
        G.add_edges_from([('A', 'B'), ('A', 'C'), ('B', 'C'), ('A', 'D')])
        G.add_node('E')

        # Validate layout of this graph with the passed GraphViz command.
        graph_layout = nx.nx_pydot.pydot_layout(G, prog=prog)
        assert_is_instance(graph_layout, dict)

        # Convert this graph into a "pydot.Dot" instance.
        P = nx.nx_pydot.to_pydot(G)

        # Convert this "pydot.Dot" instance back into a graph of the same type.
        G2 = G.fresh_copy().__class__(nx.nx_pydot.from_pydot(P))

        # Validate the original and resulting graphs to be the same.
        assert_graphs_equal(G, G2)

        # Serialize this "pydot.Dot" instance to a temporary file in dot format.
        fname = tempfile.mktemp()
        P.write_raw(fname)

        # Deserialize a list of new "pydot.Dot" instances back from this file.
        Pin_list = pydot.graph_from_dot_file(path=fname, encoding='utf-8')

        # Validate this file to contain only one graph.
        assert_equal(len(Pin_list), 1)

        # The single "pydot.Dot" instance deserialized from this file.
        Pin = Pin_list[0]

        # Sorted list of all nodes in the original "pydot.Dot" instance.
        n1 = sorted([p.get_name() for p in P.get_node_list()])

        # Sorted list of all nodes in the deserialized "pydot.Dot" instance.
        n2 = sorted([p.get_name() for p in Pin.get_node_list()])

        # Validate these instances to contain the same nodes.
        assert_equal(n1, n2)

        # Sorted list of all edges in the original "pydot.Dot" instance.
        e1 = sorted([
            (e.get_source(), e.get_destination()) for e in P.get_edge_list()])

        # Sorted list of all edges in the original "pydot.Dot" instance.
        e2 = sorted([
            (e.get_source(), e.get_destination()) for e in Pin.get_edge_list()])

        # Validate these instances to contain the same edges.
        assert_equal(e1, e2)

        # Deserialize a new graph of the same type back from this file.
        Hin = nx.nx_pydot.read_dot(fname)
        Hin = G.fresh_copy().__class__(Hin)

        # Validate the original and resulting graphs to be the same.
        assert_graphs_equal(G, Hin)

    def test_undirected(self):
        self.pydot_checks(nx.Graph(), prog='neato')

    def test_directed(self):
        self.pydot_checks(nx.DiGraph(), prog='dot')

    def test_read_write(self):
        G = nx.MultiGraph()
        G.graph['name'] = 'G'
        G.add_edge('1', '2', key='0')  # read assumes strings
        fh = StringIO()
        nx.nx_pydot.write_dot(G, fh)
        fh.seek(0)
        H = nx.nx_pydot.read_dot(fh)
        assert_graphs_equal(G, H)
