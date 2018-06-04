from nose.tools import assert_in, assert_not_in, assert_equal
from nose.tools import assert_is, assert_is_not
from nose.tools import assert_raises, assert_true, assert_false

import networkx as nx
from networkx.testing import assert_edges_equal, assert_nodes_equal

# Note: SubGraph views are not tested here. They have their own testing file


class TestReverseView(object):
    def setup(self):
        self.G = nx.path_graph(9, create_using=nx.DiGraph())
        self.rv = nx.reverse_view(self.G)

    def test_pickle(self):
        import pickle
        rv = self.rv
        prv = pickle.loads(pickle.dumps(rv, -1))
        assert_equal(rv._node, prv._node)
        assert_equal(rv._adj, prv._adj)
        assert_equal(rv.graph, prv.graph)

    def test_contains(self):
        assert_in((2, 3), self.G.edges)
        assert_not_in((3, 2), self.G.edges)
        assert_not_in((2, 3), self.rv.edges)
        assert_in((3, 2), self.rv.edges)

    def test_iter(self):
        expected = sorted(tuple(reversed(e)) for e in self.G.edges)
        assert_equal(sorted(self.rv.edges), expected)

    def test_exceptions(self):
        nxg = nx.graphviews
        assert_raises(nx.NetworkXNotImplemented, nxg.ReverseView, nx.Graph())


class TestMultiReverseView(object):
    def setup(self):
        self.G = nx.path_graph(9, create_using=nx.MultiDiGraph())
        self.G.add_edge(4, 5)
        self.rv = nx.reverse_view(self.G)

    def test_pickle(self):
        import pickle
        rv = self.rv
        prv = pickle.loads(pickle.dumps(rv, -1))
        assert_equal(rv._node, prv._node)
        assert_equal(rv._adj, prv._adj)
        assert_equal(rv.graph, prv.graph)

    def test_contains(self):
        assert_in((2, 3, 0), self.G.edges)
        assert_not_in((3, 2, 0), self.G.edges)
        assert_not_in((2, 3, 0), self.rv.edges)
        assert_in((3, 2, 0), self.rv.edges)
        assert_in((5, 4, 1), self.rv.edges)
        assert_not_in((4, 5, 1), self.rv.edges)

    def test_iter(self):
        expected = sorted((v, u, k) for u, v, k in self.G.edges)
        assert_equal(sorted(self.rv.edges), expected)

    def test_exceptions(self):
        nxg = nx.graphviews
        MG = nx.MultiGraph(self.G)
        assert_raises(nx.NetworkXNotImplemented, nxg.MultiReverseView, MG)


class TestToDirected(object):
    def setup(self):
        self.G = nx.path_graph(9)
        self.dv = nx.to_directed(self.G)
        self.MG = nx.path_graph(9, create_using=nx.MultiGraph())
        self.Mdv = nx.to_directed(self.MG)

    def test_directed(self):
        assert_false(self.G.is_directed())
        assert_true(self.dv.is_directed())

    def test_already_directed(self):
        dd = nx.to_directed(self.dv)
        Mdd = nx.to_directed(self.Mdv)
        assert_edges_equal(dd.edges, self.dv.edges)
        assert_edges_equal(Mdd.edges, self.Mdv.edges)

    def test_pickle(self):
        import pickle
        dv = self.dv
        pdv = pickle.loads(pickle.dumps(dv, -1))
        assert_equal(dv._node, pdv._node)
        assert_equal(dv._succ, pdv._succ)
        assert_equal(dv._pred, pdv._pred)
        assert_equal(dv.graph, pdv.graph)

    def test_contains(self):
        assert_in((2, 3), self.G.edges)
        assert_in((3, 2), self.G.edges)
        assert_in((2, 3), self.dv.edges)
        assert_in((3, 2), self.dv.edges)

    def test_iter(self):
        revd = [tuple(reversed(e)) for e in self.G.edges]
        expected = sorted(list(self.G.edges) + revd)
        assert_equal(sorted(self.dv.edges), expected)

    def test_exceptions(self):
        nxg = nx.graphviews
        assert_raises(nx.NetworkXError, nxg.DiGraphView, self.MG)
        assert_raises(nx.NetworkXError, nxg.MultiDiGraphView, self.G)


class TestToUndirected(object):
    def setup(self):
        self.DG = nx.path_graph(9, create_using=nx.DiGraph())
        self.uv = nx.to_undirected(self.DG)
        self.MDG = nx.path_graph(9, create_using=nx.MultiDiGraph())
        self.Muv = nx.to_undirected(self.MDG)

    def test_directed(self):
        assert_true(self.DG.is_directed())
        assert_false(self.uv.is_directed())

    def test_already_directed(self):
        uu = nx.to_undirected(self.uv)
        Muu = nx.to_undirected(self.Muv)
        assert_edges_equal(uu.edges, self.uv.edges)
        assert_edges_equal(Muu.edges, self.Muv.edges)

    def test_pickle(self):
        import pickle
        uv = self.uv
        puv = pickle.loads(pickle.dumps(uv, -1))
        assert_equal(uv._node, puv._node)
        assert_equal(uv._adj, puv._adj)
        assert_equal(uv.graph, puv.graph)
        assert_true(hasattr(uv, '_graph'))

    def test_contains(self):
        assert_in((2, 3), self.DG.edges)
        assert_not_in((3, 2), self.DG.edges)
        assert_in((2, 3), self.uv.edges)
        assert_in((3, 2), self.uv.edges)

    def test_iter(self):
        expected = sorted(self.DG.edges)
        assert_equal(sorted(self.uv.edges), expected)

    def test_exceptions(self):
        nxg = nx.graphviews
        assert_raises(nx.NetworkXError, nxg.GraphView, self.MDG)
        assert_raises(nx.NetworkXError, nxg.MultiGraphView, self.DG)


class TestChainsOfViews(object):
    def setUp(self):
        self.G = nx.path_graph(9)
        self.DG = nx.path_graph(9, create_using=nx.DiGraph())
        self.MG = nx.path_graph(9, create_using=nx.MultiGraph())
        self.MDG = nx.path_graph(9, create_using=nx.MultiDiGraph())
        self.Gv = nx.to_undirected(self.DG)
        self.DGv = nx.to_directed(self.G)
        self.MGv = nx.to_undirected(self.MDG)
        self.MDGv = nx.to_directed(self.MG)
        self.Rv = self.DG.reverse()
        self.MRv = self.MDG.reverse()
        self.graphs = [self.G, self.DG, self.MG, self.MDG,
                       self.Gv, self.DGv, self.MGv, self.MDGv,
                       self.Rv, self.MRv]
        for G in self.graphs:
            G.edges, G.nodes, G.degree

    def test_pickle(self):
        import pickle
        for G in self.graphs:
            H = pickle.loads(pickle.dumps(G, -1))
            assert_edges_equal(H.edges, G.edges)
            assert_nodes_equal(H.nodes, G.nodes)

    def test_subgraph_of_subgraph(self):
        SGv = nx.subgraph(self.G, range(3, 7))
        SDGv = nx.subgraph(self.DG, range(3, 7))
        SMGv = nx.subgraph(self.MG, range(3, 7))
        SMDGv = nx.subgraph(self.MDG, range(3, 7))
        for G in self.graphs + [SGv, SDGv, SMGv, SMDGv]:
            SG = nx.induced_subgraph(G, [4, 5, 6])
            assert_equal(list(SG), [4, 5, 6])
            SSG = SG.subgraph([6, 7])
            assert_equal(list(SSG), [6])
            # subgraph-subgraph chain is short-cut in base class method
            assert_is(SSG._graph, G)

    def test_restricted_induced_subgraph_chains(self):
        """ Test subgraph chains that both restrict and show nodes/edges.

        A restricted_view subgraph should allow induced subgraphs using
        G.subgraph that automagically without a chain (meaning the result
        is a subgraph view of the original graph not a subgraph-of-subgraph.
        """
        hide_nodes = [3, 4, 5]
        hide_edges = [(6, 7)]
        RG = nx.restricted_view(self.G, hide_nodes, hide_edges)
        nodes = [4, 5, 6, 7, 8]
        SG = nx.induced_subgraph(RG, nodes)
        SSG = RG.subgraph(nodes)
        assert_is(SSG.root_graph, SSG._graph)
        assert_is_not(SG.root_graph, SG._graph)
        assert_edges_equal(SG.edges, SSG.edges)
        # should be same as morphing the graph
        CG = self.G.copy()
        CG.remove_nodes_from(hide_nodes)
        CG.remove_edges_from(hide_edges)
        assert_edges_equal(CG.edges(nodes), SSG.edges)
        CG.remove_nodes_from([0, 1, 2, 3])
        assert_edges_equal(CG.edges, SSG.edges)
        # switch order: subgraph first, then restricted view
        SSSG = self.G.subgraph(nodes)
        RSG = nx.restricted_view(SSSG, hide_nodes, hide_edges)
        assert_is_not(RSG.root_graph, RSG._graph)
        assert_edges_equal(RSG.edges, CG.edges)

    def test_subgraph_todirected(self):
        SG = nx.induced_subgraph(self.G, [4, 5, 6])
        SSG = SG.to_directed()
        assert_equal(sorted(SSG), [4, 5, 6])
        assert_equal(sorted(SSG.edges), [(4, 5), (5, 4), (5, 6), (6, 5)])

    def test_subgraph_toundirected(self):
        SG = nx.induced_subgraph(self.G, [4, 5, 6])
        SSG = SG.to_undirected()
        assert_equal(list(SSG), [4, 5, 6])
        assert_equal(sorted(SSG.edges), [(4, 5), (5, 6)])

    def test_reverse_subgraph_toundirected(self):
        G = self.DG.reverse(copy=False)
        SG = G.subgraph([4, 5, 6])
        SSG = SG.to_undirected()
        assert_equal(list(SSG), [4, 5, 6])
        assert_equal(sorted(SSG.edges), [(4, 5), (5, 6)])

    def test_reverse_reverse_copy(self):
        G = self.DG.reverse(copy=False)
        H = G.reverse(copy=True)
        assert_equal(H.nodes, self.DG.nodes)
        assert_equal(H.edges, self.DG.edges)
        G = self.MDG.reverse(copy=False)
        H = G.reverse(copy=True)
        assert_equal(H.nodes, self.MDG.nodes)
        assert_equal(H.edges, self.MDG.edges)

    def test_subgraph_edgesubgraph_toundirected(self):
        G = self.G.copy()
        SG = G.subgraph([4, 5, 6])
        SSG = SG.edge_subgraph([(4, 5), (5, 4)])
        USSG = SSG.to_undirected()
        assert_equal(list(USSG), [4, 5])
        assert_equal(sorted(USSG.edges), [(4, 5)])

    def test_copy_subgraph(self):
        G = self.G.copy()
        SG = G.subgraph([4, 5, 6])
        CSG = SG.copy(as_view=True)
        DCSG = SG.copy(as_view=False)
        assert_equal(CSG.__class__.__name__, 'GraphView')
        assert_equal(DCSG.__class__.__name__, 'Graph')

    def test_copy_disubgraph(self):
        G = self.DG.copy()
        SG = G.subgraph([4, 5, 6])
        CSG = SG.copy(as_view=True)
        DCSG = SG.copy(as_view=False)
        assert_equal(CSG.__class__.__name__, 'DiGraphView')
        assert_equal(DCSG.__class__.__name__, 'DiGraph')

    def test_copy_multidisubgraph(self):
        G = self.MDG.copy()
        SG = G.subgraph([4, 5, 6])
        CSG = SG.copy(as_view=True)
        DCSG = SG.copy(as_view=False)
        assert_equal(CSG.__class__.__name__, 'MultiDiGraphView')
        assert_equal(DCSG.__class__.__name__, 'MultiDiGraph')

    def test_copy_multisubgraph(self):
        G = self.MGv.copy()
        SG = G.subgraph([4, 5, 6])
        CSG = SG.copy(as_view=True)
        DCSG = SG.copy(as_view=False)
        assert_equal(CSG.__class__.__name__, 'MultiGraphView')
        assert_equal(DCSG.__class__.__name__, 'MultiGraph')
