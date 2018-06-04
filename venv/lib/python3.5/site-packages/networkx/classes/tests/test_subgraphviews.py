from nose.tools import assert_equal, assert_not_equal, \
    assert_is, assert_true, assert_raises

import networkx as nx


class TestSubGraphView(object):
    gview = nx.graphviews.SubGraph
    graph = nx.Graph
    hide_edges_filter = staticmethod(nx.filters.hide_edges)
    show_edges_filter = staticmethod(nx.filters.show_edges)

    def setUp(self):
        self.G = nx.path_graph(9, create_using=self.graph())
        self.hide_edges_w_hide_nodes = {(3, 4), (4, 5), (5, 6)}

    def test_hidden_nodes(self):
        hide_nodes = [4, 5, 111]
        nodes_gone = nx.filters.hide_nodes(hide_nodes)
        G = self.gview(self.G, filter_node=nodes_gone)
        assert_equal(self.G.nodes - G.nodes, {4, 5})
        assert_equal(self.G.edges - G.edges, self.hide_edges_w_hide_nodes)
        if G.is_directed():
            assert_equal(list(G[3]), [])
            assert_equal(list(G[2]), [3])
        else:
            assert_equal(list(G[3]), [2])
            assert_equal(set(G[2]), {1, 3})
        assert_raises(KeyError, G.__getitem__, 4)
        assert_raises(KeyError, G.__getitem__, 112)
        assert_raises(KeyError, G.__getitem__, 111)
        assert_equal(G.degree(3), 3 if G.is_multigraph() else 1)
        assert_equal(G.size(), 7 if G.is_multigraph() else 5)

    def test_hidden_edges(self):
        hide_edges = [(2, 3), (8, 7), (222, 223)]
        edges_gone = self.hide_edges_filter(hide_edges)
        G = self.gview(self.G, filter_edge=edges_gone)
        assert_equal(self.G.nodes, G.nodes)
        if G.is_directed():
            assert_equal(self.G.edges - G.edges, {(2, 3)})
            assert_equal(list(G[2]), [])
            assert_equal(list(G.pred[3]), [])
            assert_equal(list(G.pred[2]), [1])
            assert_equal(G.size(), 7)
        else:
            assert_equal(self.G.edges - G.edges, {(2, 3), (7, 8)})
            assert_equal(list(G[2]), [1])
            assert_equal(G.size(), 6)
        assert_equal(list(G[3]), [4])
        assert_raises(KeyError, G.__getitem__, 221)
        assert_raises(KeyError, G.__getitem__, 222)
        assert_equal(G.degree(3), 1)

    def test_shown_node(self):
        induced_subgraph = nx.filters.show_nodes([2, 3, 111])
        G = self.gview(self.G, filter_node=induced_subgraph)
        assert_equal(set(G.nodes), {2, 3})
        if G.is_directed():
            assert_equal(list(G[3]), [])
        else:
            assert_equal(list(G[3]), [2])
        assert_equal(list(G[2]), [3])
        assert_raises(KeyError, G.__getitem__, 4)
        assert_raises(KeyError, G.__getitem__, 112)
        assert_raises(KeyError, G.__getitem__, 111)
        assert_equal(G.degree(3), 3 if G.is_multigraph() else 1)
        assert_equal(G.size(), 3 if G.is_multigraph() else 1)

    def test_shown_edges(self):
        show_edges = [(2, 3), (8, 7), (222, 223)]
        edge_subgraph = self.show_edges_filter(show_edges)
        G = self.gview(self.G, filter_edge=edge_subgraph)
        assert_equal(self.G.nodes, G.nodes)
        if G.is_directed():
            assert_equal(G.edges, {(2, 3)})
            assert_equal(list(G[3]), [])
            assert_equal(list(G[2]), [3])
            assert_equal(list(G.pred[3]), [2])
            assert_equal(list(G.pred[2]), [])
            assert_equal(G.size(), 1)
        else:
            assert_equal(G.edges, {(2, 3), (7, 8)})
            assert_equal(list(G[3]), [2])
            assert_equal(list(G[2]), [3])
            assert_equal(G.size(), 2)
        assert_raises(KeyError, G.__getitem__, 221)
        assert_raises(KeyError, G.__getitem__, 222)
        assert_equal(G.degree(3), 1)


class TestSubDiGraphView(TestSubGraphView):
    gview = nx.graphviews.SubDiGraph
    graph = nx.DiGraph
    hide_edges_filter = staticmethod(nx.filters.hide_diedges)
    show_edges_filter = staticmethod(nx.filters.show_diedges)
    hide_edges = [(2, 3), (8, 7), (222, 223)]
    excluded = {(2, 3), (3, 4), (4, 5), (5, 6)}

    def test_inoutedges(self):
        edges_gone = self.hide_edges_filter(self.hide_edges)
        hide_nodes = [4, 5, 111]
        nodes_gone = nx.filters.hide_nodes(hide_nodes)
        G = self.gview(self.G, nodes_gone, edges_gone)

        assert_equal(self.G.in_edges - G.in_edges, self.excluded)
        assert_equal(self.G.out_edges - G.out_edges, self.excluded)

    def test_pred(self):
        edges_gone = self.hide_edges_filter(self.hide_edges)
        hide_nodes = [4, 5, 111]
        nodes_gone = nx.filters.hide_nodes(hide_nodes)
        G = self.gview(self.G, nodes_gone, edges_gone)

        assert_equal(list(G.pred[2]), [1])
        assert_equal(list(G.pred[6]), [])

    def test_inout_degree(self):
        edges_gone = self.hide_edges_filter(self.hide_edges)
        hide_nodes = [4, 5, 111]
        nodes_gone = nx.filters.hide_nodes(hide_nodes)
        G = self.gview(self.G, nodes_gone, edges_gone)

        assert_equal(G.degree(2), 1)
        assert_equal(G.out_degree(2), 0)
        assert_equal(G.in_degree(2), 1)
        assert_equal(G.size(), 4)


# multigraph
class TestMultiGraphView(TestSubGraphView):
    gview = nx.graphviews.SubMultiGraph
    graph = nx.MultiGraph
    hide_edges_filter = staticmethod(nx.filters.hide_multiedges)
    show_edges_filter = staticmethod(nx.filters.show_multiedges)

    def setUp(self):
        self.G = nx.path_graph(9, create_using=self.graph())
        multiedges = {(2, 3, 4), (2, 3, 5)}
        self.G.add_edges_from(multiedges)
        self.hide_edges_w_hide_nodes = {(3, 4, 0), (4, 5, 0), (5, 6, 0)}

    def test_hidden_edges(self):
        hide_edges = [(2, 3, 4), (2, 3, 3), (8, 7, 0), (222, 223, 0)]
        edges_gone = self.hide_edges_filter(hide_edges)
        G = self.gview(self.G, filter_edge=edges_gone)
        assert_equal(self.G.nodes, G.nodes)
        if G.is_directed():
            assert_equal(self.G.edges - G.edges, {(2, 3, 4)})
            assert_equal(list(G[3]), [4])
            assert_equal(list(G[2]), [3])
            assert_equal(list(G.pred[3]), [2])  # only one 2 but two edges
            assert_equal(list(G.pred[2]), [1])
            assert_equal(G.size(), 9)
        else:
            assert_equal(self.G.edges - G.edges, {(2, 3, 4), (7, 8, 0)})
            assert_equal(list(G[3]), [2, 4])
            assert_equal(list(G[2]), [1, 3])
            assert_equal(G.size(), 8)
        assert_equal(G.degree(3), 3)
        assert_raises(KeyError, G.__getitem__, 221)
        assert_raises(KeyError, G.__getitem__, 222)

    def test_shown_edges(self):
        show_edges = [(2, 3, 4), (2, 3, 3), (8, 7, 0), (222, 223, 0)]
        edge_subgraph = self.show_edges_filter(show_edges)
        G = self.gview(self.G, filter_edge=edge_subgraph)
        assert_equal(self.G.nodes, G.nodes)
        if G.is_directed():
            assert_equal(G.edges, {(2, 3, 4)})
            assert_equal(list(G[3]), [])
            assert_equal(list(G.pred[3]), [2])
            assert_equal(list(G.pred[2]), [])
            assert_equal(G.size(), 1)
        else:
            assert_equal(G.edges, {(2, 3, 4), (7, 8, 0)})
            assert_equal(G.size(), 2)
            assert_equal(list(G[3]), [2])
        assert_equal(G.degree(3), 1)
        assert_equal(list(G[2]), [3])
        assert_raises(KeyError, G.__getitem__, 221)
        assert_raises(KeyError, G.__getitem__, 222)


# multidigraph
class TestMultiDiGraphView(TestMultiGraphView, TestSubDiGraphView):
    gview = nx.graphviews.SubMultiDiGraph
    graph = nx.MultiDiGraph
    hide_edges_filter = staticmethod(nx.filters.hide_multidiedges)
    show_edges_filter = staticmethod(nx.filters.show_multidiedges)
    hide_edges = [(2, 3, 0), (8, 7, 0), (222, 223, 0)]
    excluded = {(2, 3, 0), (3, 4, 0), (4, 5, 0), (5, 6, 0)}

    def test_inout_degree(self):
        edges_gone = self.hide_edges_filter(self.hide_edges)
        hide_nodes = [4, 5, 111]
        nodes_gone = nx.filters.hide_nodes(hide_nodes)
        G = self.gview(self.G, nodes_gone, edges_gone)

        assert_equal(G.degree(2), 3)
        assert_equal(G.out_degree(2), 2)
        assert_equal(G.in_degree(2), 1)
        assert_equal(G.size(), 6)


# induced_subgraph
class TestInducedSubGraph(object):
    def setUp(self):
        self.K3 = G = nx.complete_graph(3)
        G.graph['foo'] = []
        G.nodes[0]['foo'] = []
        G.remove_edge(1, 2)
        ll = []
        G.add_edge(1, 2, foo=ll)
        G.add_edge(2, 1, foo=ll)

    def test_full_graph(self):
        G = self.K3
        H = nx.induced_subgraph(G, [0, 1, 2, 5])
        assert_equal(H.name, G.name)
        self.graphs_equal(H, G)
        self.same_attrdict(H, G)

    def test_partial_subgraph(self):
        G = self.K3
        H = nx.induced_subgraph(G, 0)
        assert_equal(dict(H.adj), {0: {}})
        assert_not_equal(dict(G.adj), {0: {}})

        H = nx.induced_subgraph(G, [0, 1])
        assert_equal(dict(H.adj), {0: {1: {}}, 1: {0: {}}})

    def same_attrdict(self, H, G):
        old_foo = H[1][2]['foo']
        H.edges[1, 2]['foo'] = 'baz'
        assert_equal(G.edges, H.edges)
        H.edges[1, 2]['foo'] = old_foo
        assert_equal(G.edges, H.edges)
        old_foo = H.nodes[0]['foo']
        H.nodes[0]['foo'] = 'baz'
        assert_equal(G.nodes, H.nodes)
        H.nodes[0]['foo'] = old_foo
        assert_equal(G.nodes, H.nodes)

    def graphs_equal(self, H, G):
        assert_equal(G._adj, H._adj)
        assert_equal(G._node, H._node)
        assert_equal(G.graph, H.graph)
        assert_equal(G.name, H.name)
        if not G.is_directed() and not H.is_directed():
            assert_true(H._adj[1][2] is H._adj[2][1])
            assert_true(G._adj[1][2] is G._adj[2][1])
        else:  # at least one is directed
            if not G.is_directed():
                G._pred = G._adj
                G._succ = G._adj
            if not H.is_directed():
                H._pred = H._adj
                H._succ = H._adj
            assert_equal(G._pred, H._pred)
            assert_equal(G._succ, H._succ)
            assert_true(H._succ[1][2] is H._pred[2][1])
            assert_true(G._succ[1][2] is G._pred[2][1])


# edge_subgraph
class TestEdgeSubGraph(object):
    def setup(self):
        # Create a path graph on five nodes.
        self.G = G = nx.path_graph(5)
        # Add some node, edge, and graph attributes.
        for i in range(5):
            G.nodes[i]['name'] = 'node{}'.format(i)
        G.edges[0, 1]['name'] = 'edge01'
        G.edges[3, 4]['name'] = 'edge34'
        G.graph['name'] = 'graph'
        # Get the subgraph induced by the first and last edges.
        self.H = nx.edge_subgraph(G, [(0, 1), (3, 4)])

    def test_correct_nodes(self):
        """Tests that the subgraph has the correct nodes."""
        assert_equal([0, 1, 3, 4], sorted(self.H.nodes))

    def test_correct_edges(self):
        """Tests that the subgraph has the correct edges."""
        assert_equal([(0, 1, 'edge01'), (3, 4, 'edge34')],
                     sorted(self.H.edges(data='name')))

    def test_add_node(self):
        """Tests that adding a node to the original graph does not
        affect the nodes of the subgraph.

        """
        self.G.add_node(5)
        assert_equal([0, 1, 3, 4], sorted(self.H.nodes))
        self.G.remove_node(5)

    def test_remove_node(self):
        """Tests that removing a node in the original graph
        removes the nodes of the subgraph.

        """
        self.G.remove_node(0)
        assert_equal([1, 3, 4], sorted(self.H.nodes))
        self.G.add_edge(0, 1)

    def test_node_attr_dict(self):
        """Tests that the node attribute dictionary of the two graphs is
        the same object.

        """
        for v in self.H:
            assert_equal(self.G.nodes[v], self.H.nodes[v])
        # Making a change to G should make a change in H and vice versa.
        self.G.nodes[0]['name'] = 'foo'
        assert_equal(self.G.nodes[0], self.H.nodes[0])
        self.H.nodes[1]['name'] = 'bar'
        assert_equal(self.G.nodes[1], self.H.nodes[1])

    def test_edge_attr_dict(self):
        """Tests that the edge attribute dictionary of the two graphs is
        the same object.

        """
        for u, v in self.H.edges():
            assert_equal(self.G.edges[u, v], self.H.edges[u, v])
        # Making a change to G should make a change in H and vice versa.
        self.G.edges[0, 1]['name'] = 'foo'
        assert_equal(self.G.edges[0, 1]['name'],
                     self.H.edges[0, 1]['name'])
        self.H.edges[3, 4]['name'] = 'bar'
        assert_equal(self.G.edges[3, 4]['name'],
                     self.H.edges[3, 4]['name'])

    def test_graph_attr_dict(self):
        """Tests that the graph attribute dictionary of the two graphs
        is the same object.

        """
        assert_is(self.G.graph, self.H.graph)

    def test_readonly(self):
        """Tests that the subgraph cannot change the graph structure"""
        assert_raises(nx.NetworkXError, self.H.add_node, 5)
        assert_raises(nx.NetworkXError, self.H.remove_node, 0)
        assert_raises(nx.NetworkXError, self.H.add_edge, 5, 6)
        assert_raises(nx.NetworkXError, self.H.remove_edge, 0, 1)
