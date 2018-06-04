#!/usr/bin/env python
from nose.tools import assert_equal
from nose.tools import assert_is
from nose.tools import assert_not_equal
from nose.tools import assert_raises

import networkx as nx
from networkx.testing.utils import *

from test_graph import BaseAttrGraphTester, TestGraph


class BaseMultiGraphTester(BaseAttrGraphTester):
    def test_has_edge(self):
        G = self.K3
        assert_equal(G.has_edge(0, 1), True)
        assert_equal(G.has_edge(0, -1), False)
        assert_equal(G.has_edge(0, 1, 0), True)
        assert_equal(G.has_edge(0, 1, 1), False)

    def test_get_edge_data(self):
        G = self.K3
        assert_equal(G.get_edge_data(0, 1), {0: {}})
        assert_equal(G[0][1], {0: {}})
        assert_equal(G[0][1][0], {})
        assert_equal(G.get_edge_data(10, 20), None)
        assert_equal(G.get_edge_data(0, 1, 0), {})

    def test_adjacency(self):
        G = self.K3
        assert_equal(dict(G.adjacency()),
                     {0: {1: {0: {}}, 2: {0: {}}},
                      1: {0: {0: {}}, 2: {0: {}}},
                      2: {0: {0: {}}, 1: {0: {}}}})

    def deepcopy_edge_attr(self, H, G):
        assert_equal(G[1][2][0]['foo'], H[1][2][0]['foo'])
        G[1][2][0]['foo'].append(1)
        assert_not_equal(G[1][2][0]['foo'], H[1][2][0]['foo'])

    def shallow_copy_edge_attr(self, H, G):
        assert_equal(G[1][2][0]['foo'], H[1][2][0]['foo'])
        G[1][2][0]['foo'].append(1)
        assert_equal(G[1][2][0]['foo'], H[1][2][0]['foo'])

    def graphs_equal(self, H, G):
        assert_equal(G._adj, H._adj)
        assert_equal(G._node, H._node)
        assert_equal(G.graph, H.graph)
        assert_equal(G.name, H.name)
        if not G.is_directed() and not H.is_directed():
            assert_is(H._adj[1][2][0], H._adj[2][1][0])
            assert_is(G._adj[1][2][0], G._adj[2][1][0])
        else:  # at least one is directed
            if not G.is_directed():
                G._pred = G._adj
                G._succ = G._adj
            if not H.is_directed():
                H._pred = H._adj
                H._succ = H._adj
            assert_equal(G._pred, H._pred)
            assert_equal(G._succ, H._succ)
            assert_is(H._succ[1][2][0], H._pred[2][1][0])
            assert_is(G._succ[1][2][0], G._pred[2][1][0])

    def same_attrdict(self, H, G):
        # same attrdict in the edgedata
        old_foo = H[1][2][0]['foo']
        H.adj[1][2][0]['foo'] = 'baz'
        assert_equal(G._adj, H._adj)
        H.adj[1][2][0]['foo'] = old_foo
        assert_equal(G._adj, H._adj)

        old_foo = H.nodes[0]['foo']
        H.nodes[0]['foo'] = 'baz'
        assert_equal(G._node, H._node)
        H.nodes[0]['foo'] = old_foo
        assert_equal(G._node, H._node)

    def different_attrdict(self, H, G):
        # used by graph_equal_but_different
        old_foo = H[1][2][0]['foo']
        H.adj[1][2][0]['foo'] = 'baz'
        assert_not_equal(G._adj, H._adj)
        H.adj[1][2][0]['foo'] = old_foo
        assert_equal(G._adj, H._adj)

        old_foo = H.nodes[0]['foo']
        H.nodes[0]['foo'] = 'baz'
        assert_not_equal(G._node, H._node)
        H.nodes[0]['foo'] = old_foo
        assert_equal(G._node, H._node)

    def test_to_undirected(self):
        G = self.K3
        self.add_attributes(G)
        H = nx.MultiGraph(G)
        self.is_shallow_copy(H, G)
        H = G.to_undirected()
        self.is_deepcopy(H, G)

    def test_to_directed(self):
        G = self.K3
        self.add_attributes(G)
        H = nx.MultiDiGraph(G)
        self.is_shallow_copy(H, G)
        H = G.to_directed()
        self.is_deepcopy(H, G)

    def test_number_of_edges_selfloops(self):
        G = self.K3
        G.add_edge(0, 0)
        G.add_edge(0, 0)
        G.add_edge(0, 0, key='parallel edge')
        G.remove_edge(0, 0, key='parallel edge')
        assert_equal(G.number_of_edges(0, 0), 2)
        G.remove_edge(0, 0)
        assert_equal(G.number_of_edges(0, 0), 1)

    def test_edge_lookup(self):
        G = self.Graph()
        G.add_edge(1, 2, foo='bar')
        G.add_edge(1, 2, 'key', foo='biz')
        assert_edges_equal(G.edges[1, 2, 0], {'foo': 'bar'})
        assert_edges_equal(G.edges[1, 2, 'key'], {'foo': 'biz'})

    def test_edge_attr4(self):
        G = self.Graph()
        G.add_edge(1, 2, key=0, data=7, spam='bar', bar='foo')
        assert_edges_equal(G.edges(data=True),
                           [(1, 2, {'data': 7, 'spam': 'bar', 'bar': 'foo'})])
        G[1][2][0]['data'] = 10  # OK to set data like this
        assert_edges_equal(G.edges(data=True),
                           [(1, 2, {'data': 10, 'spam': 'bar', 'bar': 'foo'})])

        G.adj[1][2][0]['data'] = 20
        assert_edges_equal(G.edges(data=True),
                           [(1, 2, {'data': 20, 'spam': 'bar', 'bar': 'foo'})])
        G.edges[1, 2, 0]['data'] = 21  # another spelling, "edge"
        assert_edges_equal(G.edges(data=True),
                           [(1, 2, {'data': 21, 'spam': 'bar', 'bar': 'foo'})])
        G.adj[1][2][0]['listdata'] = [20, 200]
        G.adj[1][2][0]['weight'] = 20
        assert_edges_equal(G.edges(data=True),
                           [(1, 2, {'data': 21, 'spam': 'bar', 'bar': 'foo',
                                    'listdata': [20, 200], 'weight':20})])


class TestMultiGraph(BaseMultiGraphTester, TestGraph):
    def setUp(self):
        self.Graph = nx.MultiGraph
        # build K3
        ed1, ed2, ed3 = ({0: {}}, {0: {}}, {0: {}})
        self.k3adj = {0: {1: ed1, 2: ed2},
                      1: {0: ed1, 2: ed3},
                      2: {0: ed2, 1: ed3}}
        self.k3edges = [(0, 1), (0, 2), (1, 2)]
        self.k3nodes = [0, 1, 2]
        self.K3 = self.Graph()
        self.K3._adj = self.k3adj
        self.K3._node = {}
        self.K3._node[0] = {}
        self.K3._node[1] = {}
        self.K3._node[2] = {}

    def test_data_input(self):
        G = self.Graph({1: [2], 2: [1]}, name="test")
        assert_equal(G.name, "test")
        expected = [(1, {2: {0: {}}}), (2, {1: {0: {}}})]
        assert_equal(sorted(G.adj.items()), expected)

    def test_getitem(self):
        G = self.K3
        assert_equal(G[0], {1: {0: {}}, 2: {0: {}}})
        assert_raises(KeyError, G.__getitem__, 'j')
        assert_raises((TypeError, nx.NetworkXError), G.__getitem__, ['A'])

    def test_remove_node(self):
        G = self.K3
        G.remove_node(0)
        assert_equal(G.adj, {1: {2: {0: {}}}, 2: {1: {0: {}}}})
        assert_raises((KeyError, nx.NetworkXError), G.remove_node, -1)

    def test_add_edge(self):
        G = self.Graph()
        G.add_edge(0, 1)
        assert_equal(G.adj, {0: {1: {0: {}}}, 1: {0: {0: {}}}})
        G = self.Graph()
        G.add_edge(*(0, 1))
        assert_equal(G.adj, {0: {1: {0: {}}}, 1: {0: {0: {}}}})

    def test_add_edge_conflicting_key(self):
        G = self.Graph()
        G.add_edge(0, 1, key=1)
        G.add_edge(0, 1)
        assert_equal(G.number_of_edges(), 2)
        G = self.Graph()
        G.add_edges_from([(0, 1, 1, {})])
        G.add_edges_from([(0, 1)])
        assert_equal(G.number_of_edges(), 2)

    def test_add_edges_from(self):
        G = self.Graph()
        G.add_edges_from([(0, 1), (0, 1, {'weight': 3})])
        assert_equal(G.adj, {0: {1: {0: {}, 1: {'weight': 3}}},
                             1: {0: {0: {}, 1: {'weight': 3}}}})
        G.add_edges_from([(0, 1), (0, 1, {'weight': 3})], weight=2)
        assert_equal(G.adj, {0: {1: {0: {}, 1: {'weight': 3},
                                     2: {'weight': 2}, 3: {'weight': 3}}},
                             1: {0: {0: {}, 1: {'weight': 3},
                                     2: {'weight': 2}, 3: {'weight': 3}}}})
        G = self.Graph()
        edges = [(0, 1, {'weight': 3}), (0, 1, (('weight', 2),)),
                 (0, 1, 5), (0, 1, 's')]
        G.add_edges_from(edges)
        keydict = {0: {'weight': 3}, 1: {'weight': 2}, 5: {}, 's': {}}
        assert_equal(G._adj, {0: {1: keydict}, 1: {0: keydict}})

        # too few in tuple
        assert_raises(nx.NetworkXError, G.add_edges_from, [(0,)])
        # too many in tuple
        assert_raises(nx.NetworkXError, G.add_edges_from, [(0, 1, 2, 3, 4)])
        # not a tuple
        assert_raises(TypeError, G.add_edges_from, [0])

    def test_remove_edge(self):
        G = self.K3
        G.remove_edge(0, 1)
        assert_equal(G.adj, {0: {2: {0: {}}},
                             1: {2: {0: {}}},
                             2: {0: {0: {}},
                                 1: {0: {}}}})

        assert_raises((KeyError, nx.NetworkXError), G.remove_edge, -1, 0)
        assert_raises((KeyError, nx.NetworkXError), G.remove_edge, 0, 2,
                      key=1)

    def test_remove_edges_from(self):
        G = self.K3.copy()
        G.remove_edges_from([(0, 1)])
        kd = {0: {}}
        assert_equal(G.adj, {0: {2: kd}, 1: {2: kd}, 2: {0: kd, 1: kd}})
        G.remove_edges_from([(0, 0)])  # silent fail
        self.K3.add_edge(0, 1)
        G = self.K3.copy()
        G.remove_edges_from(list(G.edges(data=True, keys=True)))
        assert_equal(G.adj, {0: {}, 1: {}, 2: {}})
        G = self.K3.copy()
        G.remove_edges_from(list(G.edges(data=False, keys=True)))
        assert_equal(G.adj, {0: {}, 1: {}, 2: {}})
        G = self.K3.copy()
        G.remove_edges_from(list(G.edges(data=False, keys=False)))
        assert_equal(G.adj, {0: {}, 1: {}, 2: {}})
        G = self.K3.copy()
        G.remove_edges_from([(0, 1, 0), (0, 2, 0, {}), (1, 2)])
        assert_equal(G.adj, {0: {1: {1: {}}}, 1: {0: {1: {}}}, 2: {}})

    def test_remove_multiedge(self):
        G = self.K3
        G.add_edge(0, 1, key='parallel edge')
        G.remove_edge(0, 1, key='parallel edge')
        assert_equal(G.adj, {0: {1: {0: {}}, 2: {0: {}}},
                             1: {0: {0: {}}, 2: {0: {}}},
                             2: {0: {0: {}}, 1: {0: {}}}})
        G.remove_edge(0, 1)
        kd = {0: {}}
        assert_equal(G.adj, {0: {2: kd}, 1: {2: kd}, 2: {0: kd, 1: kd}})
        assert_raises((KeyError, nx.NetworkXError), G.remove_edge, -1, 0)


class TestEdgeSubgraph(object):
    """Unit tests for the :meth:`MultiGraph.edge_subgraph` method."""

    def setup(self):
        # Create a doubly-linked path graph on five nodes.
        G = nx.MultiGraph()
        nx.add_path(G, range(5))
        nx.add_path(G, range(5))
        # Add some node, edge, and graph attributes.
        for i in range(5):
            G.nodes[i]['name'] = 'node{}'.format(i)
        G.adj[0][1][0]['name'] = 'edge010'
        G.adj[0][1][1]['name'] = 'edge011'
        G.adj[3][4][0]['name'] = 'edge340'
        G.adj[3][4][1]['name'] = 'edge341'
        G.graph['name'] = 'graph'
        # Get the subgraph induced by one of the first edges and one of
        # the last edges.
        self.G = G
        self.H = G.edge_subgraph([(0, 1, 0), (3, 4, 1)])

    def test_correct_nodes(self):
        """Tests that the subgraph has the correct nodes."""
        assert_equal([0, 1, 3, 4], sorted(self.H.nodes()))

    def test_correct_edges(self):
        """Tests that the subgraph has the correct edges."""
        assert_equal([(0, 1, 0, 'edge010'), (3, 4, 1, 'edge341')],
                     sorted(self.H.edges(keys=True, data='name')))

    def test_add_node(self):
        """Tests that adding a node to the original graph does not
        affect the nodes of the subgraph.

        """
        self.G.add_node(5)
        assert_equal([0, 1, 3, 4], sorted(self.H.nodes()))

    def test_remove_node(self):
        """Tests that removing a node in the original graph does
        affect the nodes of the subgraph.

        """
        self.G.remove_node(0)
        assert_equal([1, 3, 4], sorted(self.H.nodes()))

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
        for u, v, k in self.H.edges(keys=True):
            assert_equal(self.G._adj[u][v][k], self.H._adj[u][v][k])
        # Making a change to G should make a change in H and vice versa.
        self.G._adj[0][1][0]['name'] = 'foo'
        assert_equal(self.G._adj[0][1][0]['name'],
                     self.H._adj[0][1][0]['name'])
        self.H._adj[3][4][1]['name'] = 'bar'
        assert_equal(self.G._adj[3][4][1]['name'],
                     self.H._adj[3][4][1]['name'])

    def test_graph_attr_dict(self):
        """Tests that the graph attribute dictionary of the two graphs
        is the same object.

        """
        assert_is(self.G.graph, self.H.graph)
