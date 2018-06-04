from nose.tools import assert_equal
from nose.tools import assert_is
from nose.tools import assert_not_equal
from nose.tools import assert_raises
from nose.tools import raises

import networkx as nx
from networkx.testing.utils import *


def test_deprecated():
    # for backwards compatibility with 1.x, will be removed for 3.x
    G = nx.complete_graph(3)
    assert_equal(G.node, {0: {}, 1: {}, 2: {}})

    G = nx.DiGraph()
    G.add_path([3, 4])
    assert_equal(G.adj, {3: {4: {}}, 4: {}})

    G = nx.DiGraph()
    G.add_cycle([3, 4, 5])
    assert_equal(G.adj, {3: {4: {}}, 4: {5: {}}, 5: {3: {}}})

    G = nx.DiGraph()
    G.add_star([3, 4, 5])
    assert_equal(G.adj, {3: {4: {}, 5: {}}, 4: {}, 5: {}})

    G = nx.DiGraph([(0, 0), (0, 1), (1, 2)])
    assert_equal(G.number_of_selfloops(), 1)
    assert_equal(list(G.nodes_with_selfloops()), [0])
    assert_equal(list(G.selfloop_edges()), [(0, 0)])


class BaseGraphTester(object):
    """ Tests for data-structure independent graph class features."""

    def test_contains(self):
        G = self.K3
        assert(1 in G)
        assert(4 not in G)
        assert('b' not in G)
        assert([] not in G)   # no exception for nonhashable
        assert({1: 1} not in G)  # no exception for nonhashable

    def test_order(self):
        G = self.K3
        assert_equal(len(G), 3)
        assert_equal(G.order(), 3)
        assert_equal(G.number_of_nodes(), 3)

    def test_nodes(self):
        G = self.K3
        assert_equal(sorted(G.nodes()), self.k3nodes)
        assert_equal(sorted(G.nodes(data=True)), [(0, {}), (1, {}), (2, {})])

    def test_has_node(self):
        G = self.K3
        assert(G.has_node(1))
        assert(not G.has_node(4))
        assert(not G.has_node([]))   # no exception for nonhashable
        assert(not G.has_node({1: 1}))  # no exception for nonhashable

    def test_has_edge(self):
        G = self.K3
        assert_equal(G.has_edge(0, 1), True)
        assert_equal(G.has_edge(0, -1), False)

    def test_neighbors(self):
        G = self.K3
        assert_equal(sorted(G.neighbors(0)), [1, 2])
        assert_raises((KeyError, nx.NetworkXError), G.neighbors, -1)

    def test_edges(self):
        G = self.K3
        assert_edges_equal(G.edges(), [(0, 1), (0, 2), (1, 2)])
        assert_edges_equal(G.edges(0), [(0, 1), (0, 2)])
        assert_edges_equal(G.edges([0, 1]), [(0, 1), (0, 2), (1, 2)])
        assert_raises((KeyError, nx.NetworkXError), G.edges, -1)

    def test_weighted_degree(self):
        G = self.Graph()
        G.add_edge(1, 2, weight=2)
        G.add_edge(2, 3, weight=3)
        assert_equal(sorted(d for n, d in G.degree(weight='weight')), [2, 3, 5])
        assert_equal(dict(G.degree(weight='weight')), {1: 2, 2: 5, 3: 3})
        assert_equal(G.degree(1, weight='weight'), 2)
        assert_equal(G.degree([1], weight='weight'), [(1, 2)])

    def test_degree(self):
        G = self.K3
        assert_equal(sorted(G.degree()), [(0, 2), (1, 2), (2, 2)])
        assert_equal(dict(G.degree()), {0: 2, 1: 2, 2: 2})
        assert_equal(G.degree(0), 2)
        assert_raises(nx.NetworkXError, G.degree, -1)  # node not in graph

    def test_size(self):
        G = self.K3
        assert_equal(G.size(), 3)
        assert_equal(G.number_of_edges(), 3)

    def test_nbunch_iter(self):
        G = self.K3
        assert_nodes_equal(G.nbunch_iter(), self.k3nodes)  # all nodes
        assert_nodes_equal(G.nbunch_iter(0), [0])  # single node
        assert_nodes_equal(G.nbunch_iter([0, 1]), [0, 1])  # sequence
        # sequence with none in graph
        assert_nodes_equal(G.nbunch_iter([-1]), [])
        # string sequence with none in graph
        assert_nodes_equal(G.nbunch_iter("foo"), [])
        # node not in graph doesn't get caught upon creation of iterator
        bunch = G.nbunch_iter(-1)
        # but gets caught when iterator used
        assert_raises(nx.NetworkXError, list, bunch)
        # unhashable doesn't get caught upon creation of iterator
        bunch = G.nbunch_iter([0, 1, 2, {}])
        # but gets caught when iterator hits the unhashable
        assert_raises(nx.NetworkXError, list, bunch)

    @raises(nx.NetworkXError)
    def test_nbunch_iter_node_format_raise(self):
        # Tests that a node that would have failed string formatting
        # doesn't cause an error when attempting to raise a
        # :exc:`nx.NetworkXError`.

        # For more information, see pull request #1813.
        G = self.Graph()
        nbunch = [('x', set())]
        list(G.nbunch_iter(nbunch))

    def test_selfloop_degree(self):
        G = self.Graph()
        G.add_edge(1, 1)
        assert_equal(sorted(G.degree()), [(1, 2)])
        assert_equal(dict(G.degree()), {1: 2})
        assert_equal(G.degree(1), 2)
        assert_equal(sorted(G.degree([1])), [(1, 2)])
        assert_equal(G.degree(1, weight='weight'), 2)

    def test_selfloops(self):
        G = self.K3.copy()
        G.add_edge(0, 0)
        assert_nodes_equal(nx.nodes_with_selfloops(G), [0])
        assert_edges_equal(nx.selfloop_edges(G), [(0, 0)])
        assert_equal(nx.number_of_selfloops(G), 1)
        G.remove_edge(0, 0)
        G.add_edge(0, 0)
        G.remove_edges_from([(0, 0)])
        G.add_edge(1, 1)
        G.remove_node(1)
        G.add_edge(0, 0)
        G.add_edge(1, 1)
        G.remove_nodes_from([0, 1])


class BaseAttrGraphTester(BaseGraphTester):
    """ Tests of graph class attribute features."""

    def test_weighted_degree(self):
        G = self.Graph()
        G.add_edge(1, 2, weight=2, other=3)
        G.add_edge(2, 3, weight=3, other=4)
        assert_nodes_equal((d for n, d in G.degree(weight='weight')), [2, 5, 3])
        assert_equal(dict(G.degree(weight='weight')), {1: 2, 2: 5, 3: 3})
        assert_equal(G.degree(1, weight='weight'), 2)
        assert_nodes_equal((G.degree([1], weight='weight')), [(1, 2)])

        assert_nodes_equal((d for n, d in G.degree(weight='other')), [3, 7, 4])
        assert_equal(dict(G.degree(weight='other')), {1: 3, 2: 7, 3: 4})
        assert_equal(G.degree(1, weight='other'), 3)
        assert_edges_equal((G.degree([1], weight='other')), [(1, 3)])

    def add_attributes(self, G):
        G.graph['foo'] = []
        G.nodes[0]['foo'] = []
        G.remove_edge(1, 2)
        ll = []
        G.add_edge(1, 2, foo=ll)
        G.add_edge(2, 1, foo=ll)

    def test_name(self):
        G = self.Graph(name='')
        assert_equal(G.name, "")
        G = self.Graph(name='test')
        assert_equal(G.__str__(), "test")
        assert_equal(G.name, "test")

    def test_copy(self):
        G = self.Graph()
        G.add_node(0)
        G.add_edge(1, 2)
        self.add_attributes(G)
        # copy edge datadict but any container attr are same
        H = G.copy()
        self.graphs_equal(H, G)
        self.different_attrdict(H, G)
        self.shallow_copy_attrdict(H, G)

    def test_class_copy(self):
        G = self.Graph()
        G.add_node(0)
        G.add_edge(1, 2)
        self.add_attributes(G)
        # copy edge datadict but any container attr are same
        H = G.__class__(G)
        self.graphs_equal(H, G)
        self.different_attrdict(H, G)
        self.shallow_copy_attrdict(H, G)

    def test_root_graph(self):
        G = self.Graph([(0, 1), (1, 2)])
        assert_is(G, G.root_graph)
        DG = G.to_directed(as_view=True)
        SDG = DG.subgraph([0, 1])
        RSDG = SDG.reverse(copy=False)
        assert_is(G, RSDG.root_graph)

    def test_fresh_copy(self):
        G = self.Graph()
        G.add_node(0)
        G.add_edge(1, 2)
        self.add_attributes(G)
        # copy graph structure but use fresh datadict
        H = G.fresh_copy()
        H.add_nodes_from(G)
        H.add_edges_from(G.edges())
        assert_equal(len(G.nodes[0]), 1)
        ddict = G.adj[1][2][0] if G.is_multigraph() else G.adj[1][2]
        assert_equal(len(ddict), 1)
        assert_equal(len(H.nodes[0]), 0)
        ddict = H.adj[1][2][0] if H.is_multigraph() else H.adj[1][2]
        assert_equal(len(ddict), 0)

    def is_deepcopy(self, H, G):
        self.graphs_equal(H, G)
        self.different_attrdict(H, G)
        self.deep_copy_attrdict(H, G)

    def deep_copy_attrdict(self, H, G):
        self.deepcopy_graph_attr(H, G)
        self.deepcopy_node_attr(H, G)
        self.deepcopy_edge_attr(H, G)

    def deepcopy_graph_attr(self, H, G):
        assert_equal(G.graph['foo'], H.graph['foo'])
        G.graph['foo'].append(1)
        assert_not_equal(G.graph['foo'], H.graph['foo'])

    def deepcopy_node_attr(self, H, G):
        assert_equal(G.nodes[0]['foo'], H.nodes[0]['foo'])
        G.nodes[0]['foo'].append(1)
        assert_not_equal(G.nodes[0]['foo'], H.nodes[0]['foo'])

    def deepcopy_edge_attr(self, H, G):
        assert_equal(G[1][2]['foo'], H[1][2]['foo'])
        G[1][2]['foo'].append(1)
        assert_not_equal(G[1][2]['foo'], H[1][2]['foo'])

    def is_shallow_copy(self, H, G):
        self.graphs_equal(H, G)
        self.shallow_copy_attrdict(H, G)

    def shallow_copy_attrdict(self, H, G):
        self.shallow_copy_graph_attr(H, G)
        self.shallow_copy_node_attr(H, G)
        self.shallow_copy_edge_attr(H, G)

    def shallow_copy_graph_attr(self, H, G):
        assert_equal(G.graph['foo'], H.graph['foo'])
        G.graph['foo'].append(1)
        assert_equal(G.graph['foo'], H.graph['foo'])

    def shallow_copy_node_attr(self, H, G):
        assert_equal(G.nodes[0]['foo'], H.nodes[0]['foo'])
        G.nodes[0]['foo'].append(1)
        assert_equal(G.nodes[0]['foo'], H.nodes[0]['foo'])

    def shallow_copy_edge_attr(self, H, G):
        assert_equal(G[1][2]['foo'], H[1][2]['foo'])
        G[1][2]['foo'].append(1)
        assert_equal(G[1][2]['foo'], H[1][2]['foo'])

    def same_attrdict(self, H, G):
        old_foo = H[1][2]['foo']
        H.adj[1][2]['foo'] = 'baz'
        assert_equal(G.edges, H.edges)
        H.adj[1][2]['foo'] = old_foo
        assert_equal(G.edges, H.edges)

        old_foo = H.nodes[0]['foo']
        H.nodes[0]['foo'] = 'baz'
        assert_equal(G.nodes, H.nodes)
        H.nodes[0]['foo'] = old_foo
        assert_equal(G.nodes, H.nodes)

    def different_attrdict(self, H, G):
        old_foo = H[1][2]['foo']
        H.adj[1][2]['foo'] = 'baz'
        assert_not_equal(G._adj, H._adj)
        H.adj[1][2]['foo'] = old_foo
        assert_equal(G._adj, H._adj)

        old_foo = H.nodes[0]['foo']
        H.nodes[0]['foo'] = 'baz'
        assert_not_equal(G._node, H._node)
        H.nodes[0]['foo'] = old_foo
        assert_equal(G._node, H._node)

    def graphs_equal(self, H, G):
        assert_equal(G._adj, H._adj)
        assert_equal(G._node, H._node)
        assert_equal(G.graph, H.graph)
        assert_equal(G.name, H.name)
        if not G.is_directed() and not H.is_directed():
            assert_is(H._adj[1][2], H._adj[2][1])
            assert_is(G._adj[1][2], G._adj[2][1])
        else:  # at least one is directed
            if not G.is_directed():
                G._pred = G._adj
                G._succ = G._adj
            if not H.is_directed():
                H._pred = H._adj
                H._succ = H._adj
            assert_equal(G._pred, H._pred)
            assert_equal(G._succ, H._succ)
            assert_is(H._succ[1][2], H._pred[2][1])
            assert_is(G._succ[1][2], G._pred[2][1])

    def test_graph_attr(self):
        G = self.K3
        G.graph['foo'] = 'bar'
        assert_equal(G.graph['foo'], 'bar')
        del G.graph['foo']
        assert_equal(G.graph, {})
        H = self.Graph(foo='bar')
        assert_equal(H.graph['foo'], 'bar')

    def test_node_attr(self):
        G = self.K3
        G.add_node(1, foo='bar')
        assert_nodes_equal(G.nodes(), [0, 1, 2])
        assert_nodes_equal(G.nodes(data=True),
                           [(0, {}), (1, {'foo': 'bar'}), (2, {})])
        G.nodes[1]['foo'] = 'baz'
        assert_nodes_equal(G.nodes(data=True),
                           [(0, {}), (1, {'foo': 'baz'}), (2, {})])
        assert_nodes_equal(G.nodes(data='foo'),
                           [(0, None), (1, 'baz'), (2, None)])
        assert_nodes_equal(G.nodes(data='foo', default='bar'),
                           [(0, 'bar'), (1, 'baz'), (2, 'bar')])

    def test_node_attr2(self):
        G = self.K3
        a = {'foo': 'bar'}
        G.add_node(3, **a)
        assert_nodes_equal(G.nodes(), [0, 1, 2, 3])
        assert_nodes_equal(G.nodes(data=True),
                           [(0, {}), (1, {}), (2, {}), (3, {'foo': 'bar'})])

    def test_edge_lookup(self):
        G = self.Graph()
        G.add_edge(1, 2, foo='bar')
        assert_edges_equal(G.edges[1, 2], {'foo': 'bar'})

    def test_edge_attr(self):
        G = self.Graph()
        G.add_edge(1, 2, foo='bar')
        assert_edges_equal(G.edges(data=True), [(1, 2, {'foo': 'bar'})])
        assert_edges_equal(G.edges(data='foo'), [(1, 2, 'bar')])

    def test_edge_attr2(self):
        G = self.Graph()
        G.add_edges_from([(1, 2), (3, 4)], foo='foo')
        assert_edges_equal(G.edges(data=True),
                           [(1, 2, {'foo': 'foo'}), (3, 4, {'foo': 'foo'})])
        assert_edges_equal(G.edges(data='foo'),
                           [(1, 2, 'foo'), (3, 4, 'foo')])

    def test_edge_attr3(self):
        G = self.Graph()
        G.add_edges_from([(1, 2, {'weight': 32}), (3, 4, {'weight': 64})], foo='foo')
        assert_edges_equal(G.edges(data=True),
                           [(1, 2, {'foo': 'foo', 'weight': 32}),
                            (3, 4, {'foo': 'foo', 'weight': 64})])

        G.remove_edges_from([(1, 2), (3, 4)])
        G.add_edge(1, 2, data=7, spam='bar', bar='foo')
        assert_edges_equal(G.edges(data=True),
                           [(1, 2, {'data': 7, 'spam': 'bar', 'bar': 'foo'})])

    def test_edge_attr4(self):
        G = self.Graph()
        G.add_edge(1, 2, data=7, spam='bar', bar='foo')
        assert_edges_equal(G.edges(data=True),
                           [(1, 2, {'data': 7, 'spam': 'bar', 'bar': 'foo'})])
        G[1][2]['data'] = 10  # OK to set data like this
        assert_edges_equal(G.edges(data=True),
                           [(1, 2, {'data': 10, 'spam': 'bar', 'bar': 'foo'})])

        G.adj[1][2]['data'] = 20
        assert_edges_equal(G.edges(data=True),
                           [(1, 2, {'data': 20, 'spam': 'bar', 'bar': 'foo'})])
        G.edges[1, 2]['data'] = 21  # another spelling, "edge"
        assert_edges_equal(G.edges(data=True),
                           [(1, 2, {'data': 21, 'spam': 'bar', 'bar': 'foo'})])
        G.adj[1][2]['listdata'] = [20, 200]
        G.adj[1][2]['weight'] = 20
        assert_edges_equal(G.edges(data=True),
                           [(1, 2, {'data': 21, 'spam': 'bar',
                                    'bar': 'foo', 'listdata': [20, 200], 'weight':20})])

    def test_to_undirected(self):
        G = self.K3
        self.add_attributes(G)
        H = nx.Graph(G)
        self.is_shallow_copy(H, G)
        self.different_attrdict(H, G)
        H = G.to_undirected()
        self.is_deepcopy(H, G)

    def test_to_directed(self):
        G = self.K3
        self.add_attributes(G)
        H = nx.DiGraph(G)
        self.is_shallow_copy(H, G)
        self.different_attrdict(H, G)
        H = G.to_directed()
        self.is_deepcopy(H, G)

    def test_subgraph(self):
        G = self.K3
        self.add_attributes(G)
        H = G.subgraph([0, 1, 2, 5])
        self.graphs_equal(H, G)
        self.same_attrdict(H, G)
        self.shallow_copy_attrdict(H, G)

        H = G.subgraph(0)
        assert_equal(H.adj, {0: {}})
        H = G.subgraph([])
        assert_equal(H.adj, {})
        assert_not_equal(G.adj, {})

    def test_selfloops_attr(self):
        G = self.K3.copy()
        G.add_edge(0, 0)
        G.add_edge(1, 1, weight=2)
        assert_edges_equal(nx.selfloop_edges(G, data=True),
                           [(0, 0, {}), (1, 1, {'weight': 2})])
        assert_edges_equal(nx.selfloop_edges(G, data='weight'),
                           [(0, 0, None), (1, 1, 2)])


class TestGraph(BaseAttrGraphTester):
    """Tests specific to dict-of-dict-of-dict graph data structure"""

    def setUp(self):
        self.Graph = nx.Graph
        # build dict-of-dict-of-dict K3
        ed1, ed2, ed3 = ({}, {}, {})
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
        assert_equal(sorted(G.adj.items()), [(1, {2: {}}), (2, {1: {}})])
        G = self.Graph({1: [2], 2: [1]}, name="test")
        assert_equal(G.name, "test")
        assert_equal(sorted(G.adj.items()), [(1, {2: {}}), (2, {1: {}})])

    def test_adjacency(self):
        G = self.K3
        assert_equal(dict(G.adjacency()),
                     {0: {1: {}, 2: {}}, 1: {0: {}, 2: {}}, 2: {0: {}, 1: {}}})

    def test_getitem(self):
        G = self.K3
        assert_equal(G[0], {1: {}, 2: {}})
        assert_raises(KeyError, G.__getitem__, 'j')
        assert_raises((TypeError, nx.NetworkXError), G.__getitem__, ['A'])

    def test_add_node(self):
        G = self.Graph()
        G.add_node(0)
        assert_equal(G.adj, {0: {}})
        # test add attributes
        G.add_node(1, c='red')
        G.add_node(2, c='blue')
        G.add_node(3, c='red')
        assert_equal(G.nodes[1]['c'], 'red')
        assert_equal(G.nodes[2]['c'], 'blue')
        assert_equal(G.nodes[3]['c'], 'red')
        # test updating attributes
        G.add_node(1, c='blue')
        G.add_node(2, c='red')
        G.add_node(3, c='blue')
        assert_equal(G.nodes[1]['c'], 'blue')
        assert_equal(G.nodes[2]['c'], 'red')
        assert_equal(G.nodes[3]['c'], 'blue')

    def test_add_nodes_from(self):
        G = self.Graph()
        G.add_nodes_from([0, 1, 2])
        assert_equal(G.adj, {0: {}, 1: {}, 2: {}})
        # test add attributes
        G.add_nodes_from([0, 1, 2], c='red')
        assert_equal(G.nodes[0]['c'], 'red')
        assert_equal(G.nodes[2]['c'], 'red')
        # test that attribute dicts are not the same
        assert(G.nodes[0] is not G.nodes[1])
        # test updating attributes
        G.add_nodes_from([0, 1, 2], c='blue')
        assert_equal(G.nodes[0]['c'], 'blue')
        assert_equal(G.nodes[2]['c'], 'blue')
        assert(G.nodes[0] is not G.nodes[1])
        # test tuple input
        H = self.Graph()
        H.add_nodes_from(G.nodes(data=True))
        assert_equal(H.nodes[0]['c'], 'blue')
        assert_equal(H.nodes[2]['c'], 'blue')
        assert(H.nodes[0] is not H.nodes[1])
        # specific overrides general
        H.add_nodes_from([0, (1, {'c': 'green'}), (3, {'c': 'cyan'})], c='red')
        assert_equal(H.nodes[0]['c'], 'red')
        assert_equal(H.nodes[1]['c'], 'green')
        assert_equal(H.nodes[2]['c'], 'blue')
        assert_equal(H.nodes[3]['c'], 'cyan')

    def test_remove_node(self):
        G = self.K3
        G.remove_node(0)
        assert_equal(G.adj, {1: {2: {}}, 2: {1: {}}})
        assert_raises((KeyError, nx.NetworkXError), G.remove_node, -1)

        # generator here to implement list,set,string...
    def test_remove_nodes_from(self):
        G = self.K3
        G.remove_nodes_from([0, 1])
        assert_equal(G.adj, {2: {}})
        G.remove_nodes_from([-1])  # silent fail

    def test_add_edge(self):
        G = self.Graph()
        G.add_edge(0, 1)
        assert_equal(G.adj, {0: {1: {}}, 1: {0: {}}})
        G = self.Graph()
        G.add_edge(*(0, 1))
        assert_equal(G.adj, {0: {1: {}}, 1: {0: {}}})

    def test_add_edges_from(self):
        G = self.Graph()
        G.add_edges_from([(0, 1), (0, 2, {'weight': 3})])
        assert_equal(G.adj, {0: {1: {}, 2: {'weight': 3}}, 1: {0: {}},
                             2: {0: {'weight': 3}}})
        G = self.Graph()
        G.add_edges_from([(0, 1), (0, 2, {'weight': 3}), (1, 2, {'data': 4})], data=2)
        assert_equal(G.adj, {
            0: {1: {'data': 2}, 2: {'weight': 3, 'data': 2}},
            1: {0: {'data': 2}, 2: {'data': 4}},
            2: {0: {'weight': 3, 'data': 2}, 1: {'data': 4}}
        })

        assert_raises(nx.NetworkXError,
                      G.add_edges_from, [(0,)])  # too few in tuple
        assert_raises(nx.NetworkXError,
                      G.add_edges_from, [(0, 1, 2, 3)])  # too many in tuple
        assert_raises(TypeError, G.add_edges_from, [0])  # not a tuple

    def test_remove_edge(self):
        G = self.K3
        G.remove_edge(0, 1)
        assert_equal(G.adj, {0: {2: {}}, 1: {2: {}}, 2: {0: {}, 1: {}}})
        assert_raises((KeyError, nx.NetworkXError), G.remove_edge, -1, 0)

    def test_remove_edges_from(self):
        G = self.K3
        G.remove_edges_from([(0, 1)])
        assert_equal(G.adj, {0: {2: {}}, 1: {2: {}}, 2: {0: {}, 1: {}}})
        G.remove_edges_from([(0, 0)])  # silent fail

    def test_clear(self):
        G = self.K3
        G.clear()
        assert_equal(G.adj, {})

    def test_edges_data(self):
        G = self.K3
        all_edges = [(0, 1, {}), (0, 2, {}), (1, 2, {})]
        assert_edges_equal(G.edges(data=True), all_edges)
        assert_edges_equal(G.edges(0, data=True), [(0, 1, {}), (0, 2, {})])
        assert_edges_equal(G.edges([0, 1], data=True), all_edges)
        assert_raises((KeyError, nx.NetworkXError), G.edges, -1, True)

    def test_get_edge_data(self):
        G = self.K3
        assert_equal(G.get_edge_data(0, 1), {})
        assert_equal(G[0][1], {})
        assert_equal(G.get_edge_data(10, 20), None)
        assert_equal(G.get_edge_data(-1, 0), None)
        assert_equal(G.get_edge_data(-1, 0, default=1), 1)


class TestEdgeSubgraph(object):
    """Unit tests for the :meth:`Graph.edge_subgraph` method."""

    def setup(self):
        # Create a path graph on five nodes.
        G = nx.path_graph(5)
        # Add some node, edge, and graph attributes.
        for i in range(5):
            G.nodes[i]['name'] = 'node{}'.format(i)
        G.edges[0, 1]['name'] = 'edge01'
        G.edges[3, 4]['name'] = 'edge34'
        G.graph['name'] = 'graph'
        # Get the subgraph induced by the first and last edges.
        self.G = G
        self.H = G.edge_subgraph([(0, 1), (3, 4)])

    def test_correct_nodes(self):
        """Tests that the subgraph has the correct nodes."""
        assert_equal([0, 1, 3, 4], sorted(self.H.nodes()))

    def test_correct_edges(self):
        """Tests that the subgraph has the correct edges."""
        assert_equal([(0, 1, 'edge01'), (3, 4, 'edge34')],
                     sorted(self.H.edges(data='name')))

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
