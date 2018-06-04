#!/usr/bin/env python
"""Original NetworkX graph tests"""
from nose.tools import *
import networkx as nx
from networkx import convert_node_labels_to_integers as cnlti
from networkx.testing import *


class HistoricalTests(object):

    def setUp(self):
        self.null = nx.null_graph()
        self.P1 = cnlti(nx.path_graph(1), first_label=1)
        self.P3 = cnlti(nx.path_graph(3), first_label=1)
        self.P10 = cnlti(nx.path_graph(10), first_label=1)
        self.K1 = cnlti(nx.complete_graph(1), first_label=1)
        self.K3 = cnlti(nx.complete_graph(3), first_label=1)
        self.K4 = cnlti(nx.complete_graph(4), first_label=1)
        self.K5 = cnlti(nx.complete_graph(5), first_label=1)
        self.K10 = cnlti(nx.complete_graph(10), first_label=1)
        self.G = nx.Graph

    def test_name(self):
        G = self.G(name="test")
        assert_equal(str(G), 'test')
        assert_equal(G.name, 'test')
        H = self.G()
        assert_equal(H.name, '')

    # Nodes

    def test_add_remove_node(self):
        G = self.G()
        G.add_node('A')
        assert_true(G.has_node('A'))
        G.remove_node('A')
        assert_false(G.has_node('A'))

    def test_nonhashable_node(self):
        # Test if a non-hashable object is in the Graph.  A python dict will
        # raise a TypeError, but for a Graph class a simple  False should be
        # returned (see Graph __contains__). If it cannot be a node then it is
        # not a node.
        G = self.G()
        assert_false(G.has_node(['A']))
        assert_false(G.has_node({'A': 1}))

    def test_add_nodes_from(self):
        G = self.G()
        G.add_nodes_from(list("ABCDEFGHIJKL"))
        assert_true(G.has_node("L"))
        G.remove_nodes_from(['H', 'I', 'J', 'K', 'L'])
        G.add_nodes_from([1, 2, 3, 4])
        assert_equal(sorted(G.nodes(), key=str),
                     [1, 2, 3, 4, 'A', 'B', 'C', 'D', 'E', 'F', 'G'])
        # test __iter__
        assert_equal(sorted(G, key=str),
                     [1, 2, 3, 4, 'A', 'B', 'C', 'D', 'E', 'F', 'G'])

    def test_contains(self):
        G = self.G()
        G.add_node('A')
        assert_true('A' in G)
        assert_false([] in G)  # never raise a Key or TypeError in this test
        assert_false({1: 1} in G)

    def test_add_remove(self):
        # Test add_node and remove_node acting for various nbunch
        G = self.G()
        G.add_node('m')
        assert_true(G.has_node('m'))
        G.add_node('m')   # no complaints
        assert_raises(nx.NetworkXError, G.remove_node, 'j')
        G.remove_node('m')
        assert_equal(list(G), [])

    def test_nbunch_is_list(self):
        G = self.G()
        G.add_nodes_from(list("ABCD"))
        G.add_nodes_from(self.P3)  # add nbunch of nodes (nbunch=Graph)
        assert_equal(sorted(G.nodes(), key=str),
                     [1, 2, 3, 'A', 'B', 'C', 'D'])
        G.remove_nodes_from(self.P3)  # remove nbunch of nodes (nbunch=Graph)
        assert_equal(sorted(G.nodes(), key=str),
                     ['A', 'B', 'C', 'D'])

    def test_nbunch_is_set(self):
        G = self.G()
        nbunch = set("ABCDEFGHIJKL")
        G.add_nodes_from(nbunch)
        assert_true(G.has_node("L"))

    def test_nbunch_dict(self):
        # nbunch is a dict with nodes as keys
        G = self.G()
        nbunch = set("ABCDEFGHIJKL")
        G.add_nodes_from(nbunch)
        nbunch = {'I': "foo", 'J': 2, 'K': True, 'L': "spam"}
        G.remove_nodes_from(nbunch)
        assert_true(sorted(G.nodes(), key=str),
                    ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'])

    def test_nbunch_iterator(self):
        G = self.G()
        G.add_nodes_from(['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'])
        n_iter = self.P3.nodes()
        G.add_nodes_from(n_iter)
        assert_equal(sorted(G.nodes(), key=str),
                     [1, 2, 3, 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'])
        n_iter = self.P3.nodes()  # rebuild same iterator
        G.remove_nodes_from(n_iter)  # remove nbunch of nodes (nbunch=iterator)
        assert_equal(sorted(G.nodes(), key=str),
                     ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'])

    def test_nbunch_graph(self):
        G = self.G()
        G.add_nodes_from(['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'])
        nbunch = self.K3
        G.add_nodes_from(nbunch)
        assert_true(sorted(G.nodes(), key=str),
                    [1, 2, 3, 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'])

    # Edges

    def test_add_edge(self):
        G = self.G()
        assert_raises(TypeError, G.add_edge, 'A')

        G.add_edge('A', 'B')     # testing add_edge()
        G.add_edge('A', 'B')  # should fail silently
        assert_true(G.has_edge('A', 'B'))
        assert_false(G.has_edge('A', 'C'))
        assert_true(G.has_edge(*('A', 'B')))
        if G.is_directed():
            assert_false(G.has_edge('B', 'A'))
        else:
            # G is undirected, so B->A is an edge
            assert_true(G.has_edge('B', 'A'))

        G.add_edge('A', 'C')  # test directedness
        G.add_edge('C', 'A')
        G.remove_edge('C', 'A')
        if G.is_directed():
            assert_true(G.has_edge('A', 'C'))
        else:
            assert_false(G.has_edge('A', 'C'))
        assert_false(G.has_edge('C', 'A'))

    def test_self_loop(self):
        G = self.G()
        G.add_edge('A', 'A')  # test self loops
        assert_true(G.has_edge('A', 'A'))
        G.remove_edge('A', 'A')
        G.add_edge('X', 'X')
        assert_true(G.has_node('X'))
        G.remove_node('X')
        G.add_edge('A', 'Z')  # should add the node silently
        assert_true(G.has_node('Z'))

    def test_add_edges_from(self):
        G = self.G()
        G.add_edges_from([('B', 'C')])   # test add_edges_from()
        assert_true(G.has_edge('B', 'C'))
        if G.is_directed():
            assert_false(G.has_edge('C', 'B'))
        else:
            assert_true(G.has_edge('C', 'B'))  # undirected

        G.add_edges_from([('D', 'F'), ('B', 'D')])
        assert_true(G.has_edge('D', 'F'))
        assert_true(G.has_edge('B', 'D'))

        if G.is_directed():
            assert_false(G.has_edge('D', 'B'))
        else:
            assert_true(G.has_edge('D', 'B'))  # undirected

    def test_add_edges_from2(self):
        G = self.G()
        # after failing silently, should add 2nd edge
        G.add_edges_from([tuple('IJ'), list('KK'), tuple('JK')])
        assert_true(G.has_edge(*('I', 'J')))
        assert_true(G.has_edge(*('K', 'K')))
        assert_true(G.has_edge(*('J', 'K')))
        if G.is_directed():
            assert_false(G.has_edge(*('K', 'J')))
        else:
            assert_true(G.has_edge(*('K', 'J')))

    def test_add_edges_from3(self):
        G = self.G()
        G.add_edges_from(zip(list('ACD'), list('CDE')))
        assert_true(G.has_edge('D', 'E'))
        assert_false(G.has_edge('E', 'C'))

    def test_remove_edge(self):
        G = self.G()
        G.add_nodes_from([1, 2, 3, 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'])

        G.add_edges_from(zip(list('MNOP'), list('NOPM')))
        assert_true(G.has_edge('O', 'P'))
        assert_true(G.has_edge('P', 'M'))
        G.remove_node('P')    # tests remove_node()'s handling of edges.
        assert_false(G.has_edge('P', 'M'))
        assert_raises(TypeError, G.remove_edge, 'M')

        G.add_edge('N', 'M')
        assert_true(G.has_edge('M', 'N'))
        G.remove_edge('M', 'N')
        assert_false(G.has_edge('M', 'N'))

        # self loop fails silently
        G.remove_edges_from([list('HI'), list('DF'),
                             tuple('KK'), tuple('JK')])
        assert_false(G.has_edge('H', 'I'))
        assert_false(G.has_edge('J', 'K'))
        G.remove_edges_from([list('IJ'), list('KK'), list('JK')])
        assert_false(G.has_edge('I', 'J'))
        G.remove_nodes_from(set('ZEFHIMNO'))
        G.add_edge('J', 'K')

    def test_edges_nbunch(self):
        # Test G.edges(nbunch) with various forms of nbunch
        G = self.G()
        G.add_edges_from([('A', 'B'), ('A', 'C'), ('B', 'D'),
                          ('C', 'B'), ('C', 'D')])
        # node not in nbunch should be quietly ignored
        assert_raises(nx.NetworkXError, G.edges, 6)
        assert_equals(G.edges('Z'), [])  # iterable non-node
        # nbunch can be an empty list
        assert_equals(G.edges([]), [])
        if G.is_directed():
            elist = [('A', 'B'), ('A', 'C'), ('B', 'D')]
        else:
            elist = [('A', 'B'), ('A', 'C'), ('B', 'C'), ('B', 'D')]
        # nbunch can be a list
        assert_edges_equal(G.edges(['A', 'B']), elist)
        # nbunch can be a set
        assert_edges_equal(G.edges(set(['A', 'B'])), elist)
        # nbunch can be a graph
        G1 = self.G()
        G1.add_nodes_from('AB')
        assert_edges_equal(G.edges(G1), elist)
        # nbunch can be a dict with nodes as keys
        ndict = {'A': "thing1", 'B': "thing2"}
        assert_edges_equal(G.edges(ndict), elist)
        # nbunch can be a single node
        assert_edges_equal(G.edges('A'), [('A', 'B'), ('A', 'C')])
        assert_nodes_equal(sorted(G), ['A', 'B', 'C', 'D'])

    def test_edges_nbunch(self):
        G = self.G()
        G.add_edges_from([('A', 'B'), ('A', 'C'), ('B', 'D'),
                          ('C', 'B'), ('C', 'D')])
        # Test G.edges(nbunch) with various forms of nbunch
        # node not in nbunch should be quietly ignored
        assert_equals(list(G.edges('Z')), [])
        # nbunch can be an empty list
        assert_equals(sorted(G.edges([])), [])
        if G.is_directed():
            elist = [('A', 'B'), ('A', 'C'), ('B', 'D')]
        else:
            elist = [('A', 'B'), ('A', 'C'), ('B', 'C'), ('B', 'D')]
        # nbunch can be a list
        assert_edges_equal(G.edges(['A', 'B']), elist)
        # nbunch can be a set
        assert_edges_equal(G.edges(set(['A', 'B'])), elist)
        # nbunch can be a graph
        G1 = self.G()
        G1.add_nodes_from(['A', 'B'])
        assert_edges_equal(G.edges(G1), elist)
        # nbunch can be a dict with nodes as keys
        ndict = {'A': "thing1", 'B': "thing2"}
        assert_edges_equal(G.edges(ndict), elist)
        # nbunch can be a single node
        assert_edges_equal(G.edges('A'), [('A', 'B'), ('A', 'C')])

        # nbunch can be nothing (whole graph)
        assert_edges_equal(G.edges(), [('A', 'B'), ('A', 'C'), ('B', 'D'),
                                       ('C', 'B'), ('C', 'D')])

    def test_degree(self):
        G = self.G()
        G.add_edges_from([('A', 'B'), ('A', 'C'), ('B', 'D'),
                          ('C', 'B'), ('C', 'D')])
        assert_equal(G.degree('A'), 2)

        # degree of single node in iterable container must return dict
        assert_equal(list(G.degree(['A'])), [('A', 2)])
        assert_equal(sorted(d for n, d in G.degree(['A', 'B'])), [2, 3])
        assert_equal(sorted(d for n, d in G.degree()), [2, 2, 3, 3])

    def test_degree2(self):
        H = self.G()
        H.add_edges_from([(1, 24), (1, 2)])
        assert_equal(sorted(d for n, d in H.degree([1, 24])), [1, 2])

    def test_degree_graph(self):
        P3 = nx.path_graph(3)
        P5 = nx.path_graph(5)
        # silently ignore nodes not in P3
        assert_equal(dict(d for n, d in P3.degree(['A', 'B'])), {})
        # nbunch can be a graph
        assert_equal(sorted(d for n, d in P5.degree(P3)), [1, 2, 2])
        # nbunch can be a graph thats way to big
        assert_equal(sorted(d for n, d in P3.degree(P5)), [1, 1, 2])
        assert_equal(list(P5.degree([])), [])
        assert_equal(dict(P5.degree([])), {})

    def test_null(self):
        null = nx.null_graph()
        assert_equal(list(null.degree()), [])
        assert_equal(dict(null.degree()), {})

    def test_order_size(self):
        G = self.G()
        G.add_edges_from([('A', 'B'), ('A', 'C'), ('B', 'D'),
                          ('C', 'B'), ('C', 'D')])
        assert_equal(G.order(), 4)
        assert_equal(G.size(), 5)
        assert_equal(G.number_of_edges(), 5)
        assert_equal(G.number_of_edges('A', 'B'), 1)
        assert_equal(G.number_of_edges('A', 'D'), 0)

    def test_copy(self):
        G = self.G()
        H = G.copy()      # copy
        assert_equal(H.adj, G.adj)
        assert_equal(H.name, G.name)
        assert_not_equal(H, G)

    def test_subgraph(self):
        G = self.G()
        G.add_edges_from([('A', 'B'), ('A', 'C'), ('B', 'D'),
                          ('C', 'B'), ('C', 'D')])
        SG = G.subgraph(['A', 'B', 'D'])
        assert_nodes_equal(list(SG), ['A', 'B', 'D'])
        assert_edges_equal(list(SG.edges()), [('A', 'B'), ('B', 'D')])

    def test_to_directed(self):
        G = self.G()
        if not G.is_directed():
            G.add_edges_from([('A', 'B'), ('A', 'C'), ('B', 'D'),
                              ('C', 'B'), ('C', 'D')])

            DG = G.to_directed()
            assert_not_equal(DG, G)  # directed copy or copy

            assert_true(DG.is_directed())
            assert_equal(DG.name, G.name)
            assert_equal(DG.adj, G.adj)
            assert_equal(sorted(DG.out_edges(list('AB'))),
                         [('A', 'B'), ('A', 'C'), ('B', 'A'),
                          ('B', 'C'), ('B', 'D')])
            DG.remove_edge('A', 'B')
            assert_true(DG.has_edge('B', 'A'))  # this removes B-A but not  A-B
            assert_false(DG.has_edge('A', 'B'))

    def test_to_undirected(self):
        G = self.G()
        if G.is_directed():
            G.add_edges_from([('A', 'B'), ('A', 'C'), ('B', 'D'),
                              ('C', 'B'), ('C', 'D')])
            UG = G.to_undirected()       # to_undirected
            assert_not_equal(UG, G)
            assert_false(UG.is_directed())
            assert_true(G.is_directed())
            assert_equal(UG.name, G.name)
            assert_not_equal(UG.adj, G.adj)
            assert_equal(sorted(UG.edges(list('AB'))),
                         [('A', 'B'), ('A', 'C'), ('B', 'C'), ('B', 'D')])
            assert_equal(sorted(UG.edges(['A', 'B'])),
                         [('A', 'B'), ('A', 'C'), ('B', 'C'), ('B', 'D')])
            UG.remove_edge('A', 'B')
            assert_false(UG.has_edge('B', 'A'))
            assert_false(UG.has_edge('A', 'B'))

    def test_neighbors(self):
        G = self.G()
        G.add_edges_from([('A', 'B'), ('A', 'C'), ('B', 'D'),
                          ('C', 'B'), ('C', 'D')])
        G.add_nodes_from('GJK')
        assert_equal(sorted(G['A']), ['B', 'C'])
        assert_equal(sorted(G.neighbors('A')), ['B', 'C'])
        assert_equal(sorted(G.neighbors('A')), ['B', 'C'])
        assert_equal(sorted(G.neighbors('G')), [])
        assert_raises(nx.NetworkXError, G.neighbors, 'j')

    def test_iterators(self):
        G = self.G()
        G.add_edges_from([('A', 'B'), ('A', 'C'), ('B', 'D'),
                          ('C', 'B'), ('C', 'D')])
        G.add_nodes_from('GJK')
        assert_equal(sorted(G.nodes()),
                     ['A', 'B', 'C', 'D', 'G', 'J', 'K'])
        assert_edges_equal(G.edges(),
                           [('A', 'B'), ('A', 'C'), ('B', 'D'), ('C', 'B'), ('C', 'D')])

        assert_equal(sorted([v for k, v in G.degree()]),
                     [0, 0, 0, 2, 2, 3, 3])
        assert_equal(sorted(G.degree(), key=str),
                     [('A', 2), ('B', 3), ('C', 3), ('D', 2),
                      ('G', 0), ('J', 0), ('K', 0)])
        assert_equal(sorted(G.neighbors('A')), ['B', 'C'])
        assert_raises(nx.NetworkXError, G.neighbors, 'X')
        G.clear()
        assert_equal(nx.number_of_nodes(G), 0)
        assert_equal(nx.number_of_edges(G), 0)

    def test_null_subgraph(self):
        # Subgraph of a null graph is a null graph
        nullgraph = nx.null_graph()
        G = nx.null_graph()
        H = G.subgraph([])
        assert_true(nx.is_isomorphic(H, nullgraph))

    def test_empty_subgraph(self):
        # Subgraph of an empty graph is an empty graph. test 1
        nullgraph = nx.null_graph()
        E5 = nx.empty_graph(5)
        E10 = nx.empty_graph(10)
        H = E10.subgraph([])
        assert_true(nx.is_isomorphic(H, nullgraph))
        H = E10.subgraph([1, 2, 3, 4, 5])
        assert_true(nx.is_isomorphic(H, E5))

    def test_complete_subgraph(self):
        # Subgraph of a complete graph is a complete graph
        K1 = nx.complete_graph(1)
        K3 = nx.complete_graph(3)
        K5 = nx.complete_graph(5)
        H = K5.subgraph([1, 2, 3])
        assert_true(nx.is_isomorphic(H, K3))

    def test_subgraph_nbunch(self):
        nullgraph = nx.null_graph()
        K1 = nx.complete_graph(1)
        K3 = nx.complete_graph(3)
        K5 = nx.complete_graph(5)
        # Test G.subgraph(nbunch), where nbunch is a single node
        H = K5.subgraph(1)
        assert_true(nx.is_isomorphic(H, K1))
        # Test G.subgraph(nbunch), where nbunch is a set
        H = K5.subgraph(set([1]))
        assert_true(nx.is_isomorphic(H, K1))
        # Test G.subgraph(nbunch), where nbunch is an iterator
        H = K5.subgraph(iter(K3))
        assert_true(nx.is_isomorphic(H, K3))
        # Test G.subgraph(nbunch), where nbunch is another graph
        H = K5.subgraph(K3)
        assert_true(nx.is_isomorphic(H, K3))
        H = K5.subgraph([9])
        assert_true(nx.is_isomorphic(H, nullgraph))

    def test_node_tuple_issue(self):
        H = self.G()
        # Test error handling of tuple as a node
        assert_raises(nx.NetworkXError, H.remove_node, (1, 2))
        H.remove_nodes_from([(1, 2)])  # no error
        assert_raises(nx.NetworkXError, H.neighbors, (1, 2))
