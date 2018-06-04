#!/usr/bin/env python
import random
from nose.tools import *
import networkx as nx
from networkx.testing.utils import *


class TestFunction(object):
    def setUp(self):
        self.G = nx.Graph({0: [1, 2, 3], 1: [1, 2, 0], 4: []}, name='Test')
        self.Gdegree = {0: 3, 1: 2, 2: 2, 3: 1, 4: 0}
        self.Gnodes = list(range(5))
        self.Gedges = [(0, 1), (0, 2), (0, 3), (1, 0), (1, 1), (1, 2)]
        self.DG = nx.DiGraph({0: [1, 2, 3], 1: [1, 2, 0], 4: []})
        self.DGin_degree = {0: 1, 1: 2, 2: 2, 3: 1, 4: 0}
        self.DGout_degree = {0: 3, 1: 3, 2: 0, 3: 0, 4: 0}
        self.DGnodes = list(range(5))
        self.DGedges = [(0, 1), (0, 2), (0, 3), (1, 0), (1, 1), (1, 2)]

    def test_nodes(self):
        assert_nodes_equal(self.G.nodes(), list(nx.nodes(self.G)))
        assert_nodes_equal(self.DG.nodes(), list(nx.nodes(self.DG)))

    def test_edges(self):
        assert_edges_equal(self.G.edges(), list(nx.edges(self.G)))
        assert_equal(sorted(self.DG.edges()), sorted(nx.edges(self.DG)))
        assert_edges_equal(self.G.edges(nbunch=[0, 1, 3]),
                           list(nx.edges(self.G, nbunch=[0, 1, 3])))
        assert_equal(sorted(self.DG.edges(nbunch=[0, 1, 3])),
                     sorted(nx.edges(self.DG, nbunch=[0, 1, 3])))

    def test_degree(self):
        assert_edges_equal(self.G.degree(), list(nx.degree(self.G)))
        assert_equal(sorted(self.DG.degree()), sorted(nx.degree(self.DG)))
        assert_edges_equal(self.G.degree(nbunch=[0, 1]),
                           list(nx.degree(self.G, nbunch=[0, 1])))
        assert_equal(sorted(self.DG.degree(nbunch=[0, 1])),
                     sorted(nx.degree(self.DG, nbunch=[0, 1])))
        assert_edges_equal(self.G.degree(weight='weight'),
                           list(nx.degree(self.G, weight='weight')))
        assert_equal(sorted(self.DG.degree(weight='weight')),
                     sorted(nx.degree(self.DG, weight='weight')))

    def test_neighbors(self):
        assert_equal(self.G.neighbors(1), nx.neighbors(self.G, 1))
        assert_equal(self.DG.neighbors(1), nx.neighbors(self.DG, 1))

    def test_number_of_nodes(self):
        assert_equal(self.G.number_of_nodes(), nx.number_of_nodes(self.G))
        assert_equal(self.DG.number_of_nodes(), nx.number_of_nodes(self.DG))

    def test_number_of_edges(self):
        assert_equal(self.G.number_of_edges(), nx.number_of_edges(self.G))
        assert_equal(self.DG.number_of_edges(), nx.number_of_edges(self.DG))

    def test_is_directed(self):
        assert_equal(self.G.is_directed(), nx.is_directed(self.G))
        assert_equal(self.DG.is_directed(), nx.is_directed(self.DG))

    def test_add_star(self):
        G = self.G.copy()
        nlist = [12, 13, 14, 15]
        nx.add_star(G, nlist)
        assert_edges_equal(G.edges(nlist), [(12, 13), (12, 14), (12, 15)])
        G = self.G.copy()
        nx.add_star(G, nlist, weight=2.0)
        assert_edges_equal(G.edges(nlist, data=True),
                           [(12, 13, {'weight': 2.}),
                            (12, 14, {'weight': 2.}),
                            (12, 15, {'weight': 2.})])

    def test_add_path(self):
        G = self.G.copy()
        nlist = [12, 13, 14, 15]
        nx.add_path(G, nlist)
        assert_edges_equal(G.edges(nlist), [(12, 13), (13, 14), (14, 15)])
        G = self.G.copy()
        nx.add_path(G, nlist, weight=2.0)
        assert_edges_equal(G.edges(nlist, data=True),
                           [(12, 13, {'weight': 2.}),
                            (13, 14, {'weight': 2.}),
                            (14, 15, {'weight': 2.})])

        G = self.G.copy()
        nlist = [None]
        nx.add_path(G, nlist)
        assert_edges_equal(G.edges(nlist), [])
        assert_nodes_equal(G, list(self.G) + [None])

        G = self.G.copy()
        nlist = iter([None])
        nx.add_path(G, nlist)
        assert_edges_equal(G.edges([None]), [])
        assert_nodes_equal(G, list(self.G) + [None])

        G = self.G.copy()
        nlist = [12]
        nx.add_path(G, nlist)
        assert_edges_equal(G.edges(nlist), [])
        assert_nodes_equal(G, list(self.G) + [12])

        G = self.G.copy()
        nlist = iter([12])
        nx.add_path(G, nlist)
        assert_edges_equal(G.edges([12]), [])
        assert_nodes_equal(G, list(self.G) + [12])

        G = self.G.copy()
        nlist = []
        nx.add_path(G, nlist)
        assert_edges_equal(G.edges, self.G.edges)
        assert_nodes_equal(G, list(self.G))

        G = self.G.copy()
        nlist = iter([])
        nx.add_path(G, nlist)
        assert_edges_equal(G.edges, self.G.edges)
        assert_nodes_equal(G, list(self.G))

    def test_add_cycle(self):
        G = self.G.copy()
        nlist = [12, 13, 14, 15]
        oklists = [[(12, 13), (12, 15), (13, 14), (14, 15)],
                   [(12, 13), (13, 14), (14, 15), (15, 12)]]
        nx.add_cycle(G, nlist)
        assert_true(sorted(G.edges(nlist)) in oklists)
        G = self.G.copy()
        oklists = [[(12, 13, {'weight': 1.}),
                    (12, 15, {'weight': 1.}),
                    (13, 14, {'weight': 1.}),
                    (14, 15, {'weight': 1.})],
                   [(12, 13, {'weight': 1.}),
                    (13, 14, {'weight': 1.}),
                    (14, 15, {'weight': 1.}),
                    (15, 12, {'weight': 1.})]]
        nx.add_cycle(G, nlist, weight=1.0)
        assert_true(sorted(G.edges(nlist, data=True)) in oklists)

    def test_subgraph(self):
        assert_equal(self.G.subgraph([0, 1, 2, 4]).adj,
                     nx.subgraph(self.G, [0, 1, 2, 4]).adj)
        assert_equal(self.DG.subgraph([0, 1, 2, 4]).adj,
                     nx.subgraph(self.DG, [0, 1, 2, 4]).adj)
        assert_equal(self.G.subgraph([0, 1, 2, 4]).adj,
                     nx.induced_subgraph(self.G, [0, 1, 2, 4]).adj)
        assert_equal(self.DG.subgraph([0, 1, 2, 4]).adj,
                     nx.induced_subgraph(self.DG, [0, 1, 2, 4]).adj)
        # subgraph-subgraph chain is allowed in function interface
        H = nx.induced_subgraph(self.G.subgraph([0, 1, 2, 4]), [0, 1, 4])
        assert_is_not(H._graph, self.G)
        assert_equal(H.adj, self.G.subgraph([0, 1, 4]).adj)

    def test_edge_subgraph(self):
        assert_equal(self.G.edge_subgraph([(1, 2), (0, 3)]).adj,
                     nx.edge_subgraph(self.G, [(1, 2), (0, 3)]).adj)
        assert_equal(self.DG.edge_subgraph([(1, 2), (0, 3)]).adj,
                     nx.edge_subgraph(self.DG, [(1, 2), (0, 3)]).adj)

    def test_restricted_view(self):
        H = nx.restricted_view(self.G, [0, 2, 5], [(1, 2), (3, 4)])
        assert_equal(set(H.nodes), {1, 3, 4})
        assert_equal(set(H.edges), {(1, 1)})

    def test_create_empty_copy(self):
        G = nx.create_empty_copy(self.G, with_data=False)
        assert_nodes_equal(G, list(self.G))
        assert_equal(G.graph, {})
        assert_equal(G._node, {}.fromkeys(self.G.nodes(), {}))
        assert_equal(G._adj, {}.fromkeys(self.G.nodes(), {}))
        G = nx.create_empty_copy(self.G)
        assert_nodes_equal(G, list(self.G))
        assert_equal(G.graph, self.G.graph)
        assert_equal(G._node, self.G._node)
        assert_equal(G._adj, {}.fromkeys(self.G.nodes(), {}))

    def test_degree_histogram(self):
        assert_equal(nx.degree_histogram(self.G), [1, 1, 1, 1, 1])

    def test_density(self):
        assert_equal(nx.density(self.G), 0.5)
        assert_equal(nx.density(self.DG), 0.3)
        G = nx.Graph()
        G.add_node(1)
        assert_equal(nx.density(G), 0.0)

    def test_density_selfloop(self):
        G = nx.Graph()
        G.add_edge(1, 1)
        assert_equal(nx.density(G), 0.0)
        G.add_edge(1, 2)
        assert_equal(nx.density(G), 2.0)

    def test_freeze(self):
        G = nx.freeze(self.G)
        assert_equal(G.frozen, True)
        assert_raises(nx.NetworkXError, G.add_node, 1)
        assert_raises(nx.NetworkXError, G.add_nodes_from, [1])
        assert_raises(nx.NetworkXError, G.remove_node, 1)
        assert_raises(nx.NetworkXError, G.remove_nodes_from, [1])
        assert_raises(nx.NetworkXError, G.add_edge, 1, 2)
        assert_raises(nx.NetworkXError, G.add_edges_from, [(1, 2)])
        assert_raises(nx.NetworkXError, G.remove_edge, 1, 2)
        assert_raises(nx.NetworkXError, G.remove_edges_from, [(1, 2)])
        assert_raises(nx.NetworkXError, G.clear)

    def test_is_frozen(self):
        assert_equal(nx.is_frozen(self.G), False)
        G = nx.freeze(self.G)
        assert_equal(G.frozen, nx.is_frozen(self.G))
        assert_equal(G.frozen, True)

    def test_info(self):
        G = nx.path_graph(5)
        G.name = "path_graph(5)"
        info = nx.info(G)
        expected_graph_info = '\n'.join(['Name: path_graph(5)',
                                         'Type: Graph',
                                         'Number of nodes: 5',
                                         'Number of edges: 4',
                                         'Average degree:   1.6000'])
        assert_equal(info, expected_graph_info)

        info = nx.info(G, n=1)
        expected_node_info = '\n'.join(
            ['Node 1 has the following properties:',
             'Degree: 2',
             'Neighbors: 0 2'])
        assert_equal(info, expected_node_info)

    def test_info_digraph(self):
        G = nx.DiGraph(name='path_graph(5)')
        nx.add_path(G, [0, 1, 2, 3, 4])
        info = nx.info(G)
        expected_graph_info = '\n'.join(['Name: path_graph(5)',
                                         'Type: DiGraph',
                                         'Number of nodes: 5',
                                         'Number of edges: 4',
                                         'Average in degree:   0.8000',
                                         'Average out degree:   0.8000'])
        assert_equal(info, expected_graph_info)

        info = nx.info(G, n=1)
        expected_node_info = '\n'.join(
            ['Node 1 has the following properties:',
             'Degree: 2',
             'Neighbors: 2'])
        assert_equal(info, expected_node_info)

        assert_raises(nx.NetworkXError, nx.info, G, n=-1)

    def test_neighbors(self):
        graph = nx.complete_graph(100)
        pop = random.sample(list(graph), 1)
        nbors = list(nx.neighbors(graph, pop[0]))
        # should be all the other vertices in the graph
        assert_equal(len(nbors), len(graph) - 1)

        graph = nx.path_graph(100)
        node = random.sample(list(graph), 1)[0]
        nbors = list(nx.neighbors(graph, node))
        # should be all the other vertices in the graph
        if node != 0 and node != 99:
            assert_equal(len(nbors), 2)
        else:
            assert_equal(len(nbors), 1)

        # create a star graph with 99 outer nodes
        graph = nx.star_graph(99)
        nbors = list(nx.neighbors(graph, 0))
        assert_equal(len(nbors), 99)

    def test_non_neighbors(self):
        graph = nx.complete_graph(100)
        pop = random.sample(list(graph), 1)
        nbors = list(nx.non_neighbors(graph, pop[0]))
        # should be all the other vertices in the graph
        assert_equal(len(nbors), 0)

        graph = nx.path_graph(100)
        node = random.sample(list(graph), 1)[0]
        nbors = list(nx.non_neighbors(graph, node))
        # should be all the other vertices in the graph
        if node != 0 and node != 99:
            assert_equal(len(nbors), 97)
        else:
            assert_equal(len(nbors), 98)

        # create a star graph with 99 outer nodes
        graph = nx.star_graph(99)
        nbors = list(nx.non_neighbors(graph, 0))
        assert_equal(len(nbors), 0)

        # disconnected graph
        graph = nx.Graph()
        graph.add_nodes_from(range(10))
        nbors = list(nx.non_neighbors(graph, 0))
        assert_equal(len(nbors), 9)

    def test_non_edges(self):
        # All possible edges exist
        graph = nx.complete_graph(5)
        nedges = list(nx.non_edges(graph))
        assert_equal(len(nedges), 0)

        graph = nx.path_graph(4)
        expected = [(0, 2), (0, 3), (1, 3)]
        nedges = list(nx.non_edges(graph))
        for (u, v) in expected:
            assert_true((u, v) in nedges or (v, u) in nedges)

        graph = nx.star_graph(4)
        expected = [(1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)]
        nedges = list(nx.non_edges(graph))
        for (u, v) in expected:
            assert_true((u, v) in nedges or (v, u) in nedges)

        # Directed graphs
        graph = nx.DiGraph()
        graph.add_edges_from([(0, 2), (2, 0), (2, 1)])
        expected = [(0, 1), (1, 0), (1, 2)]
        nedges = list(nx.non_edges(graph))
        for e in expected:
            assert_true(e in nedges)

    def test_is_weighted(self):
        G = nx.Graph()
        assert_false(nx.is_weighted(G))

        G = nx.path_graph(4)
        assert_false(nx.is_weighted(G))
        assert_false(nx.is_weighted(G, (2, 3)))

        G.add_node(4)
        G.add_edge(3, 4, weight=4)
        assert_false(nx.is_weighted(G))
        assert_true(nx.is_weighted(G, (3, 4)))

        G = nx.DiGraph()
        G.add_weighted_edges_from([('0', '3', 3), ('0', '1', -5),
                                   ('1', '0', -5), ('0', '2', 2),
                                   ('1', '2', 4), ('2', '3', 1)])
        assert_true(nx.is_weighted(G))
        assert_true(nx.is_weighted(G, ('1', '0')))

        G = G.to_undirected()
        assert_true(nx.is_weighted(G))
        assert_true(nx.is_weighted(G, ('1', '0')))

        assert_raises(nx.NetworkXError, nx.is_weighted, G, (1, 2))

    def test_is_negatively_weighted(self):
        G = nx.Graph()
        assert_false(nx.is_negatively_weighted(G))

        G.add_node(1)
        G.add_nodes_from([2, 3, 4, 5])
        assert_false(nx.is_negatively_weighted(G))

        G.add_edge(1, 2, weight=4)
        assert_false(nx.is_negatively_weighted(G, (1, 2)))

        G.add_edges_from([(1, 3), (2, 4), (2, 6)])
        G[1][3]['color'] = 'blue'
        assert_false(nx.is_negatively_weighted(G))
        assert_false(nx.is_negatively_weighted(G, (1, 3)))

        G[2][4]['weight'] = -2
        assert_true(nx.is_negatively_weighted(G, (2, 4)))
        assert_true(nx.is_negatively_weighted(G))

        G = nx.DiGraph()
        G.add_weighted_edges_from([('0', '3', 3), ('0', '1', -5),
                                   ('1', '0', -2), ('0', '2', 2),
                                   ('1', '2', -3), ('2', '3', 1)])
        assert_true(nx.is_negatively_weighted(G))
        assert_false(nx.is_negatively_weighted(G, ('0', '3')))
        assert_true(nx.is_negatively_weighted(G, ('1', '0')))

        assert_raises(nx.NetworkXError, nx.is_negatively_weighted, G, (1, 4))


class TestCommonNeighbors():
    def setUp(self):
        self.func = nx.common_neighbors

        def test_func(G, u, v, expected):
            result = sorted(self.func(G, u, v))
            assert_equal(result, expected)
        self.test = test_func

    def test_K5(self):
        G = nx.complete_graph(5)
        self.test(G, 0, 1, [2, 3, 4])

    def test_P3(self):
        G = nx.path_graph(3)
        self.test(G, 0, 2, [1])

    def test_S4(self):
        G = nx.star_graph(4)
        self.test(G, 1, 2, [0])

    @raises(nx.NetworkXNotImplemented)
    def test_digraph(self):
        G = nx.DiGraph()
        G.add_edges_from([(0, 1), (1, 2)])
        self.func(G, 0, 2)

    def test_nonexistent_nodes(self):
        G = nx.complete_graph(5)
        assert_raises(nx.NetworkXError, nx.common_neighbors, G, 5, 4)
        assert_raises(nx.NetworkXError, nx.common_neighbors, G, 4, 5)
        assert_raises(nx.NetworkXError, nx.common_neighbors, G, 5, 6)

    def test_custom1(self):
        """Case of no common neighbors."""
        G = nx.Graph()
        G.add_nodes_from([0, 1])
        self.test(G, 0, 1, [])

    def test_custom2(self):
        """Case of equal nodes."""
        G = nx.complete_graph(4)
        self.test(G, 0, 0, [1, 2, 3])


def test_set_node_attributes():
    graphs = [nx.Graph(), nx.DiGraph(), nx.MultiGraph(), nx.MultiDiGraph()]
    for G in graphs:
        # Test single value
        G = nx.path_graph(3, create_using=G)
        vals = 100
        attr = 'hello'
        nx.set_node_attributes(G, vals, attr)
        assert_equal(G.nodes[0][attr], vals)
        assert_equal(G.nodes[1][attr], vals)
        assert_equal(G.nodes[2][attr], vals)

        # Test dictionary
        G = nx.path_graph(3, create_using=G)
        vals = dict(zip(sorted(G.nodes()), range(len(G))))
        attr = 'hi'
        nx.set_node_attributes(G, vals, attr)
        assert_equal(G.nodes[0][attr], 0)
        assert_equal(G.nodes[1][attr], 1)
        assert_equal(G.nodes[2][attr], 2)

        # Test dictionary of dictionaries
        G = nx.path_graph(3, create_using=G)
        d = {'hi': 0, 'hello': 200}
        vals = dict.fromkeys(G.nodes(), d)
        vals.pop(0)
        nx.set_node_attributes(G, vals)
        assert_equal(G.nodes[0], {})
        assert_equal(G.nodes[1]["hi"], 0)
        assert_equal(G.nodes[2]["hello"], 200)


def test_set_edge_attributes():
    graphs = [nx.Graph(), nx.DiGraph()]
    for G in graphs:
        # Test single value
        G = nx.path_graph(3, create_using=G)
        attr = 'hello'
        vals = 3
        nx.set_edge_attributes(G, vals, attr)
        assert_equal(G[0][1][attr], vals)
        assert_equal(G[1][2][attr], vals)

        # Test multiple values
        G = nx.path_graph(3, create_using=G)
        attr = 'hi'
        edges = [(0, 1), (1, 2)]
        vals = dict(zip(edges, range(len(edges))))
        nx.set_edge_attributes(G, vals, attr)
        assert_equal(G[0][1][attr], 0)
        assert_equal(G[1][2][attr], 1)

        # Test dictionary of dictionaries
        G = nx.path_graph(3, create_using=G)
        d = {'hi': 0, 'hello': 200}
        edges = [(0, 1)]
        vals = dict.fromkeys(edges, d)
        nx.set_edge_attributes(G, vals)
        assert_equal(G[0][1]['hi'], 0)
        assert_equal(G[0][1]['hello'], 200)
        assert_equal(G[1][2], {})


def test_set_edge_attributes_multi():
    graphs = [nx.MultiGraph(), nx.MultiDiGraph()]
    for G in graphs:
        # Test single value
        G = nx.path_graph(3, create_using=G)
        attr = 'hello'
        vals = 3
        nx.set_edge_attributes(G, vals, attr)
        assert_equal(G[0][1][0][attr], vals)
        assert_equal(G[1][2][0][attr], vals)

        # Test multiple values
        G = nx.path_graph(3, create_using=G)
        attr = 'hi'
        edges = [(0, 1, 0), (1, 2, 0)]
        vals = dict(zip(edges, range(len(edges))))
        nx.set_edge_attributes(G, vals, attr)
        assert_equal(G[0][1][0][attr], 0)
        assert_equal(G[1][2][0][attr], 1)

        # Test dictionary of dictionaries
        G = nx.path_graph(3, create_using=G)
        d = {'hi': 0, 'hello': 200}
        edges = [(0, 1, 0)]
        vals = dict.fromkeys(edges, d)
        nx.set_edge_attributes(G, vals)
        assert_equal(G[0][1][0]['hi'], 0)
        assert_equal(G[0][1][0]['hello'], 200)
        assert_equal(G[1][2][0], {})


def test_get_node_attributes():
    graphs = [nx.Graph(), nx.DiGraph(), nx.MultiGraph(), nx.MultiDiGraph()]
    for G in graphs:
        G = nx.path_graph(3, create_using=G)
        attr = 'hello'
        vals = 100
        nx.set_node_attributes(G, vals, attr)
        attrs = nx.get_node_attributes(G, attr)
        assert_equal(attrs[0], vals)
        assert_equal(attrs[1], vals)
        assert_equal(attrs[2], vals)


def test_get_edge_attributes():
    graphs = [nx.Graph(), nx.DiGraph(), nx.MultiGraph(), nx.MultiDiGraph()]
    for G in graphs:
        G = nx.path_graph(3, create_using=G)
        attr = 'hello'
        vals = 100
        nx.set_edge_attributes(G, vals, attr)
        attrs = nx.get_edge_attributes(G, attr)

        assert_equal(len(attrs), 2)
        if G.is_multigraph():
            keys = [(0, 1, 0), (1, 2, 0)]
            for u, v, k in keys:
                try:
                    assert_equal(attrs[(u, v, k)], 100)
                except KeyError:
                    assert_equal(attrs[(v, u, k)], 100)
        else:
            keys = [(0, 1), (1, 2)]
            for u, v in keys:
                try:
                    assert_equal(attrs[(u, v)], 100)
                except KeyError:
                    assert_equal(attrs[(v, u)], 100)


def test_is_empty():
    graphs = [nx.Graph(), nx.DiGraph(), nx.MultiGraph(), nx.MultiDiGraph()]
    for G in graphs:
        assert_true(nx.is_empty(G))
        G.add_nodes_from(range(5))
        assert_true(nx.is_empty(G))
        G.add_edges_from([(1, 2), (3, 4)])
        assert_false(nx.is_empty(G))


def test_selfloops():
    graphs = [nx.Graph(), nx.DiGraph(), nx.MultiGraph(), nx.MultiDiGraph()]
    for graph in graphs:
        G = nx.complete_graph(3, create_using=graph)
        G.add_edge(0, 0)
        assert_nodes_equal(nx.nodes_with_selfloops(G), [0])
        assert_edges_equal(nx.selfloop_edges(G), [(0, 0)])
        assert_edges_equal(nx.selfloop_edges(G, data=True), [(0, 0, {})])
        assert_equal(nx.number_of_selfloops(G), 1)
        # test selfloop attr
        G.add_edge(1, 1, weight=2)
        assert_edges_equal(nx.selfloop_edges(G, data=True),
                           [(0, 0, {}), (1, 1, {'weight': 2})])
        assert_edges_equal(nx.selfloop_edges(G, data='weight'),
                           [(0, 0, None), (1, 1, 2)])
