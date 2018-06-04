#!/usr/bin/env python
from nose.tools import assert_equal
import networkx as nx
from networkx.algorithms import bipartite
from networkx.testing import assert_edges_equal, assert_nodes_equal


class TestBipartiteProject:

    def test_path_projected_graph(self):
        G = nx.path_graph(4)
        P = bipartite.projected_graph(G, [1, 3])
        assert_nodes_equal(list(P), [1, 3])
        assert_edges_equal(list(P.edges()), [(1, 3)])
        P = bipartite.projected_graph(G, [0, 2])
        assert_nodes_equal(list(P), [0, 2])
        assert_edges_equal(list(P.edges()), [(0, 2)])

    def test_path_projected_properties_graph(self):
        G = nx.path_graph(4)
        G.add_node(1, name='one')
        G.add_node(2, name='two')
        P = bipartite.projected_graph(G, [1, 3])
        assert_nodes_equal(list(P), [1, 3])
        assert_edges_equal(list(P.edges()), [(1, 3)])
        assert_equal(P.nodes[1]['name'], G.nodes[1]['name'])
        P = bipartite.projected_graph(G, [0, 2])
        assert_nodes_equal(list(P), [0, 2])
        assert_edges_equal(list(P.edges()), [(0, 2)])
        assert_equal(P.nodes[2]['name'], G.nodes[2]['name'])

    def test_path_collaboration_projected_graph(self):
        G = nx.path_graph(4)
        P = bipartite.collaboration_weighted_projected_graph(G, [1, 3])
        assert_nodes_equal(list(P), [1, 3])
        assert_edges_equal(list(P.edges()), [(1, 3)])
        P[1][3]['weight'] = 1
        P = bipartite.collaboration_weighted_projected_graph(G, [0, 2])
        assert_nodes_equal(list(P), [0, 2])
        assert_edges_equal(list(P.edges()), [(0, 2)])
        P[0][2]['weight'] = 1

    def test_directed_path_collaboration_projected_graph(self):
        G = nx.DiGraph()
        nx.add_path(G, range(4))
        P = bipartite.collaboration_weighted_projected_graph(G, [1, 3])
        assert_nodes_equal(list(P), [1, 3])
        assert_edges_equal(list(P.edges()), [(1, 3)])
        P[1][3]['weight'] = 1
        P = bipartite.collaboration_weighted_projected_graph(G, [0, 2])
        assert_nodes_equal(list(P), [0, 2])
        assert_edges_equal(list(P.edges()), [(0, 2)])
        P[0][2]['weight'] = 1

    def test_path_weighted_projected_graph(self):
        G = nx.path_graph(4)
        P = bipartite.weighted_projected_graph(G, [1, 3])
        assert_nodes_equal(list(P), [1, 3])
        assert_edges_equal(list(P.edges()), [(1, 3)])
        P[1][3]['weight'] = 1
        P = bipartite.weighted_projected_graph(G, [0, 2])
        assert_nodes_equal(list(P), [0, 2])
        assert_edges_equal(list(P.edges()), [(0, 2)])
        P[0][2]['weight'] = 1

    def test_path_weighted_projected_directed_graph(self):
        G = nx.DiGraph()
        nx.add_path(G, range(4))
        P = bipartite.weighted_projected_graph(G, [1, 3])
        assert_nodes_equal(list(P), [1, 3])
        assert_edges_equal(list(P.edges()), [(1, 3)])
        P[1][3]['weight'] = 1
        P = bipartite.weighted_projected_graph(G, [0, 2])
        assert_nodes_equal(list(P), [0, 2])
        assert_edges_equal(list(P.edges()), [(0, 2)])
        P[0][2]['weight'] = 1

    def test_star_projected_graph(self):
        G = nx.star_graph(3)
        P = bipartite.projected_graph(G, [1, 2, 3])
        assert_nodes_equal(list(P), [1, 2, 3])
        assert_edges_equal(list(P.edges()), [(1, 2), (1, 3), (2, 3)])
        P = bipartite.weighted_projected_graph(G, [1, 2, 3])
        assert_nodes_equal(list(P), [1, 2, 3])
        assert_edges_equal(list(P.edges()), [(1, 2), (1, 3), (2, 3)])

        P = bipartite.projected_graph(G, [0])
        assert_nodes_equal(list(P), [0])
        assert_edges_equal(list(P.edges()), [])

    def test_project_multigraph(self):
        G = nx.Graph()
        G.add_edge('a', 1)
        G.add_edge('b', 1)
        G.add_edge('a', 2)
        G.add_edge('b', 2)
        P = bipartite.projected_graph(G, 'ab')
        assert_edges_equal(list(P.edges()), [('a', 'b')])
        P = bipartite.weighted_projected_graph(G, 'ab')
        assert_edges_equal(list(P.edges()), [('a', 'b')])
        P = bipartite.projected_graph(G, 'ab', multigraph=True)
        assert_edges_equal(list(P.edges()), [('a', 'b'), ('a', 'b')])

    def test_project_collaboration(self):
        G = nx.Graph()
        G.add_edge('a', 1)
        G.add_edge('b', 1)
        G.add_edge('b', 2)
        G.add_edge('c', 2)
        G.add_edge('c', 3)
        G.add_edge('c', 4)
        G.add_edge('b', 4)
        P = bipartite.collaboration_weighted_projected_graph(G, 'abc')
        assert_equal(P['a']['b']['weight'], 1)
        assert_equal(P['b']['c']['weight'], 2)

    def test_directed_projection(self):
        G = nx.DiGraph()
        G.add_edge('A', 1)
        G.add_edge(1, 'B')
        G.add_edge('A', 2)
        G.add_edge('B', 2)
        P = bipartite.projected_graph(G, 'AB')
        assert_edges_equal(list(P.edges()), [('A', 'B')])
        P = bipartite.weighted_projected_graph(G, 'AB')
        assert_edges_equal(list(P.edges()), [('A', 'B')])
        assert_equal(P['A']['B']['weight'], 1)

        P = bipartite.projected_graph(G, 'AB', multigraph=True)
        assert_edges_equal(list(P.edges()), [('A', 'B')])

        G = nx.DiGraph()
        G.add_edge('A', 1)
        G.add_edge(1, 'B')
        G.add_edge('A', 2)
        G.add_edge(2, 'B')
        P = bipartite.projected_graph(G, 'AB')
        assert_edges_equal(list(P.edges()), [('A', 'B')])
        P = bipartite.weighted_projected_graph(G, 'AB')
        assert_edges_equal(list(P.edges()), [('A', 'B')])
        assert_equal(P['A']['B']['weight'], 2)

        P = bipartite.projected_graph(G, 'AB', multigraph=True)
        assert_edges_equal(list(P.edges()), [('A', 'B'), ('A', 'B')])


class TestBipartiteWeightedProjection:

    def setUp(self):
        # Tore Opsahl's example
        # http://toreopsahl.com/2009/05/01/projecting-two-mode-networks-onto-weighted-one-mode-networks/
        self.G = nx.Graph()
        self.G.add_edge('A', 1)
        self.G.add_edge('A', 2)
        self.G.add_edge('B', 1)
        self.G.add_edge('B', 2)
        self.G.add_edge('B', 3)
        self.G.add_edge('B', 4)
        self.G.add_edge('B', 5)
        self.G.add_edge('C', 1)
        self.G.add_edge('D', 3)
        self.G.add_edge('E', 4)
        self.G.add_edge('E', 5)
        self.G.add_edge('E', 6)
        self.G.add_edge('F', 6)
        # Graph based on figure 6 from Newman (2001)
        self.N = nx.Graph()
        self.N.add_edge('A', 1)
        self.N.add_edge('A', 2)
        self.N.add_edge('A', 3)
        self.N.add_edge('B', 1)
        self.N.add_edge('B', 2)
        self.N.add_edge('B', 3)
        self.N.add_edge('C', 1)
        self.N.add_edge('D', 1)
        self.N.add_edge('E', 3)

    def test_project_weighted_shared(self):
        edges = [('A', 'B', 2),
                 ('A', 'C', 1),
                 ('B', 'C', 1),
                 ('B', 'D', 1),
                 ('B', 'E', 2),
                 ('E', 'F', 1)]
        Panswer = nx.Graph()
        Panswer.add_weighted_edges_from(edges)
        P = bipartite.weighted_projected_graph(self.G, 'ABCDEF')
        assert_edges_equal(list(P.edges()), Panswer.edges())
        for u, v in list(P.edges()):
            assert_equal(P[u][v]['weight'], Panswer[u][v]['weight'])

        edges = [('A', 'B', 3),
                 ('A', 'E', 1),
                 ('A', 'C', 1),
                 ('A', 'D', 1),
                 ('B', 'E', 1),
                 ('B', 'C', 1),
                 ('B', 'D', 1),
                 ('C', 'D', 1)]
        Panswer = nx.Graph()
        Panswer.add_weighted_edges_from(edges)
        P = bipartite.weighted_projected_graph(self.N, 'ABCDE')
        assert_edges_equal(list(P.edges()), Panswer.edges())
        for u, v in list(P.edges()):
            assert_equal(P[u][v]['weight'], Panswer[u][v]['weight'])

    def test_project_weighted_newman(self):
        edges = [('A', 'B', 1.5),
                 ('A', 'C', 0.5),
                 ('B', 'C', 0.5),
                 ('B', 'D', 1),
                 ('B', 'E', 2),
                 ('E', 'F', 1)]
        Panswer = nx.Graph()
        Panswer.add_weighted_edges_from(edges)
        P = bipartite.collaboration_weighted_projected_graph(self.G, 'ABCDEF')
        assert_edges_equal(list(P.edges()), Panswer.edges())
        for u, v in list(P.edges()):
            assert_equal(P[u][v]['weight'], Panswer[u][v]['weight'])

        edges = [('A', 'B', 11 / 6.0),
                 ('A', 'E', 1 / 2.0),
                 ('A', 'C', 1 / 3.0),
                 ('A', 'D', 1 / 3.0),
                 ('B', 'E', 1 / 2.0),
                 ('B', 'C', 1 / 3.0),
                 ('B', 'D', 1 / 3.0),
                 ('C', 'D', 1 / 3.0)]
        Panswer = nx.Graph()
        Panswer.add_weighted_edges_from(edges)
        P = bipartite.collaboration_weighted_projected_graph(self.N, 'ABCDE')
        assert_edges_equal(list(P.edges()), Panswer.edges())
        for u, v in list(P.edges()):
            assert_equal(P[u][v]['weight'], Panswer[u][v]['weight'])

    def test_project_weighted_ratio(self):
        edges = [('A', 'B', 2 / 6.0),
                 ('A', 'C', 1 / 6.0),
                 ('B', 'C', 1 / 6.0),
                 ('B', 'D', 1 / 6.0),
                 ('B', 'E', 2 / 6.0),
                 ('E', 'F', 1 / 6.0)]
        Panswer = nx.Graph()
        Panswer.add_weighted_edges_from(edges)
        P = bipartite.weighted_projected_graph(self.G, 'ABCDEF', ratio=True)
        assert_edges_equal(list(P.edges()), Panswer.edges())
        for u, v in list(P.edges()):
            assert_equal(P[u][v]['weight'], Panswer[u][v]['weight'])

        edges = [('A', 'B', 3 / 3.0),
                 ('A', 'E', 1 / 3.0),
                 ('A', 'C', 1 / 3.0),
                 ('A', 'D', 1 / 3.0),
                 ('B', 'E', 1 / 3.0),
                 ('B', 'C', 1 / 3.0),
                 ('B', 'D', 1 / 3.0),
                 ('C', 'D', 1 / 3.0)]
        Panswer = nx.Graph()
        Panswer.add_weighted_edges_from(edges)
        P = bipartite.weighted_projected_graph(self.N, 'ABCDE', ratio=True)
        assert_edges_equal(list(P.edges()), Panswer.edges())
        for u, v in list(P.edges()):
            assert_equal(P[u][v]['weight'], Panswer[u][v]['weight'])

    def test_project_weighted_overlap(self):
        edges = [('A', 'B', 2 / 2.0),
                 ('A', 'C', 1 / 1.0),
                 ('B', 'C', 1 / 1.0),
                 ('B', 'D', 1 / 1.0),
                 ('B', 'E', 2 / 3.0),
                 ('E', 'F', 1 / 1.0)]
        Panswer = nx.Graph()
        Panswer.add_weighted_edges_from(edges)
        P = bipartite.overlap_weighted_projected_graph(self.G, 'ABCDEF', jaccard=False)
        assert_edges_equal(list(P.edges()), Panswer.edges())
        for u, v in list(P.edges()):
            assert_equal(P[u][v]['weight'], Panswer[u][v]['weight'])

        edges = [('A', 'B', 3 / 3.0),
                 ('A', 'E', 1 / 1.0),
                 ('A', 'C', 1 / 1.0),
                 ('A', 'D', 1 / 1.0),
                 ('B', 'E', 1 / 1.0),
                 ('B', 'C', 1 / 1.0),
                 ('B', 'D', 1 / 1.0),
                 ('C', 'D', 1 / 1.0)]
        Panswer = nx.Graph()
        Panswer.add_weighted_edges_from(edges)
        P = bipartite.overlap_weighted_projected_graph(self.N, 'ABCDE', jaccard=False)
        assert_edges_equal(list(P.edges()), Panswer.edges())
        for u, v in list(P.edges()):
            assert_equal(P[u][v]['weight'], Panswer[u][v]['weight'])

    def test_project_weighted_jaccard(self):
        edges = [('A', 'B', 2 / 5.0),
                 ('A', 'C', 1 / 2.0),
                 ('B', 'C', 1 / 5.0),
                 ('B', 'D', 1 / 5.0),
                 ('B', 'E', 2 / 6.0),
                 ('E', 'F', 1 / 3.0)]
        Panswer = nx.Graph()
        Panswer.add_weighted_edges_from(edges)
        P = bipartite.overlap_weighted_projected_graph(self.G, 'ABCDEF')
        assert_edges_equal(list(P.edges()), Panswer.edges())
        for u, v in list(P.edges()):
            assert_equal(P[u][v]['weight'], Panswer[u][v]['weight'])

        edges = [('A', 'B', 3 / 3.0),
                 ('A', 'E', 1 / 3.0),
                 ('A', 'C', 1 / 3.0),
                 ('A', 'D', 1 / 3.0),
                 ('B', 'E', 1 / 3.0),
                 ('B', 'C', 1 / 3.0),
                 ('B', 'D', 1 / 3.0),
                 ('C', 'D', 1 / 1.0)]
        Panswer = nx.Graph()
        Panswer.add_weighted_edges_from(edges)
        P = bipartite.overlap_weighted_projected_graph(self.N, 'ABCDE')
        assert_edges_equal(list(P.edges()), Panswer.edges())
        for u, v in P.edges():
            assert_equal(P[u][v]['weight'], Panswer[u][v]['weight'])

    def test_generic_weighted_projected_graph_simple(self):
        def shared(G, u, v):
            return len(set(G[u]) & set(G[v]))
        B = nx.path_graph(5)
        G = bipartite.generic_weighted_projected_graph(B, [0, 2, 4], weight_function=shared)
        assert_nodes_equal(list(G), [0, 2, 4])
        assert_edges_equal(list(list(G.edges(data=True))),
                           [(0, 2, {'weight': 1}), (2, 4, {'weight': 1})])

        G = bipartite.generic_weighted_projected_graph(B, [0, 2, 4])
        assert_nodes_equal(list(G), [0, 2, 4])
        assert_edges_equal(list(list(G.edges(data=True))),
                           [(0, 2, {'weight': 1}), (2, 4, {'weight': 1})])
        B = nx.DiGraph()
        nx.add_path(B, range(5))
        G = bipartite.generic_weighted_projected_graph(B, [0, 2, 4])
        assert_nodes_equal(list(G), [0, 2, 4])
        assert_edges_equal(list(G.edges(data=True)),
                           [(0, 2, {'weight': 1}), (2, 4, {'weight': 1})])

    def test_generic_weighted_projected_graph_custom(self):
        def jaccard(G, u, v):
            unbrs = set(G[u])
            vnbrs = set(G[v])
            return float(len(unbrs & vnbrs)) / len(unbrs | vnbrs)

        def my_weight(G, u, v, weight='weight'):
            w = 0
            for nbr in set(G[u]) & set(G[v]):
                w += G.edges[u, nbr].get(weight, 1) + G.edges[v, nbr].get(weight, 1)
            return w
        B = nx.bipartite.complete_bipartite_graph(2, 2)
        for i, (u, v) in enumerate(B.edges()):
            B.edges[u, v]['weight'] = i + 1
        G = bipartite.generic_weighted_projected_graph(B, [0, 1],
                                                       weight_function=jaccard)
        assert_edges_equal(list(G.edges(data=True)), [(0, 1, {'weight': 1.0})])
        G = bipartite.generic_weighted_projected_graph(B, [0, 1],
                                                       weight_function=my_weight)
        assert_edges_equal(list(G.edges(data=True)), [(0, 1, {'weight': 10})])
        G = bipartite.generic_weighted_projected_graph(B, [0, 1])
        assert_edges_equal(list(G.edges(data=True)), [(0, 1, {'weight': 2})])
