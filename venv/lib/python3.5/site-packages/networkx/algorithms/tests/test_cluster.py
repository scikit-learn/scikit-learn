#!/usr/bin/env python
from nose.tools import *
import networkx as nx


class TestTriangles:

    def test_empty(self):
        G = nx.Graph()
        assert_equal(list(nx.triangles(G).values()), [])

    def test_path(self):
        G = nx.path_graph(10)
        assert_equal(list(nx.triangles(G).values()),
                     [0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
        assert_equal(nx.triangles(G),
                     {0: 0, 1: 0, 2: 0, 3: 0, 4: 0,
                      5: 0, 6: 0, 7: 0, 8: 0, 9: 0})

    def test_cubical(self):
        G = nx.cubical_graph()
        assert_equal(list(nx.triangles(G).values()),
                     [0, 0, 0, 0, 0, 0, 0, 0])
        assert_equal(nx.triangles(G, 1), 0)
        assert_equal(list(nx.triangles(G, [1, 2]).values()), [0, 0])
        assert_equal(nx.triangles(G, 1), 0)
        assert_equal(nx.triangles(G, [1, 2]), {1: 0, 2: 0})

    def test_k5(self):
        G = nx.complete_graph(5)
        assert_equal(list(nx.triangles(G).values()), [6, 6, 6, 6, 6])
        assert_equal(sum(nx.triangles(G).values()) / 3.0, 10)
        assert_equal(nx.triangles(G, 1), 6)
        G.remove_edge(1, 2)
        assert_equal(list(nx.triangles(G).values()), [5, 3, 3, 5, 5])
        assert_equal(nx.triangles(G, 1), 3)


class TestWeightedClustering:

    def test_clustering(self):
        G = nx.Graph()
        assert_equal(list(nx.clustering(G, weight='weight').values()), [])
        assert_equal(nx.clustering(G), {})

    def test_path(self):
        G = nx.path_graph(10)
        assert_equal(list(nx.clustering(G, weight='weight').values()),
                     [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
        assert_equal(nx.clustering(G, weight='weight'),
                     {0: 0.0, 1: 0.0, 2: 0.0, 3: 0.0, 4: 0.0,
                      5: 0.0, 6: 0.0, 7: 0.0, 8: 0.0, 9: 0.0})

    def test_cubical(self):
        G = nx.cubical_graph()
        assert_equal(list(nx.clustering(G, weight='weight').values()),
                     [0, 0, 0, 0, 0, 0, 0, 0])
        assert_equal(nx.clustering(G, 1), 0)
        assert_equal(list(nx.clustering(G, [1, 2], weight='weight').values()), [0, 0])
        assert_equal(nx.clustering(G, 1, weight='weight'), 0)
        assert_equal(nx.clustering(G, [1, 2], weight='weight'), {1: 0, 2: 0})

    def test_k5(self):
        G = nx.complete_graph(5)
        assert_equal(list(nx.clustering(G, weight='weight').values()), [1, 1, 1, 1, 1])
        assert_equal(nx.average_clustering(G, weight='weight'), 1)
        G.remove_edge(1, 2)
        assert_equal(list(nx.clustering(G, weight='weight').values()),
                     [5. / 6., 1.0, 1.0, 5. / 6., 5. / 6.])
        assert_equal(nx.clustering(G, [1, 4], weight='weight'), {1: 1.0, 4: 0.83333333333333337})

    def test_triangle_and_edge(self):
        G = nx.cycle_graph(3)
        G.add_edge(0, 4, weight=2)
        assert_equal(nx.clustering(G)[0], 1.0 / 3.0)
        assert_equal(nx.clustering(G, weight='weight')[0], 1.0 / 6.0)


class TestClustering:

    def test_clustering(self):
        G = nx.Graph()
        assert_equal(list(nx.clustering(G).values()), [])
        assert_equal(nx.clustering(G), {})

    def test_path(self):
        G = nx.path_graph(10)
        assert_equal(list(nx.clustering(G).values()),
                     [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
        assert_equal(nx.clustering(G),
                     {0: 0.0, 1: 0.0, 2: 0.0, 3: 0.0, 4: 0.0,
                      5: 0.0, 6: 0.0, 7: 0.0, 8: 0.0, 9: 0.0})

    def test_cubical(self):
        G = nx.cubical_graph()
        assert_equal(list(nx.clustering(G).values()),
                     [0, 0, 0, 0, 0, 0, 0, 0])
        assert_equal(nx.clustering(G, 1), 0)
        assert_equal(list(nx.clustering(G, [1, 2]).values()), [0, 0])
        assert_equal(nx.clustering(G, 1), 0)
        assert_equal(nx.clustering(G, [1, 2]), {1: 0, 2: 0})

    def test_k5(self):
        G = nx.complete_graph(5)
        assert_equal(list(nx.clustering(G).values()), [1, 1, 1, 1, 1])
        assert_equal(nx.average_clustering(G), 1)
        G.remove_edge(1, 2)
        assert_equal(list(nx.clustering(G).values()),
                     [5. / 6., 1.0, 1.0, 5. / 6., 5. / 6.])
        assert_equal(nx.clustering(G, [1, 4]), {1: 1.0, 4: 0.83333333333333337})


class TestTransitivity:

    def test_transitivity(self):
        G = nx.Graph()
        assert_equal(nx.transitivity(G), 0.0)

    def test_path(self):
        G = nx.path_graph(10)
        assert_equal(nx.transitivity(G), 0.0)

    def test_cubical(self):
        G = nx.cubical_graph()
        assert_equal(nx.transitivity(G), 0.0)

    def test_k5(self):
        G = nx.complete_graph(5)
        assert_equal(nx.transitivity(G), 1.0)
        G.remove_edge(1, 2)
        assert_equal(nx.transitivity(G), 0.875)

    # def test_clustering_transitivity(self):
    #     # check that weighted average of clustering is transitivity
    #     G = nx.complete_graph(5)
    #     G.remove_edge(1,2)
    #     t1=nx.transitivity(G)
    #     (cluster_d2,weights)=nx.clustering(G,weights=True)
    #     trans=[]
    #     for v in G.nodes():
    #         trans.append(cluster_d2[v]*weights[v])
    #     t2=sum(trans)
    #     assert_almost_equal(abs(t1-t2),0)


class TestSquareClustering:

    def test_clustering(self):
        G = nx.Graph()
        assert_equal(list(nx.square_clustering(G).values()), [])
        assert_equal(nx.square_clustering(G), {})

    def test_path(self):
        G = nx.path_graph(10)
        assert_equal(list(nx.square_clustering(G).values()),
                     [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
        assert_equal(nx.square_clustering(G),
                     {0: 0.0, 1: 0.0, 2: 0.0, 3: 0.0, 4: 0.0,
                      5: 0.0, 6: 0.0, 7: 0.0, 8: 0.0, 9: 0.0})

    def test_cubical(self):
        G = nx.cubical_graph()
        assert_equal(list(nx.square_clustering(G).values()),
                     [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5])
        assert_equal(list(nx.square_clustering(G, [1, 2]).values()), [0.5, 0.5])
        assert_equal(nx.square_clustering(G, [1])[1], 0.5)
        assert_equal(nx.square_clustering(G, [1, 2]), {1: 0.5, 2: 0.5})

    def test_k5(self):
        G = nx.complete_graph(5)
        assert_equal(list(nx.square_clustering(G).values()), [1, 1, 1, 1, 1])

    def test_bipartite_k5(self):
        G = nx.complete_bipartite_graph(5, 5)
        assert_equal(list(nx.square_clustering(G).values()),
                     [1, 1, 1, 1, 1, 1, 1, 1, 1, 1])

    def test_lind_square_clustering(self):
        """Test C4 for figure 1 Lind et al (2005)"""
        G = nx.Graph([(1, 2), (1, 3), (1, 6), (1, 7), (2, 4), (2, 5),
                      (3, 4), (3, 5), (6, 7), (7, 8), (6, 8), (7, 9),
                      (7, 10), (6, 11), (6, 12), (2, 13), (2, 14), (3, 15), (3, 16)])
        G1 = G.subgraph([1, 2, 3, 4, 5, 13, 14, 15, 16])
        G2 = G.subgraph([1, 6, 7, 8, 9, 10, 11, 12])
        assert_equal(nx.square_clustering(G, [1])[1], 3 / 75.0)
        assert_equal(nx.square_clustering(G1, [1])[1], 2 / 6.0)
        assert_equal(nx.square_clustering(G2, [1])[1], 1 / 5.0)


def test_average_clustering():
    G = nx.cycle_graph(3)
    G.add_edge(2, 3)
    assert_equal(nx.average_clustering(G), (1 + 1 + 1 / 3.0) / 4.0)
    assert_equal(nx.average_clustering(G, count_zeros=True), (1 + 1 + 1 / 3.0) / 4.0)
    assert_equal(nx.average_clustering(G, count_zeros=False), (1 + 1 + 1 / 3.0) / 3.0)


class TestGeneralizedDegree:

    def test_generalized_degree(self):
        G = nx.Graph()
        assert_equal(nx.generalized_degree(G), {})

    def test_path(self):
        G = nx.path_graph(5)
        assert_equal(nx.generalized_degree(G, 0), {0: 1})
        assert_equal(nx.generalized_degree(G, 1), {0: 2})

    def test_cubical(self):
        G = nx.cubical_graph()
        assert_equal(nx.generalized_degree(G, 0), {0: 3})

    def test_k5(self):
        G = nx.complete_graph(5)
        assert_equal(nx.generalized_degree(G, 0), {3: 4})
        G.remove_edge(0, 1)
        assert_equal(nx.generalized_degree(G, 0), {2: 3})
