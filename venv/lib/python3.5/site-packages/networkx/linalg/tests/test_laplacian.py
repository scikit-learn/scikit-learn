from nose import SkipTest

import networkx as nx
from networkx.generators.degree_seq import havel_hakimi_graph


class TestLaplacian(object):
    numpy = 1  # nosetests attribute, use nosetests -a 'not numpy' to skip test

    @classmethod
    def setupClass(cls):
        global numpy
        global scipy
        global assert_equal
        global assert_almost_equal
        try:
            import numpy
            import scipy
            from numpy.testing import assert_equal, assert_almost_equal
        except ImportError:
            raise SkipTest('SciPy not available.')

    def setUp(self):
        deg = [3, 2, 2, 1, 0]
        self.G = havel_hakimi_graph(deg)
        self.WG = nx.Graph((u, v, {'weight': 0.5, 'other': 0.3})
                           for (u, v) in self.G.edges())
        self.WG.add_node(4)
        self.MG = nx.MultiGraph(self.G)

        # Graph with selfloops
        self.Gsl = self.G.copy()
        for node in self.Gsl.nodes():
            self.Gsl.add_edge(node, node)

    def test_laplacian(self):
        "Graph Laplacian"
        NL = numpy.array([[3, -1, -1, -1, 0],
                          [-1,  2, -1,  0, 0],
                          [-1, -1,  2,  0, 0],
                          [-1,  0,  0,  1, 0],
                          [0,  0,  0,  0, 0]])
        WL = 0.5 * NL
        OL = 0.3 * NL
        assert_equal(nx.laplacian_matrix(self.G).todense(), NL)
        assert_equal(nx.laplacian_matrix(self.MG).todense(), NL)
        assert_equal(nx.laplacian_matrix(self.G, nodelist=[0, 1]).todense(),
                     numpy.array([[1, -1], [-1, 1]]))
        assert_equal(nx.laplacian_matrix(self.WG).todense(), WL)
        assert_equal(nx.laplacian_matrix(self.WG, weight=None).todense(), NL)
        assert_equal(nx.laplacian_matrix(self.WG, weight='other').todense(), OL)

    def test_normalized_laplacian(self):
        "Generalized Graph Laplacian"
        GL = numpy.array([[1.00, -0.408, -0.408, -0.577,  0.00],
                          [-0.408,  1.00, -0.50,  0.00, 0.00],
                          [-0.408, -0.50,  1.00,  0.00,  0.00],
                          [-0.577,  0.00,  0.00,  1.00,  0.00],
                          [0.00,  0.00,  0.00,  0.00,  0.00]])
        Lsl = numpy.array([[0.75, -0.2887, -0.2887, -0.3536,  0.],
                           [-0.2887,  0.6667, -0.3333,  0.,  0.],
                           [-0.2887, -0.3333,  0.6667,  0.,  0.],
                           [-0.3536,  0.,  0.,  0.5,  0.],
                           [0.,  0.,  0.,  0.,  0.]])

        assert_almost_equal(nx.normalized_laplacian_matrix(self.G).todense(),
                            GL, decimal=3)
        assert_almost_equal(nx.normalized_laplacian_matrix(self.MG).todense(),
                            GL, decimal=3)
        assert_almost_equal(nx.normalized_laplacian_matrix(self.WG).todense(),
                            GL, decimal=3)
        assert_almost_equal(nx.normalized_laplacian_matrix(self.WG, weight='other').todense(),
                            GL, decimal=3)
        assert_almost_equal(nx.normalized_laplacian_matrix(self.Gsl).todense(),
                            Lsl, decimal=3)

    def test_directed_laplacian(self):
        "Directed Laplacian"
        # Graph used as an example in Sec. 4.1 of Langville and Meyer,
        # "Google's PageRank and Beyond". The graph contains dangling nodes, so
        # the pagerank random walk is selected by directed_laplacian
        G = nx.DiGraph()
        G.add_edges_from(((1, 2), (1, 3), (3, 1), (3, 2), (3, 5), (4, 5), (4, 6),
                          (5, 4), (5, 6), (6, 4)))
        GL = numpy.array([[0.9833, -0.2941, -0.3882, -0.0291, -0.0231, -0.0261],
                          [-0.2941,  0.8333, -0.2339, -0.0536, -0.0589, -0.0554],
                          [-0.3882, -0.2339,  0.9833, -0.0278, -0.0896, -0.0251],
                          [-0.0291, -0.0536, -0.0278,  0.9833, -0.4878, -0.6675],
                          [-0.0231, -0.0589, -0.0896, -0.4878,  0.9833, -0.2078],
                          [-0.0261, -0.0554, -0.0251, -0.6675, -0.2078,  0.9833]])
        L = nx.directed_laplacian_matrix(G, alpha=0.9, nodelist=sorted(G))
        assert_almost_equal(L, GL, decimal=3)

        # Make the graph strongly connected, so we can use a random and lazy walk
        G.add_edges_from((((2, 5), (6, 1))))
        GL = numpy.array([[1., -0.3062, -0.4714,  0.,  0., -0.3227],
                          [-0.3062,  1., -0.1443,  0., -0.3162,  0.],
                          [-0.4714, -0.1443,  1.,  0., -0.0913,  0.],
                          [0.,  0.,  0.,  1., -0.5, -0.5],
                          [0., -0.3162, -0.0913, -0.5,  1., -0.25],
                          [-0.3227,  0.,  0., -0.5, -0.25,  1.]])
        L = nx.directed_laplacian_matrix(G, alpha=0.9, nodelist=sorted(G), walk_type='random')
        assert_almost_equal(L, GL, decimal=3)

        GL = numpy.array([[0.5, -0.1531, -0.2357,  0.,  0., -0.1614],
                          [-0.1531,  0.5, -0.0722,  0., -0.1581,  0.],
                          [-0.2357, -0.0722,  0.5,  0., -0.0456,  0.],
                          [0.,  0.,  0.,  0.5, -0.25, -0.25],
                          [0., -0.1581, -0.0456, -0.25,  0.5, -0.125],
                          [-0.1614,  0.,  0., -0.25, -0.125,  0.5]])
        L = nx.directed_laplacian_matrix(G, alpha=0.9, nodelist=sorted(G), walk_type='lazy')
        assert_almost_equal(L, GL, decimal=3)
