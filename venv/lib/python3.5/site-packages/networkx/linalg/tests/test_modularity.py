from nose import SkipTest

import networkx as nx
from networkx.generators.degree_seq import havel_hakimi_graph


class TestModularity(object):
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
        # Graph used as an example in Sec. 4.1 of Langville and Meyer,
        # "Google's PageRank and Beyond". (Used for test_directed_laplacian)
        self.DG = nx.DiGraph()
        self.DG.add_edges_from(((1, 2), (1, 3), (3, 1), (3, 2), (3, 5), (4, 5), (4, 6),
                                (5, 4), (5, 6), (6, 4)))

    def test_modularity(self):
        "Modularity matrix"
        B = numpy.matrix([[-1.125,  0.25,  0.25,  0.625,  0.],
                          [0.25, -0.5,  0.5, -0.25,  0.],
                          [0.25,  0.5, -0.5, -0.25,  0.],
                          [0.625, -0.25, -0.25, -0.125,  0.],
                          [0.,  0.,  0.,  0.,  0.]])

        permutation = [4, 0, 1, 2, 3]
        assert_equal(nx.modularity_matrix(self.G), B)
        assert_equal(nx.modularity_matrix(self.G, nodelist=permutation),
                     B[numpy.ix_(permutation, permutation)])

    def test_modularity_weight(self):
        "Modularity matrix with weights"
        B = numpy.matrix([[-1.125,  0.25,  0.25,  0.625,  0.],
                          [0.25, -0.5,  0.5, -0.25,  0.],
                          [0.25,  0.5, -0.5, -0.25,  0.],
                          [0.625, -0.25, -0.25, -0.125,  0.],
                          [0.,  0.,  0.,  0.,  0.]])

        G_weighted = self.G.copy()
        for n1, n2 in G_weighted.edges():
            G_weighted.edges[n1, n2]["weight"] = 0.5
        # The following test would fail in networkx 1.1
        assert_equal(nx.modularity_matrix(G_weighted), B)
        # The following test that the modularity matrix get rescaled accordingly
        assert_equal(nx.modularity_matrix(G_weighted, weight="weight"), 0.5 * B)

    def test_directed_modularity(self):
        "Directed Modularity matrix"
        B = numpy.matrix([[-0.2,  0.6,  0.8, -0.4, -0.4, -0.4],
                          [0.,  0.,  0.,  0.,  0.,  0.],
                          [0.7,  0.4, -0.3, -0.6,  0.4, -0.6],
                          [-0.2, -0.4, -0.2, -0.4,  0.6,  0.6],
                          [-0.2, -0.4, -0.2,  0.6, -0.4,  0.6],
                          [-0.1, -0.2, -0.1,  0.8, -0.2, -0.2]])
        node_permutation = [5, 1, 2, 3, 4, 6]
        idx_permutation = [4, 0, 1, 2, 3, 5]
        mm = nx.directed_modularity_matrix(self.DG,  nodelist=sorted(self.DG))
        assert_equal(mm, B)
        assert_equal(nx.directed_modularity_matrix(self.DG,
                                                   nodelist=node_permutation),
                     B[numpy.ix_(idx_permutation, idx_permutation)])
