#!/usr/bin/env python
import random

import networkx
from nose.tools import *
from nose import SkipTest

# Example from
# A. Langville and C. Meyer, "A survey of eigenvector methods of web
# information retrieval."  http://citeseer.ist.psu.edu/713792.html


class TestPageRank(object):

    @classmethod
    def setupClass(cls):
        global numpy
        try:
            import numpy
        except ImportError:
            raise SkipTest('NumPy not available.')

    def setUp(self):
        G = networkx.DiGraph()
        edges = [(1, 2), (1, 3),
                 # 2 is a dangling node
                 (3, 1), (3, 2), (3, 5),
                 (4, 5), (4, 6),
                 (5, 4), (5, 6),
                 (6, 4)]
        G.add_edges_from(edges)
        self.G = G
        self.G.pagerank = dict(zip(sorted(G),
                                   [0.03721197, 0.05395735, 0.04150565,
                                    0.37508082, 0.20599833, 0.28624589]))
        self.dangling_node_index = 1
        self.dangling_edges = {1: 2, 2: 3,
                               3: 0, 4: 0, 5: 0, 6: 0}
        self.G.dangling_pagerank = dict(zip(sorted(G),
                                            [0.10844518, 0.18618601, 0.0710892,
                                             0.2683668, 0.15919783, 0.20671497]))

    def test_pagerank(self):
        G = self.G
        p = networkx.pagerank(G, alpha=0.9, tol=1.e-08)
        for n in G:
            assert_almost_equal(p[n], G.pagerank[n], places=4)

        nstart = dict((n, random.random()) for n in G)
        p = networkx.pagerank(G, alpha=0.9, tol=1.e-08, nstart=nstart)
        for n in G:
            assert_almost_equal(p[n], G.pagerank[n], places=4)

    @raises(networkx.PowerIterationFailedConvergence)
    def test_pagerank_max_iter(self):
        networkx.pagerank(self.G, max_iter=0)

    def test_numpy_pagerank(self):
        G = self.G
        p = networkx.pagerank_numpy(G, alpha=0.9)
        for n in G:
            assert_almost_equal(p[n], G.pagerank[n], places=4)
        personalize = dict((n, random.random()) for n in G)
        p = networkx.pagerank_numpy(G, alpha=0.9, personalization=personalize)

    def test_google_matrix(self):
        G = self.G
        M = networkx.google_matrix(G, alpha=0.9, nodelist=sorted(G))
        e, ev = numpy.linalg.eig(M.T)
        p = numpy.array(ev[:, 0] / ev[:, 0].sum())[:, 0]
        for (a, b) in zip(p, self.G.pagerank.values()):
            assert_almost_equal(a, b)

    def test_personalization(self):
        G = networkx.complete_graph(4)
        personalize = {0: 1, 1: 1, 2: 4, 3: 4}
        answer = {0: 0.23246732615667579, 1: 0.23246732615667579, 2: 0.267532673843324, 3: 0.2675326738433241}
        p = networkx.pagerank(G, alpha=0.85, personalization=personalize)
        for n in G:
            assert_almost_equal(p[n], answer[n], places=4)

    def test_zero_personalization_vector(self):
        G = networkx.complete_graph(4)
        personalize = {0: 0, 1: 0, 2: 0, 3: 0}
        assert_raises(ZeroDivisionError, networkx.pagerank, G,
                      personalization=personalize)

    def test_one_nonzero_personalization_value(self):
        G = networkx.complete_graph(4)
        personalize = {0: 0, 1: 0, 2: 0, 3: 1}
        answer = {0: 0.22077931820379187, 1: 0.22077931820379187, 2: 0.22077931820379187, 3: 0.3376620453886241}
        p = networkx.pagerank(G, alpha=0.85, personalization=personalize)
        for n in G:
            assert_almost_equal(p[n], answer[n], places=4)

    def test_incomplete_personalization(self):
        G = networkx.complete_graph(4)
        personalize = {3: 1}
        answer = {0: 0.22077931820379187, 1: 0.22077931820379187, 2: 0.22077931820379187, 3: 0.3376620453886241}
        p = networkx.pagerank(G, alpha=0.85, personalization=personalize)
        for n in G:
            assert_almost_equal(p[n], answer[n], places=4)

    def test_dangling_matrix(self):
        """
        Tests that the google_matrix doesn't change except for the dangling
        nodes.
        """
        G = self.G
        dangling = self.dangling_edges
        dangling_sum = float(sum(dangling.values()))
        M1 = networkx.google_matrix(G, personalization=dangling)
        M2 = networkx.google_matrix(G, personalization=dangling,
                                    dangling=dangling)
        for i in range(len(G)):
            for j in range(len(G)):
                if i == self.dangling_node_index and (j + 1) in dangling:
                    assert_almost_equal(M2[i, j],
                                        dangling[j + 1] / dangling_sum,
                                        places=4)
                else:
                    assert_almost_equal(M2[i, j], M1[i, j], places=4)

    def test_dangling_pagerank(self):
        pr = networkx.pagerank(self.G, dangling=self.dangling_edges)
        for n in self.G:
            assert_almost_equal(pr[n], self.G.dangling_pagerank[n], places=4)

    def test_dangling_numpy_pagerank(self):
        pr = networkx.pagerank_numpy(self.G, dangling=self.dangling_edges)
        for n in self.G:
            assert_almost_equal(pr[n], self.G.dangling_pagerank[n], places=4)

    def test_empty(self):
        G = networkx.Graph()
        assert_equal(networkx.pagerank(G), {})
        assert_equal(networkx.pagerank_numpy(G), {})
        assert_equal(networkx.google_matrix(G).shape, (0, 0))


class TestPageRankScipy(TestPageRank):

    @classmethod
    def setupClass(cls):
        global scipy
        try:
            import scipy
        except ImportError:
            raise SkipTest('SciPy not available.')

    def test_scipy_pagerank(self):
        G = self.G
        p = networkx.pagerank_scipy(G, alpha=0.9, tol=1.e-08)
        for n in G:
            assert_almost_equal(p[n], G.pagerank[n], places=4)
        personalize = dict((n, random.random()) for n in G)
        p = networkx.pagerank_scipy(G, alpha=0.9, tol=1.e-08,
                                    personalization=personalize)

    @raises(networkx.PowerIterationFailedConvergence)
    def test_scipy_pagerank_max_iter(self):
        networkx.pagerank_scipy(self.G, max_iter=0)

    def test_dangling_scipy_pagerank(self):
        pr = networkx.pagerank_scipy(self.G, dangling=self.dangling_edges)
        for n in self.G:
            assert_almost_equal(pr[n], self.G.dangling_pagerank[n], places=4)

    def test_empty_scipy(self):
        G = networkx.Graph()
        assert_equal(networkx.pagerank_scipy(G), {})
