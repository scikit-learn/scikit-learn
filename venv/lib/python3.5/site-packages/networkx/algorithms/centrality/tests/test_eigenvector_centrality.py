#!/usr/bin/env python
import math
from nose import SkipTest
from nose.tools import *
import networkx as nx


class TestEigenvectorCentrality(object):
    numpy = 1  # nosetests attribute, use nosetests -a 'not numpy' to skip test

    @classmethod
    def setupClass(cls):
        global np
        try:
            import numpy as np
            import scipy
        except ImportError:
            raise SkipTest('SciPy not available.')

    def test_K5(self):
        """Eigenvector centrality: K5"""
        G = nx.complete_graph(5)
        b = nx.eigenvector_centrality(G)
        v = math.sqrt(1 / 5.0)
        b_answer = dict.fromkeys(G, v)
        for n in sorted(G):
            assert_almost_equal(b[n], b_answer[n])
        nstart = dict([(n, 1) for n in G])
        b = nx.eigenvector_centrality(G, nstart=nstart)
        for n in sorted(G):
            assert_almost_equal(b[n], b_answer[n])

        b = nx.eigenvector_centrality_numpy(G)
        for n in sorted(G):
            assert_almost_equal(b[n], b_answer[n], places=3)

    def test_P3(self):
        """Eigenvector centrality: P3"""
        G = nx.path_graph(3)
        b_answer = {0: 0.5, 1: 0.7071, 2: 0.5}
        b = nx.eigenvector_centrality_numpy(G)
        for n in sorted(G):
            assert_almost_equal(b[n], b_answer[n], places=4)
        b = nx.eigenvector_centrality(G)
        for n in sorted(G):
            assert_almost_equal(b[n], b_answer[n], places=4)

    def test_P3_unweighted(self):
        """Eigenvector centrality: P3"""
        G = nx.path_graph(3)
        b_answer = {0: 0.5, 1: 0.7071, 2: 0.5}
        b = nx.eigenvector_centrality_numpy(G, weight=None)
        for n in sorted(G):
            assert_almost_equal(b[n], b_answer[n], places=4)

    @raises(nx.PowerIterationFailedConvergence)
    def test_maxiter(self):
        G = nx.path_graph(3)
        b = nx.eigenvector_centrality(G, max_iter=0)


class TestEigenvectorCentralityDirected(object):
    numpy = 1  # nosetests attribute, use nosetests -a 'not numpy' to skip test

    @classmethod
    def setupClass(cls):
        global np
        try:
            import numpy as np
            import scipy
        except ImportError:
            raise SkipTest('SciPy not available.')

    def setUp(self):

        G = nx.DiGraph()

        edges = [(1, 2), (1, 3), (2, 4), (3, 2), (3, 5), (4, 2), (4, 5), (4, 6),
                 (5, 6), (5, 7), (5, 8), (6, 8), (7, 1), (7, 5),
                 (7, 8), (8, 6), (8, 7)]

        G.add_edges_from(edges, weight=2.0)
        self.G = G.reverse()
        self.G.evc = [0.25368793,  0.19576478,  0.32817092,  0.40430835,
                      0.48199885, 0.15724483,  0.51346196,  0.32475403]

        H = nx.DiGraph()

        edges = [(1, 2), (1, 3), (2, 4), (3, 2), (3, 5), (4, 2), (4, 5), (4, 6),
                 (5, 6), (5, 7), (5, 8), (6, 8), (7, 1), (7, 5),
                 (7, 8), (8, 6), (8, 7)]

        G.add_edges_from(edges)
        self.H = G.reverse()
        self.H.evc = [0.25368793,  0.19576478,  0.32817092,  0.40430835,
                      0.48199885, 0.15724483,  0.51346196,  0.32475403]

    def test_eigenvector_centrality_weighted(self):
        G = self.G
        p = nx.eigenvector_centrality(G)
        for (a, b) in zip(list(p.values()), self.G.evc):
            assert_almost_equal(a, b, places=4)

    def test_eigenvector_centrality_weighted_numpy(self):
        G = self.G
        p = nx.eigenvector_centrality_numpy(G)
        for (a, b) in zip(list(p.values()), self.G.evc):
            assert_almost_equal(a, b)

    def test_eigenvector_centrality_unweighted(self):
        G = self.H
        p = nx.eigenvector_centrality(G)
        for (a, b) in zip(list(p.values()), self.G.evc):
            assert_almost_equal(a, b, places=4)

    def test_eigenvector_centrality_unweighted_numpy(self):
        G = self.H
        p = nx.eigenvector_centrality_numpy(G)
        for (a, b) in zip(list(p.values()), self.G.evc):
            assert_almost_equal(a, b)


class TestEigenvectorCentralityExceptions(object):
    numpy = 1  # nosetests attribute, use nosetests -a 'not numpy' to skip test

    @classmethod
    def setupClass(cls):
        global np
        try:
            import numpy as np
            import scipy
        except ImportError:
            raise SkipTest('SciPy not available.')
    numpy = 1  # nosetests attribute, use nosetests -a 'not numpy' to skip test

    @raises(nx.NetworkXException)
    def test_multigraph(self):
        e = nx.eigenvector_centrality(nx.MultiGraph())

    @raises(nx.NetworkXException)
    def test_multigraph_numpy(self):
        e = nx.eigenvector_centrality_numpy(nx.MultiGraph())

    @raises(nx.NetworkXException)
    def test_empty(self):
        e = nx.eigenvector_centrality(nx.Graph())

    @raises(nx.NetworkXException)
    def test_empty_numpy(self):
        e = nx.eigenvector_centrality_numpy(nx.Graph())
