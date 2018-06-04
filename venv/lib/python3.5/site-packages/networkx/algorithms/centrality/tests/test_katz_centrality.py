# -*- coding: utf-8 -*-
import math

import networkx as nx
from nose import SkipTest
from nose.tools import assert_almost_equal, assert_equal, raises


class TestKatzCentrality(object):

    def test_K5(self):
        """Katz centrality: K5"""
        G = nx.complete_graph(5)
        alpha = 0.1
        b = nx.katz_centrality(G, alpha)
        v = math.sqrt(1 / 5.0)
        b_answer = dict.fromkeys(G, v)
        for n in sorted(G):
            assert_almost_equal(b[n], b_answer[n])
        nstart = dict([(n, 1) for n in G])
        b = nx.katz_centrality(G, alpha, nstart=nstart)
        for n in sorted(G):
            assert_almost_equal(b[n], b_answer[n])

    def test_P3(self):
        """Katz centrality: P3"""
        alpha = 0.1
        G = nx.path_graph(3)
        b_answer = {0: 0.5598852584152165, 1: 0.6107839182711449,
                    2: 0.5598852584152162}
        b = nx.katz_centrality(G, alpha)
        for n in sorted(G):
            assert_almost_equal(b[n], b_answer[n], places=4)

    @raises(nx.PowerIterationFailedConvergence)
    def test_maxiter(self):
        alpha = 0.1
        G = nx.path_graph(3)
        max_iter = 0
        try:
            b = nx.katz_centrality(G, alpha, max_iter=max_iter)
        except nx.NetworkXError as e:
            assert str(max_iter) in e.args[0], "max_iter value not in error msg"
            raise  # So that the decorater sees the exception.

    def test_beta_as_scalar(self):
        alpha = 0.1
        beta = 0.1
        b_answer = {0: 0.5598852584152165, 1: 0.6107839182711449,
                    2: 0.5598852584152162}
        G = nx.path_graph(3)
        b = nx.katz_centrality(G, alpha, beta)
        for n in sorted(G):
            assert_almost_equal(b[n], b_answer[n], places=4)

    def test_beta_as_dict(self):
        alpha = 0.1
        beta = {0: 1.0, 1: 1.0, 2: 1.0}
        b_answer = {0: 0.5598852584152165, 1: 0.6107839182711449,
                    2: 0.5598852584152162}
        G = nx.path_graph(3)
        b = nx.katz_centrality(G, alpha, beta)
        for n in sorted(G):
            assert_almost_equal(b[n], b_answer[n], places=4)

    def test_multiple_alpha(self):
        alpha_list = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6]
        for alpha in alpha_list:
            b_answer = {0.1: {0: 0.5598852584152165, 1: 0.6107839182711449,
                              2: 0.5598852584152162},
                        0.2: {0: 0.5454545454545454, 1: 0.6363636363636365,
                              2: 0.5454545454545454},
                        0.3: {0: 0.5333964609104419, 1: 0.6564879518897746,
                              2: 0.5333964609104419},
                        0.4: {0: 0.5232045649263551, 1: 0.6726915834767423,
                              2: 0.5232045649263551},
                        0.5: {0: 0.5144957746691622, 1: 0.6859943117075809,
                              2: 0.5144957746691622},
                        0.6: {0: 0.5069794004195823, 1: 0.6970966755769258,
                              2: 0.5069794004195823}}
            G = nx.path_graph(3)
            b = nx.katz_centrality(G, alpha)
            for n in sorted(G):
                assert_almost_equal(b[n], b_answer[alpha][n], places=4)

    @raises(nx.NetworkXException)
    def test_multigraph(self):
        e = nx.katz_centrality(nx.MultiGraph(), 0.1)

    def test_empty(self):
        e = nx.katz_centrality(nx.Graph(), 0.1)
        assert_equal(e, {})

    @raises(nx.NetworkXException)
    def test_bad_beta(self):
        G = nx.Graph([(0, 1)])
        beta = {0: 77}
        e = nx.katz_centrality(G, 0.1, beta=beta)

    @raises(nx.NetworkXException)
    def test_bad_beta_numbe(self):
        G = nx.Graph([(0, 1)])
        e = nx.katz_centrality(G, 0.1, beta='foo')


class TestKatzCentralityNumpy(object):
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
        """Katz centrality: K5"""
        G = nx.complete_graph(5)
        alpha = 0.1
        b = nx.katz_centrality(G, alpha)
        v = math.sqrt(1 / 5.0)
        b_answer = dict.fromkeys(G, v)
        for n in sorted(G):
            assert_almost_equal(b[n], b_answer[n])
        nstart = dict([(n, 1) for n in G])
        b = nx.eigenvector_centrality_numpy(G)
        for n in sorted(G):
            assert_almost_equal(b[n], b_answer[n], places=3)

    def test_P3(self):
        """Katz centrality: P3"""
        alpha = 0.1
        G = nx.path_graph(3)
        b_answer = {0: 0.5598852584152165, 1: 0.6107839182711449,
                    2: 0.5598852584152162}
        b = nx.katz_centrality_numpy(G, alpha)
        for n in sorted(G):
            assert_almost_equal(b[n], b_answer[n], places=4)

    def test_beta_as_scalar(self):
        alpha = 0.1
        beta = 0.1
        b_answer = {0: 0.5598852584152165, 1: 0.6107839182711449,
                    2: 0.5598852584152162}
        G = nx.path_graph(3)
        b = nx.katz_centrality_numpy(G, alpha, beta)
        for n in sorted(G):
            assert_almost_equal(b[n], b_answer[n], places=4)

    def test_beta_as_dict(self):
        alpha = 0.1
        beta = {0: 1.0, 1: 1.0, 2: 1.0}
        b_answer = {0: 0.5598852584152165, 1: 0.6107839182711449,
                    2: 0.5598852584152162}
        G = nx.path_graph(3)
        b = nx.katz_centrality_numpy(G, alpha, beta)
        for n in sorted(G):
            assert_almost_equal(b[n], b_answer[n], places=4)

    def test_multiple_alpha(self):
        alpha_list = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6]
        for alpha in alpha_list:
            b_answer = {0.1: {0: 0.5598852584152165, 1: 0.6107839182711449,
                              2: 0.5598852584152162},
                        0.2: {0: 0.5454545454545454, 1: 0.6363636363636365,
                              2: 0.5454545454545454},
                        0.3: {0: 0.5333964609104419, 1: 0.6564879518897746,
                              2: 0.5333964609104419},
                        0.4: {0: 0.5232045649263551, 1: 0.6726915834767423,
                              2: 0.5232045649263551},
                        0.5: {0: 0.5144957746691622, 1: 0.6859943117075809,
                              2: 0.5144957746691622},
                        0.6: {0: 0.5069794004195823, 1: 0.6970966755769258,
                              2: 0.5069794004195823}}
            G = nx.path_graph(3)
            b = nx.katz_centrality_numpy(G, alpha)
            for n in sorted(G):
                assert_almost_equal(b[n], b_answer[alpha][n], places=4)

    @raises(nx.NetworkXException)
    def test_multigraph(self):
        e = nx.katz_centrality(nx.MultiGraph(), 0.1)

    def test_empty(self):
        e = nx.katz_centrality(nx.Graph(), 0.1)
        assert_equal(e, {})

    @raises(nx.NetworkXException)
    def test_bad_beta(self):
        G = nx.Graph([(0, 1)])
        beta = {0: 77}
        e = nx.katz_centrality_numpy(G, 0.1, beta=beta)

    @raises(nx.NetworkXException)
    def test_bad_beta_numbe(self):
        G = nx.Graph([(0, 1)])
        e = nx.katz_centrality_numpy(G, 0.1, beta='foo')

    def test_K5_unweighted(self):
        """Katz centrality: K5"""
        G = nx.complete_graph(5)
        alpha = 0.1
        b = nx.katz_centrality(G, alpha, weight=None)
        v = math.sqrt(1 / 5.0)
        b_answer = dict.fromkeys(G, v)
        for n in sorted(G):
            assert_almost_equal(b[n], b_answer[n])
        nstart = dict([(n, 1) for n in G])
        b = nx.eigenvector_centrality_numpy(G, weight=None)
        for n in sorted(G):
            assert_almost_equal(b[n], b_answer[n], places=3)

    def test_P3_unweighted(self):
        """Katz centrality: P3"""
        alpha = 0.1
        G = nx.path_graph(3)
        b_answer = {0: 0.5598852584152165, 1: 0.6107839182711449,
                    2: 0.5598852584152162}
        b = nx.katz_centrality_numpy(G, alpha, weight=None)
        for n in sorted(G):
            assert_almost_equal(b[n], b_answer[n], places=4)


class TestKatzCentralityDirected(object):
    def setUp(self):
        G = nx.DiGraph()
        edges = [(1, 2), (1, 3), (2, 4), (3, 2), (3, 5), (4, 2), (4, 5),
                 (4, 6), (5, 6), (5, 7), (5, 8), (6, 8), (7, 1), (7, 5),
                 (7, 8), (8, 6), (8, 7)]
        G.add_edges_from(edges, weight=2.0)
        self.G = G.reverse()
        self.G.alpha = 0.1
        self.G.evc = [
            0.3289589783189635,
            0.2832077296243516,
            0.3425906003685471,
            0.3970420865198392,
            0.41074871061646284,
            0.272257430756461,
            0.4201989685435462,
            0.34229059218038554,
        ]

        H = nx.DiGraph(edges)
        self.H = G.reverse()
        self.H.alpha = 0.1
        self.H.evc = [
            0.3289589783189635,
            0.2832077296243516,
            0.3425906003685471,
            0.3970420865198392,
            0.41074871061646284,
            0.272257430756461,
            0.4201989685435462,
            0.34229059218038554,
        ]

    def test_katz_centrality_weighted(self):
        G = self.G
        alpha = self.G.alpha
        p = nx.katz_centrality(G, alpha, weight='weight')
        for (a, b) in zip(list(p.values()), self.G.evc):
            assert_almost_equal(a, b)

    def test_katz_centrality_unweighted(self):
        H = self.H
        alpha = self.H.alpha
        p = nx.katz_centrality(H, alpha, weight='weight')
        for (a, b) in zip(list(p.values()), self.H.evc):
            assert_almost_equal(a, b)


class TestKatzCentralityDirectedNumpy(TestKatzCentralityDirected):
    numpy = 1  # nosetests attribute, use nosetests -a 'not numpy' to skip test

    @classmethod
    def setupClass(cls):
        global np
        try:
            import numpy as np
            import scipy
        except ImportError:
            raise SkipTest('SciPy not available.')

    def test_katz_centrality_weighted(self):
        G = self.G
        alpha = self.G.alpha
        p = nx.katz_centrality_numpy(G, alpha, weight='weight')
        for (a, b) in zip(list(p.values()), self.G.evc):
            assert_almost_equal(a, b)

    def test_katz_centrality_unweighted(self):
        H = self.H
        alpha = self.H.alpha
        p = nx.katz_centrality_numpy(H, alpha, weight='weight')
        for (a, b) in zip(list(p.values()), self.H.evc):
            assert_almost_equal(a, b)


class TestKatzEigenvectorVKatz(object):
    numpy = 1  # nosetests attribute, use nosetests -a 'not numpy' to skip test

    @classmethod
    def setupClass(cls):
        global np
        global eigvals
        try:
            import numpy as np
            import scipy
            from numpy.linalg import eigvals
        except ImportError:
            raise SkipTest('SciPy not available.')

    def test_eigenvector_v_katz_random(self):
        G = nx.gnp_random_graph(10, 0.5, seed=1234)
        l = float(max(eigvals(nx.adjacency_matrix(G).todense())))
        e = nx.eigenvector_centrality_numpy(G)
        k = nx.katz_centrality_numpy(G, 1.0 / l)
        for n in G:
            assert_almost_equal(e[n], k[n])
