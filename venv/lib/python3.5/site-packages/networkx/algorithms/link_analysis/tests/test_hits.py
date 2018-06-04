#!/usr/bin/env python
from nose.tools import *
from nose import SkipTest
from nose.plugins.attrib import attr
import networkx

# Example from
# A. Langville and C. Meyer, "A survey of eigenvector methods of web
# information retrieval."  http://citeseer.ist.psu.edu/713792.html


class TestHITS:

    def setUp(self):

        G = networkx.DiGraph()

        edges = [(1, 3), (1, 5),
                 (2, 1),
                 (3, 5),
                 (5, 4), (5, 3),
                 (6, 5)]

        G.add_edges_from(edges, weight=1)
        self.G = G
        self.G.a = dict(zip(sorted(G), [0.000000, 0.000000, 0.366025,
                                        0.133975, 0.500000, 0.000000]))
        self.G.h = dict(zip(sorted(G), [0.366025, 0.000000, 0.211325,
                                        0.000000, 0.211325, 0.211325]))

    def test_hits(self):
        G = self.G
        h, a = networkx.hits(G, tol=1.e-08)
        for n in G:
            assert_almost_equal(h[n], G.h[n], places=4)
        for n in G:
            assert_almost_equal(a[n], G.a[n], places=4)

    def test_hits_nstart(self):
        G = self.G
        nstart = dict([(i, 1. / 2) for i in G])
        h, a = networkx.hits(G, nstart=nstart)

    @attr('numpy')
    def test_hits_numpy(self):
        try:
            import numpy as np
        except ImportError:
            raise SkipTest('NumPy not available.')

        G = self.G
        h, a = networkx.hits_numpy(G)
        for n in G:
            assert_almost_equal(h[n], G.h[n], places=4)
        for n in G:
            assert_almost_equal(a[n], G.a[n], places=4)

    def test_hits_scipy(self):
        try:
            import scipy as sp
        except ImportError:
            raise SkipTest('SciPy not available.')

        G = self.G
        h, a = networkx.hits_scipy(G, tol=1.e-08)
        for n in G:
            assert_almost_equal(h[n], G.h[n], places=4)
        for n in G:
            assert_almost_equal(a[n], G.a[n], places=4)

    @attr('numpy')
    def test_empty(self):
        try:
            import numpy
        except ImportError:
            raise SkipTest('numpy not available.')
        G = networkx.Graph()
        assert_equal(networkx.hits(G), ({}, {}))
        assert_equal(networkx.hits_numpy(G), ({}, {}))
        assert_equal(networkx.authority_matrix(G).shape, (0, 0))
        assert_equal(networkx.hub_matrix(G).shape, (0, 0))

    def test_empty_scipy(self):
        try:
            import scipy
        except ImportError:
            raise SkipTest('scipy not available.')
        G = networkx.Graph()
        assert_equal(networkx.hits_scipy(G), ({}, {}))

    @raises(networkx.PowerIterationFailedConvergence)
    def test_hits_not_convergent(self):
        G = self.G
        networkx.hits(G, max_iter=0)
