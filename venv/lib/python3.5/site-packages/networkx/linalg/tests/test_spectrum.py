from nose import SkipTest

import networkx as nx
from networkx.generators.degree_seq import havel_hakimi_graph


class TestSpectrum(object):
    numpy = 1  # nosetests attribute, use nosetests -a 'not numpy' to skip test

    @classmethod
    def setupClass(cls):
        global numpy
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
        self.P = nx.path_graph(3)
        self.WG = nx.Graph((u, v, {'weight': 0.5, 'other': 0.3})
                           for (u, v) in self.G.edges())
        self.WG.add_node(4)
        self.DG = nx.DiGraph()
        nx.add_path(self.DG, [0, 1, 2])

    def test_laplacian_spectrum(self):
        "Laplacian eigenvalues"
        evals = numpy.array([0, 0, 1, 3, 4])
        e = sorted(nx.laplacian_spectrum(self.G))
        assert_almost_equal(e, evals)
        e = sorted(nx.laplacian_spectrum(self.WG, weight=None))
        assert_almost_equal(e, evals)
        e = sorted(nx.laplacian_spectrum(self.WG))
        assert_almost_equal(e, 0.5 * evals)
        e = sorted(nx.laplacian_spectrum(self.WG, weight='other'))
        assert_almost_equal(e, 0.3 * evals)

    def test_adjacency_spectrum(self):
        "Adjacency eigenvalues"
        evals = numpy.array([-numpy.sqrt(2), 0, numpy.sqrt(2)])
        e = sorted(nx.adjacency_spectrum(self.P))
        assert_almost_equal(e, evals)

    def test_modularity_spectrum(self):
        "Modularity eigenvalues"
        evals = numpy.array([-1.5, 0., 0.])
        e = sorted(nx.modularity_spectrum(self.P))
        assert_almost_equal(e, evals)
        # Directed modularity eigenvalues
        evals = numpy.array([-0.5, 0., 0.])
        e = sorted(nx.modularity_spectrum(self.DG))
        assert_almost_equal(e, evals)
