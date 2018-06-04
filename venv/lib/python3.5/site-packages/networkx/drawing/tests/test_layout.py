"""Unit tests for layout functions."""
from nose import SkipTest
from nose.tools import assert_almost_equal, assert_equal, \
    assert_false, assert_raises
import networkx as nx


class TestLayout(object):
    numpy = 1  # nosetests attribute, use nosetests -a 'not numpy' to skip test
    scipy = None

    @classmethod
    def setupClass(cls):
        global numpy, scipy
        try:
            import numpy
        except ImportError:
            raise SkipTest('NumPy not available.')
        try:
            import scipy
        except ImportError:
            pass    # Almost all tests still viable

    def setUp(self):
        self.Gi = nx.grid_2d_graph(5, 5)
        self.Gs = nx.Graph()
        nx.add_path(self.Gs, 'abcdef')
        self.bigG = nx.grid_2d_graph(25, 25)  # bigger than 500 nodes for sparse

    def test_spring_init_pos(self):
        # Tests GH #2448
        import math
        G = nx.Graph()
        G.add_edges_from([(0, 1), (1, 2), (2, 0), (2, 3)])

        init_pos = {0: (0.0, 0.0)}
        fixed_pos = [0]
        pos = nx.fruchterman_reingold_layout(G, pos=init_pos, fixed=fixed_pos)
        has_nan = any(math.isnan(c) for coords in pos.values() for c in coords)
        assert_false(has_nan, 'values should not be nan')

    def test_smoke_empty_graph(self):
        G = []
        vpos = nx.random_layout(G)
        vpos = nx.circular_layout(G)
        vpos = nx.spring_layout(G)
        vpos = nx.fruchterman_reingold_layout(G)
        vpos = nx.spectral_layout(G)
        vpos = nx.shell_layout(G)
        if self.scipy is not None:
            vpos = nx.kamada_kawai_layout(G)

    def test_smoke_int(self):
        G = self.Gi
        vpos = nx.random_layout(G)
        vpos = nx.circular_layout(G)
        vpos = nx.spring_layout(G)
        vpos = nx.fruchterman_reingold_layout(G)
        vpos = nx.fruchterman_reingold_layout(self.bigG)
        vpos = nx.spectral_layout(G)
        vpos = nx.spectral_layout(G.to_directed())
        vpos = nx.spectral_layout(self.bigG)
        vpos = nx.spectral_layout(self.bigG.to_directed())
        vpos = nx.shell_layout(G)
        if self.scipy is not None:
            vpos = nx.kamada_kawai_layout(G)

    def test_smoke_string(self):
        G = self.Gs
        vpos = nx.random_layout(G)
        vpos = nx.circular_layout(G)
        vpos = nx.spring_layout(G)
        vpos = nx.fruchterman_reingold_layout(G)
        vpos = nx.spectral_layout(G)
        vpos = nx.shell_layout(G)
        if self.scipy is not None:
            vpos = nx.kamada_kawai_layout(G)

    def check_scale_and_center(self, pos, scale, center):
        center = numpy.array(center)
        low = center - scale
        hi = center + scale
        vpos = numpy.array(list(pos.values()))
        length = vpos.max(0) - vpos.min(0)
        assert (length <= 2 * scale).all()
        assert (vpos >= low).all()
        assert (vpos <= hi).all()

    def test_scale_and_center_arg(self):
        sc = self.check_scale_and_center
        c = (4, 5)
        G = nx.complete_graph(9)
        G.add_node(9)
        sc(nx.random_layout(G, center=c), scale=0.5, center=(4.5, 5.5))
        # rest can have 2*scale length: [-scale, scale]
        sc(nx.spring_layout(G, scale=2, center=c), scale=2, center=c)
        sc(nx.spectral_layout(G, scale=2, center=c), scale=2, center=c)
        sc(nx.circular_layout(G, scale=2, center=c), scale=2, center=c)
        sc(nx.shell_layout(G, scale=2, center=c), scale=2, center=c)
        if self.scipy is not None:
            sc(nx.kamada_kawai_layout(G, scale=2, center=c), scale=2, center=c)

    def test_default_scale_and_center(self):
        sc = self.check_scale_and_center
        c = (0, 0)
        G = nx.complete_graph(9)
        G.add_node(9)
        sc(nx.random_layout(G), scale=0.5, center=(0.5, 0.5))
        sc(nx.spring_layout(G), scale=1, center=c)
        sc(nx.spectral_layout(G), scale=1, center=c)
        sc(nx.circular_layout(G), scale=1, center=c)
        sc(nx.shell_layout(G), scale=1, center=c)
        if self.scipy is not None:
            sc(nx.kamada_kawai_layout(G), scale=1, center=c)

    def test_adjacency_interface_numpy(self):
        A = nx.to_numpy_matrix(self.Gs)
        pos = nx.drawing.layout._fruchterman_reingold(A)
        assert_equal(pos.shape, (6, 2))
        pos = nx.drawing.layout._fruchterman_reingold(A, dim=3)
        assert_equal(pos.shape, (6, 3))

    def test_adjacency_interface_scipy(self):
        try:
            import scipy
        except ImportError:
            raise SkipTest('scipy not available.')
        A = nx.to_scipy_sparse_matrix(self.Gs, dtype='d')
        pos = nx.drawing.layout._sparse_fruchterman_reingold(A)
        assert_equal(pos.shape, (6, 2))
        pos = nx.drawing.layout._sparse_spectral(A)
        assert_equal(pos.shape, (6, 2))
        pos = nx.drawing.layout._sparse_fruchterman_reingold(A, dim=3)
        assert_equal(pos.shape, (6, 3))

    def test_single_nodes(self):
        G = nx.path_graph(1)
        vpos = nx.shell_layout(G)
        assert_false(vpos[0].any())
        G = nx.path_graph(3)
        vpos = nx.shell_layout(G, [[0], [1, 2]])
        assert_false(vpos[0].any())

    def test_smoke_initial_pos_fruchterman_reingold(self):
        pos = nx.circular_layout(self.Gi)
        npos = nx.fruchterman_reingold_layout(self.Gi, pos=pos)

    def test_fixed_node_fruchterman_reingold(self):
        # Dense version (numpy based)
        pos = nx.circular_layout(self.Gi)
        npos = nx.fruchterman_reingold_layout(self.Gi, pos=pos, fixed=[(0, 0)])
        assert_equal(tuple(pos[(0, 0)]), tuple(npos[(0, 0)]))
        # Sparse version (scipy based)
        pos = nx.circular_layout(self.bigG)
        npos = nx.fruchterman_reingold_layout(self.bigG, pos=pos, fixed=[(0, 0)])
        for axis in range(2):
            assert_almost_equal(pos[(0, 0)][axis], npos[(0, 0)][axis])

    def test_center_parameter(self):
        G = nx.path_graph(1)
        vpos = nx.random_layout(G, center=(1, 1))
        vpos = nx.circular_layout(G, center=(1, 1))
        assert_equal(tuple(vpos[0]), (1, 1))
        vpos = nx.spring_layout(G, center=(1, 1))
        assert_equal(tuple(vpos[0]), (1, 1))
        vpos = nx.fruchterman_reingold_layout(G, center=(1, 1))
        assert_equal(tuple(vpos[0]), (1, 1))
        vpos = nx.spectral_layout(G, center=(1, 1))
        assert_equal(tuple(vpos[0]), (1, 1))
        vpos = nx.shell_layout(G, center=(1, 1))
        assert_equal(tuple(vpos[0]), (1, 1))

    def test_center_wrong_dimensions(self):
        G = nx.path_graph(1)
        assert_raises(ValueError, nx.random_layout, G, center=(1, 1, 1))
        assert_raises(ValueError, nx.circular_layout, G, center=(1, 1, 1))
        assert_raises(ValueError, nx.spring_layout, G, center=(1, 1, 1))
        assert_raises(ValueError, nx.fruchterman_reingold_layout, G, center=(1, 1, 1))
        assert_raises(ValueError, nx.fruchterman_reingold_layout, G, dim=3, center=(1, 1))
        assert_raises(ValueError, nx.spectral_layout, G, center=(1, 1, 1))
        assert_raises(ValueError, nx.spectral_layout, G, dim=3, center=(1, 1))
        assert_raises(ValueError, nx.shell_layout, G, center=(1, 1, 1))

    def test_empty_graph(self):
        G = nx.empty_graph()
        vpos = nx.random_layout(G, center=(1, 1))
        assert_equal(vpos, {})
        vpos = nx.circular_layout(G, center=(1, 1))
        assert_equal(vpos, {})
        vpos = nx.spring_layout(G, center=(1, 1))
        assert_equal(vpos, {})
        vpos = nx.fruchterman_reingold_layout(G, center=(1, 1))
        assert_equal(vpos, {})
        vpos = nx.spectral_layout(G, center=(1, 1))
        assert_equal(vpos, {})
        vpos = nx.shell_layout(G, center=(1, 1))
        assert_equal(vpos, {})

    def test_kamada_kawai_costfn_1d(self):
        costfn = nx.drawing.layout._kamada_kawai_costfn

        pos = numpy.array([4.0, 7.0])
        invdist = 1 / numpy.array([[0.1, 2.0], [2.0, 0.3]])

        cost, grad = costfn(pos, numpy, invdist, meanweight=0, dim=1)

        assert_almost_equal(cost, ((3 / 2.0 - 1) ** 2))
        assert_almost_equal(grad[0], -0.5)
        assert_almost_equal(grad[1], 0.5)

    def test_kamada_kawai_costfn_2d(self):
        costfn = nx.drawing.layout._kamada_kawai_costfn

        pos = numpy.array([[1.3, -3.2],
                           [2.7, -0.3],
                           [5.1, 2.5]])
        invdist = 1 / numpy.array([[0.1, 2.1, 1.7],
                                   [2.1, 0.2, 0.6],
                                   [1.7, 0.6, 0.3]])
        meanwt = 0.3

        cost, grad = costfn(pos.ravel(), numpy, invdist,
                            meanweight=meanwt, dim=2)

        expected_cost = 0.5 * meanwt * numpy.sum(numpy.sum(pos, axis=0) ** 2)
        for i in range(pos.shape[0]):
            for j in range(i + 1, pos.shape[0]):
                expected_cost += (numpy.linalg.norm(pos[i] - pos[j]) * invdist[i][j] - 1.0) ** 2

        assert_almost_equal(cost, expected_cost)

        dx = 1e-4
        for nd in range(pos.shape[0]):
            for dm in range(pos.shape[1]):
                idx = nd * pos.shape[1] + dm
                pos0 = pos.flatten()

                pos0[idx] += dx
                cplus = costfn(pos0, numpy, invdist,
                               meanweight=meanwt, dim=pos.shape[1])[0]

                pos0[idx] -= 2 * dx
                cminus = costfn(pos0, numpy, invdist,
                                meanweight=meanwt, dim=pos.shape[1])[0]

                assert_almost_equal(grad[idx], (cplus - cminus) / (2 * dx),
                                    places=5)
