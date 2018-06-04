from nose import SkipTest
from nose.tools import assert_raises, assert_true, raises

import networkx as nx
from networkx.testing import assert_graphs_equal
from networkx.generators.classic import barbell_graph, cycle_graph, path_graph


class TestConvertNumpy(object):
    @classmethod
    def setupClass(cls):
        global np, sp, sparse, np_assert_equal
        try:
            import numpy as np
            import scipy as sp
            import scipy.sparse as sparse
            np_assert_equal = np.testing.assert_equal
        except ImportError:
            raise SkipTest('SciPy sparse library not available.')

    def __init__(self):
        self.G1 = barbell_graph(10, 3)
        self.G2 = cycle_graph(10, create_using=nx.DiGraph())

        self.G3 = self.create_weighted(nx.Graph())
        self.G4 = self.create_weighted(nx.DiGraph())

    def test_exceptions(self):
        class G(object):
            format = None

        assert_raises(nx.NetworkXError, nx.to_networkx_graph, G)

    def create_weighted(self, G):
        g = cycle_graph(4)
        e = list(g.edges())
        source = [u for u, v in e]
        dest = [v for u, v in e]
        weight = [s + 10 for s in source]
        ex = zip(source, dest, weight)
        G.add_weighted_edges_from(ex)
        return G

    def assert_isomorphic(self, G1, G2):
        assert_true(nx.is_isomorphic(G1, G2))

    def identity_conversion(self, G, A, create_using):
        GG = nx.from_scipy_sparse_matrix(A, create_using=create_using)
        self.assert_isomorphic(G, GG)

        GW = nx.to_networkx_graph(A, create_using=create_using)
        self.assert_isomorphic(G, GW)

        GI = create_using.__class__(A)
        self.assert_isomorphic(G, GI)

        ACSR = A.tocsr()
        GI = create_using.__class__(ACSR)
        self.assert_isomorphic(G, GI)

        ACOO = A.tocoo()
        GI = create_using.__class__(ACOO)
        self.assert_isomorphic(G, GI)

        ACSC = A.tocsc()
        GI = create_using.__class__(ACSC)
        self.assert_isomorphic(G, GI)

        AD = A.todense()
        GI = create_using.__class__(AD)
        self.assert_isomorphic(G, GI)

        AA = A.toarray()
        GI = create_using.__class__(AA)
        self.assert_isomorphic(G, GI)

    def test_shape(self):
        "Conversion from non-square sparse array."
        A = sp.sparse.lil_matrix([[1, 2, 3], [4, 5, 6]])
        assert_raises(nx.NetworkXError, nx.from_scipy_sparse_matrix, A)

    def test_identity_graph_matrix(self):
        "Conversion from graph to sparse matrix to graph."
        A = nx.to_scipy_sparse_matrix(self.G1)
        self.identity_conversion(self.G1, A, nx.Graph())

    def test_identity_digraph_matrix(self):
        "Conversion from digraph to sparse matrix to digraph."
        A = nx.to_scipy_sparse_matrix(self.G2)
        self.identity_conversion(self.G2, A, nx.DiGraph())

    def test_identity_weighted_graph_matrix(self):
        """Conversion from weighted graph to sparse matrix to weighted graph."""
        A = nx.to_scipy_sparse_matrix(self.G3)
        self.identity_conversion(self.G3, A, nx.Graph())

    def test_identity_weighted_digraph_matrix(self):
        """Conversion from weighted digraph to sparse matrix to weighted digraph."""
        A = nx.to_scipy_sparse_matrix(self.G4)
        self.identity_conversion(self.G4, A, nx.DiGraph())

    def test_nodelist(self):
        """Conversion from graph to sparse matrix to graph with nodelist."""
        P4 = path_graph(4)
        P3 = path_graph(3)
        nodelist = list(P3.nodes())
        A = nx.to_scipy_sparse_matrix(P4, nodelist=nodelist)
        GA = nx.Graph(A)
        self.assert_isomorphic(GA, P3)

        # Make nodelist ambiguous by containing duplicates.
        nodelist += [nodelist[0]]
        assert_raises(nx.NetworkXError, nx.to_numpy_matrix, P3,
                      nodelist=nodelist)

    def test_weight_keyword(self):
        WP4 = nx.Graph()
        WP4.add_edges_from((n, n + 1, dict(weight=0.5, other=0.3))
                           for n in range(3))
        P4 = path_graph(4)
        A = nx.to_scipy_sparse_matrix(P4)
        np_assert_equal(A.todense(),
                        nx.to_scipy_sparse_matrix(WP4, weight=None).todense())
        np_assert_equal(0.5 * A.todense(),
                        nx.to_scipy_sparse_matrix(WP4).todense())
        np_assert_equal(0.3 * A.todense(),
                        nx.to_scipy_sparse_matrix(WP4, weight='other').todense())

    def test_format_keyword(self):
        WP4 = nx.Graph()
        WP4.add_edges_from((n, n + 1, dict(weight=0.5, other=0.3))
                           for n in range(3))
        P4 = path_graph(4)
        A = nx.to_scipy_sparse_matrix(P4, format='csr')
        np_assert_equal(A.todense(),
                        nx.to_scipy_sparse_matrix(WP4, weight=None).todense())

        A = nx.to_scipy_sparse_matrix(P4, format='csc')
        np_assert_equal(A.todense(),
                        nx.to_scipy_sparse_matrix(WP4, weight=None).todense())

        A = nx.to_scipy_sparse_matrix(P4, format='coo')
        np_assert_equal(A.todense(),
                        nx.to_scipy_sparse_matrix(WP4, weight=None).todense())

        A = nx.to_scipy_sparse_matrix(P4, format='bsr')
        np_assert_equal(A.todense(),
                        nx.to_scipy_sparse_matrix(WP4, weight=None).todense())

        A = nx.to_scipy_sparse_matrix(P4, format='lil')
        np_assert_equal(A.todense(),
                        nx.to_scipy_sparse_matrix(WP4, weight=None).todense())

        A = nx.to_scipy_sparse_matrix(P4, format='dia')
        np_assert_equal(A.todense(),
                        nx.to_scipy_sparse_matrix(WP4, weight=None).todense())

        A = nx.to_scipy_sparse_matrix(P4, format='dok')
        np_assert_equal(A.todense(),
                        nx.to_scipy_sparse_matrix(WP4, weight=None).todense())

    @raises(nx.NetworkXError)
    def test_format_keyword_raise(self):
        WP4 = nx.Graph()
        WP4.add_edges_from((n, n + 1, dict(weight=0.5, other=0.3))
                           for n in range(3))
        P4 = path_graph(4)
        nx.to_scipy_sparse_matrix(P4, format='any_other')

    @raises(nx.NetworkXError)
    def test_null_raise(self):
        nx.to_scipy_sparse_matrix(nx.Graph())

    def test_empty(self):
        G = nx.Graph()
        G.add_node(1)
        M = nx.to_scipy_sparse_matrix(G)
        np_assert_equal(M.todense(), np.matrix([[0]]))

    def test_ordering(self):
        G = nx.DiGraph()
        G.add_edge(1, 2)
        G.add_edge(2, 3)
        G.add_edge(3, 1)
        M = nx.to_scipy_sparse_matrix(G, nodelist=[3, 2, 1])
        np_assert_equal(M.todense(), np.matrix([[0, 0, 1], [1, 0, 0], [0, 1, 0]]))

    def test_selfloop_graph(self):
        G = nx.Graph([(1, 1)])
        M = nx.to_scipy_sparse_matrix(G)
        np_assert_equal(M.todense(), np.matrix([[1]]))

    def test_selfloop_digraph(self):
        G = nx.DiGraph([(1, 1)])
        M = nx.to_scipy_sparse_matrix(G)
        np_assert_equal(M.todense(), np.matrix([[1]]))

    def test_from_scipy_sparse_matrix_parallel_edges(self):
        """Tests that the :func:`networkx.from_scipy_sparse_matrix` function
        interprets integer weights as the number of parallel edges when
        creating a multigraph.

        """
        A = sparse.csr_matrix([[1, 1], [1, 2]])
        # First, with a simple graph, each integer entry in the adjacency
        # matrix is interpreted as the weight of a single edge in the graph.
        expected = nx.DiGraph()
        edges = [(0, 0), (0, 1), (1, 0)]
        expected.add_weighted_edges_from([(u, v, 1) for (u, v) in edges])
        expected.add_edge(1, 1, weight=2)
        actual = nx.from_scipy_sparse_matrix(A, parallel_edges=True,
                                             create_using=nx.DiGraph())
        assert_graphs_equal(actual, expected)
        actual = nx.from_scipy_sparse_matrix(A, parallel_edges=False,
                                             create_using=nx.DiGraph())
        assert_graphs_equal(actual, expected)
        # Now each integer entry in the adjacency matrix is interpreted as the
        # number of parallel edges in the graph if the appropriate keyword
        # argument is specified.
        edges = [(0, 0), (0, 1), (1, 0), (1, 1), (1, 1)]
        expected = nx.MultiDiGraph()
        expected.add_weighted_edges_from([(u, v, 1) for (u, v) in edges])
        actual = nx.from_scipy_sparse_matrix(A, parallel_edges=True,
                                             create_using=nx.MultiDiGraph())
        assert_graphs_equal(actual, expected)
        expected = nx.MultiDiGraph()
        expected.add_edges_from(set(edges), weight=1)
        # The sole self-loop (edge 0) on vertex 1 should have weight 2.
        expected[1][1][0]['weight'] = 2
        actual = nx.from_scipy_sparse_matrix(A, parallel_edges=False,
                                             create_using=nx.MultiDiGraph())
        assert_graphs_equal(actual, expected)

    def test_symmetric(self):
        """Tests that a symmetric matrix has edges added only once to an
        undirected multigraph when using
        :func:`networkx.from_scipy_sparse_matrix`.

        """
        A = sparse.csr_matrix([[0, 1], [1, 0]])
        G = nx.from_scipy_sparse_matrix(A, create_using=nx.MultiGraph())
        expected = nx.MultiGraph()
        expected.add_edge(0, 1, weight=1)
        assert_graphs_equal(G, expected)
