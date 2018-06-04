from nose import SkipTest

import networkx as nx
from networkx.generators.degree_seq import havel_hakimi_graph


class TestGraphMatrix(object):
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
        self.OI = numpy.array([[-1, -1, -1, 0],
                               [1, 0, 0, -1],
                               [0, 1, 0, 1],
                               [0, 0, 1, 0],
                               [0, 0, 0, 0]])
        self.A = numpy.array([[0, 1, 1, 1, 0],
                              [1, 0, 1, 0, 0],
                              [1, 1, 0, 0, 0],
                              [1, 0, 0, 0, 0],
                              [0, 0, 0, 0, 0]])
        self.WG = havel_hakimi_graph(deg)
        self.WG.add_edges_from((u, v, {'weight': 0.5, 'other': 0.3})
                               for (u, v) in self.G.edges())
        self.WA = numpy.array([[0, 0.5, 0.5, 0.5, 0],
                               [0.5, 0, 0.5, 0, 0],
                               [0.5, 0.5, 0, 0, 0],
                               [0.5, 0, 0, 0, 0],
                               [0, 0, 0, 0, 0]])
        self.MG = nx.MultiGraph(self.G)
        self.MG2 = self.MG.copy()
        self.MG2.add_edge(0, 1)
        self.MG2A = numpy.array([[0, 2, 1, 1, 0],
                                 [2, 0, 1, 0, 0],
                                 [1, 1, 0, 0, 0],
                                 [1, 0, 0, 0, 0],
                                 [0, 0, 0, 0, 0]])
        self.MGOI = numpy.array([[-1, -1, -1, -1, 0],
                                 [1, 1, 0, 0, -1],
                                 [0, 0, 1, 0, 1],
                                 [0, 0, 0, 1, 0],
                                 [0, 0, 0, 0, 0]])
        self.no_edges_G = nx.Graph([(1, 2), (3, 2, {'weight': 8})])
        self.no_edges_A = numpy.array([[0, 0], [0, 0]])

    def test_incidence_matrix(self):
        "Conversion to incidence matrix"
        I = nx.incidence_matrix(self.G,
                                nodelist=sorted(self.G),
                                edgelist=sorted(self.G.edges()),
                                oriented=True).todense().astype(int)
        assert_equal(I, self.OI)
        I = nx.incidence_matrix(self.G,
                                nodelist=sorted(self.G),
                                edgelist=sorted(self.G.edges()),
                                oriented=False).todense().astype(int)
        assert_equal(I, numpy.abs(self.OI))

        I = nx.incidence_matrix(self.MG,
                                nodelist=sorted(self.MG),
                                edgelist=sorted(self.MG.edges()),
                                oriented=True).todense().astype(int)
        assert_equal(I, self.OI)
        I = nx.incidence_matrix(self.MG,
                                nodelist=sorted(self.MG),
                                edgelist=sorted(self.MG.edges()),
                                oriented=False).todense().astype(int)
        assert_equal(I, numpy.abs(self.OI))

        I = nx.incidence_matrix(self.MG2,
                                nodelist=sorted(self.MG2),
                                edgelist=sorted(self.MG2.edges()),
                                oriented=True).todense().astype(int)
        assert_equal(I, self.MGOI)
        I = nx.incidence_matrix(self.MG2,
                                nodelist=sorted(self.MG),
                                edgelist=sorted(self.MG2.edges()),
                                oriented=False).todense().astype(int)
        assert_equal(I, numpy.abs(self.MGOI))

    def test_weighted_incidence_matrix(self):
        I = nx.incidence_matrix(self.WG,
                                nodelist=sorted(self.WG),
                                edgelist=sorted(self.WG.edges()),
                                oriented=True).todense().astype(int)
        assert_equal(I, self.OI)
        I = nx.incidence_matrix(self.WG,
                                nodelist=sorted(self.WG),
                                edgelist=sorted(self.WG.edges()),
                                oriented=False).todense().astype(int)
        assert_equal(I, numpy.abs(self.OI))

        # assert_equal(nx.incidence_matrix(self.WG,oriented=True,
        #                                  weight='weight').todense(),0.5*self.OI)
        # assert_equal(nx.incidence_matrix(self.WG,weight='weight').todense(),
        #              numpy.abs(0.5*self.OI))
        # assert_equal(nx.incidence_matrix(self.WG,oriented=True,weight='other').todense(),
        #              0.3*self.OI)

        I = nx.incidence_matrix(self.WG,
                                nodelist=sorted(self.WG),
                                edgelist=sorted(self.WG.edges()),
                                oriented=True,
                                weight='weight').todense()
        assert_equal(I, 0.5 * self.OI)
        I = nx.incidence_matrix(self.WG,
                                nodelist=sorted(self.WG),
                                edgelist=sorted(self.WG.edges()),
                                oriented=False,
                                weight='weight').todense()
        assert_equal(I, numpy.abs(0.5 * self.OI))
        I = nx.incidence_matrix(self.WG,
                                nodelist=sorted(self.WG),
                                edgelist=sorted(self.WG.edges()),
                                oriented=True,
                                weight='other').todense()
        assert_equal(I, 0.3 * self.OI)

        # WMG=nx.MultiGraph(self.WG)
        # WMG.add_edge(0,1,weight=0.5,other=0.3)
        # assert_equal(nx.incidence_matrix(WMG,weight='weight').todense(),
        #              numpy.abs(0.5*self.MGOI))
        # assert_equal(nx.incidence_matrix(WMG,weight='weight',oriented=True).todense(),
        #              0.5*self.MGOI)
        # assert_equal(nx.incidence_matrix(WMG,weight='other',oriented=True).todense(),
        #              0.3*self.MGOI)

        WMG = nx.MultiGraph(self.WG)
        WMG.add_edge(0, 1, weight=0.5, other=0.3)
        I = nx.incidence_matrix(WMG,
                                nodelist=sorted(WMG),
                                edgelist=sorted(WMG.edges(keys=True)),
                                oriented=True,
                                weight='weight').todense()
        assert_equal(I, 0.5 * self.MGOI)
        I = nx.incidence_matrix(WMG,
                                nodelist=sorted(WMG),
                                edgelist=sorted(WMG.edges(keys=True)),
                                oriented=False,
                                weight='weight').todense()
        assert_equal(I, numpy.abs(0.5 * self.MGOI))
        I = nx.incidence_matrix(WMG,
                                nodelist=sorted(WMG),
                                edgelist=sorted(WMG.edges(keys=True)),
                                oriented=True,
                                weight='other').todense()
        assert_equal(I, 0.3 * self.MGOI)

    def test_adjacency_matrix(self):
        "Conversion to adjacency matrix"
        assert_equal(nx.adj_matrix(self.G).todense(), self.A)
        assert_equal(nx.adj_matrix(self.MG).todense(), self.A)
        assert_equal(nx.adj_matrix(self.MG2).todense(), self.MG2A)
        assert_equal(nx.adj_matrix(self.G, nodelist=[0, 1]).todense(), self.A[:2, :2])
        assert_equal(nx.adj_matrix(self.WG).todense(), self.WA)
        assert_equal(nx.adj_matrix(self.WG, weight=None).todense(), self.A)
        assert_equal(nx.adj_matrix(self.MG2, weight=None).todense(), self.MG2A)
        assert_equal(nx.adj_matrix(self.WG, weight='other').todense(), 0.6 * self.WA)
        assert_equal(nx.adj_matrix(self.no_edges_G, nodelist=[1, 3]).todense(), self.no_edges_A)
