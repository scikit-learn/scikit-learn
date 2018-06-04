from nose import SkipTest
from nose.tools import assert_raises

import networkx as nx
from networkx.testing import assert_nodes_equal, assert_edges_equal, assert_graphs_equal


class TestConvertPandas(object):
    numpy = 1  # nosetests attribute, use nosetests -a 'not numpy' to skip test

    @classmethod
    def setupClass(cls):
        try:
            import pandas as pd
        except ImportError:
            raise SkipTest('Pandas not available.')

    def __init__(self):
        global pd
        import pandas as pd

        self.r = pd.np.random.RandomState(seed=5)
        ints = self.r.random_integers(1, 10, size=(3, 2))
        a = ['A', 'B', 'C']
        b = ['D', 'A', 'E']
        df = pd.DataFrame(ints, columns=['weight', 'cost'])
        df[0] = a  # Column label 0 (int)
        df['b'] = b  # Column label 'b' (str)
        self.df = df
        mdf = pd.DataFrame([[4, 16, 'A', 'D']],
                           columns=['weight', 'cost', 0, 'b'])
        self.mdf = df.append(mdf)

    def test_exceptions(self):
        G = pd.DataFrame(["a"])  # adj
        assert_raises(nx.NetworkXError, nx.to_networkx_graph, G)
        G = pd.DataFrame(["a", 0.0])  # elist
        assert_raises(nx.NetworkXError, nx.to_networkx_graph, G)
        df = pd.DataFrame([[1, 1], [1, 0]], dtype=int, index=[1, 2], columns=["a", "b"])
        assert_raises(nx.NetworkXError, nx.from_pandas_adjacency, df)

    def test_from_edgelist_all_attr(self):
        Gtrue = nx.Graph([('E', 'C', {'cost': 9, 'weight': 10}),
                          ('B', 'A', {'cost': 1, 'weight': 7}),
                          ('A', 'D', {'cost': 7, 'weight': 4})])
        G = nx.from_pandas_edgelist(self.df, 0, 'b', True)
        assert_graphs_equal(G, Gtrue)
        # MultiGraph
        MGtrue = nx.MultiGraph(Gtrue)
        MGtrue.add_edge('A', 'D', cost=16, weight=4)
        MG = nx.from_pandas_edgelist(self.mdf, 0, 'b', True, nx.MultiGraph())
        assert_graphs_equal(MG, MGtrue)

    def test_from_edgelist_multi_attr(self):
        Gtrue = nx.Graph([('E', 'C', {'cost': 9, 'weight': 10}),
                          ('B', 'A', {'cost': 1, 'weight': 7}),
                          ('A', 'D', {'cost': 7, 'weight': 4})])
        G = nx.from_pandas_edgelist(self.df, 0, 'b', ['weight', 'cost'])
        assert_graphs_equal(G, Gtrue)

    def test_from_edgelist_multidigraph_and_edge_attr(self):
        # example from issue #2374
        Gtrue = nx.MultiDiGraph([('X1', 'X4', {'Co': 'zA', 'Mi': 0, 'St': 'X1'}),
                                 ('X1', 'X4', {'Co': 'zB', 'Mi': 54, 'St': 'X2'}),
                                 ('X1', 'X4', {'Co': 'zB', 'Mi': 49, 'St': 'X3'}),
                                 ('X1', 'X4', {'Co': 'zB', 'Mi': 44, 'St': 'X4'}),
                                 ('Y1', 'Y3', {'Co': 'zC', 'Mi': 0, 'St': 'Y1'}),
                                 ('Y1', 'Y3', {'Co': 'zC', 'Mi': 34, 'St': 'Y2'}),
                                 ('Y1', 'Y3', {'Co': 'zC', 'Mi': 29, 'St': 'X2'}),
                                 ('Y1', 'Y3', {'Co': 'zC', 'Mi': 24, 'St': 'Y3'}),
                                 ('Z1', 'Z3', {'Co': 'zD', 'Mi': 0, 'St': 'Z1'}),
                                 ('Z1', 'Z3', {'Co': 'zD', 'Mi': 14, 'St': 'X3'}),
                                 ('Z1', 'Z3', {'Co': 'zE', 'Mi': 9, 'St': 'Z2'}),
                                 ('Z1', 'Z3', {'Co': 'zE', 'Mi': 4, 'St': 'Z3'})])
        df = pd.DataFrame.from_items([
            ('O', ['X1', 'X1', 'X1', 'X1', 'Y1', 'Y1', 'Y1', 'Y1', 'Z1', 'Z1', 'Z1', 'Z1']),
            ('D', ['X4', 'X4', 'X4', 'X4', 'Y3', 'Y3', 'Y3', 'Y3', 'Z3', 'Z3', 'Z3', 'Z3']),
            ('St', ['X1', 'X2', 'X3', 'X4', 'Y1', 'Y2', 'X2', 'Y3', 'Z1', 'X3', 'Z2', 'Z3']),
            ('Co', ['zA', 'zB', 'zB', 'zB', 'zC', 'zC', 'zC', 'zC', 'zD', 'zD', 'zE', 'zE']),
            ('Mi', [0,   54,   49,   44,    0,   34,   29,   24,    0,   14,    9,   4])])
        G1 = nx.from_pandas_edgelist(df, source='O', target='D',
                                     edge_attr=True,
                                     create_using=nx.MultiDiGraph())
        G2 = nx.from_pandas_edgelist(df, source='O', target='D',
                                     edge_attr=['St', 'Co', 'Mi'],
                                     create_using=nx.MultiDiGraph())
        assert_graphs_equal(G1, Gtrue)
        assert_graphs_equal(G2, Gtrue)

    def test_from_edgelist_one_attr(self):
        Gtrue = nx.Graph([('E', 'C', {'weight': 10}),
                          ('B', 'A', {'weight': 7}),
                          ('A', 'D', {'weight': 4})])
        G = nx.from_pandas_edgelist(self.df, 0, 'b', 'weight')
        assert_graphs_equal(G, Gtrue)

    def test_from_edgelist_no_attr(self):
        Gtrue = nx.Graph([('E', 'C', {}),
                          ('B', 'A', {}),
                          ('A', 'D', {})])
        G = nx.from_pandas_edgelist(self.df, 0, 'b',)
        assert_graphs_equal(G, Gtrue)

    def test_from_edgelist(self):
        # Pandas DataFrame
        g = nx.cycle_graph(10)
        G = nx.Graph()
        G.add_nodes_from(g)
        G.add_weighted_edges_from((u, v, u) for u, v in g.edges())
        edgelist = nx.to_edgelist(G)
        source = [s for s, t, d in edgelist]
        target = [t for s, t, d in edgelist]
        weight = [d['weight'] for s, t, d in edgelist]
        edges = pd.DataFrame({'source': source,
                              'target': target,
                              'weight': weight})
        GG = nx.from_pandas_edgelist(edges, edge_attr='weight')
        assert_nodes_equal(G.nodes(), GG.nodes())
        assert_edges_equal(G.edges(), GG.edges())
        GW = nx.to_networkx_graph(edges, create_using=nx.Graph())
        assert_nodes_equal(G.nodes(), GW.nodes())
        assert_edges_equal(G.edges(), GW.edges())

    def test_from_adjacency(self):
        nodelist = [1, 2]
        dftrue = pd.DataFrame([[1, 1], [1, 0]], dtype=int, index=nodelist, columns=nodelist)
        G = nx.Graph([(1, 1), (1, 2)])
        df = nx.to_pandas_adjacency(G, dtype=int)
        pd.testing.assert_frame_equal(df, dftrue)

    def test_roundtrip(self):
        # edgelist
        Gtrue = nx.Graph([(1, 1), (1, 2)])
        df = nx.to_pandas_edgelist(Gtrue)
        G = nx.from_pandas_edgelist(df)
        assert_graphs_equal(Gtrue, G)
        # adjacency
        Gtrue = nx.Graph(({1: {1: {'weight': 1}, 2: {'weight': 1}}, 2: {1: {'weight': 1}}}))
        df = nx.to_pandas_adjacency(Gtrue, dtype=int)
        G = nx.from_pandas_adjacency(df)
        assert_graphs_equal(Gtrue, G)
