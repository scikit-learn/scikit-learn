import networkx as nx
from nose.tools import *


class TestImmediateDominators(object):

    def test_exceptions(self):
        G = nx.Graph()
        G.add_node(0)
        assert_raises(nx.NetworkXNotImplemented, nx.immediate_dominators, G, 0)
        G = nx.MultiGraph(G)
        assert_raises(nx.NetworkXNotImplemented, nx.immediate_dominators, G, 0)
        G = nx.DiGraph([[0, 0]])
        assert_raises(nx.NetworkXError, nx.immediate_dominators, G, 1)

    def test_singleton(self):
        G = nx.DiGraph()
        G.add_node(0)
        assert_equal(nx.immediate_dominators(G, 0), {0: 0})
        G.add_edge(0, 0)
        assert_equal(nx.immediate_dominators(G, 0), {0: 0})

    def test_path(self):
        n = 5
        G = nx.path_graph(n, create_using=nx.DiGraph())
        assert_equal(nx.immediate_dominators(G, 0),
                     {i: max(i - 1, 0) for i in range(n)})

    def test_cycle(self):
        n = 5
        G = nx.cycle_graph(n, create_using=nx.DiGraph())
        assert_equal(nx.immediate_dominators(G, 0),
                     {i: max(i - 1, 0) for i in range(n)})

    def test_unreachable(self):
        n = 5
        assert_greater(n, 1)
        G = nx.path_graph(n, create_using=nx.DiGraph())
        assert_equal(nx.immediate_dominators(G, n // 2),
                     {i: max(i - 1, n // 2) for i in range(n // 2, n)})

    def test_irreducible1(self):
        # Graph taken from Figure 2 of
        # K. D. Cooper, T. J. Harvey, and K. Kennedy.
        # A simple, fast dominance algorithm.
        # Software Practice & Experience, 4:110, 2001.
        edges = [(1, 2), (2, 1), (3, 2), (4, 1), (5, 3), (5, 4)]
        G = nx.DiGraph(edges)
        assert_equal(nx.immediate_dominators(G, 5),
                     {i: 5 for i in range(1, 6)})

    def test_irreducible2(self):
        # Graph taken from Figure 4 of
        # K. D. Cooper, T. J. Harvey, and K. Kennedy.
        # A simple, fast dominance algorithm.
        # Software Practice & Experience, 4:110, 2001.
        edges = [(1, 2), (2, 1), (2, 3), (3, 2), (4, 2), (4, 3), (5, 1),
                 (6, 4), (6, 5)]
        G = nx.DiGraph(edges)
        assert_equal(nx.immediate_dominators(G, 6),
                     {i: 6 for i in range(1, 7)})

    def test_domrel_png(self):
        # Graph taken from https://commons.wikipedia.org/wiki/File:Domrel.png
        edges = [(1, 2), (2, 3), (2, 4), (2, 6), (3, 5), (4, 5), (5, 2)]
        G = nx.DiGraph(edges)
        assert_equal(nx.immediate_dominators(G, 1),
                     {1: 1, 2: 1, 3: 2, 4: 2, 5: 2, 6: 2})
        # Test postdominance.
        with nx.utils.reversed(G):
            assert_equal(nx.immediate_dominators(G, 6),
                         {1: 2, 2: 6, 3: 5, 4: 5, 5: 2, 6: 6})

    def test_boost_example(self):
        # Graph taken from Figure 1 of
        # http://www.boost.org/doc/libs/1_56_0/libs/graph/doc/lengauer_tarjan_dominator.htm
        edges = [(0, 1), (1, 2), (1, 3), (2, 7), (3, 4), (4, 5), (4, 6),
                 (5, 7), (6, 4)]
        G = nx.DiGraph(edges)
        assert_equal(nx.immediate_dominators(G, 0),
                     {0: 0, 1: 0, 2: 1, 3: 1, 4: 3, 5: 4, 6: 4, 7: 1})
        # Test postdominance.
        with nx.utils.reversed(G):
            assert_equal(nx.immediate_dominators(G, 7),
                         {0: 1, 1: 7, 2: 7, 3: 4, 4: 5, 5: 7, 6: 4, 7: 7})


class TestDominanceFrontiers(object):

    def test_exceptions(self):
        G = nx.Graph()
        G.add_node(0)
        assert_raises(nx.NetworkXNotImplemented, nx.dominance_frontiers, G, 0)
        G = nx.MultiGraph(G)
        assert_raises(nx.NetworkXNotImplemented, nx.dominance_frontiers, G, 0)
        G = nx.DiGraph([[0, 0]])
        assert_raises(nx.NetworkXError, nx.dominance_frontiers, G, 1)

    def test_singleton(self):
        G = nx.DiGraph()
        G.add_node(0)
        assert_equal(nx.dominance_frontiers(G, 0), {0: set()})
        G.add_edge(0, 0)
        assert_equal(nx.dominance_frontiers(G, 0), {0: set()})

    def test_path(self):
        n = 5
        G = nx.path_graph(n, create_using=nx.DiGraph())
        assert_equal(nx.dominance_frontiers(G, 0),
                     {i: set() for i in range(n)})

    def test_cycle(self):
        n = 5
        G = nx.cycle_graph(n, create_using=nx.DiGraph())
        assert_equal(nx.dominance_frontiers(G, 0),
                     {i: set() for i in range(n)})

    def test_unreachable(self):
        n = 5
        assert_greater(n, 1)
        G = nx.path_graph(n, create_using=nx.DiGraph())
        assert_equal(nx.dominance_frontiers(G, n // 2),
                     {i: set() for i in range(n // 2, n)})

    def test_irreducible1(self):
        # Graph taken from Figure 2 of
        # K. D. Cooper, T. J. Harvey, and K. Kennedy.
        # A simple, fast dominance algorithm.
        # Software Practice & Experience, 4:110, 2001.
        edges = [(1, 2), (2, 1), (3, 2), (4, 1), (5, 3), (5, 4)]
        G = nx.DiGraph(edges)
        assert_equal({u: df
                      for u, df in nx.dominance_frontiers(G, 5).items()},
                     {1: set([2]), 2: set([1]), 3: set([2]),
                      4: set([1]), 5: set()})

    def test_irreducible2(self):
        # Graph taken from Figure 4 of
        # K. D. Cooper, T. J. Harvey, and K. Kennedy.
        # A simple, fast dominance algorithm.
        # Software Practice & Experience, 4:110, 2001.
        edges = [(1, 2), (2, 1), (2, 3), (3, 2), (4, 2), (4, 3), (5, 1),
                 (6, 4), (6, 5)]
        G = nx.DiGraph(edges)
        assert_equal(nx.dominance_frontiers(G, 6),
                     {1: set([2]), 2: set([1, 3]), 3: set([2]), 4: set([2, 3]), 5: set([1]), 6: set([])})

    def test_domrel_png(self):
        # Graph taken from https://commons.wikipedia.org/wiki/File:Domrel.png
        edges = [(1, 2), (2, 3), (2, 4), (2, 6), (3, 5), (4, 5), (5, 2)]
        G = nx.DiGraph(edges)
        assert_equal(nx.dominance_frontiers(G, 1),
                     {1: set([]), 2: set([2]), 3: set([5]), 4: set([5]),
                      5: set([2]), 6: set()})
        # Test postdominance.
        with nx.utils.reversed(G):
            assert_equal(nx.dominance_frontiers(G, 6),
                         {1: set(), 2: set([2]), 3: set([2]), 4: set([2]),
                          5: set([2]), 6: set()})

    def test_boost_example(self):
        # Graph taken from Figure 1 of
        # http://www.boost.org/doc/libs/1_56_0/libs/graph/doc/lengauer_tarjan_dominator.htm
        edges = [(0, 1), (1, 2), (1, 3), (2, 7), (3, 4), (4, 5), (4, 6),
                 (5, 7), (6, 4)]
        G = nx.DiGraph(edges)
        assert_equal(nx.dominance_frontiers(G, 0),
                     {0: set(), 1: set(), 2: set([7]), 3: set([7]),
                      4: set([4, 7]), 5: set([7]), 6: set([4]), 7: set()})
        # Test postdominance.
        with nx.utils.reversed(G):
            assert_equal(nx.dominance_frontiers(G, 7),
                         {0: set(), 1: set(), 2: set([1]), 3: set([1]),
                          4: set([1, 4]), 5: set([1]), 6: set([4]), 7: set()})

    def test_discard_issue(self):
        # https://github.com/networkx/networkx/issues/2071
        g = nx.DiGraph()
        g.add_edges_from([
            ('b0', 'b1'),
            ('b1', 'b2'),
            ('b2', 'b3'),
            ('b3', 'b1'),
            ('b1', 'b5'),
            ('b5', 'b6'),
            ('b5', 'b8'),
            ('b6', 'b7'),
            ('b8', 'b7'),
            ('b7', 'b3'),
            ('b3', 'b4')
        ]
        )
        df = nx.dominance_frontiers(g, 'b0')
        assert_equal(df, {'b4': set(), 'b5': set(['b3']), 'b6': set(['b7']),
                          'b7': set(['b3']),
                          'b0': set(), 'b1': set(['b1']), 'b2': set(['b3']),
                          'b3': set(['b1']), 'b8': set(['b7'])})

    def test_loop(self):
        g = nx.DiGraph()
        g.add_edges_from([('a', 'b'), ('b', 'c'), ('b', 'a')])
        df = nx.dominance_frontiers(g, 'a')
        assert_equal(df, {'a': set(), 'b': set(), 'c': set()})

    def test_missing_immediate_doms(self):
        # see https://github.com/networkx/networkx/issues/2070
        g = nx.DiGraph()
        edges = [
            ('entry_1', 'b1'),
            ('b1', 'b2'),
            ('b2', 'b3'),
            ('b3', 'exit'),
            ('entry_2', 'b3')
        ]

        # entry_1
        #   |
        #   b1
        #   |
        #   b2  entry_2
        #    |  /
        #    b3
        #    |
        #   exit

        g.add_edges_from(edges)
        # formerly raised KeyError on entry_2 when parsing b3
        # because entry_2 does not have immediate doms (no path)
        nx.dominance_frontiers(g, 'entry_1')

    def test_loops_larger(self):
        # from
        # http://ecee.colorado.edu/~waite/Darmstadt/motion.html
        g = nx.DiGraph()
        edges = [
            ('entry', 'exit'),
            ('entry', '1'),
            ('1', '2'),
            ('2', '3'),
            ('3', '4'),
            ('4', '5'),
            ('5', '6'),
            ('6', 'exit'),
            ('6', '2'),
            ('5', '3'),
            ('4', '4')
        ]

        g.add_edges_from(edges)
        df = nx.dominance_frontiers(g, 'entry')
        answer = {'entry': set(),
                  '1': set(['exit']),
                  '2': set(['exit', '2']),
                  '3': set(['exit', '3', '2']),
                  '4': set(['exit', '4', '3', '2']),
                  '5': set(['exit', '3', '2']),
                  '6': set(['exit', '2']),
                  'exit': set()}
        for n in df:
            assert_equal(set(df[n]), set(answer[n]))
