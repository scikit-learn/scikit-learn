#!/usr/bin/env python
from nose.tools import *
from nose import SkipTest
import networkx as nx


class TestFlowClosenessCentrality(object):
    numpy = 1  # nosetests attribute, use nosetests -a 'not numpy' to skip test

    @classmethod
    def setupClass(cls):
        global np
        try:
            import numpy as np
            import scipy
        except ImportError:
            raise SkipTest('NumPy not available.')

    def test_K4(self):
        """Closeness centrality: K4"""
        G = nx.complete_graph(4)
        b = nx.current_flow_closeness_centrality(G)
        b_answer = {0: 2.0 / 3, 1: 2.0 / 3, 2: 2.0 / 3, 3: 2.0 / 3}
        for n in sorted(G):
            assert_almost_equal(b[n], b_answer[n])

    def test_P4(self):
        """Closeness centrality: P4"""
        G = nx.path_graph(4)
        b = nx.current_flow_closeness_centrality(G)
        b_answer = {0: 1.0 / 6, 1: 1.0 / 4, 2: 1.0 / 4, 3: 1.0 / 6}
        for n in sorted(G):
            assert_almost_equal(b[n], b_answer[n])

    def test_star(self):
        """Closeness centrality: star """
        G = nx.Graph()
        nx.add_star(G, ['a', 'b', 'c', 'd'])
        b = nx.current_flow_closeness_centrality(G)
        b_answer = {'a': 1.0 / 3, 'b': 0.6 / 3, 'c': 0.6 / 3, 'd': 0.6 / 3}
        for n in sorted(G):
            assert_almost_equal(b[n], b_answer[n])


class TestWeightedFlowClosenessCentrality(object):
    pass
