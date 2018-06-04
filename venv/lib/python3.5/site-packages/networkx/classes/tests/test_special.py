#!/usr/bin/env python
from nose.tools import *
import networkx as nx
from test_graph import TestGraph
from test_digraph import TestDiGraph
from test_multigraph import TestMultiGraph
from test_multidigraph import TestMultiDiGraph
try:  # python 2.7+
    from collections import OrderedDict
except ImportError:  # python 2.6
    try:
        from ordereddict import OrderedDict
    except ImportError:
        from nose import SkipTest
        raise SkipTest('ordereddict not available')


class SpecialGraphTester(TestGraph):
    def setUp(self):
        TestGraph.setUp(self)
        self.Graph = nx.Graph


class OrderedGraphTester(TestGraph):
    def setUp(self):
        TestGraph.setUp(self)

        class MyGraph(nx.Graph):
            node_dict_factory = OrderedDict
            adjlist_outer_dict_factory = OrderedDict
            adjlist_inner_dict_factory = OrderedDict
            edge_attr_dict_factory = OrderedDict
        self.Graph = MyGraph


class ThinGraphTester(TestGraph):
    def setUp(self):
        all_edge_dict = {'weight': 1}

        class MyGraph(nx.Graph):
            def edge_attr_dict_factory(): return all_edge_dict
        self.Graph = MyGraph
        # build dict-of-dict-of-dict K3
        ed1, ed2, ed3 = (all_edge_dict, all_edge_dict, all_edge_dict)
        self.k3adj = {0: {1: ed1, 2: ed2},
                      1: {0: ed1, 2: ed3},
                      2: {0: ed2, 1: ed3}}
        self.k3edges = [(0, 1), (0, 2), (1, 2)]
        self.k3nodes = [0, 1, 2]
        self.K3 = self.Graph()
        self.K3._adj = self.k3adj
        self.K3._node = {}
        self.K3._node[0] = {}
        self.K3._node[1] = {}
        self.K3._node[2] = {}


class SpecialDiGraphTester(TestDiGraph):
    def setUp(self):
        TestDiGraph.setUp(self)
        self.Graph = nx.DiGraph


class OrderedDiGraphTester(TestDiGraph):
    def setUp(self):
        TestGraph.setUp(self)

        class MyGraph(nx.DiGraph):
            node_dict_factory = OrderedDict
            adjlist_outer_dict_factory = OrderedDict
            adjlist_inner_dict_factory = OrderedDict
            edge_attr_dict_factory = OrderedDict
        self.Graph = MyGraph


class ThinDiGraphTester(TestDiGraph):
    def setUp(self):
        all_edge_dict = {'weight': 1}

        class MyGraph(nx.DiGraph):
            def edge_attr_dict_factory(): return all_edge_dict
        self.Graph = MyGraph
        # build dict-of-dict-of-dict K3
        ed1, ed2, ed3 = (all_edge_dict, all_edge_dict, all_edge_dict)
        self.k3adj = {0: {1: ed1, 2: ed2},
                      1: {0: ed1, 2: ed3},
                      2: {0: ed2, 1: ed3}}
        self.k3edges = [(0, 1), (0, 2), (1, 2)]
        self.k3nodes = [0, 1, 2]
        self.K3 = self.Graph()
        self.K3.adj = self.k3adj
        self.K3._node = {}
        self.K3._node[0] = {}
        self.K3._node[1] = {}
        self.K3._node[2] = {}


class SpecialMultiGraphTester(TestMultiGraph):
    def setUp(self):
        TestMultiGraph.setUp(self)
        self.Graph = nx.MultiGraph


class OrderedMultiGraphTester(TestMultiGraph):
    def setUp(self):
        TestMultiGraph.setUp(self)

        class MyGraph(nx.MultiGraph):
            node_dict_factory = OrderedDict
            adjlist_outer_dict_factory = OrderedDict
            adjlist_inner_dict_factory = OrderedDict
            edge_key_dict_factory = OrderedDict
            edge_attr_dict_factory = OrderedDict
        self.Graph = MyGraph


class SpecialMultiDiGraphTester(TestMultiDiGraph):
    def setUp(self):
        TestMultiDiGraph.setUp(self)
        self.Graph = nx.MultiDiGraph


class OrderedMultiDiGraphTester(TestMultiDiGraph):
    def setUp(self):
        TestMultiDiGraph.setUp(self)

        class MyGraph(nx.MultiDiGraph):
            node_dict_factory = OrderedDict
            adjlist_outer_dict_factory = OrderedDict
            adjlist_inner_dict_factory = OrderedDict
            edge_key_dict_factory = OrderedDict
            edge_attr_dict_factory = OrderedDict
        self.Graph = MyGraph
