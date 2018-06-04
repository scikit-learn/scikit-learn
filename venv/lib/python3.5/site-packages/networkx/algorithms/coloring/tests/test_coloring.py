# -*- coding: utf-8 -*-
"""Greedy coloring test suite.

Run with nose: nosetests -v test_coloring.py
"""

__author__ = "\n".join(["Christian Olsson <chro@itu.dk>",
                        "Jan Aagaard Meier <jmei@itu.dk>",
                        "Henrik Haugbølle <hhau@itu.dk>",
                        "Jake VanderPlas <jakevdp@uw.edu>"])

import networkx as nx
from nose.tools import *

ALL_STRATEGIES = [
    'largest_first',
    'random_sequential',
    'smallest_last',
    'independent_set',
    'connected_sequential_bfs',
    'connected_sequential_dfs',
    'connected_sequential',
    'saturation_largest_first',
    'DSATUR',
]

# List of strategies where interchange=True results in an error
INTERCHANGE_INVALID = [
    'independent_set',
    'saturation_largest_first',
    'DSATUR'
]


class TestColoring:
    def test_basic_cases(self):
        def check_basic_case(graph_func, n_nodes, strategy, interchange):
            graph = graph_func()
            coloring = nx.coloring.greedy_color(graph,
                                                strategy=strategy,
                                                interchange=interchange)
            assert_true(verify_length(coloring, n_nodes))
            assert_true(verify_coloring(graph, coloring))

        for graph_func, n_nodes in BASIC_TEST_CASES.items():
            for interchange in [True, False]:
                for strategy in ALL_STRATEGIES:
                    if interchange and (strategy in INTERCHANGE_INVALID):
                        continue
                    yield (check_basic_case, graph_func,
                           n_nodes, strategy, interchange)

    def test_special_cases(self):
        def check_special_case(strategy, graph_func, interchange, colors):
            graph = graph_func()
            coloring = nx.coloring.greedy_color(graph,
                                                strategy=strategy,
                                                interchange=interchange)
            if not hasattr(colors, '__len__'):
                colors = [colors]
            assert_true(any(verify_length(coloring, n_colors)
                            for n_colors in colors))
            assert_true(verify_coloring(graph, coloring))

        for strategy, arglist in SPECIAL_TEST_CASES.items():
            for args in arglist:
                yield (check_special_case, strategy, args[0], args[1], args[2])

    def test_interchange_invalid(self):
        graph = one_node_graph()

        def check_raises(strategy):
            assert_raises(nx.NetworkXPointlessConcept,
                          nx.coloring.greedy_color,
                          graph, strategy=strategy, interchange=True)

        for strategy in INTERCHANGE_INVALID:
            yield check_raises, strategy

    def test_bad_inputs(self):
        graph = one_node_graph()
        assert_raises(nx.NetworkXError, nx.coloring.greedy_color,
                      graph, strategy='invalid strategy')

    def test_strategy_as_function(self):
        graph = lf_shc()
        colors_1 = nx.coloring.greedy_color(graph,
                                            'largest_first')
        colors_2 = nx.coloring.greedy_color(graph,
                                            nx.coloring.strategy_largest_first)
        assert_equal(colors_1, colors_2)


############################## Utility functions ##############################
def verify_coloring(graph, coloring):
    for node in graph.nodes():
        if node not in coloring:
            return False

        color = coloring[node]
        for neighbor in graph.neighbors(node):
            if coloring[neighbor] == color:
                return False

    return True


def verify_length(coloring, expected):
    coloring = dict_to_sets(coloring)
    return len(coloring) == expected


def dict_to_sets(colors):
    if len(colors) == 0:
        return []

    k = max(colors.values()) + 1
    sets = [set() for _ in range(k)]

    for (node, color) in colors.items():
        sets[color].add(node)

    return sets

############################## Graph Generation ##############################


def empty_graph():
    return nx.Graph()


def one_node_graph():
    graph = nx.Graph()
    graph.add_nodes_from([1])
    return graph


def two_node_graph():
    graph = nx.Graph()
    graph.add_nodes_from([1, 2])
    graph.add_edges_from([(1, 2)])
    return graph


def three_node_clique():
    graph = nx.Graph()
    graph.add_nodes_from([1, 2, 3])
    graph.add_edges_from([(1, 2), (1, 3), (2, 3)])
    return graph


def disconnected():
    graph = nx.Graph()
    graph.add_edges_from([
        (1, 2),
        (2, 3),
        (4, 5),
        (5, 6)
    ])
    return graph


def rs_shc():
    graph = nx.Graph()
    graph.add_nodes_from([1, 2, 3, 4])
    graph.add_edges_from([
        (1, 2),
        (2, 3),
        (3, 4)
    ])
    return graph


def slf_shc():
    graph = nx.Graph()
    graph.add_nodes_from([1, 2, 3, 4, 5, 6, 7])
    graph.add_edges_from([
        (1, 2),
        (1, 5),
        (1, 6),
        (2, 3),
        (2, 7),
        (3, 4),
        (3, 7),
        (4, 5),
        (4, 6),
        (5, 6)
    ])
    return graph


def slf_hc():
    graph = nx.Graph()
    graph.add_nodes_from([1, 2, 3, 4, 5, 6, 7, 8])
    graph.add_edges_from([
        (1, 2),
        (1, 3),
        (1, 4),
        (1, 5),
        (2, 3),
        (2, 4),
        (2, 6),
        (5, 7),
        (5, 8),
        (6, 7),
        (6, 8),
        (7, 8)
    ])
    return graph


def lf_shc():
    graph = nx.Graph()
    graph.add_nodes_from([1, 2, 3, 4, 5, 6])
    graph.add_edges_from([
        (6, 1),
        (1, 4),
        (4, 3),
        (3, 2),
        (2, 5)
    ])
    return graph


def lf_hc():
    graph = nx.Graph()
    graph.add_nodes_from([1, 2, 3, 4, 5, 6, 7])
    graph.add_edges_from([
        (1, 7),
        (1, 6),
        (1, 3),
        (1, 4),
        (7, 2),
        (2, 6),
        (2, 3),
        (2, 5),
        (5, 3),
        (5, 4),
        (4, 3)
    ])
    return graph


def sl_shc():
    graph = nx.Graph()
    graph.add_nodes_from([1, 2, 3, 4, 5, 6])
    graph.add_edges_from([
        (1, 2),
        (1, 3),
        (2, 3),
        (1, 4),
        (2, 5),
        (3, 6),
        (4, 5),
        (4, 6),
        (5, 6)
    ])
    return graph


def sl_hc():
    graph = nx.Graph()
    graph.add_nodes_from([1, 2, 3, 4, 5, 6, 7, 8])
    graph.add_edges_from([
        (1, 2),
        (1, 3),
        (1, 5),
        (1, 7),
        (2, 3),
        (2, 4),
        (2, 8),
        (8, 4),
        (8, 6),
        (8, 7),
        (7, 5),
        (7, 6),
        (3, 4),
        (4, 6),
        (6, 5),
        (5, 3)
    ])
    return graph


def gis_shc():
    graph = nx.Graph()
    graph.add_nodes_from([1, 2, 3, 4])
    graph.add_edges_from([
        (1, 2),
        (2, 3),
        (3, 4)
    ])
    return graph


def gis_hc():
    graph = nx.Graph()
    graph.add_nodes_from([1, 2, 3, 4, 5, 6])
    graph.add_edges_from([
        (1, 5),
        (2, 5),
        (3, 6),
        (4, 6),
        (5, 6)
    ])
    return graph


def cs_shc():
    graph = nx.Graph()
    graph.add_nodes_from([1, 2, 3, 4, 5])
    graph.add_edges_from([
        (1, 2),
        (1, 5),
        (2, 3),
        (2, 4),
        (2, 5),
        (3, 4),
        (4, 5)
    ])
    return graph


def rsi_shc():
    graph = nx.Graph()
    graph.add_nodes_from([1, 2, 3, 4, 5, 6])
    graph.add_edges_from([
        (1, 2),
        (1, 5),
        (1, 6),
        (2, 3),
        (3, 4),
        (4, 5),
        (4, 6),
        (5, 6)
    ])
    return graph


def lfi_shc():
    graph = nx.Graph()
    graph.add_nodes_from([1, 2, 3, 4, 5, 6, 7])
    graph.add_edges_from([
        (1, 2),
        (1, 5),
        (1, 6),
        (2, 3),
        (2, 7),
        (3, 4),
        (3, 7),
        (4, 5),
        (4, 6),
        (5, 6)
    ])
    return graph


def lfi_hc():
    graph = nx.Graph()
    graph.add_nodes_from([1, 2, 3, 4, 5, 6, 7, 8, 9])
    graph.add_edges_from([
        (1, 2),
        (1, 5),
        (1, 6),
        (1, 7),
        (2, 3),
        (2, 8),
        (2, 9),
        (3, 4),
        (3, 8),
        (3, 9),
        (4, 5),
        (4, 6),
        (4, 7),
        (5, 6)
    ])
    return graph


def sli_shc():
    graph = nx.Graph()
    graph.add_nodes_from([1, 2, 3, 4, 5, 6, 7])
    graph.add_edges_from([
        (1, 2),
        (1, 3),
        (1, 5),
        (1, 7),
        (2, 3),
        (2, 6),
        (3, 4),
        (4, 5),
        (4, 6),
        (5, 7),
        (6, 7)
    ])
    return graph


def sli_hc():
    graph = nx.Graph()
    graph.add_nodes_from([1, 2, 3, 4, 5, 6, 7, 8, 9])
    graph.add_edges_from([
        (1, 2),
        (1, 3),
        (1, 4),
        (1, 5),
        (2, 3),
        (2, 7),
        (2, 8),
        (2, 9),
        (3, 6),
        (3, 7),
        (3, 9),
        (4, 5),
        (4, 6),
        (4, 8),
        (4, 9),
        (5, 6),
        (5, 7),
        (5, 8),
        (6, 7),
        (6, 9),
        (7, 8),
        (8, 9)
    ])
    return graph


#---------------------------------------------------------------------------
# Basic tests for all strategies
# For each basic graph function, specify the number of expected colors.
BASIC_TEST_CASES = {empty_graph: 0,
                    one_node_graph: 1,
                    two_node_graph: 2,
                    disconnected: 2,
                    three_node_clique: 3}


#---------------------------------------------------------------------------
# Special test cases. Each strategy has a list of tuples of the form
# (graph function, interchange, valid # of colors)
SPECIAL_TEST_CASES = {
    'random_sequential': [
        (rs_shc, False, (2, 3)),
        (rs_shc, True, 2),
        (rsi_shc, True, (3, 4))],
    'saturation_largest_first': [
        (slf_shc, False, (3, 4)),
        (slf_hc, False, 4)],
    'largest_first': [
        (lf_shc, False, (2, 3)),
        (lf_hc, False, 4),
        (lf_shc, True, 2),
        (lf_hc, True, 3),
        (lfi_shc, True, (3, 4)),
        (lfi_hc, True, 4)],
    'smallest_last': [
        (sl_shc, False, (3, 4)),
        (sl_hc, False, 5),
        (sl_shc, True, 3),
        (sl_hc, True, 4),
        (sli_shc, True, (3, 4)),
        (sli_hc, True, 5)],
    'independent_set': [
        (gis_shc, False, (2, 3)),
        (gis_hc, False, 3)],
    'connected_sequential': [
        (cs_shc, False, (3, 4)),
        (cs_shc, True, 3)],
    'connected_sequential_dfs': [
        (cs_shc, False, (3, 4))],
}
