from __future__ import division, print_function, absolute_import

import numpy as np
from numpy.testing import assert_array_almost_equal, assert_array_equal
from pytest import raises as assert_raises
from scipy.sparse.csgraph import (shortest_path, dijkstra, johnson,
    bellman_ford, construct_dist_matrix, NegativeCycleError)


directed_G = np.array([[0, 3, 3, 0, 0],
                       [0, 0, 0, 2, 4],
                       [0, 0, 0, 0, 0],
                       [1, 0, 0, 0, 0],
                       [2, 0, 0, 2, 0]], dtype=float)

undirected_G = np.array([[0, 3, 3, 1, 2],
                         [3, 0, 0, 2, 4],
                         [3, 0, 0, 0, 0],
                         [1, 2, 0, 0, 2],
                         [2, 4, 0, 2, 0]], dtype=float)

unweighted_G = (directed_G > 0).astype(float)

directed_SP = [[0, 3, 3, 5, 7],
               [3, 0, 6, 2, 4],
               [np.inf, np.inf, 0, np.inf, np.inf],
               [1, 4, 4, 0, 8],
               [2, 5, 5, 2, 0]]

directed_pred = np.array([[-9999, 0, 0, 1, 1],
                          [3, -9999, 0, 1, 1],
                          [-9999, -9999, -9999, -9999, -9999],
                          [3, 0, 0, -9999, 1],
                          [4, 0, 0, 4, -9999]], dtype=float)

undirected_SP = np.array([[0, 3, 3, 1, 2],
                          [3, 0, 6, 2, 4],
                          [3, 6, 0, 4, 5],
                          [1, 2, 4, 0, 2],
                          [2, 4, 5, 2, 0]], dtype=float)

undirected_SP_limit_2 = np.array([[0, np.inf, np.inf, 1, 2],
                                  [np.inf, 0, np.inf, 2, np.inf],
                                  [np.inf, np.inf, 0, np.inf, np.inf],
                                  [1, 2, np.inf, 0, 2],
                                  [2, np.inf, np.inf, 2, 0]], dtype=float)

undirected_SP_limit_0 = np.ones((5, 5), dtype=float) - np.eye(5)
undirected_SP_limit_0[undirected_SP_limit_0 > 0] = np.inf

undirected_pred = np.array([[-9999, 0, 0, 0, 0],
                            [1, -9999, 0, 1, 1],
                            [2, 0, -9999, 0, 0],
                            [3, 3, 0, -9999, 3],
                            [4, 4, 0, 4, -9999]], dtype=float)

methods = ['auto', 'FW', 'D', 'BF', 'J']


def test_dijkstra_limit():
    limits = [0, 2, np.inf]
    results = [undirected_SP_limit_0,
               undirected_SP_limit_2,
               undirected_SP]

    def check(limit, result):
        SP = dijkstra(undirected_G, directed=False, limit=limit)
        assert_array_almost_equal(SP, result)

    for limit, result in zip(limits, results):
        check(limit, result)


def test_directed():
    def check(method):
        SP = shortest_path(directed_G, method=method, directed=True,
                           overwrite=False)
        assert_array_almost_equal(SP, directed_SP)

    for method in methods:
        check(method)


def test_undirected():
    def check(method, directed_in):
        if directed_in:
            SP1 = shortest_path(directed_G, method=method, directed=False,
                                overwrite=False)
            assert_array_almost_equal(SP1, undirected_SP)
        else:
            SP2 = shortest_path(undirected_G, method=method, directed=True,
                                overwrite=False)
            assert_array_almost_equal(SP2, undirected_SP)

    for method in methods:
        for directed_in in (True, False):
            check(method, directed_in)


def test_shortest_path_indices():
    indices = np.arange(4)

    def check(func, indshape):
        outshape = indshape + (5,)
        SP = func(directed_G, directed=False,
                  indices=indices.reshape(indshape))
        assert_array_almost_equal(SP, undirected_SP[indices].reshape(outshape))

    for indshape in [(4,), (4, 1), (2, 2)]:
        for func in (dijkstra, bellman_ford, johnson, shortest_path):
            check(func, indshape)

    assert_raises(ValueError, shortest_path, directed_G, method='FW',
                  indices=indices)


def test_predecessors():
    SP_res = {True: directed_SP,
              False: undirected_SP}
    pred_res = {True: directed_pred,
                False: undirected_pred}

    def check(method, directed):
        SP, pred = shortest_path(directed_G, method, directed=directed,
                                 overwrite=False,
                                 return_predecessors=True)
        assert_array_almost_equal(SP, SP_res[directed])
        assert_array_almost_equal(pred, pred_res[directed])

    for method in methods:
        for directed in (True, False):
            check(method, directed)


def test_construct_shortest_path():
    def check(method, directed):
        SP1, pred = shortest_path(directed_G,
                                  directed=directed,
                                  overwrite=False,
                                  return_predecessors=True)
        SP2 = construct_dist_matrix(directed_G, pred, directed=directed)
        assert_array_almost_equal(SP1, SP2)

    for method in methods:
        for directed in (True, False):
            check(method, directed)


def test_unweighted_path():
    def check(method, directed):
        SP1 = shortest_path(directed_G,
                            directed=directed,
                            overwrite=False,
                            unweighted=True)
        SP2 = shortest_path(unweighted_G,
                            directed=directed,
                            overwrite=False,
                            unweighted=False)
        assert_array_almost_equal(SP1, SP2)

    for method in methods:
        for directed in (True, False):
            check(method, directed)


def test_negative_cycles():
    # create a small graph with a negative cycle
    graph = np.ones([5, 5])
    graph.flat[::6] = 0
    graph[1, 2] = -2

    def check(method, directed):
        assert_raises(NegativeCycleError, shortest_path, graph, method,
                      directed)

    for method in ['FW', 'J', 'BF']:
        for directed in (True, False):
            check(method, directed)


def test_masked_input():
    G = np.ma.masked_equal(directed_G, 0)

    def check(method):
        SP = shortest_path(directed_G, method=method, directed=True,
                           overwrite=False)
        assert_array_almost_equal(SP, directed_SP)

    for method in methods:
        check(method)


def test_overwrite():
    G = np.array([[0, 3, 3, 1, 2],
                  [3, 0, 0, 2, 4],
                  [3, 0, 0, 0, 0],
                  [1, 2, 0, 0, 2],
                  [2, 4, 0, 2, 0]], dtype=float)
    foo = G.copy()
    shortest_path(foo, overwrite=False)
    assert_array_equal(foo, G)

