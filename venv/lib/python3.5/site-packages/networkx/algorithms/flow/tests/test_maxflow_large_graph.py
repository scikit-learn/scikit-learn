# -*- coding: utf-8 -*-
"""Maximum flow algorithms test suite on large graphs.
"""

__author__ = """Loïc Séguin-C. <loicseguin@gmail.com>"""
# Copyright (C) 2010 Loïc Séguin-C. <loicseguin@gmail.com>
# All rights reserved.
# BSD license.
import os
from nose.tools import *

import networkx as nx
from networkx.algorithms.flow import build_flow_dict, build_residual_network
from networkx.algorithms.flow import boykov_kolmogorov
from networkx.algorithms.flow import dinitz
from networkx.algorithms.flow import edmonds_karp
from networkx.algorithms.flow import preflow_push
from networkx.algorithms.flow import shortest_augmenting_path

flow_funcs = [
    boykov_kolmogorov,
    dinitz,
    edmonds_karp,
    preflow_push,
    shortest_augmenting_path,
]

msg = "Assertion failed in function: {0}"


def gen_pyramid(N):
    # This graph admits a flow of value 1 for which every arc is at
    # capacity (except the arcs incident to the sink which have
    # infinite capacity).
    G = nx.DiGraph()

    for i in range(N - 1):
        cap = 1. / (i + 2)
        for j in range(i + 1):
            G.add_edge((i, j), (i + 1, j),
                       capacity=cap)
            cap = 1. / (i + 1) - cap
            G.add_edge((i, j), (i + 1, j + 1),
                       capacity=cap)
            cap = 1. / (i + 2) - cap

    for j in range(N):
        G.add_edge((N - 1, j), 't')

    return G


def read_graph(name):
    dirname = os.path.dirname(__file__)
    path = os.path.join(dirname, name + '.gpickle.bz2')
    return nx.read_gpickle(path)


def validate_flows(G, s, t, soln_value, R, flow_func):
    flow_value = R.graph['flow_value']
    flow_dict = build_flow_dict(G, R)
    assert_equal(soln_value, flow_value, msg=msg.format(flow_func.__name__))
    assert_equal(set(G), set(flow_dict), msg=msg.format(flow_func.__name__))
    for u in G:
        assert_equal(set(G[u]), set(flow_dict[u]),
                     msg=msg.format(flow_func.__name__))
    excess = {u: 0 for u in flow_dict}
    for u in flow_dict:
        for v, flow in flow_dict[u].items():
            ok_(flow <= G[u][v].get('capacity', float('inf')),
                msg=msg.format(flow_func.__name__))
            ok_(flow >= 0, msg=msg.format(flow_func.__name__))
            excess[u] -= flow
            excess[v] += flow
    for u, exc in excess.items():
        if u == s:
            assert_equal(exc, -soln_value, msg=msg.format(flow_func.__name__))
        elif u == t:
            assert_equal(exc, soln_value, msg=msg.format(flow_func.__name__))
        else:
            assert_equal(exc, 0, msg=msg.format(flow_func.__name__))


class TestMaxflowLargeGraph:

    def test_complete_graph(self):
        N = 50
        G = nx.complete_graph(N)
        nx.set_edge_attributes(G, 5, 'capacity')
        R = build_residual_network(G, 'capacity')
        kwargs = dict(residual=R)

        for flow_func in flow_funcs:
            kwargs['flow_func'] = flow_func
            flow_value = nx.maximum_flow_value(G, 1, 2, **kwargs)
            assert_equal(flow_value, 5 * (N - 1),
                         msg=msg.format(flow_func.__name__))

    def test_pyramid(self):
        N = 10
        # N = 100 # this gives a graph with 5051 nodes
        G = gen_pyramid(N)
        R = build_residual_network(G, 'capacity')
        kwargs = dict(residual=R)

        for flow_func in flow_funcs:
            kwargs['flow_func'] = flow_func
            flow_value = nx.maximum_flow_value(G, (0, 0), 't', **kwargs)
            assert_almost_equal(flow_value, 1.,
                                msg=msg.format(flow_func.__name__))

    def test_gl1(self):
        G = read_graph('gl1')
        s = 1
        t = len(G)
        R = build_residual_network(G, 'capacity')
        kwargs = dict(residual=R)

        # do one flow_func to save time
        flow_func = flow_funcs[0]
        validate_flows(G, s, t, 156545, flow_func(G, s, t, **kwargs),
                       flow_func)
#        for flow_func in flow_funcs:
#            validate_flows(G, s, t, 156545, flow_func(G, s, t, **kwargs),
#                           flow_func)

    def test_gw1(self):
        G = read_graph('gw1')
        s = 1
        t = len(G)
        R = build_residual_network(G, 'capacity')
        kwargs = dict(residual=R)

        for flow_func in flow_funcs:
            validate_flows(G, s, t, 1202018, flow_func(G, s, t, **kwargs),
                           flow_func)

    def test_wlm3(self):
        G = read_graph('wlm3')
        s = 1
        t = len(G)
        R = build_residual_network(G, 'capacity')
        kwargs = dict(residual=R)

        # do one flow_func to save time
        flow_func = flow_funcs[0]
        validate_flows(G, s, t, 11875108, flow_func(G, s, t, **kwargs),
                       flow_func)
#        for flow_func in flow_funcs:
#            validate_flows(G, s, t, 11875108, flow_func(G, s, t, **kwargs),
#                           flow_func)

    def test_preflow_push_global_relabel(self):
        G = read_graph('gw1')
        R = preflow_push(G, 1, len(G), global_relabel_freq=50)
        assert_equal(R.graph['flow_value'], 1202018)
