# -*- coding: utf-8 -*-
#
# kernighan_lin.py - Kernighan–Lin bipartition algorithm
#
# Copyright 2011 Ben Edwards <bedwards@cs.unm.edu>.
# Copyright 2011 Aric Hagberg <hagberg@lanl.gov>.
# Copyright 2015 NetworkX developers.
#
# This file is part of NetworkX.
#
# NetworkX is distributed under a BSD license; see LICENSE.txt for more
# information.
"""Functions for computing the Kernighan–Lin bipartition algorithm."""
from __future__ import division

from collections import defaultdict
from itertools import islice
from operator import itemgetter
import random

import networkx as nx
from networkx.utils import not_implemented_for
from networkx.algorithms.community.community_utils import is_partition

__all__ = ['kernighan_lin_bisection']


def _compute_delta(G, A, B, weight):
    # helper to compute initial swap deltas for a pass
    delta = defaultdict(float)
    for u, v, d in G.edges(data=True):
        w = d.get(weight, 1)
        if u in A:
            if v in A:
                delta[u] -= w
                delta[v] -= w
            elif v in B:
                delta[u] += w
                delta[v] += w
        elif u in B:
            if v in A:
                delta[u] += w
                delta[v] += w
            elif v in B:
                delta[u] -= w
                delta[v] -= w
    return delta


def _update_delta(delta, G, A, B, u, v, weight):
    # helper to update swap deltas during single pass
    for _, nbr, d in G.edges(u, data=True):
        w = d.get(weight, 1)
        if nbr in A:
            delta[nbr] += 2 * w
        if nbr in B:
            delta[nbr] -= 2 * w
    for _, nbr, d in G.edges(v, data=True):
        w = d.get(weight, 1)
        if nbr in A:
            delta[nbr] -= 2 * w
        if nbr in B:
            delta[nbr] += 2 * w
    return delta


def _kernighan_lin_pass(G, A, B, weight):
    # do a single iteration of Kernighan–Lin algorithm
    # returns list of  (g_i,u_i,v_i) for i node pairs u_i,v_i
    multigraph = G.is_multigraph()
    delta = _compute_delta(G, A, B, weight)
    swapped = set()
    gains = []
    while len(swapped) < len(G):
        gain = []
        for u in A - swapped:
            for v in B - swapped:
                try:
                    if multigraph:
                        w = sum(d.get(weight, 1) for d in G[u][v].values())
                    else:
                        w = G[u][v].get(weight, 1)
                except KeyError:
                    w = 0
                gain.append((delta[u] + delta[v] - 2 * w, u, v))
        if len(gain) == 0:
            break
        maxg, u, v = max(gain, key=itemgetter(0))
        swapped |= {u, v}
        gains.append((maxg, u, v))
        delta = _update_delta(delta, G, A - swapped, B - swapped, u, v, weight)
    return gains


@not_implemented_for('directed')
def kernighan_lin_bisection(G, partition=None, max_iter=10, weight='weight'):
    """Partition a graph into two blocks using the Kernighan–Lin
    algorithm.

    This algorithm paritions a network into two sets by iteratively
    swapping pairs of nodes to reduce the edge cut between the two sets.

    Parameters
    ----------
    G : graph

    partition : tuple
        Pair of iterables containing an intial partition. If not
        specified, a random balanced partition is used.

    max_iter : int
        Maximum number of times to attempt swaps to find an
        improvemement before giving up.

    weight : key
        Edge data key to use as weight. If None, the weights are all
        set to one.

    Returns
    -------
    partition : tuple
        A pair of sets of nodes representing the bipartition.

    Raises
    -------
    NetworkXError
        If partition is not a valid partition of the nodes of the graph.

    References
    ----------
    .. [1] Kernighan, B. W.; Lin, Shen (1970).
       "An efficient heuristic procedure for partitioning graphs."
       *Bell Systems Technical Journal* 49: 291--307.
       Oxford University Press 2011.

    """
    # If no partition is provided, split the nodes randomly into a
    # balanced partition.
    if partition is None:
        nodes = list(G)
        random.shuffle(nodes)
        h = len(nodes) // 2
        partition = (nodes[:h], nodes[h:])
    # Make a copy of the partition as a pair of sets.
    try:
        A, B = set(partition[0]), set(partition[1])
    except:
        raise ValueError('partition must be two sets')
    if not is_partition(G, (A, B)):
        raise nx.NetworkXError('partition invalid')
    for i in range(max_iter):
        # `gains` is a list of triples of the form (g, u, v) for each
        # node pair (u, v), where `g` is the gain of that node pair.
        gains = _kernighan_lin_pass(G, A, B, weight)
        csum = list(nx.utils.accumulate(g for g, u, v in gains))
        max_cgain = max(csum)
        if max_cgain <= 0:
            break
        # Get the node pairs up to the index of the maximum cumulative
        # gain, and collect each `u` into `anodes` and each `v` into
        # `bnodes`, for each pair `(u, v)`.
        index = csum.index(max_cgain)
        nodesets = islice(zip(*gains[:index + 1]), 1, 3)
        anodes, bnodes = (set(s) for s in nodesets)
        A |= bnodes
        A -= anodes
        B |= anodes
        B -= bnodes
    return A, B
