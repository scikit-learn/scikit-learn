# duplication.py - functions for generating graphs by duplicating nodes
#
# Copyright 2016-2018 NetworkX developers.
# Copyright (C) 2004-2018 by
# Aric Hagberg <hagberg@lanl.gov>
# Dan Schult <dschult@colgate.edu>
# Pieter Swart <swart@lanl.gov>
#
# This file is part of NetworkX.
#
# NetworkX is distributed under a BSD license; see LICENSE.txt for more
# information.
"""Functions for generating graphs based on the "duplication" method.

These graph generators start with a small initial graph then duplicate
nodes and (partially) duplicate their edges. These functions are
generally inspired by biological networks.

"""
import random

import networkx as nx
from networkx.exception import NetworkXError

__all__ = ['partial_duplication_graph', 'duplication_divergence_graph']


def partial_duplication_graph(N, n, p, q, seed=None):
    """Return a random graph using the partial duplication model.

    Parameters
    ----------
    N : int
        The total number of nodes in the final graph.

    n : int
        The number of nodes in the initial clique.

    p : float
        The probability of joining each neighbor of a node to the
        duplicate node. Must be a number in the between zero and one,
        inclusive.

    q : float
        The probability of joining the source node to the duplicate
        node. Must be a number in the between zero and one, inclusive.

    seed : int, optional
        Seed for random number generator (default=None).

    Notes
    -----
    A graph of nodes is grown by creating a fully connected graph
    of size `n`. The following procedure is then repeated until
    a total of `N` nodes have been reached.

    1. A random node, *u*, is picked and a new node, *v*, is created.
    2. For each neighbor of *u* an edge from the neighbor to *v* is created
       with probability `p`.
    3. An edge from *u* to *v* is created with probability `q`.

    This algorithm appears in [1].

    This implementation allows the possibility of generating
    disconnected graphs.

    References
    ----------
    .. [1] Knudsen Michael, and Carsten Wiuf. "A Markov chain approach to
           randomly grown graphs." Journal of Applied Mathematics 2008.
           <https://dx.doi.org/10.1155/2008/190836>

    """
    if p < 0 or p > 1 or q < 0 or q > 1:
        msg = "partial duplication graph must have 0 <= p, q <= 1."
        raise NetworkXError(msg)
    if n > N:
        raise NetworkXError("partial duplication graph must have n <= N.")
    if seed is not None:
        random.seed(seed)

    G = nx.complete_graph(n)
    for new_node in range(n, N):
        # Add a new vertex, v, to the graph.
        G.add_node(new_node)

        # Pick a random vertex, u, already in the graph.
        src_node = random.randint(0, new_node)

        # Join v and u with probability q.
        if random.random() < q:
            G.add_edge(new_node, src_node)

        # For each neighbor of u...
        for neighbor_node in list(nx.all_neighbors(G, src_node)):
            # Add the neighbor to v with probability p.
            if random.random() < p:
                G.add_edge(new_node, neighbor_node)
    return G


def duplication_divergence_graph(n, p, seed=None):
    """Returns an undirected graph using the duplication-divergence model.

    A graph of `n` nodes is created by duplicating the initial nodes
    and retaining edges incident to the original nodes with a retention
    probability `p`.

    Parameters
    ----------
    n : int
        The desired number of nodes in the graph.
    p : float
        The probability for retaining the edge of the replicated node.
    seed : int, optional
        A seed for the random number generator of :mod:`random` (default=None).

    Returns
    -------
    G : Graph

    Raises
    ------
    NetworkXError
        If `p` is not a valid probability.
        If `n` is less than 2.

    Notes
    -----
    This algorithm appears in [1].

    This implementation disallows the possibility of generating
    disconnected graphs.

    References
    ----------
    .. [1] I. Ispolatov, P. L. Krapivsky, A. Yuryev,
       "Duplication-divergence model of protein interaction network",
       Phys. Rev. E, 71, 061911, 2005.

    """
    if p > 1 or p < 0:
        msg = "NetworkXError p={0} is not in [0,1].".format(p)
        raise nx.NetworkXError(msg)
    if n < 2:
        msg = 'n must be greater than or equal to 2'
        raise nx.NetworkXError(msg)
    if seed is not None:
        random.seed(seed)

    G = nx.Graph()

    # Initialize the graph with two connected nodes.
    G.add_edge(0, 1)
    i = 2
    while i < n:
        # Choose a random node from current graph to duplicate.
        random_node = random.choice(list(G))
        # Make the replica.
        G.add_node(i)
        # flag indicates whether at least one edge is connected on the replica.
        flag = False
        for nbr in G.neighbors(random_node):
            if random.random() < p:
                # Link retention step.
                G.add_edge(i, nbr)
                flag = True
        if not flag:
            # Delete replica if no edges retained.
            G.remove_node(i)
        else:
            # Successful duplication.
            i += 1
    return G
