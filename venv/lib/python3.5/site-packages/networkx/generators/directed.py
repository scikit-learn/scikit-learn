# -*- coding: utf-8 -*-
# Copyright (C) 2006-2018 by
#   Aric Hagberg <hagberg@lanl.gov>
#   Dan Schult <dschult@colgate.edu>
#   Pieter Swart <swart@lanl.gov>
# Copyright (C) 2009 by Willem Ligtenberg <W.P.A.Ligtenberg@tue.nl>
# All rights reserved.
# BSD license.
#
# Authors: Aric Hagberg (hagberg@lanl.gov)
#          Willem Ligtenberg (W.P.A.Ligtenberg@tue.nl)
"""
Generators for some directed graphs, including growing network (GN) graphs and
scale-free graphs.

"""
from __future__ import division

from collections import Counter
import random

import networkx as nx
from networkx.generators.classic import empty_graph
from networkx.utils import discrete_sequence
from networkx.utils import weighted_choice

__all__ = ['gn_graph', 'gnc_graph', 'gnr_graph', 'random_k_out_graph',
           'scale_free_graph']


def gn_graph(n, kernel=None, create_using=None, seed=None):
    """Return the growing network (GN) digraph with `n` nodes.

    The GN graph is built by adding nodes one at a time with a link to one
    previously added node.  The target node for the link is chosen with
    probability based on degree.  The default attachment kernel is a linear
    function of the degree of a node.

    The graph is always a (directed) tree.

    Parameters
    ----------
    n : int
        The number of nodes for the generated graph.
    kernel : function
        The attachment kernel.
    create_using : graph, optional (default DiGraph)
        Return graph of this type. The instance will be cleared.
    seed : hashable object, optional
        The seed for the random number generator.

    Examples
    --------
    To create the undirected GN graph, use the :meth:`~DiGraph.to_directed`
    method::

    >>> D = nx.gn_graph(10)  # the GN graph
    >>> G = D.to_undirected()  # the undirected version

    To specify an attachment kernel, use the `kernel` keyword argument::

    >>> D = nx.gn_graph(10, kernel=lambda x: x ** 1.5)  # A_k = k^1.5

    References
    ----------
    .. [1] P. L. Krapivsky and S. Redner,
           Organization of Growing Random Networks,
           Phys. Rev. E, 63, 066123, 2001.
    """
    if create_using is None:
        create_using = nx.DiGraph()
    elif not create_using.is_directed():
        raise nx.NetworkXError("Directed Graph required in create_using")

    if kernel is None:
        def kernel(x): return x

    if seed is not None:
        random.seed(seed)

    G = empty_graph(1, create_using)

    if n == 1:
        return G

    G.add_edge(1, 0)  # get started
    ds = [1, 1]  # degree sequence

    for source in range(2, n):
        # compute distribution from kernel and degree
        dist = [kernel(d) for d in ds]
        # choose target from discrete distribution
        target = discrete_sequence(1, distribution=dist)[0]
        G.add_edge(source, target)
        ds.append(1)  # the source has only one link (degree one)
        ds[target] += 1  # add one to the target link degree
    return G


def gnr_graph(n, p, create_using=None, seed=None):
    """Return the growing network with redirection (GNR) digraph with `n`
    nodes and redirection probability `p`.

    The GNR graph is built by adding nodes one at a time with a link to one
    previously added node.  The previous target node is chosen uniformly at
    random.  With probabiliy `p` the link is instead "redirected" to the
    successor node of the target.

    The graph is always a (directed) tree.

    Parameters
    ----------
    n : int
        The number of nodes for the generated graph.
    p : float
        The redirection probability.
    create_using : graph, optional (default DiGraph)
        Return graph of this type. The instance will be cleared.
    seed : hashable object, optional
        The seed for the random number generator.

    Examples
    --------
    To create the undirected GNR graph, use the :meth:`~DiGraph.to_directed`
    method::

    >>> D = nx.gnr_graph(10, 0.5)  # the GNR graph
    >>> G = D.to_undirected()  # the undirected version

    References
    ----------
    .. [1] P. L. Krapivsky and S. Redner,
           Organization of Growing Random Networks,
           Phys. Rev. E, 63, 066123, 2001.
    """
    if create_using is None:
        create_using = nx.DiGraph()
    elif not create_using.is_directed():
        raise nx.NetworkXError("Directed Graph required in create_using")

    if seed is not None:
        random.seed(seed)

    G = empty_graph(1, create_using)

    if n == 1:
        return G

    for source in range(1, n):
        target = random.randrange(0, source)
        if random.random() < p and target != 0:
            target = next(G.successors(target))
        G.add_edge(source, target)

    return G


def gnc_graph(n, create_using=None, seed=None):
    """Return the growing network with copying (GNC) digraph with `n` nodes.

    The GNC graph is built by adding nodes one at a time with a link to one
    previously added node (chosen uniformly at random) and to all of that
    node's successors.

    Parameters
    ----------
    n : int
        The number of nodes for the generated graph.
    create_using : graph, optional (default DiGraph)
        Return graph of this type. The instance will be cleared.
    seed : hashable object, optional
        The seed for the random number generator.

    References
    ----------
    .. [1] P. L. Krapivsky and S. Redner,
           Network Growth by Copying,
           Phys. Rev. E, 71, 036118, 2005k.},
    """
    if create_using is None:
        create_using = nx.DiGraph()
    elif not create_using.is_directed():
        raise nx.NetworkXError("Directed Graph required in create_using")

    if seed is not None:
        random.seed(seed)

    G = empty_graph(1, create_using)

    if n == 1:
        return G

    for source in range(1, n):
        target = random.randrange(0, source)
        for succ in G.successors(target):
            G.add_edge(source, succ)
        G.add_edge(source, target)

    return G


def scale_free_graph(n, alpha=0.41, beta=0.54, gamma=0.05, delta_in=0.2,
                     delta_out=0, create_using=None, seed=None):
    """Returns a scale-free directed graph.

    Parameters
    ----------
    n : integer
        Number of nodes in graph
    alpha : float
        Probability for adding a new node connected to an existing node
        chosen randomly according to the in-degree distribution.
    beta : float
        Probability for adding an edge between two existing nodes.
        One existing node is chosen randomly according the in-degree
        distribution and the other chosen randomly according to the out-degree
        distribution.
    gamma : float
        Probability for adding a new node connected to an existing node
        chosen randomly according to the out-degree distribution.
    delta_in : float
        Bias for choosing ndoes from in-degree distribution.
    delta_out : float
        Bias for choosing ndoes from out-degree distribution.
    create_using : graph, optional (default MultiDiGraph)
        Use this graph instance to start the process (default=3-cycle).
    seed : integer, optional
        Seed for random number generator

    Examples
    --------
    Create a scale-free graph on one hundred nodes::

    >>> G = nx.scale_free_graph(100)

    Notes
    -----
    The sum of `alpha`, `beta`, and `gamma` must be 1.

    References
    ----------
    .. [1] B. Bollobás, C. Borgs, J. Chayes, and O. Riordan,
           Directed scale-free graphs,
           Proceedings of the fourteenth annual ACM-SIAM Symposium on
           Discrete Algorithms, 132--139, 2003.
    """

    def _choose_node(G, distribution, delta, psum):
        cumsum = 0.0
        # normalization
        r = random.random()
        for n, d in distribution:
            cumsum += (d + delta) / psum
            if r < cumsum:
                break
        return n

    if create_using is None:
        # start with 3-cycle
        G = nx.MultiDiGraph()
        G.add_edges_from([(0, 1), (1, 2), (2, 0)])
    else:
        # keep existing graph structure?
        G = create_using
        if not (G.is_directed() and G.is_multigraph()):
            raise nx.NetworkXError("MultiDiGraph required in create_using")

    if alpha <= 0:
        raise ValueError('alpha must be >= 0.')
    if beta <= 0:
        raise ValueError('beta must be >= 0.')
    if gamma <= 0:
        raise ValueError('beta must be >= 0.')

    if alpha + beta + gamma != 1.0:
        raise ValueError('alpha+beta+gamma must equal 1.')

    # seed random number generated (uses None as default)
    random.seed(seed)
    number_of_edges = G.number_of_edges()
    while len(G) < n:
        psum_in = number_of_edges + delta_in * len(G)
        psum_out = number_of_edges + delta_out * len(G)
        r = random.random()
        # random choice in alpha,beta,gamma ranges
        if r < alpha:
            # alpha
            # add new node v
            v = len(G)
            # choose w according to in-degree and delta_in
            w = _choose_node(G, G.in_degree(), delta_in, psum_in)
        elif r < alpha + beta:
            # beta
            # choose v according to out-degree and delta_out
            v = _choose_node(G, G.out_degree(), delta_out, psum_out)
            # choose w according to in-degree and delta_in
            w = _choose_node(G, G.in_degree(), delta_in, psum_in)
        else:
            # gamma
            # choose v according to out-degree and delta_out
            v = _choose_node(G, G.out_degree(), delta_out, psum_out)
            # add new node w
            w = len(G)
        G.add_edge(v, w)
        number_of_edges += 1

    return G


def random_uniform_k_out_graph(n, k, self_loops=True, with_replacement=True,
                               seed=None):
    """Returns a random `k`-out graph with uniform attachment.

    A random `k`-out graph with uniform attachment is a multidigraph
    generated by the following algorithm. For each node *u*, choose
    `k` nodes *v* uniformly at random (with replacement). Add a
    directed edge joining *u* to *v*.

    Parameters
    ----------
    n : int
        The number of nodes in the returned graph.

    k : int
        The out-degree of each node in the returned graph.

    self_loops : bool
        If True, self-loops are allowed when generating the graph.

    with_replacement : bool
        If True, neighbors are chosen with replacement and the
        returned graph will be a directed multigraph. Otherwise,
        neighbors are chosen without replacement and the returned graph
        will be a directed graph.

    seed: int
        If provided, this is used as the seed for the random number
        generator.

    Returns
    -------
    NetworkX graph
        A `k`-out-regular directed graph generated according to the
        above algorithm. It will be a multigraph if and only if
        `with_replacement` is True.

    Raises
    ------
    ValueError
        If `with_replacement` is False and `k` is greater than
        `n`.

    See also
    --------
    random_k_out_graph

    Notes
    -----
    The return digraph or multidigraph may not be strongly connected, or
    even weakly connected.

    If `with_replacement` is True, this function is similar to
    :func:`random_k_out_graph`, if that function had parameter `alpha`
    set to positive infinity.

    """
    random.seed(seed)

    if with_replacement:
        create_using = nx.MultiDiGraph()

        def sample(v, nodes):
            if not self_loops:
                nodes = nodes - {v}
            return (random.choice(list(nodes)) for i in range(k))

    else:
        create_using = nx.DiGraph()

        def sample(v, nodes):
            if not self_loops:
                nodes = nodes - {v}
            return random.sample(nodes, k)

    G = nx.empty_graph(n, create_using=create_using)
    nodes = set(G)
    for u in G:
        G.add_edges_from((u, v) for v in sample(u, nodes))
    return G


def random_k_out_graph(n, k, alpha, self_loops=True, seed=None):
    """Returns a random `k`-out graph with preferential attachment.

    A random `k`-out graph with preferential attachment is a
    multidigraph generated by the following algorithm.

    1. Begin with an empty digraph, and initially set each node to have
       weight `alpha`.
    2. Choose a node `u` with out-degree less than `k` uniformly at
       random.
    3. Choose a node `v` from with probability proportional to its
       weight.
    4. Add a directed edge from `u` to `v`, and increase the weight
       of `v` by one.
    5. If each node has out-degree `k`, halt, otherwise repeat from
       step 2.

    For more information on this model of random graph, see [1].

    Parameters
    ----------
    n : int
        The number of nodes in the returned graph.

    k : int
        The out-degree of each node in the returned graph.

    alpha : float
        A positive :class:`float` representing the initial weight of
        each vertex. A higher number means that in step 3 above, nodes
        will be chosen more like a true uniformly random sample, and a
        lower number means that nodes are more likely to be chosen as
        their in-degree increases. If this parameter is not positive, a
        :exc:`ValueError` is raised.

    self_loops : bool
        If True, self-loops are allowed when generating the graph.

    seed: int
        If provided, this is used as the seed for the random number
        generator.

    Returns
    -------
    :class:`~networkx.classes.MultiDiGraph`
        A `k`-out-regular multidigraph generated according to the above
        algorithm.

    Raises
    ------
    ValueError
        If `alpha` is not positive.

    Notes
    -----
    The returned multidigraph may not be strongly connected, or even
    weakly connected.

    References
    ----------
    [1]: Peterson, Nicholas R., and Boris Pittel.
         "Distance between two random `k`-out digraphs, with and without
         preferential attachment."
         arXiv preprint arXiv:1311.5961 (2013).
         <https://arxiv.org/abs/1311.5961>

    """
    if alpha < 0:
        raise ValueError('alpha must be positive')
    random.seed(seed)
    G = nx.empty_graph(n, create_using=nx.MultiDiGraph())
    weights = Counter({v: alpha for v in G})
    for i in range(k * n):
        u = random.choice([v for v, d in G.out_degree() if d < k])
        # If self-loops are not allowed, make the source node `u` have
        # weight zero.
        if not self_loops:
            adjustment = Counter({u: weights[u]})
        else:
            adjustment = Counter()
        v = weighted_choice(weights - adjustment)
        G.add_edge(u, v)
        weights[v] += 1
    return G
