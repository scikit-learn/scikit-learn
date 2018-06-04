# -*- coding: utf-8 -*-
#    Copyright (C) 2004-2018 by
#    Aric Hagberg <hagberg@lanl.gov>
#    Dan Schult <dschult@colgate.edu>
#    Pieter Swart <swart@lanl.gov>
#    All rights reserved.
#    BSD license.
"""
Generators for random graphs.

"""

from __future__ import division
import itertools
import math
import random

import networkx as nx
from .classic import empty_graph, path_graph, complete_graph
from .degree_seq import degree_sequence_tree
from collections import defaultdict

__all__ = ['fast_gnp_random_graph',
           'gnp_random_graph',
           'dense_gnm_random_graph',
           'gnm_random_graph',
           'erdos_renyi_graph',
           'binomial_graph',
           'newman_watts_strogatz_graph',
           'watts_strogatz_graph',
           'connected_watts_strogatz_graph',
           'random_regular_graph',
           'barabasi_albert_graph',
           'extended_barabasi_albert_graph',
           'powerlaw_cluster_graph',
           'random_lobster',
           'random_shell_graph',
           'random_powerlaw_tree',
           'random_powerlaw_tree_sequence',
           'random_kernel_graph']


#-------------------------------------------------------------------------
#  Some Famous Random Graphs
#-------------------------------------------------------------------------


def fast_gnp_random_graph(n, p, seed=None, directed=False):
    """Returns a $G_{n,p}$ random graph, also known as an Erdős-Rényi graph or
    a binomial graph.

    Parameters
    ----------
    n : int
        The number of nodes.
    p : float
        Probability for edge creation.
    seed : int, optional
        Seed for random number generator (default=None).
    directed : bool, optional (default=False)
        If True, this function returns a directed graph.

    Notes
    -----
    The $G_{n,p}$ graph algorithm chooses each of the $[n (n - 1)] / 2$
    (undirected) or $n (n - 1)$ (directed) possible edges with probability $p$.

    This algorithm [1]_ runs in $O(n + m)$ time, where `m` is the expected number of
    edges, which equals $p n (n - 1) / 2$. This should be faster than
    :func:`gnp_random_graph` when $p$ is small and the expected number of edges
    is small (that is, the graph is sparse).

    See Also
    --------
    gnp_random_graph

    References
    ----------
    .. [1] Vladimir Batagelj and Ulrik Brandes,
       "Efficient generation of large random networks",
       Phys. Rev. E, 71, 036113, 2005.
    """
    G = empty_graph(n)

    if seed is not None:
        random.seed(seed)

    if p <= 0 or p >= 1:
        return nx.gnp_random_graph(n, p, directed=directed)

    w = -1
    lp = math.log(1.0 - p)

    if directed:
        G = nx.DiGraph(G)
        # Nodes in graph are from 0,n-1 (start with v as the first node index).
        v = 0
        while v < n:
            lr = math.log(1.0 - random.random())
            w = w + 1 + int(lr / lp)
            if v == w:  # avoid self loops
                w = w + 1
            while v < n <= w:
                w = w - n
                v = v + 1
                if v == w:  # avoid self loops
                    w = w + 1
            if v < n:
                G.add_edge(v, w)
    else:
        # Nodes in graph are from 0,n-1 (start with v as the second node index).
        v = 1
        while v < n:
            lr = math.log(1.0 - random.random())
            w = w + 1 + int(lr / lp)
            while w >= v and v < n:
                w = w - v
                v = v + 1
            if v < n:
                G.add_edge(v, w)
    return G


def gnp_random_graph(n, p, seed=None, directed=False):
    """Returns a $G_{n,p}$ random graph, also known as an Erdős-Rényi graph
    or a binomial graph.

    The $G_{n,p}$ model chooses each of the possible edges with probability $p$.

    The functions :func:`binomial_graph` and :func:`erdos_renyi_graph` are
    aliases of this function.

    Parameters
    ----------
    n : int
        The number of nodes.
    p : float
        Probability for edge creation.
    seed : int, optional
        Seed for random number generator (default=None).
    directed : bool, optional (default=False)
        If True, this function returns a directed graph.

    See Also
    --------
    fast_gnp_random_graph

    Notes
    -----
    This algorithm [2]_ runs in $O(n^2)$ time.  For sparse graphs (that is, for
    small values of $p$), :func:`fast_gnp_random_graph` is a faster algorithm.

    References
    ----------
    .. [1] P. Erdős and A. Rényi, On Random Graphs, Publ. Math. 6, 290 (1959).
    .. [2] E. N. Gilbert, Random Graphs, Ann. Math. Stat., 30, 1141 (1959).
    """
    if directed:
        G = nx.DiGraph()
    else:
        G = nx.Graph()
    G.add_nodes_from(range(n))
    if p <= 0:
        return G
    if p >= 1:
        return complete_graph(n, create_using=G)

    if seed is not None:
        random.seed(seed)

    if G.is_directed():
        edges = itertools.permutations(range(n), 2)
    else:
        edges = itertools.combinations(range(n), 2)

    for e in edges:
        if random.random() < p:
            G.add_edge(*e)
    return G


# add some aliases to common names
binomial_graph = gnp_random_graph
erdos_renyi_graph = gnp_random_graph


def dense_gnm_random_graph(n, m, seed=None):
    """Returns a $G_{n,m}$ random graph.

    In the $G_{n,m}$ model, a graph is chosen uniformly at random from the set
    of all graphs with $n$ nodes and $m$ edges.

    This algorithm should be faster than :func:`gnm_random_graph` for dense
    graphs.

    Parameters
    ----------
    n : int
        The number of nodes.
    m : int
        The number of edges.
    seed : int, optional
        Seed for random number generator (default=None).

    See Also
    --------
    gnm_random_graph()

    Notes
    -----
    Algorithm by Keith M. Briggs Mar 31, 2006.
    Inspired by Knuth's Algorithm S (Selection sampling technique),
    in section 3.4.2 of [1]_.

    References
    ----------
    .. [1] Donald E. Knuth, The Art of Computer Programming,
        Volume 2/Seminumerical algorithms, Third Edition, Addison-Wesley, 1997.
    """
    mmax = n * (n - 1) / 2
    if m >= mmax:
        G = complete_graph(n)
    else:
        G = empty_graph(n)

    if n == 1 or m >= mmax:
        return G

    if seed is not None:
        random.seed(seed)

    u = 0
    v = 1
    t = 0
    k = 0
    while True:
        if random.randrange(mmax - t) < m - k:
            G.add_edge(u, v)
            k += 1
            if k == m:
                return G
        t += 1
        v += 1
        if v == n:  # go to next row of adjacency matrix
            u += 1
            v = u + 1


def gnm_random_graph(n, m, seed=None, directed=False):
    """Returns a $G_{n,m}$ random graph.

    In the $G_{n,m}$ model, a graph is chosen uniformly at random from the set
    of all graphs with $n$ nodes and $m$ edges.

    This algorithm should be faster than :func:`dense_gnm_random_graph` for
    sparse graphs.

    Parameters
    ----------
    n : int
        The number of nodes.
    m : int
        The number of edges.
    seed : int, optional
        Seed for random number generator (default=None).
    directed : bool, optional (default=False)
        If True return a directed graph

    See also
    --------
    dense_gnm_random_graph

    """
    if directed:
        G = nx.DiGraph()
    else:
        G = nx.Graph()
    G.add_nodes_from(range(n))

    if seed is not None:
        random.seed(seed)

    if n == 1:
        return G
    max_edges = n * (n - 1)
    if not directed:
        max_edges /= 2.0
    if m >= max_edges:
        return complete_graph(n, create_using=G)

    nlist = list(G)
    edge_count = 0
    while edge_count < m:
        # generate random edge,u,v
        u = random.choice(nlist)
        v = random.choice(nlist)
        if u == v or G.has_edge(u, v):
            continue
        else:
            G.add_edge(u, v)
            edge_count = edge_count + 1
    return G


def newman_watts_strogatz_graph(n, k, p, seed=None):
    """Return a Newman–Watts–Strogatz small-world graph.

    Parameters
    ----------
    n : int
        The number of nodes.
    k : int
        Each node is joined with its `k` nearest neighbors in a ring
        topology.
    p : float
        The probability of adding a new edge for each edge.
    seed : int, optional
        The seed for the random number generator (the default is None).

    Notes
    -----
    First create a ring over $n$ nodes [1]_.  Then each node in the ring is
    connected with its $k$ nearest neighbors (or $k - 1$ neighbors if $k$
    is odd).  Then shortcuts are created by adding new edges as follows: for
    each edge $(u, v)$ in the underlying "$n$-ring with $k$ nearest
    neighbors" with probability $p$ add a new edge $(u, w)$ with
    randomly-chosen existing node $w$.  In contrast with
    :func:`watts_strogatz_graph`, no edges are removed.

    See Also
    --------
    watts_strogatz_graph()

    References
    ----------
    .. [1] M. E. J. Newman and D. J. Watts,
       Renormalization group analysis of the small-world network model,
       Physics Letters A, 263, 341, 1999.
       http://dx.doi.org/10.1016/S0375-9601(99)00757-4
    """
    if seed is not None:
        random.seed(seed)
    if k >= n:
        raise nx.NetworkXError("k>=n, choose smaller k or larger n")
    G = empty_graph(n)
    nlist = list(G.nodes())
    fromv = nlist
    # connect the k/2 neighbors
    for j in range(1, k // 2 + 1):
        tov = fromv[j:] + fromv[0:j]  # the first j are now last
        for i in range(len(fromv)):
            G.add_edge(fromv[i], tov[i])
    # for each edge u-v, with probability p, randomly select existing
    # node w and add new edge u-w
    e = list(G.edges())
    for (u, v) in e:
        if random.random() < p:
            w = random.choice(nlist)
            # no self-loops and reject if edge u-w exists
            # is that the correct NWS model?
            while w == u or G.has_edge(u, w):
                w = random.choice(nlist)
                if G.degree(u) >= n - 1:
                    break  # skip this rewiring
            else:
                G.add_edge(u, w)
    return G


def watts_strogatz_graph(n, k, p, seed=None):
    """Return a Watts–Strogatz small-world graph.

    Parameters
    ----------
    n : int
        The number of nodes
    k : int
        Each node is joined with its `k` nearest neighbors in a ring
        topology.
    p : float
        The probability of rewiring each edge
    seed : int, optional
        Seed for random number generator (default=None)

    See Also
    --------
    newman_watts_strogatz_graph()
    connected_watts_strogatz_graph()

    Notes
    -----
    First create a ring over $n$ nodes [1]_.  Then each node in the ring is joined
    to its $k$ nearest neighbors (or $k - 1$ neighbors if $k$ is odd).
    Then shortcuts are created by replacing some edges as follows: for each
    edge $(u, v)$ in the underlying "$n$-ring with $k$ nearest neighbors"
    with probability $p$ replace it with a new edge $(u, w)$ with uniformly
    random choice of existing node $w$.

    In contrast with :func:`newman_watts_strogatz_graph`, the random rewiring
    does not increase the number of edges. The rewired graph is not guaranteed
    to be connected as in :func:`connected_watts_strogatz_graph`.

    References
    ----------
    .. [1] Duncan J. Watts and Steven H. Strogatz,
       Collective dynamics of small-world networks,
       Nature, 393, pp. 440--442, 1998.
    """
    if k >= n:
        raise nx.NetworkXError("k>=n, choose smaller k or larger n")
    if seed is not None:
        random.seed(seed)

    G = nx.Graph()
    nodes = list(range(n))  # nodes are labeled 0 to n-1
    # connect each node to k/2 neighbors
    for j in range(1, k // 2 + 1):
        targets = nodes[j:] + nodes[0:j]  # first j nodes are now last in list
        G.add_edges_from(zip(nodes, targets))
    # rewire edges from each node
    # loop over all nodes in order (label) and neighbors in order (distance)
    # no self loops or multiple edges allowed
    for j in range(1, k // 2 + 1):  # outer loop is neighbors
        targets = nodes[j:] + nodes[0:j]  # first j nodes are now last in list
        # inner loop in node order
        for u, v in zip(nodes, targets):
            if random.random() < p:
                w = random.choice(nodes)
                # Enforce no self-loops or multiple edges
                while w == u or G.has_edge(u, w):
                    w = random.choice(nodes)
                    if G.degree(u) >= n - 1:
                        break  # skip this rewiring
                else:
                    G.remove_edge(u, v)
                    G.add_edge(u, w)
    return G


def connected_watts_strogatz_graph(n, k, p, tries=100, seed=None):
    """Returns a connected Watts–Strogatz small-world graph.

    Attempts to generate a connected graph by repeated generation of
    Watts–Strogatz small-world graphs.  An exception is raised if the maximum
    number of tries is exceeded.

    Parameters
    ----------
    n : int
        The number of nodes
    k : int
        Each node is joined with its `k` nearest neighbors in a ring
        topology.
    p : float
        The probability of rewiring each edge
    tries : int
        Number of attempts to generate a connected graph.
    seed : int, optional
         The seed for random number generator.

    See Also
    --------
    newman_watts_strogatz_graph()
    watts_strogatz_graph()

    """
    for i in range(tries):
        G = watts_strogatz_graph(n, k, p, seed)
        if nx.is_connected(G):
            return G
    raise nx.NetworkXError('Maximum number of tries exceeded')


def random_regular_graph(d, n, seed=None):
    r"""Returns a random $d$-regular graph on $n$ nodes.

    The resulting graph has no self-loops or parallel edges.

    Parameters
    ----------
    d : int
      The degree of each node.
    n : integer
      The number of nodes. The value of $n \times d$ must be even.
    seed : hashable object
        The seed for random number generator.

    Notes
    -----
    The nodes are numbered from $0$ to $n - 1$.

    Kim and Vu's paper [2]_ shows that this algorithm samples in an
    asymptotically uniform way from the space of random graphs when
    $d = O(n^{1 / 3 - \epsilon})$.

    Raises
    ------

    NetworkXError
        If $n \times d$ is odd or $d$ is greater than or equal to $n$.

    References
    ----------
    .. [1] A. Steger and N. Wormald,
       Generating random regular graphs quickly,
       Probability and Computing 8 (1999), 377-396, 1999.
       http://citeseer.ist.psu.edu/steger99generating.html

    .. [2] Jeong Han Kim and Van H. Vu,
       Generating random regular graphs,
       Proceedings of the thirty-fifth ACM symposium on Theory of computing,
       San Diego, CA, USA, pp 213--222, 2003.
       http://portal.acm.org/citation.cfm?id=780542.780576
    """
    if (n * d) % 2 != 0:
        raise nx.NetworkXError("n * d must be even")

    if not 0 <= d < n:
        raise nx.NetworkXError("the 0 <= d < n inequality must be satisfied")

    if d == 0:
        return empty_graph(n)

    if seed is not None:
        random.seed(seed)

    def _suitable(edges, potential_edges):
        # Helper subroutine to check if there are suitable edges remaining
        # If False, the generation of the graph has failed
        if not potential_edges:
            return True
        for s1 in potential_edges:
            for s2 in potential_edges:
                # Two iterators on the same dictionary are guaranteed
                # to visit it in the same order if there are no
                # intervening modifications.
                if s1 == s2:
                    # Only need to consider s1-s2 pair one time
                    break
                if s1 > s2:
                    s1, s2 = s2, s1
                if (s1, s2) not in edges:
                    return True
        return False

    def _try_creation():
        # Attempt to create an edge set

        edges = set()
        stubs = list(range(n)) * d

        while stubs:
            potential_edges = defaultdict(lambda: 0)
            random.shuffle(stubs)
            stubiter = iter(stubs)
            for s1, s2 in zip(stubiter, stubiter):
                if s1 > s2:
                    s1, s2 = s2, s1
                if s1 != s2 and ((s1, s2) not in edges):
                    edges.add((s1, s2))
                else:
                    potential_edges[s1] += 1
                    potential_edges[s2] += 1

            if not _suitable(edges, potential_edges):
                return None  # failed to find suitable edge set

            stubs = [node for node, potential in potential_edges.items()
                     for _ in range(potential)]
        return edges

    # Even though a suitable edge set exists,
    # the generation of such a set is not guaranteed.
    # Try repeatedly to find one.
    edges = _try_creation()
    while edges is None:
        edges = _try_creation()

    G = nx.Graph()
    G.add_edges_from(edges)

    return G


def _random_subset(seq, m):
    """ Return m unique elements from seq.

    This differs from random.sample which can return repeated
    elements if seq holds repeated elements.
    """
    targets = set()
    while len(targets) < m:
        x = random.choice(seq)
        targets.add(x)
    return targets


def barabasi_albert_graph(n, m, seed=None):
    """Returns a random graph according to the Barabási–Albert preferential
    attachment model.

    A graph of $n$ nodes is grown by attaching new nodes each with $m$
    edges that are preferentially attached to existing nodes with high degree.

    Parameters
    ----------
    n : int
        Number of nodes
    m : int
        Number of edges to attach from a new node to existing nodes
    seed : int, optional
        Seed for random number generator (default=None).

    Returns
    -------
    G : Graph

    Raises
    ------
    NetworkXError
        If `m` does not satisfy ``1 <= m < n``.

    References
    ----------
    .. [1] A. L. Barabási and R. Albert "Emergence of scaling in
       random networks", Science 286, pp 509-512, 1999.
    """

    if m < 1 or m >= n:
        raise nx.NetworkXError("Barabási–Albert network must have m >= 1"
                               " and m < n, m = %d, n = %d" % (m, n))
    if seed is not None:
        random.seed(seed)

    # Add m initial nodes (m0 in barabasi-speak)
    G = empty_graph(m)
    # Target nodes for new edges
    targets = list(range(m))
    # List of existing nodes, with nodes repeated once for each adjacent edge
    repeated_nodes = []
    # Start adding the other n-m nodes. The first node is m.
    source = m
    while source < n:
        # Add edges to m nodes from the source.
        G.add_edges_from(zip([source] * m, targets))
        # Add one node to the list for each new edge just created.
        repeated_nodes.extend(targets)
        # And the new node "source" has m edges to add to the list.
        repeated_nodes.extend([source] * m)
        # Now choose m unique nodes from the existing nodes
        # Pick uniformly from repeated_nodes (preferential attachement)
        targets = _random_subset(repeated_nodes, m)
        source += 1
    return G


def extended_barabasi_albert_graph(n, m, p, q, seed=None):
    """Returns an extended Barabási–Albert model graph.

    An extended Barabási–Albert model graph is a random graph constructed
    using preferential attachment. The extended model allows new egdes,
    rewired edges or new nodes. Based on the probabilities $p$ and $q$
    with $p + q < 1$, the growing behavior of the graph is determined as:

    1) With $p$ probability, $m$ new edges are added to the graph,
    starting from randomly chosen existing nodes and attached preferentially at the other end.

    2) With $q$ probability, $m$ existing edges are rewired
    by randomly chosing an edge and rewiring one end to a preferentially chosen node.

    3) With $(1 - p - q)$ probability, $m$ new nodes are added to the graph
    with edges attached preferentially.

    When $p = q = 0$, the model behaves just like the Barabási–Alber mo

    Parameters
    ----------
    n : int
        Number of nodes
    m : int
        Number of edges with which a new node attaches to existing nodes
    p : float
        Probability value for adding an edge between existing nodes. p + q < 1
    q : float
        Probability value of rewiring of existing edges. p + q < 1
    seed : int (optional, default: None)
        Seed for random number generator

    Returns
    -------
    G : Graph

    Raises
    ------
    NetworkXError
        If `m` does not satisfy ``1 <= m < n`` or ``1 >= p + q``

    References
    ----------
    .. [1] Albert, R., & Barabási, A. L. (2000)
       Topology of evolving networks: local events and universality
       Physical review letters, 85(24), 5234.
    """
    if m < 1 or m >= n:
        msg = "Extended Barabasi-Albert network needs m>=1 and m<n, m=%d, n=%d"
        raise nx.NetworkXError(msg % (m, n))
    if p + q >= 1:
        msg = "Extended Barabasi-Albert network needs p + q <= 1, p=%d, q=%d"
        raise nx.NetworkXError(msg % (p, q))
    if seed is not None:
        random.seed(seed)

    # Add m initial nodes (m0 in barabasi-speak)
    G = empty_graph(m)

    # List of nodes to represent the preferential attachment random selection.
    # At the creation of the graph, all nodes are added to the list
    # so that even nodes that are not connected have a chance to get selected,
    # for rewiring and adding of edges.
    # With each new edge, nodes at the ends of the edge are added to the list.
    attachment_preference = []
    attachment_preference.extend(range(m))

    # Start adding the other n-m nodes. The first node is m.
    new_node = m
    while new_node < n:
        a_probability = random.random()

        # Total number of edges of a Clique of all the nodes
        clique_degree = len(G) - 1
        clique_size = (len(G) * clique_degree) / 2

        # Adding m new edges, if there is room to add them
        if a_probability < p and G.size() <= clique_size - m:
            # Select the nodes where an edge can be added
            elligible_nodes = [nd for nd, deg in G.degree()
                               if deg < clique_degree]
            for i in range(m):
                # Choosing a random source node from elligible_nodes
                src_node = random.choice(elligible_nodes)

                # Picking a possible node that is not 'src_node' or
                # neighbor with 'src_node', with preferential attachment
                prohibited_nodes = list(G[src_node])
                prohibited_nodes.append(src_node)
                # This will raise an exception if the sequence is empty
                dest_node = random.choice([nd for nd in attachment_preference
                                           if nd not in prohibited_nodes])
                # Adding the new edge
                G.add_edge(src_node, dest_node)

                # Appending both nodes to add to their preferential attachment
                attachment_preference.append(src_node)
                attachment_preference.append(dest_node)

                # Adjusting the elligible nodes. Degree may be saturated.
                if G.degree(src_node) == clique_degree:
                    elligible_nodes.remove(src_node)
                if G.degree(dest_node) == clique_degree \
                        and dest_node in elligible_nodes:
                    elligible_nodes.remove(dest_node)

        # Rewiring m edges, if there are enough edges
        elif p <= a_probability < (p + q) and m <= G.size() < clique_size:
            # Selecting nodes that have at least 1 edge but that are not
            # fully connected to ALL other nodes (center of star).
            # These nodes are the pivot nodes of the edges to rewire
            elligible_nodes = [nd for nd, deg in G.degree()
                               if 0 < deg < clique_degree]
            for i in range(m):
                # Choosing a random source node
                node = random.choice(elligible_nodes)

                # The available nodes do have a neighbor at least.
                neighbor_nodes = list(G[node])

                # Choosing the other end that will get dettached
                src_node = random.choice(neighbor_nodes)

                # Picking a target node that is not 'node' or
                # neighbor with 'node', with preferential attachment
                neighbor_nodes.append(node)
                dest_node = random.choice([nd for nd in attachment_preference
                                           if nd not in neighbor_nodes])
                # Rewire
                G.remove_edge(node, src_node)
                G.add_edge(node, dest_node)

                # Adjusting the preferential attachment list
                attachment_preference.remove(src_node)
                attachment_preference.append(dest_node)

                # Adjusting the elligible nodes.
                # nodes may be saturated or isolated.
                if G.degree(src_node) == 0 and src_node in elligible_nodes:
                    elligible_nodes.remove(src_node)
                if dest_node in elligible_nodes:
                    if G.degree(dest_node) == clique_degree:
                        elligible_nodes.remove(dest_node)
                else:
                    if G.degree(dest_node) == 1:
                        elligible_nodes.append(dest_node)

        # Adding new node with m edges
        else:
            # Select the edges' nodes by preferential attachment
            targets = _random_subset(attachment_preference, m)
            G.add_edges_from(zip([new_node] * m, targets))

            # Add one node to the list for each new edge just created.
            attachment_preference.extend(targets)
            # The new node has m edges to it, plus itself: m + 1
            attachment_preference.extend([new_node] * (m + 1))
            new_node += 1
    return G


def powerlaw_cluster_graph(n, m, p, seed=None):
    """Holme and Kim algorithm for growing graphs with powerlaw
    degree distribution and approximate average clustering.

    Parameters
    ----------
    n : int
        the number of nodes
    m : int
        the number of random edges to add for each new node
    p : float,
        Probability of adding a triangle after adding a random edge
    seed : int, optional
        Seed for random number generator (default=None).

    Notes
    -----
    The average clustering has a hard time getting above a certain
    cutoff that depends on `m`.  This cutoff is often quite low.  The
    transitivity (fraction of triangles to possible triangles) seems to
    decrease with network size.

    It is essentially the Barabási–Albert (BA) growth model with an
    extra step that each random edge is followed by a chance of
    making an edge to one of its neighbors too (and thus a triangle).

    This algorithm improves on BA in the sense that it enables a
    higher average clustering to be attained if desired.

    It seems possible to have a disconnected graph with this algorithm
    since the initial `m` nodes may not be all linked to a new node
    on the first iteration like the BA model.

    Raises
    ------
    NetworkXError
        If `m` does not satisfy ``1 <= m <= n`` or `p` does not
        satisfy ``0 <= p <= 1``.

    References
    ----------
    .. [1] P. Holme and B. J. Kim,
       "Growing scale-free networks with tunable clustering",
       Phys. Rev. E, 65, 026107, 2002.
    """

    if m < 1 or n < m:
        raise nx.NetworkXError(
            "NetworkXError must have m>1 and m<n, m=%d,n=%d" % (m, n))

    if p > 1 or p < 0:
        raise nx.NetworkXError(
            "NetworkXError p must be in [0,1], p=%f" % (p))
    if seed is not None:
        random.seed(seed)

    G = empty_graph(m)  # add m initial nodes (m0 in barabasi-speak)
    repeated_nodes = list(G.nodes())  # list of existing nodes to sample from
    # with nodes repeated once for each adjacent edge
    source = m               # next node is m
    while source < n:        # Now add the other n-1 nodes
        possible_targets = _random_subset(repeated_nodes, m)
        # do one preferential attachment for new node
        target = possible_targets.pop()
        G.add_edge(source, target)
        repeated_nodes.append(target)  # add one node to list for each new link
        count = 1
        while count < m:  # add m-1 more new links
            if random.random() < p:  # clustering step: add triangle
                neighborhood = [nbr for nbr in G.neighbors(target)
                                if not G.has_edge(source, nbr)
                                and not nbr == source]
                if neighborhood:  # if there is a neighbor without a link
                    nbr = random.choice(neighborhood)
                    G.add_edge(source, nbr)  # add triangle
                    repeated_nodes.append(nbr)
                    count = count + 1
                    continue  # go to top of while loop
            # else do preferential attachment step if above fails
            target = possible_targets.pop()
            G.add_edge(source, target)
            repeated_nodes.append(target)
            count = count + 1

        repeated_nodes.extend([source] * m)  # add source node to list m times
        source += 1
    return G


def random_lobster(n, p1, p2, seed=None):
    """Returns a random lobster graph.

     A lobster is a tree that reduces to a caterpillar when pruning all
     leaf nodes. A caterpillar is a tree that reduces to a path graph
     when pruning all leaf nodes; setting `p2` to zero produces a caterpillar.

     Parameters
     ----------
     n : int
         The expected number of nodes in the backbone
     p1 : float
         Probability of adding an edge to the backbone
     p2 : float
         Probability of adding an edge one level beyond backbone
     seed : int, optional
         Seed for random number generator (default=None).
    """
    # a necessary ingredient in any self-respecting graph library
    if seed is not None:
        random.seed(seed)
    llen = int(2 * random.random() * n + 0.5)
    L = path_graph(llen)
    # build caterpillar: add edges to path graph with probability p1
    current_node = llen - 1
    for n in range(llen):
        if random.random() < p1:  # add fuzzy caterpillar parts
            current_node += 1
            L.add_edge(n, current_node)
            if random.random() < p2:  # add crunchy lobster bits
                current_node += 1
                L.add_edge(current_node - 1, current_node)
    return L  # voila, un lobster!


def random_shell_graph(constructor, seed=None):
    """Returns a random shell graph for the constructor given.

    Parameters
    ----------
    constructor : list of three-tuples
        Represents the parameters for a shell, starting at the center
        shell.  Each element of the list must be of the form `(n, m,
        d)`, where `n` is the number of nodes in the shell, `m` is
        the number of edges in the shell, and `d` is the ratio of
        inter-shell (next) edges to intra-shell edges. If `d` is zero,
        there will be no intra-shell edges, and if `d` is one there
        will be all possible intra-shell edges.
    seed : int, optional
        Seed for random number generator (default=None).

    Examples
    --------
    >>> constructor = [(10, 20, 0.8), (20, 40, 0.8)]
    >>> G = nx.random_shell_graph(constructor)

    """
    G = empty_graph(0)

    if seed is not None:
        random.seed(seed)

    glist = []
    intra_edges = []
    nnodes = 0
    # create gnm graphs for each shell
    for (n, m, d) in constructor:
        inter_edges = int(m * d)
        intra_edges.append(m - inter_edges)
        g = nx.convert_node_labels_to_integers(
            gnm_random_graph(n, inter_edges),
            first_label=nnodes)
        glist.append(g)
        nnodes += n
        G = nx.operators.union(G, g)

    # connect the shells randomly
    for gi in range(len(glist) - 1):
        nlist1 = list(glist[gi])
        nlist2 = list(glist[gi + 1])
        total_edges = intra_edges[gi]
        edge_count = 0
        while edge_count < total_edges:
            u = random.choice(nlist1)
            v = random.choice(nlist2)
            if u == v or G.has_edge(u, v):
                continue
            else:
                G.add_edge(u, v)
                edge_count = edge_count + 1
    return G


def random_powerlaw_tree(n, gamma=3, seed=None, tries=100):
    """Returns a tree with a power law degree distribution.

    Parameters
    ----------
    n : int
        The number of nodes.
    gamma : float
        Exponent of the power law.
    seed : int, optional
        Seed for random number generator (default=None).
    tries : int
        Number of attempts to adjust the sequence to make it a tree.

    Raises
    ------
    NetworkXError
        If no valid sequence is found within the maximum number of
        attempts.

    Notes
    -----
    A trial power law degree sequence is chosen and then elements are
    swapped with new elements from a powerlaw distribution until the
    sequence makes a tree (by checking, for example, that the number of
    edges is one smaller than the number of nodes).

    """
    # This call may raise a NetworkXError if the number of tries is succeeded.
    seq = random_powerlaw_tree_sequence(n, gamma=gamma, seed=seed, tries=tries)
    G = degree_sequence_tree(seq)
    return G


def random_powerlaw_tree_sequence(n, gamma=3, seed=None, tries=100):
    """Returns a degree sequence for a tree with a power law distribution.

    Parameters
    ----------
    n : int,
        The number of nodes.
    gamma : float
        Exponent of the power law.
    seed : int, optional
        Seed for random number generator (default=None).
    tries : int
        Number of attempts to adjust the sequence to make it a tree.

    Raises
    ------
    NetworkXError
        If no valid sequence is found within the maximum number of
        attempts.

    Notes
    -----
    A trial power law degree sequence is chosen and then elements are
    swapped with new elements from a power law distribution until
    the sequence makes a tree (by checking, for example, that the number of
    edges is one smaller than the number of nodes).

    """
    if seed is not None:
        random.seed(seed)

    # get trial sequence
    z = nx.utils.powerlaw_sequence(n, exponent=gamma)
    # round to integer values in the range [0,n]
    zseq = [min(n, max(int(round(s)), 0)) for s in z]

    # another sequence to swap values from
    z = nx.utils.powerlaw_sequence(tries, exponent=gamma)
    # round to integer values in the range [0,n]
    swap = [min(n, max(int(round(s)), 0)) for s in z]

    for deg in swap:
        # If this degree sequence can be the degree sequence of a tree, return
        # it. It can be a tree if the number of edges is one fewer than the
        # number of nodes, or in other words, `n - sum(zseq) / 2 == 1`. We
        # use an equivalent condition below that avoids floating point
        # operations.
        if 2 * n - sum(zseq) == 2:
            return zseq
        index = random.randint(0, n - 1)
        zseq[index] = swap.pop()

    raise nx.NetworkXError('Exceeded max (%d) attempts for a valid tree'
                           ' sequence.' % tries)


def random_kernel_graph(n, kernel_integral, kernel_root=None, seed=None):
    r"""Return an random graph based on the specified kernel.

    The algorithm chooses each of the $[n(n-1)]/2$ possible edges with
    probability specified by a kernel $\kappa(x,y)$ [1]_.  The kernel
    $\kappa(x,y)$ must be a symmetric (in $x,y$), non-negative,
    bounded function.

    Parameters
    ----------
    n : int
        The number of nodes
    kernal_integral : function
        Function that returns the definite integral of the kernel $\kappa(x,y)$,
        $F(y,a,b) := \int_a^b \kappa(x,y)dx$
    kernel_root: function (optional)
        Function that returns the root $b$ of the equation $F(y,a,b) = r$.
        If None, the root is found using :func:`scipy.optimize.brentq`
        (this requires SciPy).
    seed : int, optional
        Seed for random number generator (default=None)

    Notes
    -----
    The kernel is specified through its definite integral which must be
    provided as one of the arguments. If the integral and root of the
    kernel integral can be found in $O(1)$ time then this algorithm runs in
    time $O(n+m)$ where m is the expected number of edges [2]_.

    The nodes are set to integers from $0$ to $n-1$.

    Examples
    --------
    Generate an Erdős–Rényi random graph $G(n,c/n)$, with kernel
    $\kappa(x,y)=c$ where $c$ is the mean expected degree.

    >>> def integral(u, w, z):
    ...     return c * (z - w)
    >>> def root(u, w, r):
    ...     return r / c + w
    >>> c = 1
    >>> graph = nx.random_kernel_graph(1000, integral, root)

    See Also
    --------
    gnp_random_graph
    expected_degree_graph

    References
    ----------
    .. [1] Bollobás, Béla,  Janson, S. and Riordan, O.
       "The phase transition in inhomogeneous random graphs",
       *Random Structures Algorithms*, 31, 3--122, 2007.

    .. [2] Hagberg A, Lemons N (2015),
       "Fast Generation of Sparse Random Kernel Graphs".
       PLoS ONE 10(9): e0135177, 2015. doi:10.1371/journal.pone.0135177
    """
    if seed is not None:
        random.seed(seed)
    if kernel_root is None:
        import scipy.optimize as optimize

        def kernel_root(y, a, r):
            def my_function(b):
                return kernel_integral(y, a, b) - r
            return optimize.brentq(my_function, a, 1)
    graph = nx.Graph()
    graph.add_nodes_from(range(n))
    (i, j) = (1, 1)
    while i < n:
        r = -math.log(1 - random.random())  # (1-random.random()) in (0, 1]
        if kernel_integral(i / n, j / n, 1) <= r:
            i, j = i + 1, i + 1
        else:
            j = int(math.ceil(n * kernel_root(i / n, j / n, r)))
            graph.add_edge(i - 1, j - 1)
    return graph
