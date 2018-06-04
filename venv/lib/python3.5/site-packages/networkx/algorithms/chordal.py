# -*- coding: utf-8 -*-
"""
Algorithms for chordal graphs.

A graph is chordal if every cycle of length at least 4 has a chord
(an edge joining two nodes not adjacent in the cycle).
https://en.wikipedia.org/wiki/Chordal_graph
"""
import networkx as nx
import random
import sys

__authors__ = "\n".join(['Jesus Cerquides <cerquide@iiia.csic.es>'])
#    Copyright (C) 2010 by
#    Jesus Cerquides <cerquide@iiia.csic.es>
#    All rights reserved.
#    BSD license.

__all__ = ['is_chordal',
           'find_induced_nodes',
           'chordal_graph_cliques',
           'chordal_graph_treewidth',
           'NetworkXTreewidthBoundExceeded']


class NetworkXTreewidthBoundExceeded(nx.NetworkXException):
    """Exception raised when a treewidth bound has been provided and it has
    been exceeded"""


def is_chordal(G):
    """Checks whether G is a chordal graph.

    A graph is chordal if every cycle of length at least 4 has a chord
    (an edge joining two nodes not adjacent in the cycle).

    Parameters
    ----------
    G : graph
      A NetworkX graph.

    Returns
    -------
    chordal : bool
      True if G is a chordal graph and False otherwise.

    Raises
    ------
    NetworkXError
        The algorithm does not support DiGraph, MultiGraph and MultiDiGraph.
        If the input graph is an instance of one of these classes, a
        :exc:`NetworkXError` is raised.

    Examples
    --------
    >>> import networkx as nx
    >>> e=[(1,2),(1,3),(2,3),(2,4),(3,4),(3,5),(3,6),(4,5),(4,6),(5,6)]
    >>> G=nx.Graph(e)
    >>> nx.is_chordal(G)
    True

    Notes
    -----
    The routine tries to go through every node following maximum cardinality
    search. It returns False when it finds that the separator for any node
    is not a clique.  Based on the algorithms in [1]_.

    References
    ----------
    .. [1] R. E. Tarjan and M. Yannakakis, Simple linear-time algorithms
       to test chordality of graphs, test acyclicity of hypergraphs, and
       selectively reduce acyclic hypergraphs, SIAM J. Comput., 13 (1984),
       pp. 566–579.
    """
    if G.is_directed():
        raise nx.NetworkXError('Directed graphs not supported')
    if G.is_multigraph():
        raise nx.NetworkXError('Multiply connected graphs not supported.')
    if len(_find_chordality_breaker(G)) == 0:
        return True
    else:
        return False


def find_induced_nodes(G, s, t, treewidth_bound=sys.maxsize):
    """Returns the set of induced nodes in the path from s to t.

    Parameters
    ----------
    G : graph
      A chordal NetworkX graph
    s : node
        Source node to look for induced nodes
    t : node
        Destination node to look for induced nodes
    treewith_bound: float
        Maximum treewidth acceptable for the graph H. The search
        for induced nodes will end as soon as the treewidth_bound is exceeded.

    Returns
    -------
    I : Set of nodes
        The set of induced nodes in the path from s to t in G

    Raises
    ------
    NetworkXError
        The algorithm does not support DiGraph, MultiGraph and MultiDiGraph.
        If the input graph is an instance of one of these classes, a
        :exc:`NetworkXError` is raised.
        The algorithm can only be applied to chordal graphs. If
        the input graph is found to be non-chordal, a :exc:`NetworkXError` is raised.

    Examples
    --------
    >>> import networkx as nx
    >>> G=nx.Graph()
    >>> G = nx.generators.classic.path_graph(10)
    >>> I = nx.find_induced_nodes(G,1,9,2)
    >>> list(I)
    [1, 2, 3, 4, 5, 6, 7, 8, 9]

    Notes
    -----
    G must be a chordal graph and (s,t) an edge that is not in G.

    If a treewidth_bound is provided, the search for induced nodes will end
    as soon as the treewidth_bound is exceeded.

    The algorithm is inspired by Algorithm 4 in [1]_.
    A formal definition of induced node can also be found on that reference.

    References
    ----------
    .. [1] Learning Bounded Treewidth Bayesian Networks.
       Gal Elidan, Stephen Gould; JMLR, 9(Dec):2699--2731, 2008.
       http://jmlr.csail.mit.edu/papers/volume9/elidan08a/elidan08a.pdf
    """
    if not is_chordal(G):
        raise nx.NetworkXError("Input graph is not chordal.")

    H = nx.Graph(G)
    H.add_edge(s, t)
    I = set()
    triplet = _find_chordality_breaker(H, s, treewidth_bound)
    while triplet:
        (u, v, w) = triplet
        I.update(triplet)
        for n in triplet:
            if n != s:
                H.add_edge(s, n)
        triplet = _find_chordality_breaker(H, s, treewidth_bound)
    if I:
        # Add t and the second node in the induced path from s to t.
        I.add(t)
        for u in G[s]:
            if len(I & set(G[u])) == 2:
                I.add(u)
                break
    return I


def chordal_graph_cliques(G):
    """Returns the set of maximal cliques of a chordal graph.

    The algorithm breaks the graph in connected components and performs a
    maximum cardinality search in each component to get the cliques.

    Parameters
    ----------
    G : graph
      A NetworkX graph

    Returns
    -------
    cliques : A set containing the maximal cliques in G.

    Raises
    ------
    NetworkXError
        The algorithm does not support DiGraph, MultiGraph and MultiDiGraph.
        If the input graph is an instance of one of these classes, a
        :exc:`NetworkXError` is raised.
        The algorithm can only be applied to chordal graphs. If the
        input graph is found to be non-chordal, a :exc:`NetworkXError` is raised.

    Examples
    --------
    >>> import networkx as nx
    >>> e= [(1,2),(1,3),(2,3),(2,4),(3,4),(3,5),(3,6),(4,5),(4,6),(5,6),(7,8)]
    >>> G = nx.Graph(e)
    >>> G.add_node(9)
    >>> setlist = nx.chordal_graph_cliques(G)
    """
    if not is_chordal(G):
        raise nx.NetworkXError("Input graph is not chordal.")

    cliques = set()
    for C in nx.connected.connected_component_subgraphs(G):
        cliques |= _connected_chordal_graph_cliques(C)

    return cliques


def chordal_graph_treewidth(G):
    """Returns the treewidth of the chordal graph G.

    Parameters
    ----------
    G : graph
      A NetworkX graph

    Returns
    -------
    treewidth : int
        The size of the largest clique in the graph minus one.

    Raises
    ------
    NetworkXError
        The algorithm does not support DiGraph, MultiGraph and MultiDiGraph.
        If the input graph is an instance of one of these classes, a
        :exc:`NetworkXError` is raised.
        The algorithm can only be applied to chordal graphs. If
        the input graph is found to be non-chordal, a :exc:`NetworkXError` is raised.

    Examples
    --------
    >>> import networkx as nx
    >>> e = [(1,2),(1,3),(2,3),(2,4),(3,4),(3,5),(3,6),(4,5),(4,6),(5,6),(7,8)]
    >>> G = nx.Graph(e)
    >>> G.add_node(9)
    >>> nx.chordal_graph_treewidth(G)
    3

    References
    ----------
    .. [1] https://en.wikipedia.org/wiki/Tree_decomposition#Treewidth
    """
    if not is_chordal(G):
        raise nx.NetworkXError("Input graph is not chordal.")

    max_clique = -1
    for clique in nx.chordal_graph_cliques(G):
        max_clique = max(max_clique, len(clique))
    return max_clique - 1


def _is_complete_graph(G):
    """Returns True if G is a complete graph."""
    if nx.number_of_selfloops(G) > 0:
        raise nx.NetworkXError("Self loop found in _is_complete_graph()")
    n = G.number_of_nodes()
    if n < 2:
        return True
    e = G.number_of_edges()
    max_edges = ((n * (n - 1)) / 2)
    return e == max_edges


def _find_missing_edge(G):
    """ Given a non-complete graph G, returns a missing edge."""
    nodes = set(G)
    for u in G:
        missing = nodes - set(list(G[u].keys()) + [u])
        if missing:
            return (u, missing.pop())


def _max_cardinality_node(G, choices, wanna_connect):
    """Returns a the node in choices that has more connections in G
    to nodes in wanna_connect.
    """
#    max_number = None
    max_number = -1
    for x in choices:
        number = len([y for y in G[x] if y in wanna_connect])
        if number > max_number:
            max_number = number
            max_cardinality_node = x
    return max_cardinality_node


def _find_chordality_breaker(G, s=None, treewidth_bound=sys.maxsize):
    """ Given a graph G, starts a max cardinality search
    (starting from s if s is given and from a random node otherwise)
    trying to find a non-chordal cycle.

    If it does find one, it returns (u,v,w) where u,v,w are the three
    nodes that together with s are involved in the cycle.
    """

    unnumbered = set(G)
    if s is None:
        s = random.choice(list(unnumbered))
    unnumbered.remove(s)
    numbered = set([s])
#    current_treewidth = None
    current_treewidth = -1
    while unnumbered:  # and current_treewidth <= treewidth_bound:
        v = _max_cardinality_node(G, unnumbered, numbered)
        unnumbered.remove(v)
        numbered.add(v)
        clique_wanna_be = set(G[v]) & numbered
        sg = G.subgraph(clique_wanna_be)
        if _is_complete_graph(sg):
            # The graph seems to be chordal by now. We update the treewidth
            current_treewidth = max(current_treewidth, len(clique_wanna_be))
            if current_treewidth > treewidth_bound:
                raise nx.NetworkXTreewidthBoundExceeded(
                    "treewidth_bound exceeded: %s" % current_treewidth)
        else:
            # sg is not a clique,
            # look for an edge that is not included in sg
            (u, w) = _find_missing_edge(sg)
            return (u, v, w)
    return ()


def _connected_chordal_graph_cliques(G):
    """Return the set of maximal cliques of a connected chordal graph."""
    if G.number_of_nodes() == 1:
        x = frozenset(G.nodes())
        return set([x])
    else:
        cliques = set()
        unnumbered = set(G.nodes())
        v = random.choice(list(unnumbered))
        unnumbered.remove(v)
        numbered = set([v])
        clique_wanna_be = set([v])
        while unnumbered:
            v = _max_cardinality_node(G, unnumbered, numbered)
            unnumbered.remove(v)
            numbered.add(v)
            new_clique_wanna_be = set(G.neighbors(v)) & numbered
            sg = G.subgraph(clique_wanna_be)
            if _is_complete_graph(sg):
                new_clique_wanna_be.add(v)
                if not new_clique_wanna_be >= clique_wanna_be:
                    cliques.add(frozenset(clique_wanna_be))
                clique_wanna_be = new_clique_wanna_be
            else:
                raise nx.NetworkXError("Input graph is not chordal.")
        cliques.add(frozenset(clique_wanna_be))
        return cliques
