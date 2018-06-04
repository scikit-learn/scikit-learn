# -*- coding: utf-8 -*-
#    Copyright (C) 2004-2018 by
#    Aric Hagberg <hagberg@lanl.gov>
#    Dan Schult <dschult@colgate.edu>
#    Pieter Swart <swart@lanl.gov>
#    All rights reserved.
#    BSD license.
#
# Authors: Aric Hagberg <aric.hagberg@gmail.com>
#          Sérgio Nery Simões <sergionery@gmail.com>
"""
Compute the shortest paths and path lengths between nodes in the graph.

These algorithms work with undirected and directed graphs.

"""
from __future__ import division

import networkx as nx

__all__ = ['shortest_path', 'all_shortest_paths',
           'shortest_path_length', 'average_shortest_path_length',
           'has_path']


def has_path(G, source, target):
    """Return *True* if *G* has a path from *source* to *target*.

    Parameters
    ----------
    G : NetworkX graph

    source : node
       Starting node for path

    target : node
       Ending node for path
    """
    try:
        sp = nx.shortest_path(G, source, target)
    except nx.NetworkXNoPath:
        return False
    return True


def shortest_path(G, source=None, target=None, weight=None):
    """Compute shortest paths in the graph.

    Parameters
    ----------
    G : NetworkX graph

    source : node, optional
        Starting node for path. If not specified, compute shortest
        paths for each possible starting node.

    target : node, optional
        Ending node for path. If not specified, compute shortest
        paths to all possible nodes.

    weight : None or string, optional (default = None)
        If None, every edge has weight/distance/cost 1.
        If a string, use this edge attribute as the edge weight.
        Any edge attribute not present defaults to 1.

    Returns
    -------
    path: list or dictionary
        All returned paths include both the source and target in the path.

        If the source and target are both specified, return a single list
        of nodes in a shortest path from the source to the target.

        If only the source is specified, return a dictionary keyed by
        targets with a list of nodes in a shortest path from the source
        to one of the targets.

        If only the target is specified, return a dictionary keyed by
        sources with a list of nodes in a shortest path from one of the
        sources to the target.

        If neither the source nor target are specified return a dictionary
        of dictionaries with path[source][target]=[list of nodes in path].

    Examples
    --------
    >>> G = nx.path_graph(5)
    >>> print(nx.shortest_path(G, source=0, target=4))
    [0, 1, 2, 3, 4]
    >>> p = nx.shortest_path(G, source=0) # target not specified
    >>> p[4]
    [0, 1, 2, 3, 4]
    >>> p = nx.shortest_path(G, target=4) # source not specified
    >>> p[0]
    [0, 1, 2, 3, 4]
    >>> p = nx.shortest_path(G) # source, target not specified
    >>> p[0][4]
    [0, 1, 2, 3, 4]

    Notes
    -----
    There may be more than one shortest path between a source and target.
    This returns only one of them.

    See Also
    --------
    all_pairs_shortest_path()
    all_pairs_dijkstra_path()
    single_source_shortest_path()
    single_source_dijkstra_path()
    """
    if source is None:
        if target is None:
            # Find paths between all pairs.
            if weight is None:
                paths = dict(nx.all_pairs_shortest_path(G))
            else:
                paths = dict(nx.all_pairs_dijkstra_path(G, weight=weight))
        else:
            # Find paths from all nodes co-accessible to the target.
            with nx.utils.reversed(G):
                if weight is None:
                    paths = nx.single_source_shortest_path(G, target)
                else:
                    paths = nx.single_source_dijkstra_path(G, target,
                                                           weight=weight)
                # Now flip the paths so they go from a source to the target.
                for target in paths:
                    paths[target] = list(reversed(paths[target]))

    else:
        if target is None:
            # Find paths to all nodes accessible from the source.
            if weight is None:
                paths = nx.single_source_shortest_path(G, source)
            else:
                paths = nx.single_source_dijkstra_path(G, source,
                                                       weight=weight)
        else:
            # Find shortest source-target path.
            if weight is None:
                paths = nx.bidirectional_shortest_path(G, source, target)
            else:
                paths = nx.dijkstra_path(G, source, target, weight)

    return paths


def shortest_path_length(G, source=None, target=None, weight=None):
    """Compute shortest path lengths in the graph.

    Parameters
    ----------
    G : NetworkX graph

    source : node, optional
        Starting node for path.
        If not specified, compute shortest path lengths using all nodes as
        source nodes.

    target : node, optional
        Ending node for path.
        If not specified, compute shortest path lengths using all nodes as
        target nodes.

    weight : None or string, optional (default = None)
        If None, every edge has weight/distance/cost 1.
        If a string, use this edge attribute as the edge weight.
        Any edge attribute not present defaults to 1.

    Returns
    -------
    length: int or iterator
        If the source and target are both specified, return the length of
        the shortest path from the source to the target.

        If only the source is specified, return a dict keyed by target
        to the shortest path length from the source to that target.

        If only the target is specified, return a dict keyed by source
        to the shortest path length from that source to the target.

        If neither the source nor target are specified, return an iterator
        over (source, dictionary) where dictionary is keyed by target to
        shortest path length from source to that target.

    Raises
    ------
    NetworkXNoPath
        If no path exists between source and target.

    Examples
    --------
    >>> G = nx.path_graph(5)
    >>> nx.shortest_path_length(G, source=0, target=4)
    4
    >>> p = nx.shortest_path_length(G, source=0) # target not specified
    >>> p[4]
    4
    >>> p = nx.shortest_path_length(G, target=4) # source not specified
    >>> p[0]
    4
    >>> p = dict(nx.shortest_path_length(G)) # source,target not specified
    >>> p[0][4]
    4

    Notes
    -----
    The length of the path is always 1 less than the number of nodes involved
    in the path since the length measures the number of edges followed.

    For digraphs this returns the shortest directed path length. To find path
    lengths in the reverse direction use G.reverse(copy=False) first to flip
    the edge orientation.

    See Also
    --------
    all_pairs_shortest_path_length()
    all_pairs_dijkstra_path_length()
    single_source_shortest_path_length()
    single_source_dijkstra_path_length()

    """
    if source is None:
        if target is None:
            # Find paths between all pairs.
            if weight is None:
                paths = nx.all_pairs_shortest_path_length(G)
            else:
                paths = nx.all_pairs_dijkstra_path_length(G, weight=weight)
        else:
            # Find paths from all nodes co-accessible to the target.
            with nx.utils.reversed(G):
                if weight is None:
                    # We need to exhaust the iterator as Graph needs
                    # to be reversed.
                    path_length = nx.single_source_shortest_path_length
                    paths = path_length(G, target)
                else:
                    path_length = nx.single_source_dijkstra_path_length
                    paths = path_length(G, target, weight=weight)
    else:
        if source not in G:
            raise nx.NodeNotFound("Source {} not in G".format(source))

        if target is None:
            # Find paths to all nodes accessible from the source.
            if weight is None:
                paths = nx.single_source_shortest_path_length(G, source)
            else:
                paths = nx.single_source_dijkstra_path_length(G, source,
                                                              weight=weight)
        else:
            # Find shortest source-target path.
            if weight is None:
                p = nx.bidirectional_shortest_path(G, source, target)
                paths = len(p) - 1
            else:
                paths = nx.dijkstra_path_length(G, source, target, weight)
    return paths


def average_shortest_path_length(G, weight=None):
    r"""Return the average shortest path length.

    The average shortest path length is

    .. math::

       a =\sum_{s,t \in V} \frac{d(s, t)}{n(n-1)}

    where `V` is the set of nodes in `G`,
    `d(s, t)` is the shortest path from `s` to `t`,
    and `n` is the number of nodes in `G`.

    Parameters
    ----------
    G : NetworkX graph

    weight : None or string, optional (default = None)
       If None, every edge has weight/distance/cost 1.
       If a string, use this edge attribute as the edge weight.
       Any edge attribute not present defaults to 1.

    Raises
    ------
    NetworkXPointlessConcept
        If `G` is the null graph (that is, the graph on zero nodes).

    NetworkXError
        If `G` is not connected (or not weakly connected, in the case
        of a directed graph).

    Examples
    --------
    >>> G = nx.path_graph(5)
    >>> nx.average_shortest_path_length(G)
    2.0

    For disconnected graphs, you can compute the average shortest path
    length for each component

    >>> G = nx.Graph([(1, 2), (3, 4)])
    >>> for C in nx.connected_component_subgraphs(G):
    ...     print(nx.average_shortest_path_length(C))
    1.0
    1.0

    """
    n = len(G)
    # For the special case of the null graph, raise an exception, since
    # there are no paths in the null graph.
    if n == 0:
        msg = ('the null graph has no paths, thus there is no average'
               'shortest path length')
        raise nx.NetworkXPointlessConcept(msg)
    # For the special case of the trivial graph, return zero immediately.
    if n == 1:
        return 0
    # Shortest path length is undefined if the graph is disconnected.
    if G.is_directed() and not nx.is_weakly_connected(G):
        raise nx.NetworkXError("Graph is not weakly connected.")
    if not G.is_directed() and not nx.is_connected(G):
        raise nx.NetworkXError("Graph is not connected.")
    # Compute all-pairs shortest paths.
    if weight is None:
        def path_length(v): return nx.single_source_shortest_path_length(G, v)
    else:
        ssdpl = nx.single_source_dijkstra_path_length

        def path_length(v): return ssdpl(G, v, weight=weight)
    # Sum the distances for each (ordered) pair of source and target node.
    s = sum(l for u in G for l in path_length(u).values())
    return s / (n * (n - 1))


def all_shortest_paths(G, source, target, weight=None):
    """Compute all shortest paths in the graph.

    Parameters
    ----------
    G : NetworkX graph

    source : node
       Starting node for path.

    target : node
       Ending node for path.

    weight : None or string, optional (default = None)
       If None, every edge has weight/distance/cost 1.
       If a string, use this edge attribute as the edge weight.
       Any edge attribute not present defaults to 1.

    Returns
    -------
    paths : generator of lists
        A generator of all paths between source and target.

    Examples
    --------
    >>> G = nx.Graph()
    >>> nx.add_path(G, [0, 1, 2])
    >>> nx.add_path(G, [0, 10, 2])
    >>> print([p for p in nx.all_shortest_paths(G, source=0, target=2)])
    [[0, 1, 2], [0, 10, 2]]

    Notes
    -----
    There may be many shortest paths between the source and target.

    See Also
    --------
    shortest_path()
    single_source_shortest_path()
    all_pairs_shortest_path()
    """
    if weight is not None:
        pred, dist = nx.dijkstra_predecessor_and_distance(G, source,
                                                          weight=weight)
    else:
        pred = nx.predecessor(G, source)

    if source not in G:
        raise nx.NodeNotFound('Source {} is not in G'.format(source))

    if target not in pred:
        raise nx.NetworkXNoPath()

    stack = [[target, 0]]
    top = 0
    while top >= 0:
        node, i = stack[top]
        if node == source:
            yield [p for p, n in reversed(stack[:top + 1])]
        if len(pred[node]) > i:
            top += 1
            if top == len(stack):
                stack.append([pred[node][i], 0])
            else:
                stack[top] = [pred[node][i], 0]
        else:
            stack[top - 1][1] += 1
            top -= 1
