# -*- coding: utf-8 -*-
#    Copyright (C) 2004-2018 by
#    Aric Hagberg <hagberg@lanl.gov>
#    Dan Schult <dschult@colgate.edu>
#    Pieter Swart <swart@lanl.gov>
#    All rights reserved.
#    BSD license.
#
# Authors: ysitu (ysitu@users.noreply.github.com)
"""Semiconnectedness."""
import networkx as nx
from networkx.utils import not_implemented_for, pairwise

__all__ = ['is_semiconnected']


@not_implemented_for('undirected')
def is_semiconnected(G):
    """Return True if the graph is semiconnected, False otherwise.

    A graph is semiconnected if, and only if, for any pair of nodes, either one
    is reachable from the other, or they are mutually reachable.

    Parameters
    ----------
    G : NetworkX graph
        A directed graph.

    Returns
    -------
    semiconnected : bool
        True if the graph is semiconnected, False otherwise.

    Raises
    ------
    NetworkXNotImplemented :
        If the input graph is undirected.

    NetworkXPointlessConcept :
        If the graph is empty.

    Examples
    --------
    >>> G=nx.path_graph(4,create_using=nx.DiGraph())
    >>> print(nx.is_semiconnected(G))
    True
    >>> G=nx.DiGraph([(1, 2), (3, 2)])
    >>> print(nx.is_semiconnected(G))
    False

    See Also
    --------
    is_strongly_connected
    is_weakly_connected
    is_connected
    is_biconnected
    """
    if len(G) == 0:
        raise nx.NetworkXPointlessConcept(
            'Connectivity is undefined for the null graph.')

    if not nx.is_weakly_connected(G):
        return False

    G = nx.condensation(G)
    path = nx.topological_sort(G)
    return all(G.has_edge(u, v) for u, v in pairwise(path))
