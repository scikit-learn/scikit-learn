# -*- coding: utf-8 -*-
#    Copyright (C) 2004-2018 by
#    Aric Hagberg <hagberg@lanl.gov>
#    Dan Schult <dschult@colgate.edu>
#    Pieter Swart <swart@lanl.gov>
#    All rights reserved.
#    BSD license.
#
# Authors: Christopher Ellison
"""Attracting components."""
import warnings as _warnings
import networkx as nx
from networkx.utils.decorators import not_implemented_for

__all__ = ['number_attracting_components',
           'attracting_components',
           'is_attracting_component',
           'attracting_component_subgraphs',
           ]


@not_implemented_for('undirected')
def attracting_components(G):
    """Generates the attracting components in `G`.

    An attracting component in a directed graph `G` is a strongly connected
    component with the property that a random walker on the graph will never
    leave the component, once it enters the component.

    The nodes in attracting components can also be thought of as recurrent
    nodes.  If a random walker enters the attractor containing the node, then
    the node will be visited infinitely often.

    Parameters
    ----------
    G : DiGraph, MultiDiGraph
        The graph to be analyzed.

    Returns
    -------
    attractors : generator of sets
        A generator of sets of nodes, one for each attracting component of G.

    Raises
    ------
    NetworkXNotImplemented :
        If the input graph is undirected.

    See Also
    --------
    number_attracting_components
    is_attracting_component

    """
    scc = list(nx.strongly_connected_components(G))
    cG = nx.condensation(G, scc)
    for n in cG:
        if cG.out_degree(n) == 0:
            yield scc[n]


@not_implemented_for('undirected')
def number_attracting_components(G):
    """Returns the number of attracting components in `G`.

    Parameters
    ----------
    G : DiGraph, MultiDiGraph
        The graph to be analyzed.

    Returns
    -------
    n : int
        The number of attracting components in G.

    Raises
    ------
    NetworkXNotImplemented :
        If the input graph is undirected.

    See Also
    --------
    attracting_components
    is_attracting_component

    """
    return sum(1 for ac in attracting_components(G))


@not_implemented_for('undirected')
def is_attracting_component(G):
    """Returns True if `G` consists of a single attracting component.

    Parameters
    ----------
    G : DiGraph, MultiDiGraph
        The graph to be analyzed.

    Returns
    -------
    attracting : bool
        True if `G` has a single attracting component. Otherwise, False.

    Raises
    ------
    NetworkXNotImplemented :
        If the input graph is undirected.

    See Also
    --------
    attracting_components
    number_attracting_components

    """
    ac = list(attracting_components(G))
    if len(ac) == 1:
        return len(ac[0]) == len(G)
    return False


@not_implemented_for('undirected')
def attracting_component_subgraphs(G, copy=True):
    """DEPRECATED: Use ``(G.subgraph(c) for c in attracting_components(G))``

           Or ``(G.subgraph(c).copy() for c in attracting_components(G))``
    """
    msg = "attracting_component_subgraphs is deprecated and will be removed" \
        "in 2.2. Use (G.subgraph(c).copy() for c in attracting_components(G))"
    _warnings.warn(msg, DeprecationWarning)
    for c in attracting_components(G):
        if copy:
            yield G.subgraph(c).copy()
        else:
            yield G.subgraph(c)
