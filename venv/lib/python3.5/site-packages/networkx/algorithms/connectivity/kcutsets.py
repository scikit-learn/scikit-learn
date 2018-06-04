# -*- coding: utf-8 -*-
"""
Kanevsky all minimum node k cutsets algorithm.
"""
from operator import itemgetter
from itertools import combinations

import networkx as nx
from .utils import build_auxiliary_node_connectivity
from networkx.algorithms.flow import (
    build_residual_network,
    edmonds_karp,
    shortest_augmenting_path,
)
default_flow_func = edmonds_karp


__author__ = '\n'.join(['Jordi Torrents <jtorrents@milnou.net>'])

__all__ = ['all_node_cuts']


def all_node_cuts(G, k=None, flow_func=None):
    r"""Returns all minimum k cutsets of an undirected graph G. 

    This implementation is based on Kanevsky's algorithm [1]_ for finding all
    minimum-size node cut-sets of an undirected graph G; ie the set (or sets) 
    of nodes of cardinality equal to the node connectivity of G. Thus if 
    removed, would break G into two or more connected components.

    Parameters
    ----------
    G : NetworkX graph
        Undirected graph

    k : Integer
        Node connectivity of the input graph. If k is None, then it is 
        computed. Default value: None.

    flow_func : function
        Function to perform the underlying flow computations. Default value
        edmonds_karp. This function performs better in sparse graphs with
        right tailed degree distributions. shortest_augmenting_path will
        perform better in denser graphs.


    Returns
    -------
    cuts : a generator of node cutsets
        Each node cutset has cardinality equal to the node connectivity of
        the input graph.

    Examples
    --------
    >>> # A two-dimensional grid graph has 4 cutsets of cardinality 2
    >>> G = nx.grid_2d_graph(5, 5)
    >>> cutsets = list(nx.all_node_cuts(G))
    >>> len(cutsets)
    4
    >>> all(2 == len(cutset) for cutset in cutsets)
    True
    >>> nx.node_connectivity(G)
    2

    Notes
    -----
    This implementation is based on the sequential algorithm for finding all
    minimum-size separating vertex sets in a graph [1]_. The main idea is to
    compute minimum cuts using local maximum flow computations among a set 
    of nodes of highest degree and all other non-adjacent nodes in the Graph.
    Once we find a minimum cut, we add an edge between the high degree
    node and the target node of the local maximum flow computation to make 
    sure that we will not find that minimum cut again.

    See also
    --------
    node_connectivity
    edmonds_karp
    shortest_augmenting_path

    References
    ----------
    .. [1]  Kanevsky, A. (1993). Finding all minimum-size separating vertex 
            sets in a graph. Networks 23(6), 533--541.
            http://onlinelibrary.wiley.com/doi/10.1002/net.3230230604/abstract

    """
    if not nx.is_connected(G):
        raise nx.NetworkXError('Input graph is disconnected.')

    # Addess some corner cases first.
    # For cycle graphs
    if G.order() == G.size():
        if all(2 == d for n, d in G.degree()):
            seen = set()
            for u in G:
                for v in nx.non_neighbors(G, u):
                    if (u, v) not in seen and (v, u) not in seen:
                        yield {v, u}
                        seen.add((v, u))
            return
    # For complete Graphs
    if nx.density(G) == 1:
        for cut_set in combinations(G, len(G) - 1):
            yield set(cut_set)
        return
    # Initialize data structures.
    # Keep track of the cuts already computed so we do not repeat them.
    seen = []
    # Even-Tarjan reduction is what we call auxiliary digraph
    # for node connectivity.
    H = build_auxiliary_node_connectivity(G)
    mapping = H.graph['mapping']
    R = build_residual_network(H, 'capacity')
    kwargs = dict(capacity='capacity', residual=R)
    # Define default flow function
    if flow_func is None:
        flow_func = default_flow_func
    if flow_func is shortest_augmenting_path:
        kwargs['two_phase'] = True
    # Begin the actual algorithm
    # step 1: Find node connectivity k of G
    if k is None:
        k = nx.node_connectivity(G, flow_func=flow_func)
    # step 2:
    # Find k nodes with top degree, call it X:
    X = {n for n, d in sorted(G.degree(), key=itemgetter(1), reverse=True)[:k]}
    # Check if X is a k-node-cutset
    if _is_separating_set(G, X):
        seen.append(X)
        yield X

    for x in X:
        # step 3: Compute local connectivity flow of x with all other
        # non adjacent nodes in G
        non_adjacent = set(G) - X - set(G[x])
        for v in non_adjacent:
            # step 4: compute maximum flow in an Even-Tarjan reduction H of G
            # and step:5 build the associated residual network R
            R = flow_func(H, '%sB' % mapping[x], '%sA' % mapping[v], **kwargs)
            flow_value = R.graph['flow_value']

            if flow_value == k:
                # Remove saturated edges form the residual network
                saturated_edges = [(u, w, d) for (u, w, d) in
                                   R.edges(data=True)
                                   if d['capacity'] == d['flow']]
                R.remove_edges_from(saturated_edges)
                # step 6: shrink the strongly connected components of
                # residual flow network R and call it L
                L = nx.condensation(R)
                cmap = L.graph['mapping']
                # step 7: Compute antichains of L; they map to closed sets in H
                # Any edge in H that links a closed set is part of a cutset
                for antichain in nx.antichains(L):
                    # Nodes in an antichain of the condensation graph of
                    # the residual network map to a closed set of nodes that
                    # define a node partition of the auxiliary digraph H.
                    S = {n for n, scc in cmap.items() if scc in antichain}
                    # Find the cutset that links the node partition (S,~S) in H
                    cutset = set()
                    for u in S:
                        cutset.update((u, w) for w in H[u] if w not in S)
                    # The edges in H that form the cutset are internal edges
                    # (ie edges that represent a node of the original graph G)
                    node_cut = {H.nodes[n]['id'] for edge in cutset for n in edge}

                    if len(node_cut) == k:
                        if node_cut not in seen:
                            yield node_cut
                            seen.append(node_cut)
                        # Add an edge (x, v) to make sure that we do not
                        # find this cutset again. This is equivalent
                        # of adding the edge in the input graph
                        # G.add_edge(x, v) and then regenerate H and R:
                        # Add edges to the auxiliary digraph.
                        H.add_edge('%sB' % mapping[x], '%sA' % mapping[v],
                                   capacity=1)
                        H.add_edge('%sB' % mapping[v], '%sA' % mapping[x],
                                   capacity=1)
                        # Add edges to the residual network.
                        R.add_edge('%sB' % mapping[x], '%sA' % mapping[v],
                                   capacity=1)
                        R.add_edge('%sA' % mapping[v], '%sB' % mapping[x],
                                   capacity=1)
                        break
                # Add again the saturated edges to reuse the residual network
                R.add_edges_from(saturated_edges)


def _is_separating_set(G, cut):
    """Assumes that the input graph is connected"""
    if len(cut) == len(G) - 1:
        return True

    H = nx.restricted_view(G, cut, [])
    if nx.is_connected(H):
        return False
    return True
