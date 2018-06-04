#    Copyright (C) 2015 by
#    Alessandro Luongo
#    BSD license.
#
# Authors:
#    Alessandro Luongo <alessandro.luongo@studenti.unimi.it>
#
"""Functions for computing the harmonic centrality of a graph."""
from __future__ import division
from functools import partial

import networkx as nx

__all__ = ['harmonic_centrality']


def harmonic_centrality(G, nbunch=None, distance=None):
    r"""Compute harmonic centrality for nodes.

    Harmonic centrality [1]_ of a node `u` is the sum of the reciprocal
    of the shortest path distances from all other nodes to `u`

    .. math::

        C(u) = \sum_{v \neq u} \frac{1}{d(v, u)}

    where `d(v, u)` is the shortest-path distance between `v` and `u`.

    Notice that higher values indicate higher centrality.

    Parameters
    ----------
    G : graph
      A NetworkX graph

    nbunch : container
      Container of nodes. If provided harmonic centrality will be computed
      only over the nodes in nbunch.

    distance : edge attribute key, optional (default=None)
      Use the specified edge attribute as the edge distance in shortest
      path calculations.  If `None`, then each edge will have distance equal to 1.

    Returns
    -------
    nodes : dictionary
      Dictionary of nodes with harmonic centrality as the value.

    See Also
    --------
    betweenness_centrality, load_centrality, eigenvector_centrality,
    degree_centrality, closeness_centrality

    Notes
    -----
    If the 'distance' keyword is set to an edge attribute key then the
    shortest-path length will be computed using Dijkstra's algorithm with
    that edge attribute as the edge weight.

    References
    ----------
    .. [1] Boldi, Paolo, and Sebastiano Vigna. "Axioms for centrality."
           Internet Mathematics 10.3-4 (2014): 222-262.
    """
    if G.is_directed():
        G = G.reverse()
    spl = partial(nx.shortest_path_length, G, weight=distance)
    return {u: sum(1 / d if d > 0 else 0 for v, d in spl(source=u).items())
            for u in G.nbunch_iter(nbunch)}
