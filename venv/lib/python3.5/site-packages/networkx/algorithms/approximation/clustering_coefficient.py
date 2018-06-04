# -*- coding: utf-8 -*-
#   Copyright (C) 2013 by
#   Fred Morstatter <fred.morstatter@asu.edu>
#   Jordi Torrents <jtorrents@milnou.net>
#   All rights reserved.
#   BSD license.
import random
from networkx.utils import not_implemented_for

__all__ = ['average_clustering']
__author__ = """\n""".join(['Fred Morstatter <fred.morstatter@asu.edu>',
                            'Jordi Torrents <jtorrents@milnou.net>'])


@not_implemented_for('directed')
def average_clustering(G, trials=1000):
    r"""Estimates the average clustering coefficient of G.

    The local clustering of each node in `G` is the fraction of triangles
    that actually exist over all possible triangles in its neighborhood.
    The average clustering coefficient of a graph `G` is the mean of
    local clusterings.

    This function finds an approximate average clustering coefficient
    for G by repeating `n` times (defined in `trials`) the following
    experiment: choose a node at random, choose two of its neighbors
    at random, and check if they are connected. The approximate
    coefficient is the fraction of triangles found over the number
    of trials [1]_.

    Parameters
    ----------
    G : NetworkX graph

    trials : integer
        Number of trials to perform (default 1000).

    Returns
    -------
    c : float
        Approximated average clustering coefficient.

    References
    ----------
    .. [1] Schank, Thomas, and Dorothea Wagner. Approximating clustering
       coefficient and transitivity. Universität Karlsruhe, Fakultät für
       Informatik, 2004.
       http://www.emis.ams.org/journals/JGAA/accepted/2005/SchankWagner2005.9.2.pdf

    """
    n = len(G)
    triangles = 0
    nodes = list(G)
    for i in [int(random.random() * n) for i in range(trials)]:
        nbrs = list(G[nodes[i]])
        if len(nbrs) < 2:
            continue
        u, v = random.sample(nbrs, 2)
        if u in G[v]:
            triangles += 1
    return triangles / float(trials)
