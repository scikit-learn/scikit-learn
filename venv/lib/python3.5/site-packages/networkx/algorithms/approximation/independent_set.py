# -*- coding: utf-8 -*-
"""
Independent Set

Independent set or stable set is a set of vertices in a graph, no two of
which are adjacent. That is, it is a set I of vertices such that for every
two vertices in I, there is no edge connecting the two. Equivalently, each
edge in the graph has at most one endpoint in I. The size of an independent
set is the number of vertices it contains.

A maximum independent set is a largest independent set for a given graph G
and its size is denoted α(G). The problem of finding such a set is called
the maximum independent set problem and is an NP-hard optimization problem.
As such, it is unlikely that there exists an efficient algorithm for finding
a maximum independent set of a graph.

`Wikipedia: Independent set <https://en.wikipedia.org/wiki/Independent_set_(graph_theory)>`_

Independent set algorithm is based on the following paper:

$O(|V|/(log|V|)^2)$ apx of maximum clique/independent set.

Boppana, R., & Halldórsson, M. M. (1992).
Approximating maximum independent sets by excluding subgraphs.
BIT Numerical Mathematics, 32(2), 180–196. Springer.
doi:10.1007/BF01994876

"""
#   Copyright (C) 2011-2012 by
#   Nicholas Mancuso <nick.mancuso@gmail.com>
#   All rights reserved.
#   BSD license.
from networkx.algorithms.approximation import clique_removal
__all__ = ["maximum_independent_set"]
__author__ = """Nicholas Mancuso (nick.mancuso@gmail.com)"""


def maximum_independent_set(G):
    """Return an approximate maximum independent set.

    Parameters
    ----------
    G : NetworkX graph
        Undirected graph

    Returns
    -------
    iset : Set
        The apx-maximum independent set

    Notes
    -----
    Finds the $O(|V|/(log|V|)^2)$ apx of independent set in the worst case.


    References
    ----------
    .. [1] Boppana, R., & Halldórsson, M. M. (1992).
       Approximating maximum independent sets by excluding subgraphs.
       BIT Numerical Mathematics, 32(2), 180–196. Springer.
    """
    iset, _ = clique_removal(G)
    return iset
