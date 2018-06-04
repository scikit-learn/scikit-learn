# -*- coding: utf-8 -*-
"""
Communicability.
"""
#    Copyright (C) 2011 by
#    Aric Hagberg <hagberg@lanl.gov>
#    Dan Schult <dschult@colgate.edu>
#    Pieter Swart <swart@lanl.gov>
#    Previously coded as communicability centrality
#    All rights reserved.
#    BSD license.
import networkx as nx
from networkx.utils import *
__author__ = "\n".join(['Aric Hagberg (hagberg@lanl.gov)',
                        'Franck Kalala (franckkalala@yahoo.fr'])
__all__ = ['communicability',
           'communicability_exp',
           ]


@not_implemented_for('directed')
@not_implemented_for('multigraph')
def communicability(G):
    r"""Return communicability between all pairs of nodes in G.

    The communicability between pairs of nodes in G is the sum of
    closed walks of different lengths starting at node u and ending at node v.

    Parameters
    ----------
    G: graph

    Returns
    -------
    comm: dictionary of dictionaries
        Dictionary of dictionaries keyed by nodes with communicability
        as the value.

    Raises
    ------
    NetworkXError
       If the graph is not undirected and simple.

    See Also
    --------
    communicability_exp:
       Communicability between all pairs of nodes in G  using spectral
       decomposition.
    communicability_betweenness_centrality:
       Communicability betweeness centrality for each node in G.

    Notes
    -----
    This algorithm uses a spectral decomposition of the adjacency matrix.
    Let G=(V,E) be a simple undirected graph.  Using the connection between
    the powers  of the adjacency matrix and the number of walks in the graph,
    the communicability  between nodes `u` and `v` based on the graph spectrum
    is [1]_

    .. math::
        C(u,v)=\sum_{j=1}^{n}\phi_{j}(u)\phi_{j}(v)e^{\lambda_{j}},

    where `\phi_{j}(u)` is the `u\rm{th}` element of the `j\rm{th}` orthonormal
    eigenvector of the adjacency matrix associated with the eigenvalue
    `\lambda_{j}`.

    References
    ----------
    .. [1] Ernesto Estrada, Naomichi Hatano,
       "Communicability in complex networks",
       Phys. Rev. E 77, 036111 (2008).
       https://arxiv.org/abs/0707.0756

    Examples
    --------
    >>> G = nx.Graph([(0,1),(1,2),(1,5),(5,4),(2,4),(2,3),(4,3),(3,6)])
    >>> c = nx.communicability(G)
    """
    import numpy
    import scipy.linalg
    nodelist = list(G)  # ordering of nodes in matrix
    A = nx.to_numpy_matrix(G, nodelist)
    # convert to 0-1 matrix
    A[A != 0.0] = 1
    w, vec = numpy.linalg.eigh(A)
    expw = numpy.exp(w)
    mapping = dict(zip(nodelist, range(len(nodelist))))
    c = {}
    # computing communicabilities
    for u in G:
        c[u] = {}
        for v in G:
            s = 0
            p = mapping[u]
            q = mapping[v]
            for j in range(len(nodelist)):
                s += vec[:, j][p, 0] * vec[:, j][q, 0] * expw[j]
            c[u][v] = float(s)
    return c


@not_implemented_for('directed')
@not_implemented_for('multigraph')
def communicability_exp(G):
    r"""Return communicability between all pairs of nodes in G.

    Communicability between pair of node (u,v) of node in G is the sum of
    closed walks of different lengths starting at node u and ending at node v.

    Parameters
    ----------
    G: graph

    Returns
    -------
    comm: dictionary of dictionaries
        Dictionary of dictionaries keyed by nodes with communicability
        as the value.

    Raises
    ------
    NetworkXError
        If the graph is not undirected and simple.

    See Also
    --------
    communicability:
       Communicability between pairs of nodes in G.
    communicability_betweenness_centrality:
       Communicability betweeness centrality for each node in G.

    Notes
    -----
    This algorithm uses matrix exponentiation of the adjacency matrix.

    Let G=(V,E) be a simple undirected graph.  Using the connection between
    the powers  of the adjacency matrix and the number of walks in the graph,
    the communicability between nodes u and v is [1]_,

    .. math::
        C(u,v) = (e^A)_{uv},

    where `A` is the adjacency matrix of G.

    References
    ----------
    .. [1] Ernesto Estrada, Naomichi Hatano,
       "Communicability in complex networks",
       Phys. Rev. E 77, 036111 (2008).
       https://arxiv.org/abs/0707.0756

    Examples
    --------
    >>> G = nx.Graph([(0,1),(1,2),(1,5),(5,4),(2,4),(2,3),(4,3),(3,6)])
    >>> c = nx.communicability_exp(G)
    """
    import scipy.linalg
    nodelist = list(G)  # ordering of nodes in matrix
    A = nx.to_numpy_matrix(G, nodelist)
    # convert to 0-1 matrix
    A[A != 0.0] = 1
    # communicability matrix
    expA = scipy.linalg.expm(A.A)
    mapping = dict(zip(nodelist, range(len(nodelist))))
    c = {}
    for u in G:
        c[u] = {}
        for v in G:
            c[u][v] = float(expA[mapping[u], mapping[v]])
    return c

# fixture for nose tests


def setup_module(module):
    from nose import SkipTest
    try:
        import numpy
    except:
        raise SkipTest("NumPy not available")
    try:
        import scipy
    except:
        raise SkipTest("SciPy not available")
