# -*- coding: utf-8 -*-
"""
Spectral bipartivity measure.
"""
import networkx as nx
__author__ = """Aric Hagberg (hagberg@lanl.gov)"""
#    Copyright (C) 2011 by
#    Aric Hagberg <hagberg@lanl.gov>
#    Dan Schult <dschult@colgate.edu>
#    Pieter Swart <swart@lanl.gov>
#    All rights reserved.
#    BSD license.
__all__ = ['spectral_bipartivity']


def spectral_bipartivity(G, nodes=None, weight='weight'):
    """Returns the spectral bipartivity.

    Parameters
    ----------
    G : NetworkX graph 

    nodes : list or container  optional(default is all nodes)
      Nodes to return value of spectral bipartivity contribution.

    weight : string or None  optional (default = 'weight')
      Edge data key to use for edge weights. If None, weights set to 1.

    Returns
    -------
    sb : float or dict
       A single number if the keyword nodes is not specified, or
       a dictionary keyed by node with the spectral bipartivity contribution
       of that node as the value.

    Examples
    --------
    >>> from networkx.algorithms import bipartite
    >>> G = nx.path_graph(4)
    >>> bipartite.spectral_bipartivity(G)
    1.0

    Notes
    -----
    This implementation uses Numpy (dense) matrices which are not efficient
    for storing large sparse graphs.  

    See Also
    --------
    color

    References
    ----------
    .. [1] E. Estrada and J. A. Rodríguez-Velázquez, "Spectral measures of
       bipartivity in complex networks", PhysRev E 72, 046105 (2005)
    """
    try:
        import scipy.linalg
    except ImportError:
        raise ImportError('spectral_bipartivity() requires SciPy: ',
                          'http://scipy.org/')
    nodelist = list(G)  # ordering of nodes in matrix
    A = nx.to_numpy_matrix(G, nodelist, weight=weight)
    expA = scipy.linalg.expm(A)
    expmA = scipy.linalg.expm(-A)
    coshA = 0.5 * (expA + expmA)
    if nodes is None:
        # return single number for entire graph
        return coshA.diagonal().sum() / expA.diagonal().sum()
    else:
        # contribution for individual nodes
        index = dict(zip(nodelist, range(len(nodelist))))
        sb = {}
        for n in nodes:
            i = index[n]
            sb[n] = coshA[i, i] / expA[i, i]
        return sb


def setup_module(module):
    """Fixture for nose tests."""
    from nose import SkipTest
    try:
        import numpy
    except:
        raise SkipTest("NumPy not available")
    try:
        import scipy
    except:
        raise SkipTest("SciPy not available")
