"""
Eigenvalue spectrum of graphs.
"""
#    Copyright (C) 2004-2018 by
#    Aric Hagberg <hagberg@lanl.gov>
#    Dan Schult <dschult@colgate.edu>
#    Pieter Swart <swart@lanl.gov>
#    All rights reserved.
#    BSD license.
import networkx as nx
__author__ = "\n".join(['Aric Hagberg <aric.hagberg@gmail.com>',
                        'Pieter Swart (swart@lanl.gov)',
                        'Dan Schult(dschult@colgate.edu)',
                        'Jean-Gabriel Young (jean.gabriel.young@gmail.com)'])

__all__ = ['laplacian_spectrum', 'adjacency_spectrum', 'modularity_spectrum']


def laplacian_spectrum(G, weight='weight'):
    """Return eigenvalues of the Laplacian of G

    Parameters
    ----------
    G : graph
       A NetworkX graph

    weight : string or None, optional (default='weight')
       The edge data key used to compute each value in the matrix.
       If None, then each edge has weight 1.

    Returns
    -------
    evals : NumPy array
      Eigenvalues

    Notes
    -----
    For MultiGraph/MultiDiGraph, the edges weights are summed.
    See to_numpy_matrix for other options.

    See Also
    --------
    laplacian_matrix
    """
    from scipy.linalg import eigvalsh
    return eigvalsh(nx.laplacian_matrix(G, weight=weight).todense())


def adjacency_spectrum(G, weight='weight'):
    """Return eigenvalues of the adjacency matrix of G.

    Parameters
    ----------
    G : graph
       A NetworkX graph

    weight : string or None, optional (default='weight')
       The edge data key used to compute each value in the matrix.
       If None, then each edge has weight 1.

    Returns
    -------
    evals : NumPy array
      Eigenvalues

    Notes
    -----
    For MultiGraph/MultiDiGraph, the edges weights are summed.
    See to_numpy_matrix for other options.

    See Also
    --------
    adjacency_matrix
    """
    from scipy.linalg import eigvals
    return eigvals(nx.adjacency_matrix(G, weight=weight).todense())


def modularity_spectrum(G):
    """Return eigenvalues of the modularity matrix of G.

    Parameters
    ----------
    G : Graph
       A NetworkX Graph or DiGraph

    Returns
    -------
    evals : NumPy array
      Eigenvalues

    See Also
    --------
    modularity_matrix

    References
    ----------
    .. [1] M. E. J. Newman, "Modularity and community structure in networks",
       Proc. Natl. Acad. Sci. USA, vol. 103, pp. 8577-8582, 2006.
    """
    from scipy.linalg import eigvals
    if G.is_directed():
        return eigvals(nx.directed_modularity_matrix(G))
    else:
        return eigvals(nx.modularity_matrix(G))

# fixture for nose tests


def setup_module(module):
    from nose import SkipTest
    try:
        import scipy.linalg
    except:
        raise SkipTest("scipy.linalg not available")
