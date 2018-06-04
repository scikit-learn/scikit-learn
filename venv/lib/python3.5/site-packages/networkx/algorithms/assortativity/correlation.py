#-*- coding: utf-8 -*-
"""Node assortativity coefficients and correlation measures.
"""
import networkx as nx
from networkx.algorithms.assortativity.mixing import degree_mixing_matrix, \
    attribute_mixing_matrix, numeric_mixing_matrix
from networkx.algorithms.assortativity.pairs import node_degree_xy, \
    node_attribute_xy
__author__ = ' '.join(['Aric Hagberg <aric.hagberg@gmail.com>',
                       'Oleguer Sagarra <oleguer.sagarra@gmail.com>'])
__all__ = ['degree_pearson_correlation_coefficient',
           'degree_assortativity_coefficient',
           'attribute_assortativity_coefficient',
           'numeric_assortativity_coefficient']


def degree_assortativity_coefficient(G, x='out', y='in', weight=None,
                                     nodes=None):
    """Compute degree assortativity of graph.

    Assortativity measures the similarity of connections
    in the graph with respect to the node degree.

    Parameters
    ----------
    G : NetworkX graph

    x: string ('in','out')
       The degree type for source node (directed graphs only).

    y: string ('in','out')
       The degree type for target node (directed graphs only).

    weight: string or None, optional (default=None)
       The edge attribute that holds the numerical value used 
       as a weight.  If None, then each edge has weight 1.
       The degree is the sum of the edge weights adjacent to the node.

    nodes: list or iterable (optional)
        Compute degree assortativity only for nodes in container. 
        The default is all nodes.

    Returns
    -------
    r : float
       Assortativity of graph by degree.

    Examples
    --------
    >>> G=nx.path_graph(4)
    >>> r=nx.degree_assortativity_coefficient(G)
    >>> print("%3.1f"%r)
    -0.5

    See Also
    --------
    attribute_assortativity_coefficient
    numeric_assortativity_coefficient
    neighbor_connectivity
    degree_mixing_dict
    degree_mixing_matrix

    Notes
    -----
    This computes Eq. (21) in Ref. [1]_ , where e is the joint
    probability distribution (mixing matrix) of the degrees.  If G is
    directed than the matrix e is the joint probability of the 
    user-specified degree type for the source and target.

    References
    ----------
    .. [1] M. E. J. Newman, Mixing patterns in networks,
       Physical Review E, 67 026126, 2003
    .. [2] Foster, J.G., Foster, D.V., Grassberger, P. & Paczuski, M. 
       Edge direction and the structure of networks, PNAS 107, 10815-20 (2010).
    """
    M = degree_mixing_matrix(G, x=x, y=y, nodes=nodes, weight=weight)
    return numeric_ac(M)


def degree_pearson_correlation_coefficient(G, x='out', y='in',
                                           weight=None, nodes=None):
    """Compute degree assortativity of graph. 

    Assortativity measures the similarity of connections
    in the graph with respect to the node degree.

    This is the same as degree_assortativity_coefficient but uses the 
    potentially faster scipy.stats.pearsonr function.

    Parameters
    ----------
    G : NetworkX graph

    x: string ('in','out')
       The degree type for source node (directed graphs only).

    y: string ('in','out')
       The degree type for target node (directed graphs only).

    weight: string or None, optional (default=None)
       The edge attribute that holds the numerical value used 
       as a weight.  If None, then each edge has weight 1.
       The degree is the sum of the edge weights adjacent to the node.

    nodes: list or iterable (optional)
        Compute pearson correlation of degrees only for specified nodes.
        The default is all nodes.

    Returns
    -------
    r : float
       Assortativity of graph by degree.

    Examples
    --------
    >>> G=nx.path_graph(4)
    >>> r=nx.degree_pearson_correlation_coefficient(G) 
    >>> print("%3.1f"%r)
    -0.5

    Notes
    -----
    This calls scipy.stats.pearsonr.

    References
    ----------
    .. [1] M. E. J. Newman, Mixing patterns in networks
           Physical Review E, 67 026126, 2003
    .. [2] Foster, J.G., Foster, D.V., Grassberger, P. & Paczuski, M. 
       Edge direction and the structure of networks, PNAS 107, 10815-20 (2010).
    """
    try:
        import scipy.stats as stats
    except ImportError:
        raise ImportError(
            "Assortativity requires SciPy: http://scipy.org/ ")
    xy = node_degree_xy(G, x=x, y=y, nodes=nodes, weight=weight)
    x, y = zip(*xy)
    return stats.pearsonr(x, y)[0]


def attribute_assortativity_coefficient(G, attribute, nodes=None):
    """Compute assortativity for node attributes.

    Assortativity measures the similarity of connections
    in the graph with respect to the given attribute.

    Parameters
    ----------
    G : NetworkX graph

    attribute : string 
        Node attribute key

    nodes: list or iterable (optional)
        Compute attribute assortativity for nodes in container. 
        The default is all nodes. 

    Returns
    -------
    r: float
       Assortativity of graph for given attribute

    Examples
    --------
    >>> G=nx.Graph()
    >>> G.add_nodes_from([0,1],color='red')
    >>> G.add_nodes_from([2,3],color='blue')
    >>> G.add_edges_from([(0,1),(2,3)])
    >>> print(nx.attribute_assortativity_coefficient(G,'color'))
    1.0

    Notes
    -----
    This computes Eq. (2) in Ref. [1]_ , trace(M)-sum(M))/(1-sum(M),
    where M is the joint probability distribution (mixing matrix)
    of the specified attribute.

    References
    ----------
    .. [1] M. E. J. Newman, Mixing patterns in networks,
       Physical Review E, 67 026126, 2003
    """
    M = attribute_mixing_matrix(G, attribute, nodes)
    return attribute_ac(M)


def numeric_assortativity_coefficient(G, attribute, nodes=None):
    """Compute assortativity for numerical node attributes.

    Assortativity measures the similarity of connections
    in the graph with respect to the given numeric attribute.
    The numeric attribute must be an integer.

    Parameters
    ----------
    G : NetworkX graph

    attribute : string 
        Node attribute key.  The corresponding attribute value must be an
        integer.

    nodes: list or iterable (optional)
        Compute numeric assortativity only for attributes of nodes in 
        container. The default is all nodes.

    Returns
    -------
    r: float
       Assortativity of graph for given attribute

    Examples
    --------
    >>> G=nx.Graph()
    >>> G.add_nodes_from([0,1],size=2)
    >>> G.add_nodes_from([2,3],size=3)
    >>> G.add_edges_from([(0,1),(2,3)])
    >>> print(nx.numeric_assortativity_coefficient(G,'size'))
    1.0

    Notes
    -----
    This computes Eq. (21) in Ref. [1]_ , for the mixing matrix of 
    of the specified attribute.

    References
    ----------
    .. [1] M. E. J. Newman, Mixing patterns in networks
           Physical Review E, 67 026126, 2003
    """
    a = numeric_mixing_matrix(G, attribute, nodes)
    return numeric_ac(a)


def attribute_ac(M):
    """Compute assortativity for attribute matrix M.

    Parameters
    ----------
    M : numpy array or matrix
        Attribute mixing matrix.

    Notes
    -----
    This computes Eq. (2) in Ref. [1]_ , (trace(e)-sum(e))/(1-sum(e)),
    where e is the joint probability distribution (mixing matrix)
    of the specified attribute.

    References
    ----------
    .. [1] M. E. J. Newman, Mixing patterns in networks,
       Physical Review E, 67 026126, 2003
    """
    try:
        import numpy
    except ImportError:
        raise ImportError(
            "attribute_assortativity requires NumPy: http://scipy.org/ ")
    if M.sum() != 1.0:
        M = M / float(M.sum())
    M = numpy.asmatrix(M)
    s = (M * M).sum()
    t = M.trace()
    r = (t - s) / (1 - s)
    return float(r)


def numeric_ac(M):
    # M is a numpy matrix or array
    # numeric assortativity coefficient, pearsonr
    try:
        import numpy
    except ImportError:
        raise ImportError('numeric_assortativity requires ',
                          'NumPy: http://scipy.org/')
    if M.sum() != 1.0:
        M = M / float(M.sum())
    nx, ny = M.shape  # nx=ny
    x = numpy.arange(nx)
    y = numpy.arange(ny)
    a = M.sum(axis=0)
    b = M.sum(axis=1)
    vara = (a * x**2).sum() - ((a * x).sum())**2
    varb = (b * x**2).sum() - ((b * x).sum())**2
    xy = numpy.outer(x, y)
    ab = numpy.outer(a, b)
    return (xy * (M - ab)).sum() / numpy.sqrt(vara * varb)


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
