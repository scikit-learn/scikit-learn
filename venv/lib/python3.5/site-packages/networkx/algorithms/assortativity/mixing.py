#-*- coding: utf-8 -*-
"""
Mixing matrices for node attributes and degree.
"""
import networkx as nx
from networkx.utils import dict_to_numpy_array
from networkx.algorithms.assortativity.pairs import node_degree_xy, \
    node_attribute_xy
__author__ = ' '.join(['Aric Hagberg <aric.hagberg@gmail.com>'])
__all__ = ['attribute_mixing_matrix',
           'attribute_mixing_dict',
           'degree_mixing_matrix',
           'degree_mixing_dict',
           'numeric_mixing_matrix',
           'mixing_dict']


def attribute_mixing_dict(G, attribute, nodes=None, normalized=False):
    """Return dictionary representation of mixing matrix for attribute.

    Parameters
    ----------
    G : graph
       NetworkX graph object.

    attribute : string
       Node attribute key.

    nodes: list or iterable (optional)
        Unse nodes in container to build the dict. The default is all nodes.

    normalized : bool (default=False)
       Return counts if False or probabilities if True.

    Examples
    --------
    >>> G=nx.Graph()
    >>> G.add_nodes_from([0,1],color='red')
    >>> G.add_nodes_from([2,3],color='blue')
    >>> G.add_edge(1,3)
    >>> d=nx.attribute_mixing_dict(G,'color')
    >>> print(d['red']['blue'])
    1
    >>> print(d['blue']['red']) # d symmetric for undirected graphs
    1

    Returns
    -------
    d : dictionary
       Counts or joint probability of occurrence of attribute pairs.
    """
    xy_iter = node_attribute_xy(G, attribute, nodes)
    return mixing_dict(xy_iter, normalized=normalized)


def attribute_mixing_matrix(G, attribute, nodes=None, mapping=None,
                            normalized=True):
    """Return mixing matrix for attribute.

    Parameters
    ----------
    G : graph
       NetworkX graph object.

    attribute : string
       Node attribute key.

    nodes: list or iterable (optional)
        Use only nodes in container to build the matrix. The default is
        all nodes.

    mapping : dictionary, optional
       Mapping from node attribute to integer index in matrix.
       If not specified, an arbitrary ordering will be used.

    normalized : bool (default=False)
       Return counts if False or probabilities if True.

    Returns
    -------
    m: numpy array
       Counts or joint probability of occurrence of attribute pairs.
    """
    d = attribute_mixing_dict(G, attribute, nodes)
    a = dict_to_numpy_array(d, mapping=mapping)
    if normalized:
        a = a / a.sum()
    return a


def degree_mixing_dict(G, x='out', y='in', weight=None,
                       nodes=None, normalized=False):
    """Return dictionary representation of mixing matrix for degree.

    Parameters
    ----------
    G : graph
        NetworkX graph object.

    x: string ('in','out')
       The degree type for source node (directed graphs only).

    y: string ('in','out')
       The degree type for target node (directed graphs only).

    weight: string or None, optional (default=None)
       The edge attribute that holds the numerical value used
       as a weight.  If None, then each edge has weight 1.
       The degree is the sum of the edge weights adjacent to the node.

    normalized : bool (default=False)
        Return counts if False or probabilities if True.

    Returns
    -------
    d: dictionary
       Counts or joint probability of occurrence of degree pairs.
    """
    xy_iter = node_degree_xy(G, x=x, y=y, nodes=nodes, weight=weight)
    return mixing_dict(xy_iter, normalized=normalized)


def degree_mixing_matrix(G, x='out', y='in', weight=None,
                         nodes=None, normalized=True):
    """Return mixing matrix for attribute.

    Parameters
    ----------
    G : graph
       NetworkX graph object.

    x: string ('in','out')
       The degree type for source node (directed graphs only).

    y: string ('in','out')
       The degree type for target node (directed graphs only).

    nodes: list or iterable (optional)
        Build the matrix using only nodes in container.
        The default is all nodes.

    weight: string or None, optional (default=None)
       The edge attribute that holds the numerical value used
       as a weight.  If None, then each edge has weight 1.
       The degree is the sum of the edge weights adjacent to the node.

    normalized : bool (default=False)
       Return counts if False or probabilities if True.

    Returns
    -------
    m: numpy array
       Counts, or joint probability, of occurrence of node degree.
    """
    d = degree_mixing_dict(G, x=x, y=y, nodes=nodes, weight=weight)
    s = set(d.keys())
    for k, v in d.items():
        s.update(v.keys())
    m = max(s)
    mapping = {x: x for x in range(m + 1)}
    a = dict_to_numpy_array(d, mapping=mapping)
    if normalized:
        a = a / a.sum()
    return a


def numeric_mixing_matrix(G, attribute, nodes=None, normalized=True):
    """Return numeric mixing matrix for attribute.

    The attribute must be an integer.

    Parameters
    ----------
    G : graph
       NetworkX graph object.

    attribute : string
       Node attribute key.  The corresponding attribute must be an integer.

    nodes: list or iterable (optional)
        Build the matrix only with nodes in container. The default is all nodes.

    normalized : bool (default=False)
       Return counts if False or probabilities if True.

    Returns
    -------
    m: numpy array
       Counts, or joint, probability of occurrence of node attribute pairs.
    """
    d = attribute_mixing_dict(G, attribute, nodes)
    s = set(d.keys())
    for k, v in d.items():
        s.update(v.keys())
    m = max(s)
    mapping = {x: x for x in range(m + 1)}
    a = dict_to_numpy_array(d, mapping=mapping)
    if normalized:
        a = a / a.sum()
    return a


def mixing_dict(xy, normalized=False):
    """Return a dictionary representation of mixing matrix.

    Parameters
    ----------
    xy : list or container of two-tuples
       Pairs of (x,y) items.

    attribute : string
       Node attribute key

    normalized : bool (default=False)
       Return counts if False or probabilities if True.

    Returns
    -------
    d: dictionary
       Counts or Joint probability of occurrence of values in xy.
    """
    d = {}
    psum = 0.0
    for x, y in xy:
        if x not in d:
            d[x] = {}
        if y not in d:
            d[y] = {}
        v = d[x].get(y, 0)
        d[x][y] = v + 1
        psum += 1

    if normalized:
        for k, jdict in d.items():
            for j in jdict:
                jdict[j] /= psum
    return d


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
