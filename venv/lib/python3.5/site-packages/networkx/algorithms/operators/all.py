"""Operations on many graphs.
"""
#    Copyright (C) 2013 by
#    Aric Hagberg <hagberg@lanl.gov>
#    Dan Schult <dschult@colgate.edu>
#    Pieter Swart <swart@lanl.gov>
#    All rights reserved.
#    BSD license.
try:
    from itertools import izip_longest as zip_longest
except ImportError:  # Python3 has zip_longest
    from itertools import zip_longest
import networkx as nx

__author__ = """\n""".join(['Robert King <kingrobertking@gmail.com>',
                            'Aric Hagberg <aric.hagberg@gmail.com>'])

__all__ = ['union_all', 'compose_all', 'disjoint_union_all',
           'intersection_all']


def union_all(graphs, rename=(None,)):
    """Return the union of all graphs.

    The graphs must be disjoint, otherwise an exception is raised.

    Parameters
    ----------
    graphs : list of graphs
       List of NetworkX graphs

    rename : bool , default=(None, None)
       Node names of G and H can be changed by specifying the tuple
       rename=('G-','H-') (for example).  Node "u" in G is then renamed
       "G-u" and "v" in H is renamed "H-v".

    Returns
    -------
    U : a graph with the same type as the first graph in list

    Notes
    -----
    To force a disjoint union with node relabeling, use
    disjoint_union_all(G,H) or convert_node_labels_to integers().

    Graph, edge, and node attributes are propagated to the union graph.
    If a graph attribute is present in multiple graphs, then the value
    from the last graph in the list with that attribute is used.

    See Also
    --------
    union
    disjoint_union_all
    """
    graphs_names = zip_longest(graphs, rename)
    U, gname = next(graphs_names)
    for H, hname in graphs_names:
        U = nx.union(U, H, (gname, hname))
        gname = None
    return U


def disjoint_union_all(graphs):
    """Return the disjoint union of all graphs.

    This operation forces distinct integer node labels starting with 0
    for the first graph in the list and numbering consecutively.

    Parameters
    ----------
    graphs : list
       List of NetworkX graphs

    Returns
    -------
    U : A graph with the same type as the first graph in list

    Notes
    -----
    It is recommended that the graphs be either all directed or all undirected.

    Graph, edge, and node attributes are propagated to the union graph.
    If a graph attribute is present in multiple graphs, then the value
    from the last graph in the list with that attribute is used.
    """
    graphs = iter(graphs)
    U = next(graphs)
    for H in graphs:
        U = nx.disjoint_union(U, H)
    return U


def compose_all(graphs):
    """Return the composition of all graphs.

    Composition is the simple union of the node sets and edge sets.
    The node sets of the supplied graphs need not be disjoint.

    Parameters
    ----------
    graphs : list
       List of NetworkX graphs

    Returns
    -------
    C : A graph with the same type as the first graph in list

    Notes
    -----
    It is recommended that the supplied graphs be either all directed or all
    undirected.

    Graph, edge, and node attributes are propagated to the union graph.
    If a graph attribute is present in multiple graphs, then the value
    from the last graph in the list with that attribute is used.
    """
    graphs = iter(graphs)
    C = next(graphs)
    for H in graphs:
        C = nx.compose(C, H)
    return C


def intersection_all(graphs):
    """Return a new graph that contains only the edges that exist in
    all graphs.

    All supplied graphs must have the same node set.

    Parameters
    ----------
    graphs_list : list
       List of NetworkX graphs

    Returns
    -------
    R : A new graph with the same type as the first graph in list

    Notes
    -----
    Attributes from the graph, nodes, and edges are not copied to the new
    graph.
    """
    graphs = iter(graphs)
    R = next(graphs)
    for H in graphs:
        R = nx.intersection(R, H)
    return R
