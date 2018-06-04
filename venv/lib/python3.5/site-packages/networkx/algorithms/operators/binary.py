"""
Operations on graphs including union, intersection, difference.
"""
#    Copyright (C) 2004-2018 by
#    Aric Hagberg <hagberg@lanl.gov>
#    Dan Schult <dschult@colgate.edu>
#    Pieter Swart <swart@lanl.gov>
#    All rights reserved.
#    BSD license.
import networkx as nx
from networkx.utils import is_string_like
__author__ = """\n""".join(['Aric Hagberg <aric.hagberg@gmail.com>',
                            'Pieter Swart (swart@lanl.gov)',
                            'Dan Schult(dschult@colgate.edu)'])
__all__ = ['union', 'compose', 'disjoint_union', 'intersection',
           'difference', 'symmetric_difference']


def union(G, H, rename=(None, None), name=None):
    """ Return the union of graphs G and H.

    Graphs G and H must be disjoint, otherwise an exception is raised.

    Parameters
    ----------
    G,H : graph
       A NetworkX graph

    rename : bool , default=(None, None)
       Node names of G and H can be changed by specifying the tuple
       rename=('G-','H-') (for example).  Node "u" in G is then renamed
       "G-u" and "v" in H is renamed "H-v".

    name : string
       Specify the name for the union graph

    Returns
    -------
    U : A union graph with the same type as G.

    Notes
    -----
    To force a disjoint union with node relabeling, use
    disjoint_union(G,H) or convert_node_labels_to integers().

    Graph, edge, and node attributes are propagated from G and H
    to the union graph.  If a graph attribute is present in both
    G and H the value from H is used.

    See Also
    --------
    disjoint_union
    """
    if not G.is_multigraph() == H.is_multigraph():
        raise nx.NetworkXError('G and H must both be graphs or multigraphs.')
    # Union is the same type as G
    R = G.fresh_copy()
    # add graph attributes, H attributes take precedent over G attributes
    R.graph.update(G.graph)
    R.graph.update(H.graph)

    # rename graph to obtain disjoint node labels
    def add_prefix(graph, prefix):
        if prefix is None:
            return graph

        def label(x):
            if is_string_like(x):
                name = prefix + x
            else:
                name = prefix + repr(x)
            return name
        return nx.relabel_nodes(graph, label)
    G = add_prefix(G, rename[0])
    H = add_prefix(H, rename[1])
    if set(G) & set(H):
        raise nx.NetworkXError('The node sets of G and H are not disjoint.',
                               'Use appropriate rename=(Gprefix,Hprefix)'
                               'or use disjoint_union(G,H).')
    if G.is_multigraph():
        G_edges = G.edges(keys=True, data=True)
    else:
        G_edges = G.edges(data=True)
    if H.is_multigraph():
        H_edges = H.edges(keys=True, data=True)
    else:
        H_edges = H.edges(data=True)

    # add nodes
    R.add_nodes_from(G)
    R.add_edges_from(G_edges)
    # add edges
    R.add_nodes_from(H)
    R.add_edges_from(H_edges)
    # add node attributes
    for n in G:
        R.nodes[n].update(G.nodes[n])
    for n in H:
        R.nodes[n].update(H.nodes[n])

    return R


def disjoint_union(G, H):
    """ Return the disjoint union of graphs G and H.

    This algorithm forces distinct integer node labels.

    Parameters
    ----------
    G,H : graph
       A NetworkX graph

    Returns
    -------
    U : A union graph with the same type as G.

    Notes
    -----
    A new graph is created, of the same class as G.  It is recommended
    that G and H be either both directed or both undirected.

    The nodes of G are relabeled 0 to len(G)-1, and the nodes of H are
    relabeled len(G) to len(G)+len(H)-1.

    Graph, edge, and node attributes are propagated from G and H
    to the union graph.  If a graph attribute is present in both
    G and H the value from H is used.
    """
    R1 = nx.convert_node_labels_to_integers(G)
    R2 = nx.convert_node_labels_to_integers(H, first_label=len(R1))
    R = union(R1, R2)
    R.graph.update(G.graph)
    R.graph.update(H.graph)
    return R


def intersection(G, H):
    """Return a new graph that contains only the edges that exist in
    both G and H.

    The node sets of H and G must be the same.

    Parameters
    ----------
    G,H : graph
       A NetworkX graph.  G and H must have the same node sets.

    Returns
    -------
    GH : A new graph with the same type as G.

    Notes
    -----
    Attributes from the graph, nodes, and edges are not copied to the new
    graph.  If you want a new graph of the intersection of G and H
    with the attributes (including edge data) from G use remove_nodes_from()
    as follows

    >>> G=nx.path_graph(3)
    >>> H=nx.path_graph(5)
    >>> R=G.copy()
    >>> R.remove_nodes_from(n for n in G if n not in H)
    """
    # create new graph
    R = nx.create_empty_copy(G)

    if not G.is_multigraph() == H.is_multigraph():
        raise nx.NetworkXError('G and H must both be graphs or multigraphs.')
    if set(G) != set(H):
        raise nx.NetworkXError("Node sets of graphs are not equal")

    if G.number_of_edges() <= H.number_of_edges():
        if G.is_multigraph():
            edges = G.edges(keys=True)
        else:
            edges = G.edges()
        for e in edges:
            if H.has_edge(*e):
                R.add_edge(*e)
    else:
        if H.is_multigraph():
            edges = H.edges(keys=True)
        else:
            edges = H.edges()
        for e in edges:
            if G.has_edge(*e):
                R.add_edge(*e)
    return R


def difference(G, H):
    """Return a new graph that contains the edges that exist in G but not in H.

    The node sets of H and G must be the same.

    Parameters
    ----------
    G,H : graph
       A NetworkX graph.  G and H must have the same node sets.

    Returns
    -------
    D : A new graph with the same type as G.

    Notes
    -----
    Attributes from the graph, nodes, and edges are not copied to the new
    graph.  If you want a new graph of the difference of G and H with
    with the attributes (including edge data) from G use remove_nodes_from()
    as follows:

    >>> G = nx.path_graph(3)
    >>> H = nx.path_graph(5)
    >>> R = G.copy()
    >>> R.remove_nodes_from(n for n in G if n in H)
    """
    # create new graph
    if not G.is_multigraph() == H.is_multigraph():
        raise nx.NetworkXError('G and H must both be graphs or multigraphs.')
    R = nx.create_empty_copy(G)

    if set(G) != set(H):
        raise nx.NetworkXError("Node sets of graphs not equal")

    if G.is_multigraph():
        edges = G.edges(keys=True)
    else:
        edges = G.edges()
    for e in edges:
        if not H.has_edge(*e):
            R.add_edge(*e)
    return R


def symmetric_difference(G, H):
    """Return new graph with edges that exist in either G or H but not both.

    The node sets of H and G must be the same.

    Parameters
    ----------
    G,H : graph
       A NetworkX graph.  G and H must have the same node sets.

    Returns
    -------
    D : A new graph with the same type as G.

    Notes
    -----
    Attributes from the graph, nodes, and edges are not copied to the new
    graph.
    """
    # create new graph
    if not G.is_multigraph() == H.is_multigraph():
        raise nx.NetworkXError('G and H must both be graphs or multigraphs.')
    R = nx.create_empty_copy(G)

    if set(G) != set(H):
        raise nx.NetworkXError("Node sets of graphs not equal")

    gnodes = set(G)  # set of nodes in G
    hnodes = set(H)  # set of nodes in H
    nodes = gnodes.symmetric_difference(hnodes)
    R.add_nodes_from(nodes)

    if G.is_multigraph():
        edges = G.edges(keys=True)
    else:
        edges = G.edges()
    # we could copy the data here but then this function doesn't
    # match intersection and difference
    for e in edges:
        if not H.has_edge(*e):
            R.add_edge(*e)

    if H.is_multigraph():
        edges = H.edges(keys=True)
    else:
        edges = H.edges()
    for e in edges:
        if not G.has_edge(*e):
            R.add_edge(*e)
    return R


def compose(G, H):
    """Return a new graph of G composed with H.

    Composition is the simple union of the node sets and edge sets.
    The node sets of G and H do not need to be disjoint.

    Parameters
    ----------
    G, H : graph
       A NetworkX graph

    Returns
    -------
    C: A new graph  with the same type as G

    Notes
    -----
    It is recommended that G and H be either both directed or both undirected.
    Attributes from H take precedent over attributes from G.

    For MultiGraphs, the edges are identified by incident nodes AND edge-key.
    This can cause surprises (i.e., edge `(1, 2)` may or may not be the same
    in two graphs) if you use MultiGraph without keeping track of edge keys.
    """
    if not G.is_multigraph() == H.is_multigraph():
        raise nx.NetworkXError('G and H must both be graphs or multigraphs.')

    R = G.fresh_copy()
    # add graph attributes, H attributes take precedent over G attributes
    R.graph.update(G.graph)
    R.graph.update(H.graph)

    R.add_nodes_from(G.nodes(data=True))
    R.add_nodes_from(H.nodes(data=True))

    if G.is_multigraph():
        R.add_edges_from(G.edges(keys=True, data=True))
    else:
        R.add_edges_from(G.edges(data=True))
    if H.is_multigraph():
        R.add_edges_from(H.edges(keys=True, data=True))
    else:
        R.add_edges_from(H.edges(data=True))
    return R
