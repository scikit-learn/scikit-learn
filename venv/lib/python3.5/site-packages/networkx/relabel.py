#    Copyright (C) 2006-2018 by
#    Aric Hagberg <hagberg@lanl.gov>
#    Dan Schult <dschult@colgate.edu>
#    Pieter Swart <swart@lanl.gov>
#    All rights reserved.
#    BSD license.
import networkx as nx

__all__ = ['convert_node_labels_to_integers', 'relabel_nodes']


def relabel_nodes(G, mapping, copy=True):
    """Relabel the nodes of the graph G.

    Parameters
    ----------
    G : graph
       A NetworkX graph

    mapping : dictionary
       A dictionary with the old labels as keys and new labels as values.
       A partial mapping is allowed.

    copy : bool (optional, default=True)
       If True return a copy, or if False relabel the nodes in place.

    Examples
    --------
    To create a new graph with nodes relabeled according to a given
    dictionary:

    >>> G = nx.path_graph(3)
    >>> sorted(G)
    [0, 1, 2]
    >>> mapping = {0: 'a', 1: 'b', 2: 'c'}
    >>> H = nx.relabel_nodes(G, mapping)
    >>> sorted(H)
    ['a', 'b', 'c']

    Nodes can be relabeled with any hashable object, including numbers
    and strings:

    >>> import string
    >>> G = nx.path_graph(26)  # nodes are integers 0 through 25
    >>> sorted(G)[:3]
    [0, 1, 2]
    >>> mapping = dict(zip(G, string.ascii_lowercase))
    >>> G = nx.relabel_nodes(G, mapping) # nodes are characters a through z
    >>> sorted(G)[:3]
    ['a', 'b', 'c']
    >>> mapping = dict(zip(G, range(1, 27)))
    >>> G = nx.relabel_nodes(G, mapping)  # nodes are integers 1 through 26
    >>> sorted(G)[:3]
    [1, 2, 3]

    To perform a partial in-place relabeling, provide a dictionary
    mapping only a subset of the nodes, and set the `copy` keyword
    argument to False:

    >>> G = nx.path_graph(3)  # nodes 0-1-2
    >>> mapping = {0: 'a', 1: 'b'} # 0->'a' and 1->'b'
    >>> G = nx.relabel_nodes(G, mapping, copy=False)
    >>> sorted(G, key=str)
    [2, 'a', 'b']

    A mapping can also be given as a function:

    >>> G = nx.path_graph(3)
    >>> H = nx.relabel_nodes(G, lambda x: x ** 2)
    >>> list(H)
    [0, 1, 4]

    Notes
    -----
    Only the nodes specified in the mapping will be relabeled.

    The keyword setting copy=False modifies the graph in place.
    Relabel_nodes avoids naming collisions by building a
    directed graph from ``mapping`` which specifies the order of
    relabelings. Naming collisions, such as a->b, b->c, are ordered
    such that "b" gets renamed to "c" before "a" gets renamed "b".
    In cases of circular mappings (e.g. a->b, b->a), modifying the
    graph is not possible in-place and an exception is raised.
    In that case, use copy=True.

    See Also
    --------
    convert_node_labels_to_integers
    """
    # you can pass a function f(old_label)->new_label
    # but we'll just make a dictionary here regardless
    if not hasattr(mapping, "__getitem__"):
        m = {n: mapping(n) for n in G}
    else:
        m = mapping
    if copy:
        return _relabel_copy(G, m)
    else:
        return _relabel_inplace(G, m)


def _relabel_inplace(G, mapping):
    old_labels = set(mapping.keys())
    new_labels = set(mapping.values())
    if len(old_labels & new_labels) > 0:
        # labels sets overlap
        # can we topological sort and still do the relabeling?
        D = nx.DiGraph(list(mapping.items()))
        D.remove_edges_from(nx.selfloop_edges(D))
        try:
            nodes = reversed(list(nx.topological_sort(D)))
        except nx.NetworkXUnfeasible:
            raise nx.NetworkXUnfeasible('The node label sets are overlapping '
                                        'and no ordering can resolve the '
                                        'mapping. Use copy=True.')
    else:
        # non-overlapping label sets
        nodes = old_labels

    multigraph = G.is_multigraph()
    directed = G.is_directed()

    for old in nodes:
        try:
            new = mapping[old]
        except KeyError:
            continue
        if new == old:
            continue
        try:
            G.add_node(new, **G.nodes[old])
        except KeyError:
            raise KeyError("Node %s is not in the graph" % old)
        if multigraph:
            new_edges = [(new, new if old == target else target, key, data)
                         for (_, target, key, data)
                         in G.edges(old, data=True, keys=True)]
            if directed:
                new_edges += [(new if old == source else source, new, key, data)
                              for (source, _, key, data)
                              in G.in_edges(old, data=True, keys=True)]
        else:
            new_edges = [(new, new if old == target else target, data)
                         for (_, target, data) in G.edges(old, data=True)]
            if directed:
                new_edges += [(new if old == source else source, new, data)
                              for (source, _, data) in G.in_edges(old, data=True)]
        G.remove_node(old)
        G.add_edges_from(new_edges)
    return G


def _relabel_copy(G, mapping):
    H = G.fresh_copy()
    H.add_nodes_from(mapping.get(n, n) for n in G)
    H._node.update((mapping.get(n, n), d.copy()) for n, d in G.nodes.items())
    if G.is_multigraph():
        H.add_edges_from((mapping.get(n1, n1), mapping.get(n2, n2), k, d.copy())
                         for (n1, n2, k, d) in G.edges(keys=True, data=True))
    else:
        H.add_edges_from((mapping.get(n1, n1), mapping.get(n2, n2), d.copy())
                         for (n1, n2, d) in G.edges(data=True))
    H.graph.update(G.graph)
    return H


def convert_node_labels_to_integers(G, first_label=0, ordering="default",
                                    label_attribute=None):
    """Return a copy of the graph G with the nodes relabeled using
    consecutive integers.

    Parameters
    ----------
    G : graph
       A NetworkX graph

    first_label : int, optional (default=0)
       An integer specifying the starting offset in numbering nodes.
       The new integer labels are numbered first_label, ..., n-1+first_label.

    ordering : string
       "default" : inherit node ordering from G.nodes()
       "sorted"  : inherit node ordering from sorted(G.nodes())
       "increasing degree" : nodes are sorted by increasing degree
       "decreasing degree" : nodes are sorted by decreasing degree

    label_attribute : string, optional (default=None)
       Name of node attribute to store old label.  If None no attribute
       is created.

    Notes
    -----
    Node and edge attribute data are copied to the new (relabeled) graph.

    See Also
    --------
    relabel_nodes
    """
    N = G.number_of_nodes() + first_label
    if ordering == "default":
        mapping = dict(zip(G.nodes(), range(first_label, N)))
    elif ordering == "sorted":
        nlist = sorted(G.nodes())
        mapping = dict(zip(nlist, range(first_label, N)))
    elif ordering == "increasing degree":
        dv_pairs = [(d, n) for (n, d) in G.degree()]
        dv_pairs.sort()  # in-place sort from lowest to highest degree
        mapping = dict(zip([n for d, n in dv_pairs], range(first_label, N)))
    elif ordering == "decreasing degree":
        dv_pairs = [(d, n) for (n, d) in G.degree()]
        dv_pairs.sort()  # in-place sort from lowest to highest degree
        dv_pairs.reverse()
        mapping = dict(zip([n for d, n in dv_pairs], range(first_label, N)))
    else:
        raise nx.NetworkXError('Unknown node ordering: %s' % ordering)
    H = relabel_nodes(G, mapping)
    # create node attribute with the old label
    if label_attribute is not None:
        nx.set_node_attributes(H, {v: k for k, v in mapping.items()},
                               label_attribute)
    return H
