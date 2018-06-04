#    Copyright (C) 2011 by
#    Aric Hagberg <hagberg@lanl.gov>
#    Dan Schult <dschult@colgate.edu>
#    Pieter Swart <swart@lanl.gov>
#    All rights reserved.
#    BSD license.
from itertools import chain
import networkx as nx
from networkx.utils import make_str
__author__ = """Aric Hagberg (hagberg@lanl.gov))"""
__all__ = ['tree_data', 'tree_graph']

_attrs = dict(id='id', children='children')


def tree_data(G, root, attrs=_attrs):
    """Return data in tree format that is suitable for JSON serialization
    and use in Javascript documents.

    Parameters
    ----------
    G : NetworkX graph
       G must be an oriented tree

    root : node
       The root of the tree

    attrs : dict
        A dictionary that contains two keys 'id' and 'children'. The
        corresponding values provide the attribute names for storing
        NetworkX-internal graph data. The values should be unique. Default
        value: :samp:`dict(id='id', children='children')`.

        If some user-defined graph data use these attribute names as data keys,
        they may be silently dropped.

    Returns
    -------
    data : dict
       A dictionary with node-link formatted data.

    Raises
    ------
    NetworkXError
        If values in attrs are not unique.

    Examples
    --------
    >>> from networkx.readwrite import json_graph
    >>> G = nx.DiGraph([(1,2)])
    >>> data = json_graph.tree_data(G,root=1)

    To serialize with json

    >>> import json
    >>> s = json.dumps(data)

    Notes
    -----
    Node attributes are stored in this format but keys
    for attributes must be strings if you want to serialize with JSON.

    Graph and edge attributes are not stored.

    The default value of attrs will be changed in a future release of NetworkX.

    See Also
    --------
    tree_graph, node_link_data, node_link_data
    """
    if G.number_of_nodes() != G.number_of_edges() + 1:
        raise TypeError("G is not a tree.")
    if not G.is_directed():
        raise TypeError("G is not directed.")

    id_ = attrs['id']
    children = attrs['children']
    if id_ == children:
        raise nx.NetworkXError('Attribute names are not unique.')

    def add_children(n, G):
        nbrs = G[n]
        if len(nbrs) == 0:
            return []
        children_ = []
        for child in nbrs:
            d = dict(chain(G.nodes[child].items(), [(id_, child)]))
            c = add_children(child, G)
            if c:
                d[children] = c
            children_.append(d)
        return children_

    data = dict(chain(G.nodes[root].items(), [(id_, root)]))
    data[children] = add_children(root, G)
    return data


def tree_graph(data, attrs=_attrs):
    """Return graph from tree data format.

    Parameters
    ----------
    data : dict
        Tree formatted graph data

    Returns
    -------
    G : NetworkX DiGraph

    attrs : dict
        A dictionary that contains two keys 'id' and 'children'. The
        corresponding values provide the attribute names for storing
        NetworkX-internal graph data. The values should be unique. Default
        value: :samp:`dict(id='id', children='children')`.

    Examples
    --------
    >>> from networkx.readwrite import json_graph
    >>> G = nx.DiGraph([(1,2)])
    >>> data = json_graph.tree_data(G,root=1)
    >>> H = json_graph.tree_graph(data)

    Notes
    -----
    The default value of attrs will be changed in a future release of NetworkX.

    See Also
    --------
    tree_graph, node_link_data, adjacency_data
    """
    graph = nx.DiGraph()
    id_ = attrs['id']
    children = attrs['children']

    def add_children(parent, children_):
        for data in children_:
            child = data[id_]
            graph.add_edge(parent, child)
            grandchildren = data.get(children, [])
            if grandchildren:
                add_children(child, grandchildren)
            nodedata = dict((make_str(k), v) for k, v in data.items()
                            if k != id_ and k != children)
            graph.add_node(child, **nodedata)

    root = data[id_]
    children_ = data.get(children, [])
    nodedata = dict((make_str(k), v) for k, v in data.items()
                    if k != id_ and k != children)
    graph.add_node(root, **nodedata)
    add_children(root, children_)
    return graph
