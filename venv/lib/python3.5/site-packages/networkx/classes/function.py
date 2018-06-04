#    Copyright (C) 2004-2018 by
#    Aric Hagberg <hagberg@lanl.gov>
#    Dan Schult <dschult@colgate.edu>
#    Pieter Swart <swart@lanl.gov>
#    All rights reserved.
#    BSD license.
#
# Authors: Aric Hagberg <hagberg@lanl.gov>
#          Pieter Swart <swart@lanl.gov>
#          Dan Schult <dschult@colgate.edu>
"""Functional interface to graph methods and assorted utilities.
"""
from __future__ import division

from collections import Counter
from itertools import chain
try:
    from itertools import zip_longest
except ImportError:
    from itertools import izip_longest as zip_longest

import networkx as nx
from networkx.utils import pairwise, not_implemented_for

__all__ = ['nodes', 'edges', 'degree', 'degree_histogram', 'neighbors',
           'number_of_nodes', 'number_of_edges', 'density',
           'is_directed', 'info', 'freeze', 'is_frozen', 'subgraph',
           'induced_subgraph', 'edge_subgraph', 'restricted_view',
           'reverse_view', 'to_directed', 'to_undirected',
           'add_star', 'add_path', 'add_cycle',
           'create_empty_copy', 'set_node_attributes',
           'get_node_attributes', 'set_edge_attributes',
           'get_edge_attributes', 'all_neighbors', 'non_neighbors',
           'non_edges', 'common_neighbors', 'is_weighted',
           'is_negatively_weighted', 'is_empty',
           'selfloop_edges', 'nodes_with_selfloops', 'number_of_selfloops',
           ]


def nodes(G):
    """Return an iterator over the graph nodes."""
    return G.nodes()


def edges(G, nbunch=None):
    """Return an edge view of edges incident to nodes in nbunch.

    Return all edges if nbunch is unspecified or nbunch=None.

    For digraphs, edges=out_edges
    """
    return G.edges(nbunch)


def degree(G, nbunch=None, weight=None):
    """Return a degree view of single node or of nbunch of nodes.
    If nbunch is ommitted, then return degrees of *all* nodes.
    """
    return G.degree(nbunch, weight)


def neighbors(G, n):
    """Return a list of nodes connected to node n. """
    return G.neighbors(n)


def number_of_nodes(G):
    """Return the number of nodes in the graph."""
    return G.number_of_nodes()


def number_of_edges(G):
    """Return the number of edges in the graph. """
    return G.number_of_edges()


def density(G):
    r"""Return the density of a graph.

    The density for undirected graphs is

    .. math::

       d = \frac{2m}{n(n-1)},

    and for directed graphs is

    .. math::

       d = \frac{m}{n(n-1)},

    where `n` is the number of nodes and `m`  is the number of edges in `G`.

    Notes
    -----
    The density is 0 for a graph without edges and 1 for a complete graph.
    The density of multigraphs can be higher than 1.

    Self loops are counted in the total number of edges so graphs with self
    loops can have density higher than 1.
    """
    n = number_of_nodes(G)
    m = number_of_edges(G)
    if m == 0 or n <= 1:
        return 0
    d = m / (n * (n - 1))
    if not G.is_directed():
        d *= 2
    return d


def degree_histogram(G):
    """Return a list of the frequency of each degree value.

    Parameters
    ----------
    G : Networkx graph
       A graph

    Returns
    -------
    hist : list
       A list of frequencies of degrees.
       The degree values are the index in the list.

    Notes
    -----
    Note: the bins are width one, hence len(list) can be large
    (Order(number_of_edges))
    """
    counts = Counter(d for n, d in G.degree())
    return [counts.get(i, 0) for i in range(max(counts) + 1)]


def is_directed(G):
    """ Return True if graph is directed."""
    return G.is_directed()


def frozen(*args):
    """Dummy method for raising errors when trying to modify frozen graphs"""
    raise nx.NetworkXError("Frozen graph can't be modified")


def freeze(G):
    """Modify graph to prevent further change by adding or removing
    nodes or edges.

    Node and edge data can still be modified.

    Parameters
    ----------
    G : graph
      A NetworkX graph

    Examples
    --------
    >>> G = nx.path_graph(4)
    >>> G = nx.freeze(G)
    >>> try:
    ...    G.add_edge(4, 5)
    ... except nx.NetworkXError as e:
    ...    print(str(e))
    Frozen graph can't be modified

    Notes
    -----
    To "unfreeze" a graph you must make a copy by creating a new graph object:

    >>> graph = nx.path_graph(4)
    >>> frozen_graph = nx.freeze(graph)
    >>> unfrozen_graph = nx.Graph(frozen_graph)
    >>> nx.is_frozen(unfrozen_graph)
    False

    See Also
    --------
    is_frozen
    """
    G.add_node = frozen
    G.add_nodes_from = frozen
    G.remove_node = frozen
    G.remove_nodes_from = frozen
    G.add_edge = frozen
    G.add_edges_from = frozen
    G.remove_edge = frozen
    G.remove_edges_from = frozen
    G.clear = frozen
    G.frozen = True
    return G


def is_frozen(G):
    """Return True if graph is frozen.

    Parameters
    ----------
    G : graph
      A NetworkX graph

    See Also
    --------
    freeze
    """
    try:
        return G.frozen
    except AttributeError:
        return False


def add_star(G_to_add_to, nodes_for_star, **attr):
    """Add a star to Graph G_to_add_to.

    The first node in `nodes_for_star` is the middle of the star.
    It is connected to all other nodes.

    Parameters
    ----------
    nodes_for_star : iterable container
        A container of nodes.
    attr : keyword arguments, optional (default= no attributes)
        Attributes to add to every edge in star.

    See Also
    --------
    add_path, add_cycle

    Examples
    --------
    >>> G = nx.Graph()
    >>> nx.add_star(G, [0, 1, 2, 3])
    >>> nx.add_star(G, [10, 11, 12], weight=2)
    """
    nlist = iter(nodes_for_star)
    v = next(nlist)
    edges = ((v, n) for n in nlist)
    G_to_add_to.add_edges_from(edges, **attr)


def add_path(G_to_add_to, nodes_for_path, **attr):
    """Add a path to the Graph G_to_add_to.

    Parameters
    ----------
    nodes_for_path : iterable container
        A container of nodes.  A path will be constructed from
        the nodes (in order) and added to the graph.
    attr : keyword arguments, optional (default= no attributes)
        Attributes to add to every edge in path.

    See Also
    --------
    add_star, add_cycle

    Examples
    --------
    >>> G = nx.Graph()
    >>> nx.add_path(G, [0, 1, 2, 3])
    >>> nx.add_path(G, [10, 11, 12], weight=7)
    """
    nlist = iter(nodes_for_path)
    try:
        first_node = next(nlist)
    except StopIteration:
        return
    G_to_add_to.add_node(first_node)
    G_to_add_to.add_edges_from(pairwise(chain((first_node,), nlist)), **attr)


def add_cycle(G_to_add_to, nodes_for_cycle, **attr):
    """Add a cycle to the Graph G_to_add_to.

    Parameters
    ----------
    nodes_for_cycle: iterable container
        A container of nodes.  A cycle will be constructed from
        the nodes (in order) and added to the graph.
    attr : keyword arguments, optional (default= no attributes)
        Attributes to add to every edge in cycle.

    See Also
    --------
    add_path, add_star

    Examples
    --------
    >>> G = nx.Graph()   # or DiGraph, MultiGraph, MultiDiGraph, etc
    >>> nx.add_cycle(G, [0, 1, 2, 3])
    >>> nx.add_cycle(G, [10, 11, 12], weight=7)
    """
    G_to_add_to.add_edges_from(pairwise(nodes_for_cycle, cyclic=True), **attr)


def subgraph(G, nbunch):
    """Return the subgraph induced on nodes in nbunch.

    Parameters
    ----------
    G : graph
       A NetworkX graph

    nbunch : list, iterable
       A container of nodes that will be iterated through once (thus
       it should be an iterator or be iterable).  Each element of the
       container should be a valid node type: any hashable type except
       None.  If nbunch is None, return all edges data in the graph.
       Nodes in nbunch that are not in the graph will be (quietly)
       ignored.

    Notes
    -----
    subgraph(G) calls G.subgraph()
    """
    return G.subgraph(nbunch)


def induced_subgraph(G, nbunch):
    """Return a SubGraph view of `G` showing only nodes in nbunch.

    The induced subgraph of a graph on a set of nodes N is the
    graph with nodes N and edges from G which have both ends in N.

    Parameters
    ----------
    G : NetworkX Graph
    nbunch : node, container of nodes or None (for all nodes)

    Returns
    -------
    subgraph : SubGraph View
        A read-only view of the subgraph in `G` induced by the nodes.
        Changes to the graph `G` will be reflected in the view.

    Notes
    -----
    To create a mutable subgraph with its own copies of nodes
    edges and attributes use `subgraph.copy()` or `Graph(subgraph)`

    For an inplace reduction of a graph to a subgraph you can remove nodes:
    `G.remove_nodes_from(n in G if n not in set(nbunch))`

    If you are going to compute subgraphs of your subgraphs you could
    end up with a chain of views that can be very slow once the chain
    has about 15 views in it. If they are all induced subgraphs, you
    can short-cut the chain by making them all subgraphs of the original
    graph. The graph class method `G.subgraph` does this when `G` is
    a subgraph. In contrast, this function allows you to choose to build
    chains or not, as you wish. The returned subgraph is a view on `G`.

    Examples
    --------
    >>> import networkx as nx
    >>> G = nx.path_graph(4)  # or DiGraph, MultiGraph, MultiDiGraph, etc
    >>> H = G.subgraph([0, 1, 2])
    >>> list(H.edges)
    [(0, 1), (1, 2)]
    """
    induced_nodes = nx.filters.show_nodes(G.nbunch_iter(nbunch))
    if G.is_multigraph():
        if G.is_directed():
            return nx.graphviews.SubMultiDiGraph(G, induced_nodes)
        return nx.graphviews.SubMultiGraph(G, induced_nodes)
    if G.is_directed():
        return nx.graphviews.SubDiGraph(G, induced_nodes)
    return nx.graphviews.SubGraph(G, induced_nodes)


def edge_subgraph(G, edges):
    """Returns a view of the subgraph induced by the specified edges.

    The induced subgraph contains each edge in `edges` and each
    node incident to any of those edges.

    Parameters
    ----------
    G : NetworkX Graph
    edges : iterable
        An iterable of edges. Edges not present in `G` are ignored.

    Returns
    -------
    subgraph : SubGraph View
        A read-only edge-induced subgraph of `G`.
        Changes to `G` are reflected in the view.

    Notes
    -----
    To create a mutable subgraph with its own copies of nodes
    edges and attributes use `subgraph.copy()` or `Graph(subgraph)`

    If you create a subgraph of a subgraph recursively you can end up
    with a chain of subgraphs that becomes very slow with about 15
    nested subgraph views. Luckily the edge_subgraph filter nests
    nicely so you can use the original graph (`subgraph.root_graph`)
    as G in this function to avoid chains. We do not rule out chains
    programmatically so that odd cases like an `edge_subgraph` of a
    `restricted_view` can be created.

    Examples
    --------
    >>> import networkx as nx
    >>> G = nx.path_graph(5)
    >>> H = G.edge_subgraph([(0, 1), (3, 4)])
    >>> list(H.nodes)
    [0, 1, 3, 4]
    >>> list(H.edges)
    [(0, 1), (3, 4)]
    """
    nxf = nx.filters
    nxg = nx.graphviews
    edges = set(edges)
    nodes = set()
    for e in edges:
        nodes.update(e[:2])
    induced_nodes = nxf.show_nodes(nodes)
    if G.is_multigraph():
        if G.is_directed():
            induced_edges = nxf.show_multidiedges(edges)
            return nxg.SubMultiDiGraph(G, induced_nodes, induced_edges)
        induced_edges = nxf.show_multiedges(edges)
        return nxg.SubMultiGraph(G, induced_nodes, induced_edges)
    if G.is_directed():
        induced_edges = nxf.show_diedges(edges)
        return nxg.SubDiGraph(G, induced_nodes, induced_edges)
    induced_edges = nxf.show_edges(edges)
    return nxg.SubGraph(G, induced_nodes, induced_edges)


def restricted_view(G, nodes, edges):
    """Returns a view of `G` with hidden nodes and edges.

    The resulting subgraph filters out node `nodes` and edges `edges`.
    Filtered out nodes also filter out any of their edges.

    Parameters
    ----------
    G : NetworkX Graph
    nodes : iterable
        An iterable of nodes. Nodes not present in `G` are ignored.
    edges : iterable
        An iterable of edges. Edges not present in `G` are ignored.

    Returns
    -------
    subgraph : SubGraph View
        A read-only restricted view of `G` filtering out nodes and edges.
        Changes to `G` are reflected in the view.

    Notes
    -----
    To create a mutable subgraph with its own copies of nodes
    edges and attributes use `subgraph.copy()` or `Graph(subgraph)`

    If you create a subgraph of a subgraph recursively you may end up
    with a chain of subgraph views. Such chains can get quite slow
    for lengths near 15. To avoid long chains, try to make your subgraph
    based on the original graph (`subgraph.root_graph`). We do not
    rule out chains programatically so that odd cases like an
    `edge_subgraph` of a `restricted_view` can be created.

    Examples
    --------
    >>> import networkx as nx
    >>> G = nx.path_graph(5)
    >>> H = nx.restricted_view(G, [0], [(1, 2), (3, 4)])
    >>> list(H.nodes)
    [1, 2, 3, 4]
    >>> list(H.edges)
    [(2, 3)]
    """
    nxf = nx.filters
    nxg = nx.graphviews
    h_nodes = nxf.hide_nodes(nodes)
    if G.is_multigraph():
        if G.is_directed():
            h_edges = nxf.hide_multidiedges(edges)
            return nxg.SubMultiDiGraph(G, h_nodes, h_edges)
        h_edges = nxf.hide_multiedges(edges)
        return nxg.SubMultiGraph(G, h_nodes, h_edges)
    if G.is_directed():
        h_edges = nxf.hide_diedges(edges)
        return nxg.SubDiGraph(G, h_nodes, h_edges)
    h_edges = nxf.hide_edges(edges)
    return nxg.SubGraph(G, h_nodes, h_edges)


@not_implemented_for('undirected')
def reverse_view(digraph):
    """Provide a reverse view of the digraph with edges reversed.

    Identical to digraph.reverse(copy=False)
    """
    if digraph.is_multigraph():
        return nx.graphviews.MultiReverseView(digraph)
    return nx.graphviews.ReverseView(digraph)


def to_directed(graph):
    """Return a directed view of the graph `graph`.

    Identical to graph.to_directed(as_view=True)
    """
    if graph.is_multigraph():
        return nx.graphviews.MultiDiGraphView(graph)
    return nx.graphviews.DiGraphView(graph)


def to_undirected(graph):
    """Return an undirected view of the graph `graph`.

    Identical to graph.to_undirected(as_view=True)
    """
    if graph.is_multigraph():
        return nx.graphviews.MultiGraphView(graph)
    return nx.graphviews.GraphView(graph)


def create_empty_copy(G, with_data=True):
    """Return a copy of the graph G with all of the edges removed.

    Parameters
    ----------
    G : graph
       A NetworkX graph

    with_data :  bool (default=True)
       Propagate Graph and Nodes data to the new graph.

    See Also
    -----
    empty_graph

    """
    H = G.fresh_copy()
    H.add_nodes_from(G.nodes(data=with_data))
    if with_data:
        H.graph.update(G.graph)
    return H


def info(G, n=None):
    """Print short summary of information for the graph G or the node n.

    Parameters
    ----------
    G : Networkx graph
       A graph
    n : node (any hashable)
       A node in the graph G
    """
    info = ''  # append this all to a string
    if n is None:
        info += "Name: %s\n" % G.name
        type_name = [type(G).__name__]
        info += "Type: %s\n" % ",".join(type_name)
        info += "Number of nodes: %d\n" % G.number_of_nodes()
        info += "Number of edges: %d\n" % G.number_of_edges()
        nnodes = G.number_of_nodes()
        if len(G) > 0:
            if G.is_directed():
                deg = sum(d for n, d in G.in_degree()) / float(nnodes)
                info += "Average in degree: %8.4f\n" % deg
                deg = sum(d for n, d in G.out_degree()) / float(nnodes)
                info += "Average out degree: %8.4f" % deg
            else:
                s = sum(dict(G.degree()).values())
                info += "Average degree: %8.4f" % (float(s) / float(nnodes))

    else:
        if n not in G:
            raise nx.NetworkXError("node %s not in graph" % (n,))
        info += "Node % s has the following properties:\n" % n
        info += "Degree: %d\n" % G.degree(n)
        info += "Neighbors: "
        info += ' '.join(str(nbr) for nbr in G.neighbors(n))
    return info


def set_node_attributes(G, values, name=None):
    """Sets node attributes from a given value or dictionary of values.

    Parameters
    ----------
    G : NetworkX Graph

    values : scalar value, dict-like
        What the node attribute should be set to.  If `values` is
        not a dictionary, then it is treated as a single attribute value
        that is then applied to every node in `G`.  This means that if
        you provide a mutable object, like a list, updates to that object
        will be reflected in the node attribute for each edge.  The attribute
        name will be `name`.

        If `values` is a dict or a dict of dict, the corresponding node's
        attributes will be updated to `values`.

    name : string (optional, default=None)
        Name of the node attribute to set if values is a scalar.

    Examples
    --------
    After computing some property of the nodes of a graph, you may want
    to assign a node attribute to store the value of that property for
    each node::

        >>> G = nx.path_graph(3)
        >>> bb = nx.betweenness_centrality(G)
        >>> isinstance(bb, dict)
        True
        >>> nx.set_node_attributes(G, bb, 'betweenness')
        >>> G.nodes[1]['betweenness']
        1.0

    If you provide a list as the second argument, updates to the list
    will be reflected in the node attribute for each node::

        >>> G = nx.path_graph(3)
        >>> labels = []
        >>> nx.set_node_attributes(G, labels, 'labels')
        >>> labels.append('foo')
        >>> G.nodes[0]['labels']
        ['foo']
        >>> G.nodes[1]['labels']
        ['foo']
        >>> G.nodes[2]['labels']
        ['foo']

    If you provide a dictionary of dictionaries as the second argument,
    the entire dictionary will be used to update node attributes::

        >>> G = nx.path_graph(3)
        >>> attrs = {0: {'attr1': 20, 'attr2': 'nothing'}, 1: {'attr2': 3}}
        >>> nx.set_node_attributes(G, attrs)
        >>> G.nodes[0]['attr1']
        20
        >>> G.nodes[0]['attr2']
        'nothing'
        >>> G.nodes[1]['attr2']
        3
        >>> G.nodes[2]
        {}

    """
    # Set node attributes based on type of `values`
    if name is not None:  # `values` must not be a dict of dict
        try:  # `values` is a dict
            for n, v in values.items():
                try:
                    G.nodes[n][name] = values[n]
                except KeyError:
                    pass
        except AttributeError:  # `values` is a constant
            for n in G:
                G.nodes[n][name] = values
    else:  # `values` must be dict of dict
        for n, d in values.items():
            try:
                G.nodes[n].update(d)
            except KeyError:
                pass


def get_node_attributes(G, name):
    """Get node attributes from graph

    Parameters
    ----------
    G : NetworkX Graph

    name : string
       Attribute name

    Returns
    -------
    Dictionary of attributes keyed by node.

    Examples
    --------
    >>> G = nx.Graph()
    >>> G.add_nodes_from([1, 2, 3], color='red')
    >>> color = nx.get_node_attributes(G, 'color')
    >>> color[1]
    'red'
    """
    return {n: d[name] for n, d in G.nodes.items() if name in d}


def set_edge_attributes(G, values, name=None):
    """Sets edge attributes from a given value or dictionary of values.

    Parameters
    ----------
    G : NetworkX Graph

    values : scalar value, dict-like
        What the edge attribute should be set to.  If `values` is
        not a dictionary, then it is treated as a single attribute value
        that is then applied to every edge in `G`.  This means that if
        you provide a mutable object, like a list, updates to that object
        will be reflected in the edge attribute for each edge.  The attribute
        name will be `name`.

        If `values` is a dict or a dict of dict, the corresponding edge'
        attributes will be updated to `values`.  For multigraphs, the tuples
        must be of the form ``(u, v, key)``, where `u` and `v` are nodes
        and `key` is the key corresponding to the edge.  For non-multigraphs,
        the keys must be tuples of the form ``(u, v)``.

    name : string (optional, default=None)
        Name of the edge attribute to set if values is a scalar.

    Examples
    --------
    After computing some property of the edges of a graph, you may want
    to assign a edge attribute to store the value of that property for
    each edge::

        >>> G = nx.path_graph(3)
        >>> bb = nx.edge_betweenness_centrality(G, normalized=False)
        >>> nx.set_edge_attributes(G, bb, 'betweenness')
        >>> G.edges[1, 2]['betweenness']
        2.0

    If you provide a list as the second argument, updates to the list
    will be reflected in the edge attribute for each edge::

        >>> labels = []
        >>> nx.set_edge_attributes(G, labels, 'labels')
        >>> labels.append('foo')
        >>> G.edges[0, 1]['labels']
        ['foo']
        >>> G.edges[1, 2]['labels']
        ['foo']

    If you provide a dictionary of dictionaries as the second argument,
    the entire dictionary will be used to update edge attributes::

        >>> G = nx.path_graph(3)
        >>> attrs = {(0, 1): {'attr1': 20, 'attr2': 'nothing'},
        ...          (1, 2): {'attr2': 3}}
        >>> nx.set_edge_attributes(G, attrs)
        >>> G[0][1]['attr1']
        20
        >>> G[0][1]['attr2']
        'nothing'
        >>> G[1][2]['attr2']
        3

    """
    if name is not None:
        # `values` does not contain attribute names
        try:
            # if `values` is a dict using `.items()` => {edge: value}
            if G.is_multigraph():
                for (u, v, key), value in values.items():
                    try:
                        G[u][v][key][name] = value
                    except KeyError:
                        pass
            else:
                for (u, v), value in values.items():
                    try:
                        G[u][v][name] = value
                    except KeyError:
                        pass
        except AttributeError:
            # treat `values` as a constant
            for u, v, data in G.edges(data=True):
                data[name] = values
    else:
        # `values` consists of doct-of-dict {edge: {attr: value}} shape
        if G.is_multigraph():
            for (u, v, key), d in values.items():
                try:
                    G[u][v][key].update(d)
                except KeyError:
                    pass
        else:
            for (u, v), d in values.items():
                try:
                    G[u][v].update(d)
                except KeyError:
                    pass


def get_edge_attributes(G, name):
    """Get edge attributes from graph

    Parameters
    ----------
    G : NetworkX Graph

    name : string
       Attribute name

    Returns
    -------
    Dictionary of attributes keyed by edge. For (di)graphs, the keys are
    2-tuples of the form: (u, v). For multi(di)graphs, the keys are 3-tuples of
    the form: (u, v, key).

    Examples
    --------
    >>> G = nx.Graph()
    >>> nx.add_path(G, [1, 2, 3], color='red')
    >>> color = nx.get_edge_attributes(G, 'color')
    >>> color[(1, 2)]
    'red'
    """
    if G.is_multigraph():
        edges = G.edges(keys=True, data=True)
    else:
        edges = G.edges(data=True)
    return {x[:-1]: x[-1][name] for x in edges if name in x[-1]}


def all_neighbors(graph, node):
    """Returns all of the neighbors of a node in the graph.

    If the graph is directed returns predecessors as well as successors.

    Parameters
    ----------
    graph : NetworkX graph
        Graph to find neighbors.

    node : node
        The node whose neighbors will be returned.

    Returns
    -------
    neighbors : iterator
        Iterator of neighbors
    """
    if graph.is_directed():
        values = chain(graph.predecessors(node), graph.successors(node))
    else:
        values = graph.neighbors(node)
    return values


def non_neighbors(graph, node):
    """Returns the non-neighbors of the node in the graph.

    Parameters
    ----------
    graph : NetworkX graph
        Graph to find neighbors.

    node : node
        The node whose neighbors will be returned.

    Returns
    -------
    non_neighbors : iterator
        Iterator of nodes in the graph that are not neighbors of the node.
    """
    nbors = set(neighbors(graph, node)) | {node}
    return (nnode for nnode in graph if nnode not in nbors)


def non_edges(graph):
    """Returns the non-existent edges in the graph.

    Parameters
    ----------
    graph : NetworkX graph.
        Graph to find non-existent edges.

    Returns
    -------
    non_edges : iterator
        Iterator of edges that are not in the graph.
    """
    if graph.is_directed():
        for u in graph:
            for v in non_neighbors(graph, u):
                yield (u, v)
    else:
        nodes = set(graph)
        while nodes:
            u = nodes.pop()
            for v in nodes - set(graph[u]):
                yield (u, v)


@not_implemented_for('directed')
def common_neighbors(G, u, v):
    """Return the common neighbors of two nodes in a graph.

    Parameters
    ----------
    G : graph
        A NetworkX undirected graph.

    u, v : nodes
        Nodes in the graph.

    Returns
    -------
    cnbors : iterator
        Iterator of common neighbors of u and v in the graph.

    Raises
    ------
    NetworkXError
        If u or v is not a node in the graph.

    Examples
    --------
    >>> G = nx.complete_graph(5)
    >>> sorted(nx.common_neighbors(G, 0, 1))
    [2, 3, 4]
    """
    if u not in G:
        raise nx.NetworkXError('u is not in the graph.')
    if v not in G:
        raise nx.NetworkXError('v is not in the graph.')

    # Return a generator explicitly instead of yielding so that the above
    # checks are executed eagerly.
    return (w for w in G[u] if w in G[v] and w not in (u, v))


def is_weighted(G, edge=None, weight='weight'):
    """Returns True if `G` has weighted edges.

    Parameters
    ----------
    G : graph
        A NetworkX graph.

    edge : tuple, optional
        A 2-tuple specifying the only edge in `G` that will be tested. If
        None, then every edge in `G` is tested.

    weight: string, optional
        The attribute name used to query for edge weights.

    Returns
    -------
    bool
        A boolean signifying if `G`, or the specified edge, is weighted.

    Raises
    ------
    NetworkXError
        If the specified edge does not exist.

    Examples
    --------
    >>> G = nx.path_graph(4)
    >>> nx.is_weighted(G)
    False
    >>> nx.is_weighted(G, (2, 3))
    False

    >>> G = nx.DiGraph()
    >>> G.add_edge(1, 2, weight=1)
    >>> nx.is_weighted(G)
    True

    """
    if edge is not None:
        data = G.get_edge_data(*edge)
        if data is None:
            msg = 'Edge {!r} does not exist.'.format(edge)
            raise nx.NetworkXError(msg)
        return weight in data

    if is_empty(G):
        # Special handling required since: all([]) == True
        return False

    return all(weight in data for u, v, data in G.edges(data=True))


def is_negatively_weighted(G, edge=None, weight='weight'):
    """Returns True if `G` has negatively weighted edges.

    Parameters
    ----------
    G : graph
        A NetworkX graph.

    edge : tuple, optional
        A 2-tuple specifying the only edge in `G` that will be tested. If
        None, then every edge in `G` is tested.

    weight: string, optional
        The attribute name used to query for edge weights.

    Returns
    -------
    bool
        A boolean signifying if `G`, or the specified edge, is negatively
        weighted.

    Raises
    ------
    NetworkXError
        If the specified edge does not exist.

    Examples
    --------
    >>> G = nx.Graph()
    >>> G.add_edges_from([(1, 3), (2, 4), (2, 6)])
    >>> G.add_edge(1, 2, weight=4)
    >>> nx.is_negatively_weighted(G, (1, 2))
    False
    >>> G[2][4]['weight'] = -2
    >>> nx.is_negatively_weighted(G)
    True
    >>> G = nx.DiGraph()
    >>> edges = [('0', '3', 3), ('0', '1', -5), ('1', '0', -2)]
    >>> G.add_weighted_edges_from(edges)
    >>> nx.is_negatively_weighted(G)
    True

    """
    if edge is not None:
        data = G.get_edge_data(*edge)
        if data is None:
            msg = 'Edge {!r} does not exist.'.format(edge)
            raise nx.NetworkXError(msg)
        return weight in data and data[weight] < 0

    return any(weight in data and data[weight] < 0
               for u, v, data in G.edges(data=True))


def is_empty(G):
    """Returns True if `G` has no edges.

    Parameters
    ----------
    G : graph
        A NetworkX graph.

    Returns
    -------
    bool
        True if `G` has no edges, and False otherwise.

    Notes
    -----
    An empty graph can have nodes but not edges. The empty graph with zero
    nodes is known as the null graph. This is an $O(n)$ operation where n
    is the number of nodes in the graph.

    """
    return not any(G.adj.values())


def nodes_with_selfloops(G):
    """Returns an iterator over nodes with self loops.

    A node with a self loop has an edge with both ends adjacent
    to that node.

    Returns
    -------
    nodelist : iterator
        A iterator over nodes with self loops.

    See Also
    --------
    selfloop_edges, number_of_selfloops

    Examples
    --------
    >>> G = nx.Graph()   # or DiGraph, MultiGraph, MultiDiGraph, etc
    >>> G.add_edge(1, 1)
    >>> G.add_edge(1, 2)
    >>> list(nx.nodes_with_selfloops(G))
    [1]

    """
    return (n for n, nbrs in G.adj.items() if n in nbrs)


def selfloop_edges(G, data=False, keys=False, default=None):
    """Returns an iterator over selfloop edges.

    A selfloop edge has the same node at both ends.

    Parameters
    ----------
    data : string or bool, optional (default=False)
        Return selfloop edges as two tuples (u, v) (data=False)
        or three-tuples (u, v, datadict) (data=True)
        or three-tuples (u, v, datavalue) (data='attrname')
    keys : bool, optional (default=False)
        If True, return edge keys with each edge.
    default : value, optional (default=None)
        Value used for edges that dont have the requested attribute.
        Only relevant if data is not True or False.

    Returns
    -------
    edgeiter : iterator over edge tuples
        An iterator over all selfloop edges.

    See Also
    --------
    nodes_with_selfloops, number_of_selfloops

    Examples
    --------
    >>> G = nx.MultiGraph()   # or Graph, DiGraph, MultiDiGraph, etc
    >>> ekey = G.add_edge(1, 1)
    >>> ekey = G.add_edge(1, 2)
    >>> list(nx.selfloop_edges(G))
    [(1, 1)]
    >>> list(nx.selfloop_edges(G, data=True))
    [(1, 1, {})]
    >>> list(nx.selfloop_edges(G, keys=True))
    [(1, 1, 0)]
    >>> list(nx.selfloop_edges(G, keys=True, data=True))
    [(1, 1, 0, {})]
    """
    if data is True:
        if G.is_multigraph():
            if keys is True:
                return ((n, n, k, d)
                        for n, nbrs in G.adj.items()
                        if n in nbrs for k, d in nbrs[n].items())
            else:
                return ((n, n, d)
                        for n, nbrs in G.adj.items()
                        if n in nbrs for d in nbrs[n].values())
        else:
            return ((n, n, nbrs[n]) for n, nbrs in G.adj.items() if n in nbrs)
    elif data is not False:
        if G.is_multigraph():
            if keys is True:
                return ((n, n, k, d.get(data, default))
                        for n, nbrs in G.adj.items()
                        if n in nbrs for k, d in nbrs[n].items())
            else:
                return ((n, n, d.get(data, default))
                        for n, nbrs in G.adj.items()
                        if n in nbrs for d in nbrs[n].values())
        else:
            return ((n, n, nbrs[n].get(data, default))
                    for n, nbrs in G.adj.items() if n in nbrs)
    else:
        if G.is_multigraph():
            if keys is True:
                return ((n, n, k)
                        for n, nbrs in G.adj.items()
                        if n in nbrs for k in nbrs[n])
            else:
                return ((n, n)
                        for n, nbrs in G.adj.items()
                        if n in nbrs for d in nbrs[n].values())
        else:
            return ((n, n) for n, nbrs in G.adj.items() if n in nbrs)


def number_of_selfloops(G):
    """Return the number of selfloop edges.

    A selfloop edge has the same node at both ends.

    Returns
    -------
    nloops : int
        The number of selfloops.

    See Also
    --------
    nodes_with_selfloops, selfloop_edges

    Examples
    --------
    >>> G = nx.Graph()   # or DiGraph, MultiGraph, MultiDiGraph, etc
    >>> G.add_edge(1, 1)
    >>> G.add_edge(1, 2)
    >>> nx.number_of_selfloops(G)
    1
    """
    return sum(1 for _ in nx.selfloop_edges(G))
