"""
===========================
Depth First Search on Edges
===========================

Algorithms for a depth-first traversal of edges in a graph.

"""

FORWARD = 'forward'
REVERSE = 'reverse'

__all__ = ['edge_dfs']


def helper_funcs(G, orientation):
    """
    These are various G-specific functions that help us implement the algorithm
    for all graph types: graph, multigraph, directed or not.

    """
    ignore_orientation = G.is_directed() and orientation == 'ignore'
    reverse_orientation = G.is_directed() and orientation == 'reverse'

    if ignore_orientation:
        # When we ignore the orientation, we still need to know how the edge
        # was traversed, so we add an object representing the direction.
        def out_edges(u_for_edges, **kwds):
            for edge in G.out_edges(u_for_edges, **kwds):
                yield edge + (FORWARD,)
            for edge in G.in_edges(u_for_edges, **kwds):
                yield edge + (REVERSE,)
    elif reverse_orientation:
        def out_edges(u_for_edges, **kwds):
            for edge in G.in_edges(u_for_edges, **kwds):
                yield edge + (REVERSE,)
    else:
        # If "yield from" were an option, we could pass kwds automatically.
        out_edges = G.edges

    # If every edge had a unique key, then it would be easier to track which
    # edges had been visited. Since that is not available, we will form a
    # unique identifier from the edge and key (if present). If the graph
    # is undirected, then the head and tail need to be stored as a frozenset.
    if ignore_orientation or reverse_orientation:
        # edge is a 4-tuple: (u, v, key, direction)
        # u and v always represent the true tail and head of the edge.
        def key(edge):
            # We want everything but the direction.
            return edge[:-1]
    else:
        if G.is_directed():
            def key(edge):
                return edge
        else:
            # edge is a 3-tuple:  (u, v, key)
            def key(edge):
                new_edge = (frozenset(edge[:2]),) + edge[2:]
                return new_edge

    def traversed_tailhead(edge):
        """
        Returns the tail and head of an edge, as it was traversed.

        So in general, this is different from the true tail and head.
        (Also, undirected edges have no true tail or head.)

        """
        if (ignore_orientation or reverse_orientation) and edge[-1] == REVERSE:
            tail, head = edge[1], edge[0]
        else:
            tail, head = edge[0], edge[1]
        return tail, head

    return out_edges, key, traversed_tailhead


def edge_dfs(G, source=None, orientation='original'):
    """
    A directed, depth-first traversal of edges in `G`, beginning at `source`.

    Parameters
    ----------
    G : graph
        A directed/undirected graph/multigraph.

    source : node, list of nodes
        The node from which the traversal begins. If None, then a source
        is chosen arbitrarily and repeatedly until all edges from each node in
        the graph are searched.

    orientation : 'original' | 'reverse' | 'ignore'
        For directed graphs and directed multigraphs, edge traversals need not
        respect the original orientation of the edges. When set to 'reverse',
        then every edge will be traversed in the reverse direction. When set to
        'ignore', then each directed edge is treated as a single undirected
        edge that can be traversed in either direction. For undirected graphs
        and undirected multigraphs, this parameter is meaningless and is not
        consulted by the algorithm.

    Yields
    ------
    edge : directed edge
        A directed edge indicating the path taken by the depth-first traversal.
        For graphs, `edge` is of the form `(u, v)` where `u` and `v`
        are the tail and head of the edge as determined by the traversal. For
        multigraphs, `edge` is of the form `(u, v, key)`, where `key` is
        the key of the edge. When the graph is directed, then `u` and `v`
        are always in the order of the actual directed edge. If orientation is
        'reverse' or 'ignore', then `edge` takes the form
        `(u, v, key, direction)` where direction is a string, 'forward' or
        'reverse', that indicates if the edge was traversed in the forward
        (tail to head) or reverse (head to tail) direction, respectively.

    Examples
    --------
    >>> import networkx as nx
    >>> nodes = [0, 1, 2, 3]
    >>> edges = [(0, 1), (1, 0), (1, 0), (2, 1), (3, 1)]

    >>> list(nx.edge_dfs(nx.Graph(edges), nodes))
    [(0, 1), (1, 2), (1, 3)]

    >>> list(nx.edge_dfs(nx.DiGraph(edges), nodes))
    [(0, 1), (1, 0), (2, 1), (3, 1)]

    >>> list(nx.edge_dfs(nx.MultiGraph(edges), nodes))
    [(0, 1, 0), (1, 0, 1), (0, 1, 2), (1, 2, 0), (1, 3, 0)]

    >>> list(nx.edge_dfs(nx.MultiDiGraph(edges), nodes))
    [(0, 1, 0), (1, 0, 0), (1, 0, 1), (2, 1, 0), (3, 1, 0)]

    >>> list(nx.edge_dfs(nx.DiGraph(edges), nodes, orientation='ignore'))
    [(0, 1, 'forward'), (1, 0, 'forward'), (2, 1, 'reverse'), (3, 1, 'reverse')]

    >>> list(nx.edge_dfs(nx.MultiDiGraph(edges), nodes, orientation='ignore'))
    [(0, 1, 0, 'forward'), (1, 0, 0, 'forward'), (1, 0, 1, 'reverse'), (2, 1, 0, 'reverse'), (3, 1, 0, 'reverse')]

    Notes
    -----
    The goal of this function is to visit edges. It differs from the more
    familiar depth-first traversal of nodes, as provided by
    :func:`networkx.algorithms.traversal.depth_first_search.dfs_edges`, in
    that it does not stop once every node has been visited. In a directed graph
    with edges [(0, 1), (1, 2), (2, 1)], the edge (2, 1) would not be visited
    if not for the functionality provided by this function.

    See Also
    --------
    dfs_edges

    """
    nodes = list(G.nbunch_iter(source))
    if not nodes:
        raise StopIteration

    kwds = {'data': False}
    if G.is_multigraph():
        kwds['keys'] = True

    out_edges, key, tailhead = helper_funcs(G, orientation)

    visited_edges = set()
    visited_nodes = set()
    edges = {}

    for start_node in nodes:
        stack = [start_node]
        while stack:
            current_node = stack[-1]
            if current_node not in visited_nodes:
                edges[current_node] = iter(out_edges(current_node, **kwds))
                visited_nodes.add(current_node)

            try:
                edge = next(edges[current_node])
            except StopIteration:
                # No more edges from the current node.
                stack.pop()
            else:
                edge_key = key(edge)
                if edge_key not in visited_edges:
                    visited_edges.add(edge_key)
                    # Mark the traversed "to" node as to-be-explored.
                    stack.append(tailhead(edge)[1])
                    yield edge
