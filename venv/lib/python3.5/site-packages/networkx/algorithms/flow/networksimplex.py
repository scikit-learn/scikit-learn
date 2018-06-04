# -*- coding: utf-8 -*-
"""
Minimum cost flow algorithms on directed connected graphs.
"""

__author__ = """Loïc Séguin-C. <loicseguin@gmail.com>"""
# Copyright (C) 2010 Loïc Séguin-C. <loicseguin@gmail.com>
# All rights reserved.
# BSD license.

__all__ = ['network_simplex']

from itertools import chain, islice, repeat
from math import ceil, sqrt
import networkx as nx
from networkx.utils import not_implemented_for

try:
    from itertools import izip as zip
except ImportError:
    pass
try:
    range = xrange
except NameError:
    pass


@not_implemented_for('undirected')
def network_simplex(G, demand='demand', capacity='capacity', weight='weight'):
    r"""Find a minimum cost flow satisfying all demands in digraph G.

    This is a primal network simplex algorithm that uses the leaving
    arc rule to prevent cycling.

    G is a digraph with edge costs and capacities and in which nodes
    have demand, i.e., they want to send or receive some amount of
    flow. A negative demand means that the node wants to send flow, a
    positive demand means that the node want to receive flow. A flow on
    the digraph G satisfies all demand if the net flow into each node
    is equal to the demand of that node.

    Parameters
    ----------
    G : NetworkX graph
        DiGraph on which a minimum cost flow satisfying all demands is
        to be found.

    demand : string
        Nodes of the graph G are expected to have an attribute demand
        that indicates how much flow a node wants to send (negative
        demand) or receive (positive demand). Note that the sum of the
        demands should be 0 otherwise the problem in not feasible. If
        this attribute is not present, a node is considered to have 0
        demand. Default value: 'demand'.

    capacity : string
        Edges of the graph G are expected to have an attribute capacity
        that indicates how much flow the edge can support. If this
        attribute is not present, the edge is considered to have
        infinite capacity. Default value: 'capacity'.

    weight : string
        Edges of the graph G are expected to have an attribute weight
        that indicates the cost incurred by sending one unit of flow on
        that edge. If not present, the weight is considered to be 0.
        Default value: 'weight'.

    Returns
    -------
    flowCost : integer, float
        Cost of a minimum cost flow satisfying all demands.

    flowDict : dictionary
        Dictionary of dictionaries keyed by nodes such that
        flowDict[u][v] is the flow edge (u, v).

    Raises
    ------
    NetworkXError
        This exception is raised if the input graph is not directed,
        not connected or is a multigraph.

    NetworkXUnfeasible
        This exception is raised in the following situations:

            * The sum of the demands is not zero. Then, there is no
              flow satisfying all demands.
            * There is no flow satisfying all demand.

    NetworkXUnbounded
        This exception is raised if the digraph G has a cycle of
        negative cost and infinite capacity. Then, the cost of a flow
        satisfying all demands is unbounded below.

    Notes
    -----
    This algorithm is not guaranteed to work if edge weights or demands
    are floating point numbers (overflows and roundoff errors can
    cause problems). As a workaround you can use integer numbers by
    multiplying the relevant edge attributes by a convenient
    constant factor (eg 100).

    See also
    --------
    cost_of_flow, max_flow_min_cost, min_cost_flow, min_cost_flow_cost

    Examples
    --------
    A simple example of a min cost flow problem.

    >>> import networkx as nx
    >>> G = nx.DiGraph()
    >>> G.add_node('a', demand=-5)
    >>> G.add_node('d', demand=5)
    >>> G.add_edge('a', 'b', weight=3, capacity=4)
    >>> G.add_edge('a', 'c', weight=6, capacity=10)
    >>> G.add_edge('b', 'd', weight=1, capacity=9)
    >>> G.add_edge('c', 'd', weight=2, capacity=5)
    >>> flowCost, flowDict = nx.network_simplex(G)
    >>> flowCost
    24
    >>> flowDict # doctest: +SKIP
    {'a': {'c': 1, 'b': 4}, 'c': {'d': 1}, 'b': {'d': 4}, 'd': {}}

    The mincost flow algorithm can also be used to solve shortest path
    problems. To find the shortest path between two nodes u and v,
    give all edges an infinite capacity, give node u a demand of -1 and
    node v a demand a 1. Then run the network simplex. The value of a
    min cost flow will be the distance between u and v and edges
    carrying positive flow will indicate the path.

    >>> G=nx.DiGraph()
    >>> G.add_weighted_edges_from([('s', 'u' ,10), ('s' ,'x' ,5),
    ...                            ('u', 'v' ,1), ('u' ,'x' ,2),
    ...                            ('v', 'y' ,1), ('x' ,'u' ,3),
    ...                            ('x', 'v' ,5), ('x' ,'y' ,2),
    ...                            ('y', 's' ,7), ('y' ,'v' ,6)])
    >>> G.add_node('s', demand = -1)
    >>> G.add_node('v', demand = 1)
    >>> flowCost, flowDict = nx.network_simplex(G)
    >>> flowCost == nx.shortest_path_length(G, 's', 'v', weight='weight')
    True
    >>> sorted([(u, v) for u in flowDict for v in flowDict[u] if flowDict[u][v] > 0])
    [('s', 'x'), ('u', 'v'), ('x', 'u')]
    >>> nx.shortest_path(G, 's', 'v', weight = 'weight')
    ['s', 'x', 'u', 'v']

    It is possible to change the name of the attributes used for the
    algorithm.

    >>> G = nx.DiGraph()
    >>> G.add_node('p', spam=-4)
    >>> G.add_node('q', spam=2)
    >>> G.add_node('a', spam=-2)
    >>> G.add_node('d', spam=-1)
    >>> G.add_node('t', spam=2)
    >>> G.add_node('w', spam=3)
    >>> G.add_edge('p', 'q', cost=7, vacancies=5)
    >>> G.add_edge('p', 'a', cost=1, vacancies=4)
    >>> G.add_edge('q', 'd', cost=2, vacancies=3)
    >>> G.add_edge('t', 'q', cost=1, vacancies=2)
    >>> G.add_edge('a', 't', cost=2, vacancies=4)
    >>> G.add_edge('d', 'w', cost=3, vacancies=4)
    >>> G.add_edge('t', 'w', cost=4, vacancies=1)
    >>> flowCost, flowDict = nx.network_simplex(G, demand='spam',
    ...                                         capacity='vacancies',
    ...                                         weight='cost')
    >>> flowCost
    37
    >>> flowDict  # doctest: +SKIP
    {'a': {'t': 4}, 'd': {'w': 2}, 'q': {'d': 1}, 'p': {'q': 2, 'a': 2}, 't': {'q': 1, 'w': 1}, 'w': {}}

    References
    ----------
    .. [1] Z. Kiraly, P. Kovacs.
           Efficient implementation of minimum-cost flow algorithms.
           Acta Universitatis Sapientiae, Informatica 4(1):67--118. 2012.
    .. [2] R. Barr, F. Glover, D. Klingman.
           Enhancement of spanning tree labeling procedures for network
           optimization.
           INFOR 17(1):16--34. 1979.
    """
    ###########################################################################
    # Problem essentials extraction and sanity check
    ###########################################################################

    if len(G) == 0:
        raise nx.NetworkXError('graph has no nodes')

    # Number all nodes and edges and hereafter reference them using ONLY their
    # numbers

    N = list(G)                                # nodes
    I = {u: i for i, u in enumerate(N)}        # node indices
    D = [G.nodes[u].get(demand, 0) for u in N]  # node demands

    inf = float('inf')
    for p, b in zip(N, D):
        if abs(b) == inf:
            raise nx.NetworkXError('node %r has infinite demand' % (p,))

    multigraph = G.is_multigraph()
    S = []  # edge sources
    T = []  # edge targets
    if multigraph:
        K = []  # edge keys
    E = {}  # edge indices
    U = []  # edge capacities
    C = []  # edge weights

    if not multigraph:
        edges = G.edges(data=True)
    else:
        edges = G.edges(data=True, keys=True)
    edges = (e for e in edges
             if e[0] != e[1] and e[-1].get(capacity, inf) != 0)
    for i, e in enumerate(edges):
        S.append(I[e[0]])
        T.append(I[e[1]])
        if multigraph:
            K.append(e[2])
        E[e[:-1]] = i
        U.append(e[-1].get(capacity, inf))
        C.append(e[-1].get(weight, 0))

    for e, c in zip(E, C):
        if abs(c) == inf:
            raise nx.NetworkXError('edge %r has infinite weight' % (e,))
    if not multigraph:
        edges = nx.selfloop_edges(G, data=True)
    else:
        edges = nx.selfloop_edges(G, data=True, keys=True)
    for e in edges:
        if abs(e[-1].get(weight, 0)) == inf:
            raise nx.NetworkXError('edge %r has infinite weight' % (e[:-1],))

    ###########################################################################
    # Quick infeasibility detection
    ###########################################################################

    if sum(D) != 0:
        raise nx.NetworkXUnfeasible('total node demand is not zero')
    for e, u in zip(E, U):
        if u < 0:
            raise nx.NetworkXUnfeasible('edge %r has negative capacity' % (e,))
    if not multigraph:
        edges = nx.selfloop_edges(G, data=True)
    else:
        edges = nx.selfloop_edges(G, data=True, keys=True)
    for e in edges:
        if e[-1].get(capacity, inf) < 0:
            raise nx.NetworkXUnfeasible(
                'edge %r has negative capacity' % (e[:-1],))

    ###########################################################################
    # Initialization
    ###########################################################################

    # Add a dummy node -1 and connect all existing nodes to it with infinite-
    # capacity dummy edges. Node -1 will serve as the root of the
    # spanning tree of the network simplex method. The new edges will used to
    # trivially satisfy the node demands and create an initial strongly
    # feasible spanning tree.
    n = len(N)  # number of nodes
    for p, d in enumerate(D):
        if d > 0:  # Must be greater-than here. Zero-demand nodes must have
                   # edges pointing towards the root to ensure strong
                   # feasibility.
            S.append(-1)
            T.append(p)
        else:
            S.append(p)
            T.append(-1)
    faux_inf = 3 * max(chain([sum(u for u in U if u < inf),
                              sum(abs(c) for c in C)],
                             (abs(d) for d in D))) or 1
    C.extend(repeat(faux_inf, n))
    U.extend(repeat(faux_inf, n))

    # Construct the initial spanning tree.
    e = len(E)                                           # number of edges
    x = list(chain(repeat(0, e), (abs(d) for d in D)))   # edge flows
    pi = [faux_inf if d <= 0 else -faux_inf for d in D]  # node potentials
    parent = list(chain(repeat(-1, n), [None]))  # parent nodes
    edge = list(range(e, e + n))                 # edges to parents
    size = list(chain(repeat(1, n), [n + 1]))    # subtree sizes
    next = list(chain(range(1, n), [-1, 0]))     # next nodes in depth-first thread
    prev = list(range(-1, n))                    # previous nodes in depth-first thread
    last = list(chain(range(n), [n - 1]))        # last descendants in depth-first thread

    ###########################################################################
    # Pivot loop
    ###########################################################################

    def reduced_cost(i):
        """Return the reduced cost of an edge i.
        """
        c = C[i] - pi[S[i]] + pi[T[i]]
        return c if x[i] == 0 else -c

    def find_entering_edges():
        """Yield entering edges until none can be found.
        """
        if e == 0:
            return

        # Entering edges are found by combining Dantzig's rule and Bland's
        # rule. The edges are cyclically grouped into blocks of size B. Within
        # each block, Dantzig's rule is applied to find an entering edge. The
        # blocks to search is determined following Bland's rule.
        B = int(ceil(sqrt(e)))  # pivot block size
        M = (e + B - 1) // B    # number of blocks needed to cover all edges
        m = 0                   # number of consecutive blocks without eligible
        # entering edges
        f = 0                   # first edge in block
        while m < M:
            # Determine the next block of edges.
            l = f + B
            if l <= e:
                edges = range(f, l)
            else:
                l -= e
                edges = chain(range(f, e), range(l))
            f = l
            # Find the first edge with the lowest reduced cost.
            i = min(edges, key=reduced_cost)
            c = reduced_cost(i)
            if c >= 0:
                # No entering edge found in the current block.
                m += 1
            else:
                # Entering edge found.
                if x[i] == 0:
                    p = S[i]
                    q = T[i]
                else:
                    p = T[i]
                    q = S[i]
                yield i, p, q
                m = 0
        # All edges have nonnegative reduced costs. The current flow is
        # optimal.

    def find_apex(p, q):
        """Find the lowest common ancestor of nodes p and q in the spanning
        tree.
        """
        size_p = size[p]
        size_q = size[q]
        while True:
            while size_p < size_q:
                p = parent[p]
                size_p = size[p]
            while size_p > size_q:
                q = parent[q]
                size_q = size[q]
            if size_p == size_q:
                if p != q:
                    p = parent[p]
                    size_p = size[p]
                    q = parent[q]
                    size_q = size[q]
                else:
                    return p

    def trace_path(p, w):
        """Return the nodes and edges on the path from node p to its ancestor
        w.
        """
        Wn = [p]
        We = []
        while p != w:
            We.append(edge[p])
            p = parent[p]
            Wn.append(p)
        return Wn, We

    def find_cycle(i, p, q):
        """Return the nodes and edges on the cycle containing edge i == (p, q)
        when the latter is added to the spanning tree.

        The cycle is oriented in the direction from p to q.
        """
        w = find_apex(p, q)
        Wn, We = trace_path(p, w)
        Wn.reverse()
        We.reverse()
        We.append(i)
        WnR, WeR = trace_path(q, w)
        del WnR[-1]
        Wn += WnR
        We += WeR
        return Wn, We

    def residual_capacity(i, p):
        """Return the residual capacity of an edge i in the direction away
        from its endpoint p.
        """
        return U[i] - x[i] if S[i] == p else x[i]

    def find_leaving_edge(Wn, We):
        """Return the leaving edge in a cycle represented by Wn and We.
        """
        j, s = min(zip(reversed(We), reversed(Wn)),
                   key=lambda i_p: residual_capacity(*i_p))
        t = T[j] if S[j] == s else S[j]
        return j, s, t

    def augment_flow(Wn, We, f):
        """Augment f units of flow along a cycle represented by Wn and We.
        """
        for i, p in zip(We, Wn):
            if S[i] == p:
                x[i] += f
            else:
                x[i] -= f

    def trace_subtree(p):
        """Yield the nodes in the subtree rooted at a node p.
        """
        yield p
        l = last[p]
        while p != l:
            p = next[p]
            yield p

    def remove_edge(s, t):
        """Remove an edge (s, t) where parent[t] == s from the spanning tree.
        """
        size_t = size[t]
        prev_t = prev[t]
        last_t = last[t]
        next_last_t = next[last_t]
        # Remove (s, t).
        parent[t] = None
        edge[t] = None
        # Remove the subtree rooted at t from the depth-first thread.
        next[prev_t] = next_last_t
        prev[next_last_t] = prev_t
        next[last_t] = t
        prev[t] = last_t
        # Update the subtree sizes and last descendants of the (old) acenstors
        # of t.
        while s is not None:
            size[s] -= size_t
            if last[s] == last_t:
                last[s] = prev_t
            s = parent[s]

    def make_root(q):
        """Make a node q the root of its containing subtree.
        """
        ancestors = []
        while q is not None:
            ancestors.append(q)
            q = parent[q]
        ancestors.reverse()
        for p, q in zip(ancestors, islice(ancestors, 1, None)):
            size_p = size[p]
            last_p = last[p]
            prev_q = prev[q]
            last_q = last[q]
            next_last_q = next[last_q]
            # Make p a child of q.
            parent[p] = q
            parent[q] = None
            edge[p] = edge[q]
            edge[q] = None
            size[p] = size_p - size[q]
            size[q] = size_p
            # Remove the subtree rooted at q from the depth-first thread.
            next[prev_q] = next_last_q
            prev[next_last_q] = prev_q
            next[last_q] = q
            prev[q] = last_q
            if last_p == last_q:
                last[p] = prev_q
                last_p = prev_q
            # Add the remaining parts of the subtree rooted at p as a subtree
            # of q in the depth-first thread.
            prev[p] = last_q
            next[last_q] = p
            next[last_p] = q
            prev[q] = last_p
            last[q] = last_p

    def add_edge(i, p, q):
        """Add an edge (p, q) to the spanning tree where q is the root of a
        subtree.
        """
        last_p = last[p]
        next_last_p = next[last_p]
        size_q = size[q]
        last_q = last[q]
        # Make q a child of p.
        parent[q] = p
        edge[q] = i
        # Insert the subtree rooted at q into the depth-first thread.
        next[last_p] = q
        prev[q] = last_p
        prev[next_last_p] = last_q
        next[last_q] = next_last_p
        # Update the subtree sizes and last descendants of the (new) ancestors
        # of q.
        while p is not None:
            size[p] += size_q
            if last[p] == last_p:
                last[p] = last_q
            p = parent[p]

    def update_potentials(i, p, q):
        """Update the potentials of the nodes in the subtree rooted at a node
        q connected to its parent p by an edge i.
        """
        if q == T[i]:
            d = pi[p] - C[i] - pi[q]
        else:
            d = pi[p] + C[i] - pi[q]
        for q in trace_subtree(q):
            pi[q] += d

    # Pivot loop
    for i, p, q in find_entering_edges():
        Wn, We = find_cycle(i, p, q)
        j, s, t = find_leaving_edge(Wn, We)
        augment_flow(Wn, We, residual_capacity(j, s))
        if i != j:  # Do nothing more if the entering edge is the same as the
                    # the leaving edge.
            if parent[t] != s:
                # Ensure that s is the parent of t.
                s, t = t, s
            if We.index(i) > We.index(j):
                # Ensure that q is in the subtree rooted at t.
                p, q = q, p
            remove_edge(s, t)
            make_root(q)
            add_edge(i, p, q)
            update_potentials(i, p, q)

    ###########################################################################
    # Infeasibility and unboundedness detection
    ###########################################################################

    if any(x[i] != 0 for i in range(-n, 0)):
        raise nx.NetworkXUnfeasible('no flow satisfies all node demands')

    if (any(x[i] * 2 >= faux_inf for i in range(e)) or
        any(e[-1].get(capacity, inf) == inf and e[-1].get(weight, 0) < 0
            for e in nx.selfloop_edges(G, data=True))):
        raise nx.NetworkXUnbounded(
            'negative cycle with infinite capacity found')

    ###########################################################################
    # Flow cost calculation and flow dict construction
    ###########################################################################

    del x[e:]
    flow_cost = sum(c * x for c, x in zip(C, x))
    flow_dict = {n: {} for n in N}

    def add_entry(e):
        """Add a flow dict entry.
        """
        d = flow_dict[e[0]]
        for k in e[1:-2]:
            try:
                d = d[k]
            except KeyError:
                t = {}
                d[k] = t
                d = t
        d[e[-2]] = e[-1]

    S = (N[s] for s in S)  # Use original nodes.
    T = (N[t] for t in T)  # Use original nodes.
    if not multigraph:
        for e in zip(S, T, x):
            add_entry(e)
        edges = G.edges(data=True)
    else:
        for e in zip(S, T, K, x):
            add_entry(e)
        edges = G.edges(data=True, keys=True)
    for e in edges:
        if e[0] != e[1]:
            if e[-1].get(capacity, inf) == 0:
                add_entry(e[:-1] + (0,))
        else:
            c = e[-1].get(weight, 0)
            if c >= 0:
                add_entry(e[:-1] + (0,))
            else:
                u = e[-1][capacity]
                flow_cost += c * u
                add_entry(e[:-1] + (u,))

    return flow_cost, flow_dict
