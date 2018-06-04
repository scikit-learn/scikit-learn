# -*- coding: utf-8 -*-
"""
Generators and functions for bipartite graphs.

"""
#    Copyright (C) 2006-2011 by
#    Aric Hagberg <hagberg@lanl.gov>
#    Dan Schult <dschult@colgate.edu>
#    Pieter Swart <swart@lanl.gov>
#    All rights reserved.
#    BSD license.
import math
import numbers
import random
import networkx
from functools import reduce
import networkx as nx
from networkx.utils import nodes_or_number

__author__ = """\n""".join(['Aric Hagberg (hagberg@lanl.gov)',
                            'Pieter Swart (swart@lanl.gov)',
                            'Dan Schult(dschult@colgate.edu)'])
__all__ = ['configuration_model',
           'havel_hakimi_graph',
           'reverse_havel_hakimi_graph',
           'alternating_havel_hakimi_graph',
           'preferential_attachment_graph',
           'random_graph',
           'gnmk_random_graph',
           'complete_bipartite_graph',
           ]


@nodes_or_number([0, 1])
def complete_bipartite_graph(n1, n2, create_using=None):
    """Return the complete bipartite graph `K_{n_1,n_2}`.

    Composed of two partitions with `n_1` nodes in the first
    and `n_2` nodes in the second. Each node in the first is
    connected to each node in the second.

    Parameters
    ----------
    n1 : integer
       Number of nodes for node set A.
    n2 : integer
       Number of nodes for node set B.
    create_using : NetworkX graph instance, optional
       Return graph of this type.

    Notes
    -----
    Node labels are the integers 0 to `n_1 + n_2 - 1`.

    The nodes are assigned the attribute 'bipartite' with the value 0 or 1
    to indicate which bipartite set the node belongs to.
    """
    if create_using is None:
        G = nx.Graph()
    else:
        if create_using.is_directed():
            raise nx.NetworkXError("Directed Graph not supported")
        G = create_using
        G.clear()

    n1, top = n1
    n2, bottom = n2
    if isinstance(n2, numbers.Integral):
        bottom = [n1 + i for i in bottom]
    G.add_nodes_from(top, bipartite=0)
    G.add_nodes_from(bottom, bipartite=1)
    G.add_edges_from((u, v) for u in top for v in bottom)
    G.graph['name'] = "complete_bipartite_graph(%s,%s)" % (n1, n2)
    return G


def configuration_model(aseq, bseq, create_using=None, seed=None):
    """Return a random bipartite graph from two given degree sequences.

    Parameters
    ----------
    aseq : list
       Degree sequence for node set A.
    bseq : list
       Degree sequence for node set B.
    create_using : NetworkX graph instance, optional
       Return graph of this type.
    seed : integer, optional
       Seed for random number generator.

    Nodes from the set A are connected to nodes in the set B by
    choosing randomly from the possible free stubs, one in A and
    one in B.

    Notes
    -----
    The sum of the two sequences must be equal: sum(aseq)=sum(bseq)
    If no graph type is specified use MultiGraph with parallel edges.
    If you want a graph with no parallel edges use create_using=Graph()
    but then the resulting degree sequences might not be exact.

    The nodes are assigned the attribute 'bipartite' with the value 0 or 1
    to indicate which bipartite set the node belongs to.

    This function is not imported in the main namespace.
    To use it you have to explicitly import the bipartite package.
    """
    if create_using is None:
        create_using = networkx.MultiGraph()
    elif create_using.is_directed():
        raise networkx.NetworkXError(
            "Directed Graph not supported")

    G = networkx.empty_graph(0, create_using)

    if not seed is None:
        random.seed(seed)

    # length and sum of each sequence
    lena = len(aseq)
    lenb = len(bseq)
    suma = sum(aseq)
    sumb = sum(bseq)

    if not suma == sumb:
        raise networkx.NetworkXError(
            'invalid degree sequences, sum(aseq)!=sum(bseq),%s,%s'
            % (suma, sumb))

    G = _add_nodes_with_bipartite_label(G, lena, lenb)

    if max(aseq) == 0:
        return G  # done if no edges

    # build lists of degree-repeated vertex numbers
    stubs = []
    stubs.extend([[v] * aseq[v] for v in range(0, lena)])
    astubs = []
    astubs = [x for subseq in stubs for x in subseq]

    stubs = []
    stubs.extend([[v] * bseq[v - lena] for v in range(lena, lena + lenb)])
    bstubs = []
    bstubs = [x for subseq in stubs for x in subseq]

    # shuffle lists
    random.shuffle(astubs)
    random.shuffle(bstubs)

    G.add_edges_from([[astubs[i], bstubs[i]] for i in range(suma)])

    G.name = "bipartite_configuration_model"
    return G


def havel_hakimi_graph(aseq, bseq, create_using=None):
    """Return a bipartite graph from two given degree sequences using a
    Havel-Hakimi style construction.

    Nodes from the set A are connected to nodes in the set B by
    connecting the highest degree nodes in set A to the highest degree
    nodes in set B until all stubs are connected.

    Parameters
    ----------
    aseq : list
       Degree sequence for node set A.
    bseq : list
       Degree sequence for node set B.
    create_using : NetworkX graph instance, optional
       Return graph of this type.

    Notes
    -----
    This function is not imported in the main namespace.
    To use it you have to explicitly import the bipartite package.

    The sum of the two sequences must be equal: sum(aseq)=sum(bseq)
    If no graph type is specified use MultiGraph with parallel edges.
    If you want a graph with no parallel edges use create_using=Graph()
    but then the resulting degree sequences might not be exact.

    The nodes are assigned the attribute 'bipartite' with the value 0 or 1
    to indicate which bipartite set the node belongs to.
    """
    if create_using is None:
        create_using = networkx.MultiGraph()
    elif create_using.is_directed():
        raise networkx.NetworkXError(
            "Directed Graph not supported")

    G = networkx.empty_graph(0, create_using)

    # length of the each sequence
    naseq = len(aseq)
    nbseq = len(bseq)

    suma = sum(aseq)
    sumb = sum(bseq)

    if not suma == sumb:
        raise networkx.NetworkXError(
            'invalid degree sequences, sum(aseq)!=sum(bseq),%s,%s'
            % (suma, sumb))

    G = _add_nodes_with_bipartite_label(G, naseq, nbseq)

    if max(aseq) == 0:
        return G  # done if no edges

    # build list of degree-repeated vertex numbers
    astubs = [[aseq[v], v] for v in range(0, naseq)]
    bstubs = [[bseq[v - naseq], v] for v in range(naseq, naseq + nbseq)]
    astubs.sort()
    while astubs:
        (degree, u) = astubs.pop()  # take of largest degree node in the a set
        if degree == 0:
            break  # done, all are zero
        # connect the source to largest degree nodes in the b set
        bstubs.sort()
        for target in bstubs[-degree:]:
            v = target[1]
            G.add_edge(u, v)
            target[0] -= 1  # note this updates bstubs too.
            if target[0] == 0:
                bstubs.remove(target)

    G.name = "bipartite_havel_hakimi_graph"
    return G


def reverse_havel_hakimi_graph(aseq, bseq, create_using=None):
    """Return a bipartite graph from two given degree sequences using a
    Havel-Hakimi style construction.

    Nodes from set A are connected to nodes in the set B by connecting
    the highest degree nodes in set A to the lowest degree nodes in
    set B until all stubs are connected.

    Parameters
    ----------
    aseq : list
       Degree sequence for node set A.
    bseq : list
       Degree sequence for node set B.
    create_using : NetworkX graph instance, optional
       Return graph of this type.

    Notes
    -----
    This function is not imported in the main namespace.
    To use it you have to explicitly import the bipartite package.

    The sum of the two sequences must be equal: sum(aseq)=sum(bseq)
    If no graph type is specified use MultiGraph with parallel edges.
    If you want a graph with no parallel edges use create_using=Graph()
    but then the resulting degree sequences might not be exact.

    The nodes are assigned the attribute 'bipartite' with the value 0 or 1
    to indicate which bipartite set the node belongs to.
    """
    if create_using is None:
        create_using = networkx.MultiGraph()
    elif create_using.is_directed():
        raise networkx.NetworkXError(
            "Directed Graph not supported")

    G = networkx.empty_graph(0, create_using)

    # length of the each sequence
    lena = len(aseq)
    lenb = len(bseq)
    suma = sum(aseq)
    sumb = sum(bseq)

    if not suma == sumb:
        raise networkx.NetworkXError(
            'invalid degree sequences, sum(aseq)!=sum(bseq),%s,%s'
            % (suma, sumb))

    G = _add_nodes_with_bipartite_label(G, lena, lenb)

    if max(aseq) == 0:
        return G  # done if no edges

    # build list of degree-repeated vertex numbers
    astubs = [[aseq[v], v] for v in range(0, lena)]
    bstubs = [[bseq[v - lena], v] for v in range(lena, lena + lenb)]
    astubs.sort()
    bstubs.sort()
    while astubs:
        (degree, u) = astubs.pop()  # take of largest degree node in the a set
        if degree == 0:
            break  # done, all are zero
        # connect the source to the smallest degree nodes in the b set
        for target in bstubs[0:degree]:
            v = target[1]
            G.add_edge(u, v)
            target[0] -= 1  # note this updates bstubs too.
            if target[0] == 0:
                bstubs.remove(target)

    G.name = "bipartite_reverse_havel_hakimi_graph"
    return G


def alternating_havel_hakimi_graph(aseq, bseq, create_using=None):
    """Return a bipartite graph from two given degree sequences using
    an alternating Havel-Hakimi style construction.

    Nodes from the set A are connected to nodes in the set B by
    connecting the highest degree nodes in set A to alternatively the
    highest and the lowest degree nodes in set B until all stubs are
    connected.

    Parameters
    ----------
    aseq : list
       Degree sequence for node set A.
    bseq : list
       Degree sequence for node set B.
    create_using : NetworkX graph instance, optional
       Return graph of this type.

    Notes
    -----
    This function is not imported in the main namespace.
    To use it you have to explicitly import the bipartite package.

    The sum of the two sequences must be equal: sum(aseq)=sum(bseq)
    If no graph type is specified use MultiGraph with parallel edges.
    If you want a graph with no parallel edges use create_using=Graph()
    but then the resulting degree sequences might not be exact.

    The nodes are assigned the attribute 'bipartite' with the value 0 or 1
    to indicate which bipartite set the node belongs to.
    """
    if create_using is None:
        create_using = networkx.MultiGraph()
    elif create_using.is_directed():
        raise networkx.NetworkXError(
            "Directed Graph not supported")

    G = networkx.empty_graph(0, create_using)

    # length of the each sequence
    naseq = len(aseq)
    nbseq = len(bseq)
    suma = sum(aseq)
    sumb = sum(bseq)

    if not suma == sumb:
        raise networkx.NetworkXError(
            'invalid degree sequences, sum(aseq)!=sum(bseq),%s,%s'
            % (suma, sumb))

    G = _add_nodes_with_bipartite_label(G, naseq, nbseq)

    if max(aseq) == 0:
        return G  # done if no edges
    # build list of degree-repeated vertex numbers
    astubs = [[aseq[v], v] for v in range(0, naseq)]
    bstubs = [[bseq[v - naseq], v] for v in range(naseq, naseq + nbseq)]
    while astubs:
        astubs.sort()
        (degree, u) = astubs.pop()  # take of largest degree node in the a set
        if degree == 0:
            break  # done, all are zero
        bstubs.sort()
        small = bstubs[0:degree // 2]  # add these low degree targets
        large = bstubs[(-degree + degree // 2):]  # and these high degree targets
        stubs = [x for z in zip(large, small) for x in z]  # combine, sorry
        if len(stubs) < len(small) + len(large):  # check for zip truncation
            stubs.append(large.pop())
        for target in stubs:
            v = target[1]
            G.add_edge(u, v)
            target[0] -= 1  # note this updates bstubs too.
            if target[0] == 0:
                bstubs.remove(target)

    G.name = "bipartite_alternating_havel_hakimi_graph"
    return G


def preferential_attachment_graph(aseq, p, create_using=None, seed=None):
    """Create a bipartite graph with a preferential attachment model from
    a given single degree sequence.

    Parameters
    ----------
    aseq : list
       Degree sequence for node set A.
    p :  float
       Probability that a new bottom node is added.
    create_using : NetworkX graph instance, optional
       Return graph of this type.
    seed : integer, optional
       Seed for random number generator.

    References
    ----------
    .. [1] Jean-Loup Guillaume and Matthieu Latapy,
       Bipartite structure of all complex networks,
       Inf. Process. Lett. 90, 2004, pg. 215-221
       http://dx.doi.org/10.1016/j.ipl.2004.03.007

    Notes
    -----

    This function is not imported in the main namespace.
    To use it you have to explicitly import the bipartite package.
    """
    if create_using is None:
        create_using = networkx.MultiGraph()
    elif create_using.is_directed():
        raise networkx.NetworkXError(
            "Directed Graph not supported")

    if p > 1:
        raise networkx.NetworkXError("probability %s > 1" % (p))

    G = networkx.empty_graph(0, create_using)

    if not seed is None:
        random.seed(seed)

    naseq = len(aseq)
    G = _add_nodes_with_bipartite_label(G, naseq, 0)
    vv = [[v] * aseq[v] for v in range(0, naseq)]
    while vv:
        while vv[0]:
            source = vv[0][0]
            vv[0].remove(source)
            if random.random() < p or G.number_of_nodes() == naseq:
                target = G.number_of_nodes()
                G.add_node(target, bipartite=1)
                G.add_edge(source, target)
            else:
                bb = [[b] * G.degree(b) for b in range(naseq, G.number_of_nodes())]
                # flatten the list of lists into a list.
                bbstubs = reduce(lambda x, y: x + y, bb)
                # choose preferentially a bottom node.
                target = random.choice(bbstubs)
                G.add_node(target, bipartite=1)
                G.add_edge(source, target)
        vv.remove(vv[0])
    G.name = "bipartite_preferential_attachment_model"
    return G


def random_graph(n, m, p, seed=None, directed=False):
    """Return a bipartite random graph.

    This is a bipartite version of the binomial (Erdős-Rényi) graph.

    Parameters
    ----------
    n : int
        The number of nodes in the first bipartite set.
    m : int
        The number of nodes in the second bipartite set.
    p : float
        Probability for edge creation.
    seed : int, optional
        Seed for random number generator (default=None).
    directed : bool, optional (default=False)
        If True return a directed graph

    Notes
    -----
    This function is not imported in the main namespace.
    To use it you have to explicitly import the bipartite package.

    The bipartite random graph algorithm chooses each of the n*m (undirected)
    or 2*nm (directed) possible edges with probability p.

    This algorithm is $O(n+m)$ where $m$ is the expected number of edges.

    The nodes are assigned the attribute 'bipartite' with the value 0 or 1
    to indicate which bipartite set the node belongs to.

    See Also
    --------
    gnp_random_graph, configuration_model

    References
    ----------
    .. [1] Vladimir Batagelj and Ulrik Brandes,
       "Efficient generation of large random networks",
       Phys. Rev. E, 71, 036113, 2005.
    """
    G = nx.Graph()
    G = _add_nodes_with_bipartite_label(G, n, m)
    if directed:
        G = nx.DiGraph(G)
    G.name = "fast_gnp_random_graph(%s,%s,%s)" % (n, m, p)

    if not seed is None:
        random.seed(seed)

    if p <= 0:
        return G
    if p >= 1:
        return nx.complete_bipartite_graph(n, m)

    lp = math.log(1.0 - p)

    v = 0
    w = -1
    while v < n:
        lr = math.log(1.0 - random.random())
        w = w + 1 + int(lr / lp)
        while w >= m and v < n:
            w = w - m
            v = v + 1
        if v < n:
            G.add_edge(v, n + w)

    if directed:
        # use the same algorithm to
        # add edges from the "m" to "n" set
        v = 0
        w = -1
        while v < n:
            lr = math.log(1.0 - random.random())
            w = w + 1 + int(lr / lp)
            while w >= m and v < n:
                w = w - m
                v = v + 1
            if v < n:
                G.add_edge(n + w, v)

    return G


def gnmk_random_graph(n, m, k, seed=None, directed=False):
    """Return a random bipartite graph G_{n,m,k}.

    Produces a bipartite graph chosen randomly out of the set of all graphs
    with n top nodes, m bottom nodes, and k edges.

    Parameters
    ----------
    n : int
        The number of nodes in the first bipartite set.
    m : int
        The number of nodes in the second bipartite set.
    k : int
        The number of edges
    seed : int, optional
        Seed for random number generator (default=None).
    directed : bool, optional (default=False)
        If True return a directed graph

    Examples
    --------
    from networkx.algorithms import bipartite
    G = bipartite.gnmk_random_graph(10,20,50)

    See Also
    --------
    gnm_random_graph

    Notes
    -----
    This function is not imported in the main namespace.
    To use it you have to explicitly import the bipartite package.

    If k > m * n then a complete bipartite graph is returned.

    This graph is a bipartite version of the `G_{nm}` random graph model.
    """
    G = networkx.Graph()
    G = _add_nodes_with_bipartite_label(G, n, m)
    if directed:
        G = nx.DiGraph(G)
    G.name = "bipartite_gnm_random_graph(%s,%s,%s)" % (n, m, k)
    if seed is not None:
        random.seed(seed)
    if n == 1 or m == 1:
        return G
    max_edges = n * m  # max_edges for bipartite networks
    if k >= max_edges:  # Maybe we should raise an exception here
        return networkx.complete_bipartite_graph(n, m, create_using=G)

    top = [n for n, d in G.nodes(data=True) if d['bipartite'] == 0]
    bottom = list(set(G) - set(top))
    edge_count = 0
    while edge_count < k:
        # generate random edge,u,v
        u = random.choice(top)
        v = random.choice(bottom)
        if v in G[u]:
            continue
        else:
            G.add_edge(u, v)
            edge_count += 1
    return G


def _add_nodes_with_bipartite_label(G, lena, lenb):
    G.add_nodes_from(range(0, lena + lenb))
    b = dict(zip(range(0, lena), [0] * lena))
    b.update(dict(zip(range(lena, lena + lenb), [1] * lenb)))
    nx.set_node_attributes(G, b, 'bipartite')
    return G
