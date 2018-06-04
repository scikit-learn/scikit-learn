"""Functions to convert NetworkX graphs to and from other formats.

The preferred way of converting data to a NetworkX graph is through the
graph constuctor.  The constructor calls the to_networkx_graph() function
which attempts to guess the input type and convert it automatically.

Examples
--------
Create a graph with a single edge from a dictionary of dictionaries

>>> d={0: {1: 1}} # dict-of-dicts single edge (0,1)
>>> G=nx.Graph(d)

See Also
--------
nx_agraph, nx_pydot
"""
#    Copyright (C) 2006-2013 by
#    Aric Hagberg <hagberg@lanl.gov>
#    Dan Schult <dschult@colgate.edu>
#    Pieter Swart <swart@lanl.gov>
#    All rights reserved.
#    BSD license.
import warnings
import networkx as nx
__author__ = """\n""".join(['Aric Hagberg <aric.hagberg@gmail.com>',
                            'Pieter Swart (swart@lanl.gov)',
                            'Dan Schult(dschult@colgate.edu)'])
__all__ = ['to_networkx_graph',
           'from_dict_of_dicts', 'to_dict_of_dicts',
           'from_dict_of_lists', 'to_dict_of_lists',
           'from_edgelist', 'to_edgelist']


def _prep_create_using(create_using):
    """Return a graph object ready to be populated.

    If create_using is None return the default (just networkx.Graph())
    If create_using.clear() works, assume it returns a graph object.
    Otherwise raise an exception because create_using is not a networkx graph.

    """
    if create_using is None:
        return nx.Graph()
    try:
        create_using.clear()
    except:
        raise TypeError("Input graph is not a networkx graph type")
    return create_using


def to_networkx_graph(data, create_using=None, multigraph_input=False):
    """Make a NetworkX graph from a known data structure.

    The preferred way to call this is automatically
    from the class constructor

    >>> d = {0: {1: {'weight':1}}} # dict-of-dicts single edge (0,1)
    >>> G = nx.Graph(d)

    instead of the equivalent

    >>> G = nx.from_dict_of_dicts(d)

    Parameters
    ----------
    data : object to be converted

       Current known types are:
         any NetworkX graph
         dict-of-dicts
         dict-of-lists
         list of edges
         Pandas DataFrame (row per edge)
         numpy matrix
         numpy ndarray
         scipy sparse matrix
         pygraphviz agraph

    create_using : NetworkX graph
       Use specified graph for result.  Otherwise a new graph is created.

    multigraph_input : bool (default False)
      If True and  data is a dict_of_dicts,
      try to create a multigraph assuming dict_of_dict_of_lists.
      If data and create_using are both multigraphs then create
      a multigraph from a multigraph.

    """
    # NX graph
    if hasattr(data, "adj"):
        try:
            result = from_dict_of_dicts(data.adj,
                                        create_using=create_using,
                                        multigraph_input=data.is_multigraph())
            if hasattr(data, 'graph'):  # data.graph should be dict-like
                result.graph.update(data.graph)
            if hasattr(data, 'nodes'):  # data.nodes should be dict-like
                result._node.update((n, dd.copy()) for n, dd in data.nodes.items())
            return result
        except:
            raise nx.NetworkXError("Input is not a correct NetworkX graph.")

    # pygraphviz  agraph
    if hasattr(data, "is_strict"):
        try:
            return nx.nx_agraph.from_agraph(data, create_using=create_using)
        except:
            raise nx.NetworkXError("Input is not a correct pygraphviz graph.")

    # dict of dicts/lists
    if isinstance(data, dict):
        try:
            return from_dict_of_dicts(data, create_using=create_using,
                                      multigraph_input=multigraph_input)
        except:
            try:
                return from_dict_of_lists(data, create_using=create_using)
            except:
                raise TypeError("Input is not known type.")

    # list or generator of edges

    if (isinstance(data, (list, tuple)) or
            any(hasattr(data, attr) for attr in ['_adjdict', 'next', '__next__'])):
        try:
            return from_edgelist(data, create_using=create_using)
        except:
            raise nx.NetworkXError("Input is not a valid edge list")

    # Pandas DataFrame
    try:
        import pandas as pd
        if isinstance(data, pd.DataFrame):
            if data.shape[0] == data.shape[1]:
                try:
                    return nx.from_pandas_adjacency(data, create_using=create_using)
                except:
                    msg = "Input is not a correct Pandas DataFrame adjacency matrix."
                    raise nx.NetworkXError(msg)
            else:
                try:
                    return nx.from_pandas_edgelist(data, edge_attr=True, create_using=create_using)
                except:
                    msg = "Input is not a correct Pandas DataFrame edge-list."
                    raise nx.NetworkXError(msg)
    except ImportError:
        msg = 'pandas not found, skipping conversion test.'
        warnings.warn(msg, ImportWarning)

    # numpy matrix or ndarray
    try:
        import numpy
        if isinstance(data, (numpy.matrix, numpy.ndarray)):
            try:
                return nx.from_numpy_matrix(data, create_using=create_using)
            except:
                raise nx.NetworkXError(
                    "Input is not a correct numpy matrix or array.")
    except ImportError:
        warnings.warn('numpy not found, skipping conversion test.',
                      ImportWarning)

    # scipy sparse matrix - any format
    try:
        import scipy
        if hasattr(data, "format"):
            try:
                return nx.from_scipy_sparse_matrix(data, create_using=create_using)
            except:
                raise nx.NetworkXError(
                    "Input is not a correct scipy sparse matrix type.")
    except ImportError:
        warnings.warn('scipy not found, skipping conversion test.',
                      ImportWarning)

    raise nx.NetworkXError(
        "Input is not a known data type for conversion.")


def to_dict_of_lists(G, nodelist=None):
    """Return adjacency representation of graph as a dictionary of lists.

    Parameters
    ----------
    G : graph
       A NetworkX graph

    nodelist : list
       Use only nodes specified in nodelist

    Notes
    -----
    Completely ignores edge data for MultiGraph and MultiDiGraph.

    """
    if nodelist is None:
        nodelist = G

    d = {}
    for n in nodelist:
        d[n] = [nbr for nbr in G.neighbors(n) if nbr in nodelist]
    return d


def from_dict_of_lists(d, create_using=None):
    """Return a graph from a dictionary of lists.

    Parameters
    ----------
    d : dictionary of lists
      A dictionary of lists adjacency representation.

    create_using : NetworkX graph
       Use specified graph for result.  Otherwise a new graph is created.

    Examples
    --------
    >>> dol = {0: [1]} # single edge (0,1)
    >>> G = nx.from_dict_of_lists(dol)

    or

    >>> G = nx.Graph(dol) # use Graph constructor

    """
    G = _prep_create_using(create_using)
    G.add_nodes_from(d)
    if G.is_multigraph() and not G.is_directed():
        # a dict_of_lists can't show multiedges.  BUT for undirected graphs,
        # each edge shows up twice in the dict_of_lists.
        # So we need to treat this case separately.
        seen = {}
        for node, nbrlist in d.items():
            for nbr in nbrlist:
                if nbr not in seen:
                    G.add_edge(node, nbr)
            seen[node] = 1  # don't allow reverse edge to show up
    else:
        G.add_edges_from(((node, nbr) for node, nbrlist in d.items()
                          for nbr in nbrlist))
    return G


def to_dict_of_dicts(G, nodelist=None, edge_data=None):
    """Return adjacency representation of graph as a dictionary of dictionaries.

    Parameters
    ----------
    G : graph
       A NetworkX graph

    nodelist : list
       Use only nodes specified in nodelist

    edge_data : list, optional
       If provided,  the value of the dictionary will be
       set to edge_data for all edges.  This is useful to make
       an adjacency matrix type representation with 1 as the edge data.
       If edgedata is None, the edgedata in G is used to fill the values.
       If G is a multigraph, the edgedata is a dict for each pair (u,v).
    """
    dod = {}
    if nodelist is None:
        if edge_data is None:
            for u, nbrdict in G.adjacency():
                dod[u] = nbrdict.copy()
        else:  # edge_data is not None
            for u, nbrdict in G.adjacency():
                dod[u] = dod.fromkeys(nbrdict, edge_data)
    else:  # nodelist is not None
        if edge_data is None:
            for u in nodelist:
                dod[u] = {}
                for v, data in ((v, data) for v, data in G[u].items() if v in nodelist):
                    dod[u][v] = data
        else:  # nodelist and edge_data are not None
            for u in nodelist:
                dod[u] = {}
                for v in (v for v in G[u] if v in nodelist):
                    dod[u][v] = edge_data
    return dod


def from_dict_of_dicts(d, create_using=None, multigraph_input=False):
    """Return a graph from a dictionary of dictionaries.

    Parameters
    ----------
    d : dictionary of dictionaries
      A dictionary of dictionaries adjacency representation.

    create_using : NetworkX graph
       Use specified graph for result.  Otherwise a new graph is created.

    multigraph_input : bool (default False)
       When True, the values of the inner dict are assumed
       to be containers of edge data for multiple edges.
       Otherwise this routine assumes the edge data are singletons.

    Examples
    --------
    >>> dod = {0: {1: {'weight': 1}}} # single edge (0,1)
    >>> G = nx.from_dict_of_dicts(dod)

    or

    >>> G = nx.Graph(dod) # use Graph constructor

    """
    G = _prep_create_using(create_using)
    G.add_nodes_from(d)
    # is dict a MultiGraph or MultiDiGraph?
    if multigraph_input:
        # make a copy of the list of edge data (but not the edge data)
        if G.is_directed():
            if G.is_multigraph():
                G.add_edges_from((u, v, key, data)
                                 for u, nbrs in d.items()
                                 for v, datadict in nbrs.items()
                                 for key, data in datadict.items())
            else:
                G.add_edges_from((u, v, data)
                                 for u, nbrs in d.items()
                                 for v, datadict in nbrs.items()
                                 for key, data in datadict.items())
        else:  # Undirected
            if G.is_multigraph():
                seen = set()   # don't add both directions of undirected graph
                for u, nbrs in d.items():
                    for v, datadict in nbrs.items():
                        if (u, v) not in seen:
                            G.add_edges_from((u, v, key, data)
                                             for key, data in datadict.items())
                            seen.add((v, u))
            else:
                seen = set()   # don't add both directions of undirected graph
                for u, nbrs in d.items():
                    for v, datadict in nbrs.items():
                        if (u, v) not in seen:
                            G.add_edges_from((u, v, data)
                                             for key, data in datadict.items())
                            seen.add((v, u))

    else:  # not a multigraph to multigraph transfer
        if G.is_multigraph() and not G.is_directed():
            # d can have both representations u-v, v-u in dict.  Only add one.
            # We don't need this check for digraphs since we add both directions,
            # or for Graph() since it is done implicitly (parallel edges not allowed)
            seen = set()
            for u, nbrs in d.items():
                for v, data in nbrs.items():
                    if (u, v) not in seen:
                        G.add_edge(u, v, key=0)
                        G[u][v][0].update(data)
                    seen.add((v, u))
        else:
            G.add_edges_from(((u, v, data)
                              for u, nbrs in d.items()
                              for v, data in nbrs.items()))
    return G


def to_edgelist(G, nodelist=None):
    """Return a list of edges in the graph.

    Parameters
    ----------
    G : graph
       A NetworkX graph

    nodelist : list
       Use only nodes specified in nodelist

    """
    if nodelist is None:
        return G.edges(data=True)
    return G.edges(nodelist, data=True)


def from_edgelist(edgelist, create_using=None):
    """Return a graph from a list of edges.

    Parameters
    ----------
    edgelist : list or iterator
      Edge tuples

    create_using : NetworkX graph
       Use specified graph for result.  Otherwise a new graph is created.

    Examples
    --------
    >>> edgelist = [(0, 1)] # single edge (0,1)
    >>> G = nx.from_edgelist(edgelist)

    or

    >>> G = nx.Graph(edgelist) # use Graph constructor

    """
    G = _prep_create_using(create_using)
    G.add_edges_from(edgelist)
    return G
