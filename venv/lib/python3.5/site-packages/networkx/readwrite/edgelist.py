"""
**********
Edge Lists
**********
Read and write NetworkX graphs as edge lists.

The multi-line adjacency list format is useful for graphs with nodes
that can be meaningfully represented as strings.  With the edgelist
format simple edge data can be stored but node or graph data is not.
There is no way of representing isolated nodes unless the node has a
self-loop edge.

Format
------
You can read or write three formats of edge lists with these functions.

Node pairs with no data::

 1 2

Python dictionary as data::

 1 2 {'weight':7, 'color':'green'}

Arbitrary data::

 1 2 7 green
"""
__author__ = """Aric Hagberg (hagberg@lanl.gov)\nDan Schult (dschult@colgate.edu)"""
#    Copyright (C) 2004-2018 by
#    Aric Hagberg <hagberg@lanl.gov>
#    Dan Schult <dschult@colgate.edu>
#    Pieter Swart <swart@lanl.gov>
#    All rights reserved.
#    BSD license.

__all__ = ['generate_edgelist',
           'write_edgelist',
           'parse_edgelist',
           'read_edgelist',
           'read_weighted_edgelist',
           'write_weighted_edgelist']

from networkx.utils import open_file, make_str
import networkx as nx


def generate_edgelist(G, delimiter=' ', data=True):
    """Generate a single line of the graph G in edge list format.

    Parameters
    ----------
    G : NetworkX graph

    delimiter : string, optional
       Separator for node labels

    data : bool or list of keys
       If False generate no edge data.  If True use a dictionary
       representation of edge data.  If a list of keys use a list of data
       values corresponding to the keys.

    Returns
    -------
    lines : string
        Lines of data in adjlist format.

    Examples
    --------
    >>> G = nx.lollipop_graph(4, 3)
    >>> G[1][2]['weight'] = 3
    >>> G[3][4]['capacity'] = 12
    >>> for line in nx.generate_edgelist(G, data=False):
    ...     print(line)
    0 1
    0 2
    0 3
    1 2
    1 3
    2 3
    3 4
    4 5
    5 6

    >>> for line in nx.generate_edgelist(G):
    ...     print(line)
    0 1 {}
    0 2 {}
    0 3 {}
    1 2 {'weight': 3}
    1 3 {}
    2 3 {}
    3 4 {'capacity': 12}
    4 5 {}
    5 6 {}

    >>> for line in nx.generate_edgelist(G,data=['weight']):
    ...     print(line)
    0 1
    0 2
    0 3
    1 2 3
    1 3
    2 3
    3 4
    4 5
    5 6

    See Also
    --------
    write_adjlist, read_adjlist
    """
    if data is True:
        for u, v, d in G.edges(data=True):
            e = u, v, dict(d)
            yield delimiter.join(map(make_str, e))
    elif data is False:
        for u, v in G.edges(data=False):
            e = u, v
            yield delimiter.join(map(make_str, e))
    else:
        for u, v, d in G.edges(data=True):
            e = [u, v]
            try:
                e.extend(d[k] for k in data)
            except KeyError:
                pass  # missing data for this edge, should warn?
            yield delimiter.join(map(make_str, e))


@open_file(1, mode='wb')
def write_edgelist(G, path, comments="#", delimiter=' ', data=True,
                   encoding='utf-8'):
    """Write graph as a list of edges.

    Parameters
    ----------
    G : graph
       A NetworkX graph
    path : file or string
       File or filename to write. If a file is provided, it must be
       opened in 'wb' mode. Filenames ending in .gz or .bz2 will be compressed.
    comments : string, optional
       The character used to indicate the start of a comment
    delimiter : string, optional
       The string used to separate values.  The default is whitespace.
    data : bool or list, optional
       If False write no edge data.
       If True write a string representation of the edge data dictionary..
       If a list (or other iterable) is provided, write the  keys specified
       in the list.
    encoding: string, optional
       Specify which encoding to use when writing file.

    Examples
    --------
    >>> G=nx.path_graph(4)
    >>> nx.write_edgelist(G, "test.edgelist")
    >>> G=nx.path_graph(4)
    >>> fh=open("test.edgelist",'wb')
    >>> nx.write_edgelist(G, fh)
    >>> nx.write_edgelist(G, "test.edgelist.gz")
    >>> nx.write_edgelist(G, "test.edgelist.gz", data=False)

    >>> G=nx.Graph()
    >>> G.add_edge(1,2,weight=7,color='red')
    >>> nx.write_edgelist(G,'test.edgelist',data=False)
    >>> nx.write_edgelist(G,'test.edgelist',data=['color'])
    >>> nx.write_edgelist(G,'test.edgelist',data=['color','weight'])

    See Also
    --------
    write_edgelist()
    write_weighted_edgelist()
    """

    for line in generate_edgelist(G, delimiter, data):
        line += '\n'
        path.write(line.encode(encoding))


def parse_edgelist(lines, comments='#', delimiter=None,
                   create_using=None, nodetype=None, data=True):
    """Parse lines of an edge list representation of a graph.

    Parameters
    ----------
    lines : list or iterator of strings
        Input data in edgelist format
    comments : string, optional
       Marker for comment lines
    delimiter : string, optional
       Separator for node labels
    create_using: NetworkX graph container, optional
       Use given NetworkX graph for holding nodes or edges.
    nodetype : Python type, optional
       Convert nodes to this type.
    data : bool or list of (label,type) tuples
       If False generate no edge data or if True use a dictionary
       representation of edge data or a list tuples specifying dictionary
       key names and types for edge data.

    Returns
    -------
    G: NetworkX Graph
        The graph corresponding to lines

    Examples
    --------
    Edgelist with no data:

    >>> lines = ["1 2",
    ...          "2 3",
    ...          "3 4"]
    >>> G = nx.parse_edgelist(lines, nodetype = int)
    >>> list(G)
    [1, 2, 3, 4]
    >>> list(G.edges())
    [(1, 2), (2, 3), (3, 4)]

    Edgelist with data in Python dictionary representation:

    >>> lines = ["1 2 {'weight':3}",
    ...          "2 3 {'weight':27}",
    ...          "3 4 {'weight':3.0}"]
    >>> G = nx.parse_edgelist(lines, nodetype = int)
    >>> list(G)
    [1, 2, 3, 4]
    >>> list(G.edges(data=True))
    [(1, 2, {'weight': 3}), (2, 3, {'weight': 27}), (3, 4, {'weight': 3.0})]

    Edgelist with data in a list:

    >>> lines = ["1 2 3",
    ...          "2 3 27",
    ...          "3 4 3.0"]
    >>> G = nx.parse_edgelist(lines, nodetype = int, data=(('weight',float),))
    >>> list(G)
    [1, 2, 3, 4]
    >>> list(G.edges(data=True))
    [(1, 2, {'weight': 3.0}), (2, 3, {'weight': 27.0}), (3, 4, {'weight': 3.0})]

    See Also
    --------
    read_weighted_edgelist

    """
    from ast import literal_eval
    if create_using is None:
        G = nx.Graph()
    else:
        try:
            G = create_using
            G.clear()
        except:
            raise TypeError("create_using input is not a NetworkX graph type")

    for line in lines:
        p = line.find(comments)
        if p >= 0:
            line = line[:p]
        if not len(line):
            continue
        # split line, should have 2 or more
        s = line.strip().split(delimiter)
        if len(s) < 2:
            continue
        u = s.pop(0)
        v = s.pop(0)
        d = s
        if nodetype is not None:
            try:
                u = nodetype(u)
                v = nodetype(v)
            except:
                raise TypeError("Failed to convert nodes %s,%s to type %s."
                                % (u, v, nodetype))

        if len(d) == 0 or data is False:
            # no data or data type specified
            edgedata = {}
        elif data is True:
            # no edge types specified
            try:  # try to evaluate as dictionary
                edgedata = dict(literal_eval(' '.join(d)))
            except:
                raise TypeError(
                    "Failed to convert edge data (%s) to dictionary." % (d))
        else:
            # convert edge data to dictionary with specified keys and type
            if len(d) != len(data):
                raise IndexError(
                    "Edge data %s and data_keys %s are not the same length" %
                    (d, data))
            edgedata = {}
            for (edge_key, edge_type), edge_value in zip(data, d):
                try:
                    edge_value = edge_type(edge_value)
                except:
                    raise TypeError(
                        "Failed to convert %s data %s to type %s."
                        % (edge_key, edge_value, edge_type))
                edgedata.update({edge_key: edge_value})
        G.add_edge(u, v, **edgedata)
    return G


@open_file(0, mode='rb')
def read_edgelist(path, comments="#", delimiter=None, create_using=None,
                  nodetype=None, data=True, edgetype=None, encoding='utf-8'):
    """Read a graph from a list of edges.

    Parameters
    ----------
    path : file or string
       File or filename to read. If a file is provided, it must be
       opened in 'rb' mode.
       Filenames ending in .gz or .bz2 will be uncompressed.
    comments : string, optional
       The character used to indicate the start of a comment.
    delimiter : string, optional
       The string used to separate values.  The default is whitespace.
    create_using : Graph container, optional,
       Use specified container to build graph.  The default is networkx.Graph,
       an undirected graph.
    nodetype : int, float, str, Python type, optional
       Convert node data from strings to specified type
    data : bool or list of (label,type) tuples
       Tuples specifying dictionary key names and types for edge data
    edgetype : int, float, str, Python type, optional OBSOLETE
       Convert edge data from strings to specified type and use as 'weight'
    encoding: string, optional
       Specify which encoding to use when reading file.

    Returns
    -------
    G : graph
       A networkx Graph or other type specified with create_using

    Examples
    --------
    >>> nx.write_edgelist(nx.path_graph(4), "test.edgelist")
    >>> G=nx.read_edgelist("test.edgelist")

    >>> fh=open("test.edgelist", 'rb')
    >>> G=nx.read_edgelist(fh)
    >>> fh.close()

    >>> G=nx.read_edgelist("test.edgelist", nodetype=int)
    >>> G=nx.read_edgelist("test.edgelist",create_using=nx.DiGraph())

    Edgelist with data in a list:

    >>> textline = '1 2 3'
    >>> fh = open('test.edgelist','w')
    >>> d = fh.write(textline)
    >>> fh.close()
    >>> G = nx.read_edgelist('test.edgelist', nodetype=int, data=(('weight',float),))
    >>> list(G)
    [1, 2]
    >>> list(G.edges(data=True))
    [(1, 2, {'weight': 3.0})]

    See parse_edgelist() for more examples of formatting.

    See Also
    --------
    parse_edgelist

    Notes
    -----
    Since nodes must be hashable, the function nodetype must return hashable
    types (e.g. int, float, str, frozenset - or tuples of those, etc.)
    """
    lines = (line.decode(encoding) for line in path)
    return parse_edgelist(lines, comments=comments, delimiter=delimiter,
                          create_using=create_using, nodetype=nodetype,
                          data=data)


def write_weighted_edgelist(G, path, comments="#",
                            delimiter=' ', encoding='utf-8'):
    """Write graph G as a list of edges with numeric weights.

    Parameters
    ----------
    G : graph
       A NetworkX graph
    path : file or string
       File or filename to write. If a file is provided, it must be
       opened in 'wb' mode.
       Filenames ending in .gz or .bz2 will be compressed.
    comments : string, optional
       The character used to indicate the start of a comment
    delimiter : string, optional
       The string used to separate values.  The default is whitespace.
    encoding: string, optional
       Specify which encoding to use when writing file.

    Examples
    --------
    >>> G=nx.Graph()
    >>> G.add_edge(1,2,weight=7)
    >>> nx.write_weighted_edgelist(G, 'test.weighted.edgelist')

    See Also
    --------
    read_edgelist()
    write_edgelist()
    write_weighted_edgelist()

    """
    write_edgelist(G, path, comments=comments, delimiter=delimiter,
                   data=('weight',), encoding=encoding)


def read_weighted_edgelist(path, comments="#", delimiter=None,
                           create_using=None, nodetype=None, encoding='utf-8'):
    """Read a graph as list of edges with numeric weights.

    Parameters
    ----------
    path : file or string
       File or filename to read. If a file is provided, it must be
       opened in 'rb' mode.
       Filenames ending in .gz or .bz2 will be uncompressed.
    comments : string, optional
       The character used to indicate the start of a comment.
    delimiter : string, optional
       The string used to separate values.  The default is whitespace.
    create_using : Graph container, optional,
       Use specified container to build graph.  The default is networkx.Graph,
       an undirected graph.
    nodetype : int, float, str, Python type, optional
       Convert node data from strings to specified type
    encoding: string, optional
       Specify which encoding to use when reading file.

    Returns
    -------
    G : graph
       A networkx Graph or other type specified with create_using

    Notes
    -----
    Since nodes must be hashable, the function nodetype must return hashable
    types (e.g. int, float, str, frozenset - or tuples of those, etc.)

    Example edgelist file format.

    With numeric edge data::

     # read with
     # >>> G=nx.read_weighted_edgelist(fh)
     # source target data
     a b 1
     a c 3.14159
     d e 42
    """
    return read_edgelist(path,
                         comments=comments,
                         delimiter=delimiter,
                         create_using=create_using,
                         nodetype=nodetype,
                         data=(('weight', float),),
                         encoding=encoding
                         )


# fixture for nose tests
def teardown_module(module):
    import os
    for fname in ['test.edgelist', 'test.edgelist.gz',
                  'test.weighted.edgelist']:
        if os.path.isfile(fname):
            os.unlink(fname)
