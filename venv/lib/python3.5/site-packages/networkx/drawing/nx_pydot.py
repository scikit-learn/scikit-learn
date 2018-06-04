"""
*****
Pydot
*****

Import and export NetworkX graphs in Graphviz dot format using pydot.

Either this module or nx_agraph can be used to interface with graphviz.

See Also
--------
pydot:         https://github.com/erocarrera/pydot
Graphviz:      http://www.research.att.com/sw/tools/graphviz/
DOT Language:  http://www.graphviz.org/doc/info/lang.html
"""
# Author: Aric Hagberg (aric.hagberg@gmail.com)

#    Copyright (C) 2004-2018 by
#    Aric Hagberg <hagberg@lanl.gov>
#    Dan Schult <dschult@colgate.edu>
#    Pieter Swart <swart@lanl.gov>
#    Cecil Curry <leycec@gmail.com>
#    All rights reserved.
#    BSD license.
from locale import getpreferredencoding
from networkx.utils import open_file, make_str
from pkg_resources import parse_version
import networkx as nx

__all__ = ['write_dot', 'read_dot', 'graphviz_layout', 'pydot_layout',
           'to_pydot', 'from_pydot']

# Minimum required version of pydot, which broke backwards API compatibility in
# non-trivial ways and is thus a hard NetworkX requirement. Note that, although
# pydot 1.2.0 was the first to do so, pydot 1.2.3 resolves a critical long-
# standing Python 2.x issue required for sane NetworkX operation. See also:
#     https://github.com/erocarrera/pydot/blob/master/ChangeLog
PYDOT_VERSION_MIN = '1.2.3'

# 2.x/3.x compatibility
try:
    basestring
except NameError:
    basestring = str
    unicode = str


@open_file(1, mode='w')
def write_dot(G, path):
    """Write NetworkX graph G to Graphviz dot format on path.

    Path can be a string or a file handle.
    """
    P = to_pydot(G)
    path.write(P.to_string())
    return


@open_file(0, mode='r')
def read_dot(path):
    """Return a NetworkX :class:`MultiGraph` or :class:`MultiDiGraph` from the
    dot file with the passed path.

    If this file contains multiple graphs, only the first such graph is
    returned. All graphs _except_ the first are silently ignored.

    Parameters
    ----------
    path : str or file
        Filename or file handle.

    Returns
    -------
    G : MultiGraph or MultiDiGraph
        A :class:`MultiGraph` or :class:`MultiDiGraph`.

    Notes
    -----
    Use `G = nx.Graph(read_dot(path))` to return a :class:`Graph` instead of a
    :class:`MultiGraph`.
    """
    pydot = _import_pydot()
    data = path.read()

    # List of one or more "pydot.Dot" instances deserialized from this file.
    P_list = pydot.graph_from_dot_data(data)

    # Convert only the first such instance into a NetworkX graph.
    return from_pydot(P_list[0])


def from_pydot(P):
    """Return a NetworkX graph from a Pydot graph.

    Parameters
    ----------
    P : Pydot graph
      A graph created with Pydot

    Returns
    -------
    G : NetworkX multigraph
        A MultiGraph or MultiDiGraph.

    Examples
    --------
    >>> K5 = nx.complete_graph(5)
    >>> A = nx.nx_pydot.to_pydot(K5)
    >>> G = nx.nx_pydot.from_pydot(A) # return MultiGraph

    # make a Graph instead of MultiGraph
    >>> G = nx.Graph(nx.nx_pydot.from_pydot(A))

    """
    if P.get_strict(None):  # pydot bug: get_strict() shouldn't take argument
        multiedges = False
    else:
        multiedges = True

    if P.get_type() == 'graph':  # undirected
        if multiedges:
            N = nx.MultiGraph()
        else:
            N = nx.Graph()
    else:
        if multiedges:
            N = nx.MultiDiGraph()
        else:
            N = nx.DiGraph()

    # assign defaults
    name = P.get_name().strip('"')
    if name != '':
        N.name = name

    # add nodes, attributes to N.node_attr
    for p in P.get_node_list():
        n = p.get_name().strip('"')
        if n in ('node', 'graph', 'edge'):
            continue
        N.add_node(n, **p.get_attributes())

    # add edges
    for e in P.get_edge_list():
        u = e.get_source()
        v = e.get_destination()
        attr = e.get_attributes()
        s = []
        d = []

        if isinstance(u, basestring):
            s.append(u.strip('"'))
        else:
            for unodes in u['nodes']:
                s.append(unodes.strip('"'))

        if isinstance(v, basestring):
            d.append(v.strip('"'))
        else:
            for vnodes in v['nodes']:
                d.append(vnodes.strip('"'))

        for source_node in s:
            for destination_node in d:
                N.add_edge(source_node, destination_node, **attr)

    # add default attributes for graph, nodes, edges
    pattr = P.get_attributes()
    if pattr:
        N.graph['graph'] = pattr
    try:
        N.graph['node'] = P.get_node_defaults()[0]
    except:  # IndexError,TypeError:
        pass  # N.graph['node']={}
    try:
        N.graph['edge'] = P.get_edge_defaults()[0]
    except:  # IndexError,TypeError:
        pass  # N.graph['edge']={}
    return N


def to_pydot(N):
    """Return a pydot graph from a NetworkX graph N.

    Parameters
    ----------
    N : NetworkX graph
      A graph created with NetworkX

    Examples
    --------
    >>> K5 = nx.complete_graph(5)
    >>> P = nx.nx_pydot.to_pydot(K5)

    Notes
    -----

    """
    pydot = _import_pydot()

    # set Graphviz graph type
    if N.is_directed():
        graph_type = 'digraph'
    else:
        graph_type = 'graph'
    strict = nx.number_of_selfloops(N) == 0 and not N.is_multigraph()

    name = N.name
    graph_defaults = N.graph.get('graph', {})
    if name is '':
        P = pydot.Dot('', graph_type=graph_type, strict=strict,
                      **graph_defaults)
    else:
        P = pydot.Dot('"%s"' % name, graph_type=graph_type, strict=strict,
                      **graph_defaults)
    try:
        P.set_node_defaults(**N.graph['node'])
    except KeyError:
        pass
    try:
        P.set_edge_defaults(**N.graph['edge'])
    except KeyError:
        pass

    for n, nodedata in N.nodes(data=True):
        str_nodedata = dict((k, make_str(v)) for k, v in nodedata.items())
        p = pydot.Node(make_str(n), **str_nodedata)
        P.add_node(p)

    if N.is_multigraph():
        for u, v, key, edgedata in N.edges(data=True, keys=True):
            str_edgedata = dict((k, make_str(v)) for k, v in edgedata.items() if k != 'key')
            edge = pydot.Edge(make_str(u), make_str(v),
                              key=make_str(key), **str_edgedata)
            P.add_edge(edge)

    else:
        for u, v, edgedata in N.edges(data=True):
            str_edgedata = dict((k, make_str(v)) for k, v in edgedata.items())
            edge = pydot.Edge(make_str(u), make_str(v), **str_edgedata)
            P.add_edge(edge)
    return P


def graphviz_layout(G, prog='neato', root=None, **kwds):
    """Create node positions using Pydot and Graphviz.

    Returns a dictionary of positions keyed by node.

    Examples
    --------
    >>> G = nx.complete_graph(4)
    >>> pos = nx.nx_pydot.graphviz_layout(G)
    >>> pos = nx.nx_pydot.graphviz_layout(G, prog='dot')

    Notes
    -----
    This is a wrapper for pydot_layout.
    """
    return pydot_layout(G=G, prog=prog, root=root, **kwds)


# FIXME: Document the "root" parameter.
# FIXME: Why does this function accept a variadic dictionary of keyword arguments
# (i.e., "**kwds") but fail to do anything with those arguments? This is probably
# wrong, as unrecognized keyword arguments will be silently ignored.
def pydot_layout(G, prog='neato', root=None, **kwds):
    """Create node positions using :mod:`pydot` and Graphviz.

    Parameters
    --------
    G : Graph
        NetworkX graph to be laid out.
    prog : optional[str]
        Basename of the GraphViz command with which to layout this graph.
        Defaults to `neato`, the default GraphViz command for undirected graphs.

    Returns
    --------
    dict
        Dictionary of positions keyed by node.

    Examples
    --------
    >>> G = nx.complete_graph(4)
    >>> pos = nx.nx_pydot.pydot_layout(G)
    >>> pos = nx.nx_pydot.pydot_layout(G, prog='dot')
    """
    pydot = _import_pydot()
    P = to_pydot(G)
    if root is not None:
        P.set("root", make_str(root))

    # List of low-level bytes comprising a string in the dot language converted
    # from the passed graph with the passed external GraphViz command.
    D_bytes = P.create_dot(prog=prog)

    # Unique string decoded from these bytes with the preferred locale encoding.
    D = unicode(D_bytes, encoding=getpreferredencoding())

    if D == "":  # no data returned
        print("Graphviz layout with %s failed" % (prog))
        print()
        print("To debug what happened try:")
        print("P = nx.nx_pydot.to_pydot(G)")
        print("P.write_dot(\"file.dot\")")
        print("And then run %s on file.dot" % (prog))
        return

    # List of one or more "pydot.Dot" instances deserialized from this string.
    Q_list = pydot.graph_from_dot_data(D)
    assert len(Q_list) == 1

    # The first and only such instance, as guaranteed by the above assertion.
    Q = Q_list[0]

    node_pos = {}
    for n in G.nodes():
        pydot_node = pydot.Node(make_str(n)).get_name()
        node = Q.get_node(pydot_node)

        if isinstance(node, list):
            node = node[0]
        pos = node.get_pos()[1:-1]  # strip leading and trailing double quotes
        if pos is not None:
            xx, yy = pos.split(",")
            node_pos[n] = (float(xx), float(yy))
    return node_pos


def _import_pydot():
    '''
    Import and return the `pydot` module if the currently installed version of
    this module satisfies NetworkX requirements _or_ raise an exception.

    Returns
    --------
    :mod:`pydot`
        Imported `pydot` module object.

    Raises
    --------
    ImportError
        If the `pydot` module is either unimportable _or_ importable but of
        insufficient version.
    '''

    import pydot

    # If the currently installed version of pydot is older than this minimum,
    # raise an exception. The pkg_resources.parse_version() function bundled
    # with setuptools is commonly regarded to be the most robust means of
    # comparing version strings. (Your mileage may vary.)
    if parse_version(pydot.__version__) < parse_version(PYDOT_VERSION_MIN):
        raise ImportError(
            'pydot %s < %s' % (pydot.__version__, PYDOT_VERSION_MIN))

    return pydot

# fixture for nose tests


def setup_module(module):
    from nose import SkipTest
    try:
        return _import_pydot()
    except ImportError:
        raise SkipTest("pydot not available")
