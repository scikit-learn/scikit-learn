#    Copyright (C) 2008-2014 by
#    Aric Hagberg <hagberg@lanl.gov>
#    Dan Schult <dschult@colgate.edu>
#    Pieter Swart <swart@lanl.gov>
#    All rights reserved.
#    BSD license.
#
# Authors: Aric Hagberg (hagberg@lanl.gov)
"""
*****
Pajek
*****
Read graphs in Pajek format.

This implementation handles directed and undirected graphs including
those with self loops and parallel edges.

Format
------
See http://vlado.fmf.uni-lj.si/pub/networks/pajek/doc/draweps.htm
for format information.

"""

import networkx as nx
from networkx.utils import is_string_like, open_file, make_str

__all__ = ['read_pajek', 'parse_pajek', 'generate_pajek', 'write_pajek']


def generate_pajek(G):
    """Generate lines in Pajek graph format.

    Parameters
    ----------
    G : graph
       A Networkx graph

    References
    ----------
    See http://vlado.fmf.uni-lj.si/pub/networks/pajek/doc/draweps.htm
    for format information.
    """
    if G.name == '':
        name = 'NetworkX'
    else:
        name = G.name
    # Apparently many Pajek format readers can't process this line
    # So we'll leave it out for now.
    # yield '*network %s'%name

    # write nodes with attributes
    yield '*vertices %s' % (G.order())
    nodes = list(G)
    # make dictionary mapping nodes to integers
    nodenumber = dict(zip(nodes, range(1, len(nodes) + 1)))
    for n in nodes:
        na = G.nodes.get(n, {})
        x = na.get('x', 0.0)
        y = na.get('y', 0.0)
        id = int(na.get('id', nodenumber[n]))
        nodenumber[n] = id
        shape = na.get('shape', 'ellipse')
        s = ' '.join(map(make_qstr, (id, n, x, y, shape)))
        for k, v in na.items():
            if v.strip() != '':
                s += ' %s %s' % (make_qstr(k), make_qstr(v))
        yield s

    # write edges with attributes
    if G.is_directed():
        yield '*arcs'
    else:
        yield '*edges'
    for u, v, edgedata in G.edges(data=True):
        d = edgedata.copy()
        value = d.pop('weight', 1.0)  # use 1 as default edge value
        s = ' '.join(map(make_qstr, (nodenumber[u], nodenumber[v], value)))
        for k, v in d.items():
            if v.strip() != '':
                s += ' %s %s' % (make_qstr(k), make_qstr(v))
        yield s


@open_file(1, mode='wb')
def write_pajek(G, path, encoding='UTF-8'):
    """Write graph in Pajek format to path.

    Parameters
    ----------
    G : graph
       A Networkx graph
    path : file or string
       File or filename to write.
       Filenames ending in .gz or .bz2 will be compressed.

    Examples
    --------
    >>> G=nx.path_graph(4)
    >>> nx.write_pajek(G, "test.net")

    References
    ----------
    See http://vlado.fmf.uni-lj.si/pub/networks/pajek/doc/draweps.htm
    for format information.
    """
    for line in generate_pajek(G):
        line += '\n'
        path.write(line.encode(encoding))


@open_file(0, mode='rb')
def read_pajek(path, encoding='UTF-8'):
    """Read graph in Pajek format from path.

    Parameters
    ----------
    path : file or string
       File or filename to write.
       Filenames ending in .gz or .bz2 will be uncompressed.

    Returns
    -------
    G : NetworkX MultiGraph or MultiDiGraph.

    Examples
    --------
    >>> G=nx.path_graph(4)
    >>> nx.write_pajek(G, "test.net")
    >>> G=nx.read_pajek("test.net")

    To create a Graph instead of a MultiGraph use

    >>> G1=nx.Graph(G)

    References
    ----------
    See http://vlado.fmf.uni-lj.si/pub/networks/pajek/doc/draweps.htm
    for format information.
    """
    lines = (line.decode(encoding) for line in path)
    return parse_pajek(lines)


def parse_pajek(lines):
    """Parse Pajek format graph from string or iterable.

    Parameters
    ----------
    lines : string or iterable
       Data in Pajek format.

    Returns
    -------
    G : NetworkX graph

    See Also
    --------
    read_pajek()

    """
    import shlex
    # multigraph=False
    if is_string_like(lines):
        lines = iter(lines.split('\n'))
    lines = iter([line.rstrip('\n') for line in lines])
    G = nx.MultiDiGraph()  # are multiedges allowed in Pajek? assume yes
    labels = []  # in the order of the file, needed for matrix
    while lines:
        try:
            l = next(lines)
        except:  # EOF
            break
        if l.lower().startswith("*network"):
            try:
                label, name = l.split(None, 1)
            except ValueError:
                # Line was not of the form:  *network NAME
                pass
            else:
                G.graph['name'] = name
        elif l.lower().startswith("*vertices"):
            nodelabels = {}
            l, nnodes = l.split()
            for i in range(int(nnodes)):
                l = next(lines)
                try:
                    splitline = [x.decode('utf-8') for x in
                                 shlex.split(make_str(l).encode('utf-8'))]
                except AttributeError:
                    splitline = shlex.split(str(l))
                id, label = splitline[0:2]
                labels.append(label)
                G.add_node(label)
                nodelabels[id] = label
                G.nodes[label]['id'] = id
                try:
                    x, y, shape = splitline[2:5]
                    G.nodes[label].update({'x': float(x),
                                           'y': float(y),
                                           'shape': shape})
                except:
                    pass
                extra_attr = zip(splitline[5::2], splitline[6::2])
                G.nodes[label].update(extra_attr)
        elif l.lower().startswith("*edges") or l.lower().startswith("*arcs"):
            if l.lower().startswith("*edge"):
                # switch from multidigraph to multigraph
                G = nx.MultiGraph(G)
            if l.lower().startswith("*arcs"):
                # switch to directed with multiple arcs for each existing edge
                G = G.to_directed()
            for l in lines:
                try:
                    splitline = [x.decode('utf-8') for x in
                                 shlex.split(make_str(l).encode('utf-8'))]
                except AttributeError:
                    splitline = shlex.split(str(l))

                if len(splitline) < 2:
                    continue
                ui, vi = splitline[0:2]
                u = nodelabels.get(ui, ui)
                v = nodelabels.get(vi, vi)
                # parse the data attached to this edge and put in a dictionary
                edge_data = {}
                try:
                    # there should always be a single value on the edge?
                    w = splitline[2:3]
                    edge_data.update({'weight': float(w[0])})
                except:
                    pass
                    # if there isn't, just assign a 1
#                    edge_data.update({'value':1})
                extra_attr = zip(splitline[3::2], splitline[4::2])
                edge_data.update(extra_attr)
                # if G.has_edge(u,v):
                #     multigraph=True
                G.add_edge(u, v, **edge_data)
        elif l.lower().startswith("*matrix"):
            G = nx.DiGraph(G)
            adj_list = ((labels[row], labels[col], {'weight': int(data)})
                        for (row, line) in enumerate(lines)
                        for (col, data) in enumerate(line.split())
                        if int(data) != 0)
            G.add_edges_from(adj_list)

    return G


def make_qstr(t):
    """Return the string representation of t.
    Add outer double-quotes if the string has a space.
    """
    if not is_string_like(t):
        t = str(t)
    if " " in t:
        t = r'"%s"' % t
    return t


# fixture for nose tests
def teardown_module(module):
    import os
    os.unlink('test.net')
