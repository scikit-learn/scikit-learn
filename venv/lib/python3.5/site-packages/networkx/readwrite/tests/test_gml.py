#!/usr/bin/env python
# encoding: utf-8

from ast import literal_eval
import codecs
import io
from nose.tools import *
from nose import SkipTest
import networkx as nx
from networkx.readwrite.gml import literal_stringizer, literal_destringizer
import os
import tempfile

try:
    unicode
except NameError:
    unicode = str
try:
    unichr
except NameError:
    unichr = chr


class TestGraph(object):

    def setUp(self):
        self.simple_data = """Creator "me"
Version "xx"
graph [
 comment "This is a sample graph"
 directed 1
 IsPlanar 1
 pos  [ x 0 y 1 ]
 node [
   id 1
   label "Node 1"
   pos [ x 1 y 1 ]
 ]
 node [
    id 2
    pos [ x 1 y 2 ]
    label "Node 2"
    ]
  node [
    id 3
    label "Node 3"
    pos [ x 1 y 3 ]
  ]
  edge [
    source 1
    target 2
    label "Edge from node 1 to node 2"
    color [line "blue" thickness 3]

  ]
  edge [
    source 2
    target 3
    label "Edge from node 2 to node 3"
  ]
  edge [
    source 3
    target 1
    label "Edge from node 3 to node 1"
  ]
]
"""

    def test_parse_gml_cytoscape_bug(self):
        # example from issue #321, originally #324 in trac
        cytoscape_example = """
Creator "Cytoscape"
Version 1.0
graph   [
    node    [
        root_index  -3
        id  -3
        graphics    [
            x   -96.0
            y   -67.0
            w   40.0
            h   40.0
            fill    "#ff9999"
            type    "ellipse"
            outline "#666666"
            outline_width   1.5
        ]
        label   "node2"
    ]
    node    [
        root_index  -2
        id  -2
        graphics    [
            x   63.0
            y   37.0
            w   40.0
            h   40.0
            fill    "#ff9999"
            type    "ellipse"
            outline "#666666"
            outline_width   1.5
        ]
        label   "node1"
    ]
    node    [
        root_index  -1
        id  -1
        graphics    [
            x   -31.0
            y   -17.0
            w   40.0
            h   40.0
            fill    "#ff9999"
            type    "ellipse"
            outline "#666666"
            outline_width   1.5
        ]
        label   "node0"
    ]
    edge    [
        root_index  -2
        target  -2
        source  -1
        graphics    [
            width   1.5
            fill    "#0000ff"
            type    "line"
            Line    [
            ]
            source_arrow    0
            target_arrow    3
        ]
        label   "DirectedEdge"
    ]
    edge    [
        root_index  -1
        target  -1
        source  -3
        graphics    [
            width   1.5
            fill    "#0000ff"
            type    "line"
            Line    [
            ]
            source_arrow    0
            target_arrow    3
        ]
        label   "DirectedEdge"
    ]
]
"""
        nx.parse_gml(cytoscape_example)

    def test_parse_gml(self):
        G = nx.parse_gml(self.simple_data, label='label')
        assert_equals(sorted(G.nodes()),
                      ['Node 1', 'Node 2', 'Node 3'])
        assert_equals([e for e in sorted(G.edges())],
                      [('Node 1', 'Node 2'),
                       ('Node 2', 'Node 3'),
                       ('Node 3', 'Node 1')])

        assert_equals([e for e in sorted(G.edges(data=True))],
                      [('Node 1', 'Node 2',
                        {'color': {'line': 'blue', 'thickness': 3},
                         'label': 'Edge from node 1 to node 2'}),
                       ('Node 2', 'Node 3',
                        {'label': 'Edge from node 2 to node 3'}),
                       ('Node 3', 'Node 1',
                        {'label': 'Edge from node 3 to node 1'})])

    def test_read_gml(self):
        (fd, fname) = tempfile.mkstemp()
        fh = open(fname, 'w')
        fh.write(self.simple_data)
        fh.close()
        Gin = nx.read_gml(fname, label='label')
        G = nx.parse_gml(self.simple_data, label='label')
        assert_equals(sorted(G.nodes(data=True)), sorted(Gin.nodes(data=True)))
        assert_equals(sorted(G.edges(data=True)), sorted(Gin.edges(data=True)))
        os.close(fd)
        os.unlink(fname)

    def test_labels_are_strings(self):
        # GML requires labels to be strings (i.e., in quotes)
        answer = """graph [
  node [
    id 0
    label "1203"
  ]
]"""
        G = nx.Graph()
        G.add_node(1203)
        data = '\n'.join(nx.generate_gml(G, stringizer=literal_stringizer))
        assert_equal(data, answer)

    def test_relabel_duplicate(self):
        data = """
graph
[
        label   ""
        directed        1
        node
        [
                id      0
                label   "same"
        ]
        node
        [
                id      1
                label   "same"
        ]
]
"""
        fh = io.BytesIO(data.encode('UTF-8'))
        fh.seek(0)
        assert_raises(
            nx.NetworkXError, nx.read_gml, fh, label='label')

    def test_tuplelabels(self):
        # https://github.com/networkx/networkx/pull/1048
        # Writing tuple labels to GML failed.
        G = nx.OrderedGraph()
        G.add_edge((0, 1), (1, 0))
        data = '\n'.join(nx.generate_gml(G, stringizer=literal_stringizer))
        answer = """graph [
  node [
    id 0
    label "(0,1)"
  ]
  node [
    id 1
    label "(1,0)"
  ]
  edge [
    source 0
    target 1
  ]
]"""
        assert_equal(data, answer)

    def test_quotes(self):
        # https://github.com/networkx/networkx/issues/1061
        # Encoding quotes as HTML entities.
        G = nx.path_graph(1)
        G.name = "path_graph(1)"
        attr = 'This is "quoted" and this is a copyright: ' + unichr(169)
        G.nodes[0]['demo'] = attr
        fobj = tempfile.NamedTemporaryFile()
        nx.write_gml(G, fobj)
        fobj.seek(0)
        # Should be bytes in 2.x and 3.x
        data = fobj.read().strip().decode('ascii')
        answer = """graph [
  name "path_graph(1)"
  node [
    id 0
    label "0"
    demo "This is &#34;quoted&#34; and this is a copyright: &#169;"
  ]
]"""
        assert_equal(data, answer)

    def test_unicode_node(self):
        node = 'node' + unichr(169)
        G = nx.Graph()
        G.add_node(node)
        fobj = tempfile.NamedTemporaryFile()
        nx.write_gml(G, fobj)
        fobj.seek(0)
        # Should be bytes in 2.x and 3.x
        data = fobj.read().strip().decode('ascii')
        answer = """graph [
  node [
    id 0
    label "node&#169;"
  ]
]"""
        assert_equal(data, answer)

    def test_float_label(self):
        node = 1.0
        G = nx.Graph()
        G.add_node(node)
        fobj = tempfile.NamedTemporaryFile()
        nx.write_gml(G, fobj)
        fobj.seek(0)
        # Should be bytes in 2.x and 3.x
        data = fobj.read().strip().decode('ascii')
        answer = """graph [
  node [
    id 0
    label "1.0"
  ]
]"""
        assert_equal(data, answer)

    def test_name(self):
        G = nx.parse_gml('graph [ name "x" node [ id 0 label "x" ] ]')
        assert_equal('x', G.graph['name'])
        G = nx.parse_gml('graph [ node [ id 0 label "x" ] ]')
        assert_equal('', G.name)
        assert_not_in('name', G.graph)

    def test_graph_types(self):
        for directed in [None, False, True]:
            for multigraph in [None, False, True]:
                gml = 'graph ['
                if directed is not None:
                    gml += ' directed ' + str(int(directed))
                if multigraph is not None:
                    gml += ' multigraph ' + str(int(multigraph))
                gml += ' node [ id 0 label "0" ]'
                gml += ' edge [ source 0 target 0 ]'
                gml += ' ]'
                G = nx.parse_gml(gml)
                assert_equal(bool(directed), G.is_directed())
                assert_equal(bool(multigraph), G.is_multigraph())
                gml = 'graph [\n'
                if directed is True:
                    gml += '  directed 1\n'
                if multigraph is True:
                    gml += '  multigraph 1\n'
                gml += """  node [
    id 0
    label "0"
  ]
  edge [
    source 0
    target 0
"""
                if multigraph:
                    gml += '    key 0\n'
                gml += '  ]\n]'
                assert_equal(gml, '\n'.join(nx.generate_gml(G)))

    def test_data_types(self):
        data = [True, False, 10 ** 20, -2e33, "'", '"&&amp;&&#34;"',
                [{(b'\xfd',): '\x7f', unichr(0x4444): (1, 2)}, (2, "3")]]
        try:
            data.append(unichr(0x14444))  # fails under IronPython
        except ValueError:
            data.append(unichr(0x1444))
        try:
            data.append(literal_eval('{2.3j, 1 - 2.3j, ()}'))  # fails under Python 2.7
        except ValueError:
            data.append([2.3j, 1 - 2.3j, ()])
        G = nx.Graph()
        G.name = data
        G.graph['data'] = data
        G.add_node(0, int=-1, data=dict(data=data))
        G.add_edge(0, 0, float=-2.5, data=data)
        gml = '\n'.join(nx.generate_gml(G, stringizer=literal_stringizer))
        G = nx.parse_gml(gml, destringizer=literal_destringizer)
        assert_equal(data, G.name)
        assert_equal({'name': data, unicode('data'): data}, G.graph)
        assert_equal(list(G.nodes(data=True)),
                     [(0, dict(int=-1, data=dict(data=data)))])
        assert_equal(list(G.edges(data=True)), [(0, 0, dict(float=-2.5, data=data))])
        G = nx.Graph()
        G.graph['data'] = 'frozenset([1, 2, 3])'
        G = nx.parse_gml(nx.generate_gml(G), destringizer=literal_eval)
        assert_equal(G.graph['data'], 'frozenset([1, 2, 3])')

    def test_escape_unescape(self):
        gml = """graph [
  name "&amp;&#34;&#xf;&#x4444;&#1234567890;&#x1234567890abcdef;&unknown;"
]"""
        G = nx.parse_gml(gml)
        assert_equal(
            '&"\x0f' + unichr(0x4444) + '&#1234567890;&#x1234567890abcdef;&unknown;',
            G.name)
        gml = '\n'.join(nx.generate_gml(G))
        assert_equal("""graph [
  name "&#38;&#34;&#15;&#17476;&#38;#1234567890;&#38;#x1234567890abcdef;&#38;unknown;"
]""", gml)

    def test_exceptions(self):
        assert_raises(ValueError, literal_destringizer, '(')
        assert_raises(ValueError, literal_destringizer, 'frozenset([1, 2, 3])')
        assert_raises(ValueError, literal_destringizer, literal_destringizer)
        assert_raises(ValueError, literal_stringizer, frozenset([1, 2, 3]))
        assert_raises(ValueError, literal_stringizer, literal_stringizer)
        with tempfile.TemporaryFile() as f:
            f.write(codecs.BOM_UTF8 + 'graph[]'.encode('ascii'))
            f.seek(0)
            assert_raises(nx.NetworkXError, nx.read_gml, f)

        def assert_parse_error(gml):
            assert_raises(nx.NetworkXError, nx.parse_gml, gml)

        assert_parse_error(['graph [\n\n', unicode(']')])
        assert_parse_error('')
        assert_parse_error('Creator ""')
        assert_parse_error('0')
        assert_parse_error('graph ]')
        assert_parse_error('graph [ 1 ]')
        assert_parse_error('graph [ 1.E+2 ]')
        assert_parse_error('graph [ "A" ]')
        assert_parse_error('graph [ ] graph ]')
        assert_parse_error('graph [ ] graph [ ]')
        assert_parse_error('graph [ data [1, 2, 3] ]')
        assert_parse_error('graph [ node [ ] ]')
        assert_parse_error('graph [ node [ id 0 ] ]')
        nx.parse_gml('graph [ node [ id "a" ] ]', label='id')
        assert_parse_error(
            'graph [ node [ id 0 label 0 ] node [ id 0 label 1 ] ]')
        assert_parse_error(
            'graph [ node [ id 0 label 0 ] node [ id 1 label 0 ] ]')
        assert_parse_error('graph [ node [ id 0 label 0 ] edge [ ] ]')
        assert_parse_error('graph [ node [ id 0 label 0 ] edge [ source 0 ] ]')
        nx.parse_gml(
            'graph [edge [ source 0 target 0 ] node [ id 0 label 0 ] ]')
        assert_parse_error(
            'graph [ node [ id 0 label 0 ] edge [ source 1 target 0 ] ]')
        assert_parse_error(
            'graph [ node [ id 0 label 0 ] edge [ source 0 target 1 ] ]')
        assert_parse_error(
            'graph [ node [ id 0 label 0 ] node [ id 1 label 1 ] '
            'edge [ source 0 target 1 ] edge [ source 1 target 0 ] ]')
        nx.parse_gml(
            'graph [ node [ id 0 label 0 ] node [ id 1 label 1 ] '
            'edge [ source 0 target 1 ] edge [ source 1 target 0 ] '
            'directed 1 ]')
        nx.parse_gml(
            'graph [ node [ id 0 label 0 ] node [ id 1 label 1 ] '
            'edge [ source 0 target 1 ] edge [ source 0 target 1 ]'
            'multigraph 1 ]')
        nx.parse_gml(
            'graph [ node [ id 0 label 0 ] node [ id 1 label 1 ] '
            'edge [ source 0 target 1 key 0 ] edge [ source 0 target 1 ]'
            'multigraph 1 ]')
        assert_parse_error(
            'graph [ node [ id 0 label 0 ] node [ id 1 label 1 ] '
            'edge [ source 0 target 1 key 0 ] edge [ source 0 target 1 key 0 ]'
            'multigraph 1 ]')
        nx.parse_gml(
            'graph [ node [ id 0 label 0 ] node [ id 1 label 1 ] '
            'edge [ source 0 target 1 key 0 ] edge [ source 1 target 0 key 0 ]'
            'directed 1 multigraph 1 ]')

        def assert_generate_error(*args, **kwargs):
            assert_raises(nx.NetworkXError,
                          lambda: list(nx.generate_gml(*args, **kwargs)))

        G = nx.Graph()
        G.graph[3] = 3
        assert_generate_error(G)
        G = nx.Graph()
        G.graph['3'] = 3
        assert_generate_error(G)
        G = nx.Graph()
        G.graph['data'] = frozenset([1, 2, 3])
        assert_generate_error(G, stringizer=literal_stringizer)
        G = nx.Graph()
        G.graph['data'] = []
        assert_generate_error(G)
        assert_generate_error(G, stringizer=len)
