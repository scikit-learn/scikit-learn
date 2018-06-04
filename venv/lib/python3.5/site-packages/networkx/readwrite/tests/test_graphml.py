#!/usr/bin/env python
from nose.tools import *
from nose import SkipTest
import networkx as nx
from networkx.testing.utils import *
import io
import tempfile
import os


class BaseGraphML(object):
    def setUp(self):
        self.simple_directed_data = """<?xml version="1.0" encoding="UTF-8"?>
<!-- This file was written by the JAVA GraphML Library.-->
<graphml xmlns="http://graphml.graphdrawing.org/xmlns"
         xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
         xsi:schemaLocation="http://graphml.graphdrawing.org/xmlns
         http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd">
  <graph id="G" edgedefault="directed">
    <node id="n0"/>
    <node id="n1"/>
    <node id="n2"/>
    <node id="n3"/>
    <node id="n4"/>
    <node id="n5"/>
    <node id="n6"/>
    <node id="n7"/>
    <node id="n8"/>
    <node id="n9"/>
    <node id="n10"/>
    <edge id="foo" source="n0" target="n2"/>
    <edge source="n1" target="n2"/>
    <edge source="n2" target="n3"/>
    <edge source="n3" target="n5"/>
    <edge source="n3" target="n4"/>
    <edge source="n4" target="n6"/>
    <edge source="n6" target="n5"/>
    <edge source="n5" target="n7"/>
    <edge source="n6" target="n8"/>
    <edge source="n8" target="n7"/>
    <edge source="n8" target="n9"/>
  </graph>
</graphml>"""
        self.simple_directed_graph = nx.DiGraph()
        self.simple_directed_graph.add_node('n10')
        self.simple_directed_graph.add_edge('n0', 'n2', id='foo')
        self.simple_directed_graph.add_edges_from([('n1', 'n2'),
                                                   ('n2', 'n3'),
                                                   ('n3', 'n5'),
                                                   ('n3', 'n4'),
                                                   ('n4', 'n6'),
                                                   ('n6', 'n5'),
                                                   ('n5', 'n7'),
                                                   ('n6', 'n8'),
                                                   ('n8', 'n7'),
                                                   ('n8', 'n9'),
                                                   ])
        self.simple_directed_fh = \
            io.BytesIO(self.simple_directed_data.encode('UTF-8'))

        self.attribute_data = """<?xml version="1.0" encoding="UTF-8"?>
<graphml xmlns="http://graphml.graphdrawing.org/xmlns"
      xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
      xsi:schemaLocation="http://graphml.graphdrawing.org/xmlns
        http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd">
  <key id="d0" for="node" attr.name="color" attr.type="string">
    <default>yellow</default>
  </key>
  <key id="d1" for="edge" attr.name="weight" attr.type="double"/>
  <graph id="G" edgedefault="directed">
    <node id="n0">
      <data key="d0">green</data>
    </node>
    <node id="n1"/>
    <node id="n2">
      <data key="d0">blue</data>
    </node>
    <node id="n3">
      <data key="d0">red</data>
    </node>
    <node id="n4"/>
    <node id="n5">
      <data key="d0">turquoise</data>
    </node>
    <edge id="e0" source="n0" target="n2">
      <data key="d1">1.0</data>
    </edge>
    <edge id="e1" source="n0" target="n1">
      <data key="d1">1.0</data>
    </edge>
    <edge id="e2" source="n1" target="n3">
      <data key="d1">2.0</data>
    </edge>
    <edge id="e3" source="n3" target="n2"/>
    <edge id="e4" source="n2" target="n4"/>
    <edge id="e5" source="n3" target="n5"/>
    <edge id="e6" source="n5" target="n4">
      <data key="d1">1.1</data>
    </edge>
  </graph>
</graphml>
"""
        self.attribute_graph = nx.DiGraph(id='G')
        self.attribute_graph.graph['node_default'] = {'color': 'yellow'}
        self.attribute_graph.add_node('n0', color='green')
        self.attribute_graph.add_node('n2', color='blue')
        self.attribute_graph.add_node('n3', color='red')
        self.attribute_graph.add_node('n4')
        self.attribute_graph.add_node('n5', color='turquoise')
        self.attribute_graph.add_edge('n0', 'n2', id='e0', weight=1.0)
        self.attribute_graph.add_edge('n0', 'n1', id='e1', weight=1.0)
        self.attribute_graph.add_edge('n1', 'n3', id='e2', weight=2.0)
        self.attribute_graph.add_edge('n3', 'n2', id='e3')
        self.attribute_graph.add_edge('n2', 'n4', id='e4')
        self.attribute_graph.add_edge('n3', 'n5', id='e5')
        self.attribute_graph.add_edge('n5', 'n4', id='e6', weight=1.1)
        self.attribute_fh = io.BytesIO(self.attribute_data.encode('UTF-8'))

        self.attribute_numeric_type_data = """<?xml version='1.0' encoding='utf-8'?>
<graphml xmlns="http://graphml.graphdrawing.org/xmlns"
         xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
         xsi:schemaLocation="http://graphml.graphdrawing.org/xmlns
         http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd">
  <key attr.name="weight" attr.type="double" for="node" id="d1" />
  <key attr.name="weight" attr.type="double" for="edge" id="d0" />
  <graph edgedefault="directed">
    <node id="n0">
      <data key="d1">1</data>
    </node>
    <node id="n1">
      <data key="d1">2.0</data>
    </node>
    <edge source="n0" target="n1">
      <data key="d0">1</data>
    </edge>
    <edge source="n1" target="n0">
      <data key="d0">k</data>
    </edge>
    <edge source="n1" target="n1">
      <data key="d0">1.0</data>
    </edge>
  </graph>
</graphml>
"""
        self.attribute_numeric_type_graph = nx.DiGraph()
        self.attribute_numeric_type_graph.add_node('n0', weight=1)
        self.attribute_numeric_type_graph.add_node('n1', weight=2.0)
        self.attribute_numeric_type_graph.add_edge('n0', 'n1', weight=1)
        self.attribute_numeric_type_graph.add_edge('n1', 'n1', weight=1.0)
        fh = io.BytesIO(self.attribute_numeric_type_data.encode('UTF-8'))
        self.attribute_numeric_type_fh = fh

        self.simple_undirected_data = """<?xml version="1.0" encoding="UTF-8"?>
<graphml xmlns="http://graphml.graphdrawing.org/xmlns"
         xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
         xsi:schemaLocation="http://graphml.graphdrawing.org/xmlns
         http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd">
  <graph id="G">
    <node id="n0"/>
    <node id="n1"/>
    <node id="n2"/>
    <node id="n10"/>
    <edge id="foo" source="n0" target="n2"/>
    <edge source="n1" target="n2"/>
    <edge source="n2" target="n3"/>
  </graph>
</graphml>"""
#    <edge source="n8" target="n10" directed="false"/>
        self.simple_undirected_graph = nx.Graph()
        self.simple_undirected_graph.add_node('n10')
        self.simple_undirected_graph.add_edge('n0', 'n2', id='foo')
        self.simple_undirected_graph.add_edges_from([('n1', 'n2'),
                                                     ('n2', 'n3'),
                                                     ])
        fh = io.BytesIO(self.simple_undirected_data.encode('UTF-8'))
        self.simple_undirected_fh = fh


class TestReadGraphML(BaseGraphML):
    def test_read_simple_directed_graphml(self):
        G = self.simple_directed_graph
        H = nx.read_graphml(self.simple_directed_fh)
        assert_equal(sorted(G.nodes()), sorted(H.nodes()))
        assert_equal(sorted(G.edges()), sorted(H.edges()))
        assert_equal(sorted(G.edges(data=True)),
                     sorted(H.edges(data=True)))
        self.simple_directed_fh.seek(0)

        I = nx.parse_graphml(self.simple_directed_data)
        assert_equal(sorted(G.nodes()), sorted(I.nodes()))
        assert_equal(sorted(G.edges()), sorted(I.edges()))
        assert_equal(sorted(G.edges(data=True)),
                     sorted(I.edges(data=True)))

    def test_read_simple_undirected_graphml(self):
        G = self.simple_undirected_graph
        H = nx.read_graphml(self.simple_undirected_fh)
        assert_nodes_equal(G.nodes(), H.nodes())
        assert_edges_equal(G.edges(), H.edges())
        self.simple_undirected_fh.seek(0)

        I = nx.parse_graphml(self.simple_undirected_data)
        assert_nodes_equal(G.nodes(), I.nodes())
        assert_edges_equal(G.edges(), I.edges())

    def test_read_attribute_graphml(self):
        G = self.attribute_graph
        H = nx.read_graphml(self.attribute_fh)
        assert_nodes_equal(G.nodes(True), sorted(H.nodes(data=True)))
        ge = sorted(G.edges(data=True))
        he = sorted(H.edges(data=True))
        for a, b in zip(ge, he):
            assert_equal(a, b)
        self.attribute_fh.seek(0)

        I = nx.parse_graphml(self.attribute_data)
        assert_equal(sorted(G.nodes(True)), sorted(I.nodes(data=True)))
        ge = sorted(G.edges(data=True))
        he = sorted(I.edges(data=True))
        for a, b in zip(ge, he):
            assert_equal(a, b)

    def test_directed_edge_in_undirected(self):
        s = """<?xml version="1.0" encoding="UTF-8"?>
<graphml xmlns="http://graphml.graphdrawing.org/xmlns"
         xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
         xsi:schemaLocation="http://graphml.graphdrawing.org/xmlns
         http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd">
  <graph id="G">
    <node id="n0"/>
    <node id="n1"/>
    <node id="n2"/>
    <edge source="n0" target="n1"/>
    <edge source="n1" target="n2" directed='true'/>
  </graph>
</graphml>"""
        fh = io.BytesIO(s.encode('UTF-8'))
        assert_raises(nx.NetworkXError, nx.read_graphml, fh)
        assert_raises(nx.NetworkXError, nx.parse_graphml, s)

    def test_undirected_edge_in_directed(self):
        s = """<?xml version="1.0" encoding="UTF-8"?>
<graphml xmlns="http://graphml.graphdrawing.org/xmlns"
         xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
         xsi:schemaLocation="http://graphml.graphdrawing.org/xmlns
         http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd">
  <graph id="G" edgedefault='directed'>
    <node id="n0"/>
    <node id="n1"/>
    <node id="n2"/>
    <edge source="n0" target="n1"/>
    <edge source="n1" target="n2" directed='false'/>
  </graph>
</graphml>"""
        fh = io.BytesIO(s.encode('UTF-8'))
        assert_raises(nx.NetworkXError, nx.read_graphml, fh)
        assert_raises(nx.NetworkXError, nx.parse_graphml, s)

    def test_key_raise(self):
        s = """<?xml version="1.0" encoding="UTF-8"?>
<graphml xmlns="http://graphml.graphdrawing.org/xmlns"
         xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
         xsi:schemaLocation="http://graphml.graphdrawing.org/xmlns
         http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd">
  <key id="d0" for="node" attr.name="color" attr.type="string">
    <default>yellow</default>
  </key>
  <key id="d1" for="edge" attr.name="weight" attr.type="double"/>
  <graph id="G" edgedefault="directed">
    <node id="n0">
      <data key="d0">green</data>
    </node>
    <node id="n1"/>
    <node id="n2">
      <data key="d0">blue</data>
    </node>
    <edge id="e0" source="n0" target="n2">
      <data key="d2">1.0</data>
    </edge>
  </graph>
</graphml>
"""
        fh = io.BytesIO(s.encode('UTF-8'))
        assert_raises(nx.NetworkXError, nx.read_graphml, fh)
        assert_raises(nx.NetworkXError, nx.parse_graphml, s)

    def test_hyperedge_raise(self):
        s = """<?xml version="1.0" encoding="UTF-8"?>
<graphml xmlns="http://graphml.graphdrawing.org/xmlns"
         xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
         xsi:schemaLocation="http://graphml.graphdrawing.org/xmlns
         http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd">
  <key id="d0" for="node" attr.name="color" attr.type="string">
    <default>yellow</default>
  </key>
  <key id="d1" for="edge" attr.name="weight" attr.type="double"/>
  <graph id="G" edgedefault="directed">
    <node id="n0">
      <data key="d0">green</data>
    </node>
    <node id="n1"/>
    <node id="n2">
      <data key="d0">blue</data>
    </node>
    <hyperedge id="e0" source="n0" target="n2">
       <endpoint node="n0"/>
       <endpoint node="n1"/>
       <endpoint node="n2"/>
    </hyperedge>
  </graph>
</graphml>
"""
        fh = io.BytesIO(s.encode('UTF-8'))
        assert_raises(nx.NetworkXError, nx.read_graphml, fh)
        assert_raises(nx.NetworkXError, nx.parse_graphml, s)

    def test_multigraph_keys(self):
        # Test that reading multigraphs uses edge id attributes as keys
        s = """<?xml version="1.0" encoding="UTF-8"?>
<graphml xmlns="http://graphml.graphdrawing.org/xmlns"
         xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
         xsi:schemaLocation="http://graphml.graphdrawing.org/xmlns
         http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd">
  <graph id="G" edgedefault="directed">
    <node id="n0"/>
    <node id="n1"/>
    <edge id="e0" source="n0" target="n1"/>
    <edge id="e1" source="n0" target="n1"/>
  </graph>
</graphml>
"""
        fh = io.BytesIO(s.encode('UTF-8'))
        G = nx.read_graphml(fh)
        expected = [("n0", "n1", "e0"), ("n0", "n1", "e1")]
        assert_equal(sorted(G.edges(keys=True)), expected)
        fh.seek(0)
        H = nx.parse_graphml(s)
        assert_equal(sorted(H.edges(keys=True)), expected)

    def test_preserve_multi_edge_data(self):
        """
        Test that data and keys of edges are preserved on consequent
        write and reads
        """
        G = nx.MultiGraph()
        G.add_node(1)
        G.add_node(2)
        G.add_edges_from([
            # edges with no data, no keys:
            (1, 2),
            # edges with only data:
            (1, 2, dict(key='data_key1')),
            (1, 2, dict(id='data_id2')),
            (1, 2, dict(key='data_key3', id='data_id3')),
            # edges with both data and keys:
            (1, 2, 103, dict(key='data_key4')),
            (1, 2, 104, dict(id='data_id5')),
            (1, 2, 105, dict(key='data_key6', id='data_id7')),
        ])
        fh = io.BytesIO()
        nx.write_graphml(G, fh)
        fh.seek(0)
        H = nx.read_graphml(fh, node_type=int)
        assert_edges_equal(
            G.edges(data=True, keys=True), H.edges(data=True, keys=True)
        )
        assert_equal(G._adj, H._adj)

    def test_yfiles_extension(self):
        data = """<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<graphml xmlns="http://graphml.graphdrawing.org/xmlns"
         xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
         xmlns:y="http://www.yworks.com/xml/graphml"
         xmlns:yed="http://www.yworks.com/xml/yed/3"
         xsi:schemaLocation="http://graphml.graphdrawing.org/xmlns
         http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd">
  <!--Created by yFiles for Java 2.7-->
  <key for="graphml" id="d0" yfiles.type="resources"/>
  <key attr.name="url" attr.type="string" for="node" id="d1"/>
  <key attr.name="description" attr.type="string" for="node" id="d2"/>
  <key for="node" id="d3" yfiles.type="nodegraphics"/>
  <key attr.name="Description" attr.type="string" for="graph" id="d4">
    <default/>
  </key>
  <key attr.name="url" attr.type="string" for="edge" id="d5"/>
  <key attr.name="description" attr.type="string" for="edge" id="d6"/>
  <key for="edge" id="d7" yfiles.type="edgegraphics"/>
  <graph edgedefault="directed" id="G">
    <node id="n0">
      <data key="d3">
        <y:ShapeNode>
          <y:Geometry height="30.0" width="30.0" x="125.0" y="100.0"/>
          <y:Fill color="#FFCC00" transparent="false"/>
          <y:BorderStyle color="#000000" type="line" width="1.0"/>
          <y:NodeLabel alignment="center" autoSizePolicy="content"
           borderDistance="0.0" fontFamily="Dialog" fontSize="13"
           fontStyle="plain" hasBackgroundColor="false" hasLineColor="false"
           height="19.1328125" modelName="internal" modelPosition="c"
           textColor="#000000" visible="true" width="12.27099609375"
           x="8.864501953125" y="5.43359375">1</y:NodeLabel>
          <y:Shape type="rectangle"/>
        </y:ShapeNode>
      </data>
    </node>
    <node id="n1">
      <data key="d3">
        <y:ShapeNode>
          <y:Geometry height="30.0" width="30.0" x="183.0" y="205.0"/>
          <y:Fill color="#FFCC00" transparent="false"/>
          <y:BorderStyle color="#000000" type="line" width="1.0"/>
          <y:NodeLabel alignment="center" autoSizePolicy="content"
          borderDistance="0.0" fontFamily="Dialog" fontSize="13"
          fontStyle="plain" hasBackgroundColor="false" hasLineColor="false"
          height="19.1328125" modelName="internal" modelPosition="c"
          textColor="#000000" visible="true" width="12.27099609375"
          x="8.864501953125" y="5.43359375">2</y:NodeLabel>
          <y:Shape type="rectangle"/>
        </y:ShapeNode>
      </data>
    </node>
    <edge id="e0" source="n0" target="n1">
      <data key="d7">
        <y:PolyLineEdge>
          <y:Path sx="0.0" sy="0.0" tx="0.0" ty="0.0"/>
          <y:LineStyle color="#000000" type="line" width="1.0"/>
          <y:Arrows source="none" target="standard"/>
          <y:BendStyle smoothed="false"/>
        </y:PolyLineEdge>
      </data>
    </edge>
  </graph>
  <data key="d0">
    <y:Resources/>
  </data>
</graphml>
"""
        fh = io.BytesIO(data.encode('UTF-8'))
        G = nx.read_graphml(fh)
        assert_equal(list(G.edges()), [('n0', 'n1')])
        assert_equal(G['n0']['n1']['id'], 'e0')
        assert_equal(G.nodes['n0']['label'], '1')
        assert_equal(G.nodes['n1']['label'], '2')

        H = nx.parse_graphml(data)
        assert_equal(list(H.edges()), [('n0', 'n1')])
        assert_equal(H['n0']['n1']['id'], 'e0')
        assert_equal(H.nodes['n0']['label'], '1')
        assert_equal(H.nodes['n1']['label'], '2')

    def test_bool(self):
        s = """<?xml version="1.0" encoding="UTF-8"?>
<graphml xmlns="http://graphml.graphdrawing.org/xmlns"
         xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
         xsi:schemaLocation="http://graphml.graphdrawing.org/xmlns
         http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd">
  <key id="d0" for="node" attr.name="test" attr.type="boolean">
    <default>false</default>
  </key>
  <graph id="G" edgedefault="directed">
    <node id="n0">
      <data key="d0">true</data>
    </node>
    <node id="n1"/>
    <node id="n2">
      <data key="d0">false</data>
    </node>
    <node id="n3">
      <data key="d0">FaLsE</data>
    </node>
    <node id="n4">
      <data key="d0">True</data>
    </node>
    <node id="n5">
      <data key="d0">0</data>
    </node>
    <node id="n6">
      <data key="d0">1</data>
    </node>
  </graph>
</graphml>
"""
        fh = io.BytesIO(s.encode('UTF-8'))
        G = nx.read_graphml(fh)
        H = nx.parse_graphml(s)
        for graph in [G, H]:
            assert_equal(graph.nodes['n0']['test'], True)
            assert_equal(graph.nodes['n2']['test'], False)
            assert_equal(graph.nodes['n3']['test'], False)
            assert_equal(graph.nodes['n4']['test'], True)
            assert_equal(graph.nodes['n5']['test'], False)
            assert_equal(graph.nodes['n6']['test'], True)

    def test_graphml_header_line(self):
        good = """<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<graphml xmlns="http://graphml.graphdrawing.org/xmlns"
         xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
         xsi:schemaLocation="http://graphml.graphdrawing.org/xmlns
         http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd">
  <key id="d0" for="node" attr.name="test" attr.type="boolean">
    <default>false</default>
  </key>
  <graph id="G">
    <node id="n0">
      <data key="d0">true</data>
    </node>
  </graph>
</graphml>
"""
        bad = """<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<graphml>
  <key id="d0" for="node" attr.name="test" attr.type="boolean">
    <default>false</default>
  </key>
  <graph id="G">
    <node id="n0">
      <data key="d0">true</data>
    </node>
  </graph>
</graphml>
"""
        ugly = """<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<graphml xmlns="https://ghghgh">
  <key id="d0" for="node" attr.name="test" attr.type="boolean">
    <default>false</default>
  </key>
  <graph id="G">
    <node id="n0">
      <data key="d0">true</data>
    </node>
  </graph>
</graphml>
"""
        for s in (good, bad):
            fh = io.BytesIO(s.encode('UTF-8'))
            G = nx.read_graphml(fh)
            H = nx.parse_graphml(s)
            for graph in [G, H]:
                assert_equal(graph.nodes['n0']['test'], True)

        fh = io.BytesIO(ugly.encode('UTF-8'))
        assert_raises(nx.NetworkXError, nx.read_graphml, fh)
        assert_raises(nx.NetworkXError, nx.parse_graphml, ugly)

    def test_read_attributes_with_groups(self):
        data = """\
<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<graphml xmlns="http://graphml.graphdrawing.org/xmlns" xmlns:java="http://www.yworks.com/xml/yfiles-common/1.0/java" xmlns:sys="http://www.yworks.com/xml/yfiles-common/markup/primitives/2.0" xmlns:x="http://www.yworks.com/xml/yfiles-common/markup/2.0" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xmlns:y="http://www.yworks.com/xml/graphml" xmlns:yed="http://www.yworks.com/xml/yed/3" xsi:schemaLocation="http://graphml.graphdrawing.org/xmlns http://www.yworks.com/xml/schema/graphml/1.1/ygraphml.xsd">
  <!--Created by yEd 3.17-->
  <key attr.name="Description" attr.type="string" for="graph" id="d0"/>
  <key for="port" id="d1" yfiles.type="portgraphics"/>
  <key for="port" id="d2" yfiles.type="portgeometry"/>
  <key for="port" id="d3" yfiles.type="portuserdata"/>
  <key attr.name="CustomProperty" attr.type="string" for="node" id="d4">
    <default/>
  </key>
  <key attr.name="url" attr.type="string" for="node" id="d5"/>
  <key attr.name="description" attr.type="string" for="node" id="d6"/>
  <key for="node" id="d7" yfiles.type="nodegraphics"/>
  <key for="graphml" id="d8" yfiles.type="resources"/>
  <key attr.name="url" attr.type="string" for="edge" id="d9"/>
  <key attr.name="description" attr.type="string" for="edge" id="d10"/>
  <key for="edge" id="d11" yfiles.type="edgegraphics"/>
  <graph edgedefault="directed" id="G">
    <data key="d0"/>
    <node id="n0">
      <data key="d4"><![CDATA[CustomPropertyValue]]></data>
      <data key="d6"/>
      <data key="d7">
        <y:ShapeNode>
          <y:Geometry height="30.0" width="30.0" x="125.0" y="-255.4611111111111"/>
          <y:Fill color="#FFCC00" transparent="false"/>
          <y:BorderStyle color="#000000" raised="false" type="line" width="1.0"/>
          <y:NodeLabel alignment="center" autoSizePolicy="content" fontFamily="Dialog" fontSize="12" fontStyle="plain" hasBackgroundColor="false" hasLineColor="false" height="17.96875" horizontalTextPosition="center" iconTextGap="4" modelName="custom" textColor="#000000" verticalTextPosition="bottom" visible="true" width="11.634765625" x="9.1826171875" y="6.015625">2<y:LabelModel>
              <y:SmartNodeLabelModel distance="4.0"/>
            </y:LabelModel>
            <y:ModelParameter>
              <y:SmartNodeLabelModelParameter labelRatioX="0.0" labelRatioY="0.0" nodeRatioX="0.0" nodeRatioY="0.0" offsetX="0.0" offsetY="0.0" upX="0.0" upY="-1.0"/>
            </y:ModelParameter>
          </y:NodeLabel>
          <y:Shape type="rectangle"/>
        </y:ShapeNode>
      </data>
    </node>
    <node id="n1" yfiles.foldertype="group">
      <data key="d4"><![CDATA[CustomPropertyValue]]></data>
      <data key="d5"/>
      <data key="d6"/>
      <data key="d7">
        <y:ProxyAutoBoundsNode>
          <y:Realizers active="0">
            <y:GroupNode>
              <y:Geometry height="250.38333333333333" width="140.0" x="-30.0" y="-330.3833333333333"/>
              <y:Fill color="#F5F5F5" transparent="false"/>
              <y:BorderStyle color="#000000" type="dashed" width="1.0"/>
              <y:NodeLabel alignment="right" autoSizePolicy="node_width" backgroundColor="#EBEBEB" borderDistance="0.0" fontFamily="Dialog" fontSize="15" fontStyle="plain" hasLineColor="false" height="21.4609375" horizontalTextPosition="center" iconTextGap="4" modelName="internal" modelPosition="t" textColor="#000000" verticalTextPosition="bottom" visible="true" width="140.0" x="0.0" y="0.0">Group 3</y:NodeLabel>
              <y:Shape type="roundrectangle"/>
              <y:State closed="false" closedHeight="50.0" closedWidth="50.0" innerGraphDisplayEnabled="false"/>
              <y:Insets bottom="15" bottomF="15.0" left="15" leftF="15.0" right="15" rightF="15.0" top="15" topF="15.0"/>
              <y:BorderInsets bottom="1" bottomF="1.0" left="0" leftF="0.0" right="0" rightF="0.0" top="1" topF="1.0001736111111086"/>
            </y:GroupNode>
            <y:GroupNode>
              <y:Geometry height="50.0" width="50.0" x="0.0" y="60.0"/>
              <y:Fill color="#F5F5F5" transparent="false"/>
              <y:BorderStyle color="#000000" type="dashed" width="1.0"/>
              <y:NodeLabel alignment="right" autoSizePolicy="node_width" backgroundColor="#EBEBEB" borderDistance="0.0" fontFamily="Dialog" fontSize="15" fontStyle="plain" hasLineColor="false" height="21.4609375" horizontalTextPosition="center" iconTextGap="4" modelName="internal" modelPosition="t" textColor="#000000" verticalTextPosition="bottom" visible="true" width="65.201171875" x="-7.6005859375" y="0.0">Folder 3</y:NodeLabel>
              <y:Shape type="roundrectangle"/>
              <y:State closed="true" closedHeight="50.0" closedWidth="50.0" innerGraphDisplayEnabled="false"/>
              <y:Insets bottom="5" bottomF="5.0" left="5" leftF="5.0" right="5" rightF="5.0" top="5" topF="5.0"/>
              <y:BorderInsets bottom="0" bottomF="0.0" left="0" leftF="0.0" right="0" rightF="0.0" top="0" topF="0.0"/>
            </y:GroupNode>
          </y:Realizers>
        </y:ProxyAutoBoundsNode>
      </data>
      <graph edgedefault="directed" id="n1:">
        <node id="n1::n0" yfiles.foldertype="group">
          <data key="d4"><![CDATA[CustomPropertyValue]]></data>
          <data key="d5"/>
          <data key="d6"/>
          <data key="d7">
            <y:ProxyAutoBoundsNode>
              <y:Realizers active="0">
                <y:GroupNode>
                  <y:Geometry height="83.46111111111111" width="110.0" x="-15.0" y="-292.9222222222222"/>
                  <y:Fill color="#F5F5F5" transparent="false"/>
                  <y:BorderStyle color="#000000" type="dashed" width="1.0"/>
                  <y:NodeLabel alignment="right" autoSizePolicy="node_width" backgroundColor="#EBEBEB" borderDistance="0.0" fontFamily="Dialog" fontSize="15" fontStyle="plain" hasLineColor="false" height="21.4609375" horizontalTextPosition="center" iconTextGap="4" modelName="internal" modelPosition="t" textColor="#000000" verticalTextPosition="bottom" visible="true" width="110.0" x="0.0" y="0.0">Group 1</y:NodeLabel>
                  <y:Shape type="roundrectangle"/>
                  <y:State closed="false" closedHeight="50.0" closedWidth="50.0" innerGraphDisplayEnabled="false"/>
                  <y:Insets bottom="15" bottomF="15.0" left="15" leftF="15.0" right="15" rightF="15.0" top="15" topF="15.0"/>
                  <y:BorderInsets bottom="1" bottomF="1.0" left="0" leftF="0.0" right="0" rightF="0.0" top="1" topF="1.0001736111111086"/>
                </y:GroupNode>
                <y:GroupNode>
                  <y:Geometry height="50.0" width="50.0" x="0.0" y="60.0"/>
                  <y:Fill color="#F5F5F5" transparent="false"/>
                  <y:BorderStyle color="#000000" type="dashed" width="1.0"/>
                  <y:NodeLabel alignment="right" autoSizePolicy="node_width" backgroundColor="#EBEBEB" borderDistance="0.0" fontFamily="Dialog" fontSize="15" fontStyle="plain" hasLineColor="false" height="21.4609375" horizontalTextPosition="center" iconTextGap="4" modelName="internal" modelPosition="t" textColor="#000000" verticalTextPosition="bottom" visible="true" width="65.201171875" x="-7.6005859375" y="0.0">Folder 1</y:NodeLabel>
                  <y:Shape type="roundrectangle"/>
                  <y:State closed="true" closedHeight="50.0" closedWidth="50.0" innerGraphDisplayEnabled="false"/>
                  <y:Insets bottom="5" bottomF="5.0" left="5" leftF="5.0" right="5" rightF="5.0" top="5" topF="5.0"/>
                  <y:BorderInsets bottom="0" bottomF="0.0" left="0" leftF="0.0" right="0" rightF="0.0" top="0" topF="0.0"/>
                </y:GroupNode>
              </y:Realizers>
            </y:ProxyAutoBoundsNode>
          </data>
          <graph edgedefault="directed" id="n1::n0:">
            <node id="n1::n0::n0">
              <data key="d4"><![CDATA[CustomPropertyValue]]></data>
              <data key="d6"/>
              <data key="d7">
                <y:ShapeNode>
                  <y:Geometry height="30.0" width="30.0" x="50.0" y="-255.4611111111111"/>
                  <y:Fill color="#FFCC00" transparent="false"/>
                  <y:BorderStyle color="#000000" raised="false" type="line" width="1.0"/>
                  <y:NodeLabel alignment="center" autoSizePolicy="content" fontFamily="Dialog" fontSize="12" fontStyle="plain" hasBackgroundColor="false" hasLineColor="false" height="17.96875" horizontalTextPosition="center" iconTextGap="4" modelName="custom" textColor="#000000" verticalTextPosition="bottom" visible="true" width="11.634765625" x="9.1826171875" y="6.015625">1<y:LabelModel>
                      <y:SmartNodeLabelModel distance="4.0"/>
                    </y:LabelModel>
                    <y:ModelParameter>
                      <y:SmartNodeLabelModelParameter labelRatioX="0.0" labelRatioY="0.0" nodeRatioX="0.0" nodeRatioY="0.0" offsetX="0.0" offsetY="0.0" upX="0.0" upY="-1.0"/>
                    </y:ModelParameter>
                  </y:NodeLabel>
                  <y:Shape type="rectangle"/>
                </y:ShapeNode>
              </data>
            </node>
            <node id="n1::n0::n1">
              <data key="d4"><![CDATA[CustomPropertyValue]]></data>
              <data key="d6"/>
              <data key="d7">
                <y:ShapeNode>
                  <y:Geometry height="30.0" width="30.0" x="0.0" y="-255.4611111111111"/>
                  <y:Fill color="#FFCC00" transparent="false"/>
                  <y:BorderStyle color="#000000" raised="false" type="line" width="1.0"/>
                  <y:NodeLabel alignment="center" autoSizePolicy="content" fontFamily="Dialog" fontSize="12" fontStyle="plain" hasBackgroundColor="false" hasLineColor="false" height="17.96875" horizontalTextPosition="center" iconTextGap="4" modelName="custom" textColor="#000000" verticalTextPosition="bottom" visible="true" width="11.634765625" x="9.1826171875" y="6.015625">3<y:LabelModel>
                      <y:SmartNodeLabelModel distance="4.0"/>
                    </y:LabelModel>
                    <y:ModelParameter>
                      <y:SmartNodeLabelModelParameter labelRatioX="0.0" labelRatioY="0.0" nodeRatioX="0.0" nodeRatioY="0.0" offsetX="0.0" offsetY="0.0" upX="0.0" upY="-1.0"/>
                    </y:ModelParameter>
                  </y:NodeLabel>
                  <y:Shape type="rectangle"/>
                </y:ShapeNode>
              </data>
            </node>
          </graph>
        </node>
        <node id="n1::n1" yfiles.foldertype="group">
          <data key="d4"><![CDATA[CustomPropertyValue]]></data>
          <data key="d5"/>
          <data key="d6"/>
          <data key="d7">
            <y:ProxyAutoBoundsNode>
              <y:Realizers active="0">
                <y:GroupNode>
                  <y:Geometry height="83.46111111111111" width="110.0" x="-15.0" y="-179.4611111111111"/>
                  <y:Fill color="#F5F5F5" transparent="false"/>
                  <y:BorderStyle color="#000000" type="dashed" width="1.0"/>
                  <y:NodeLabel alignment="right" autoSizePolicy="node_width" backgroundColor="#EBEBEB" borderDistance="0.0" fontFamily="Dialog" fontSize="15" fontStyle="plain" hasLineColor="false" height="21.4609375" horizontalTextPosition="center" iconTextGap="4" modelName="internal" modelPosition="t" textColor="#000000" verticalTextPosition="bottom" visible="true" width="110.0" x="0.0" y="0.0">Group 2</y:NodeLabel>
                  <y:Shape type="roundrectangle"/>
                  <y:State closed="false" closedHeight="50.0" closedWidth="50.0" innerGraphDisplayEnabled="false"/>
                  <y:Insets bottom="15" bottomF="15.0" left="15" leftF="15.0" right="15" rightF="15.0" top="15" topF="15.0"/>
                  <y:BorderInsets bottom="1" bottomF="1.0" left="0" leftF="0.0" right="0" rightF="0.0" top="1" topF="1.0001736111111086"/>
                </y:GroupNode>
                <y:GroupNode>
                  <y:Geometry height="50.0" width="50.0" x="0.0" y="60.0"/>
                  <y:Fill color="#F5F5F5" transparent="false"/>
                  <y:BorderStyle color="#000000" type="dashed" width="1.0"/>
                  <y:NodeLabel alignment="right" autoSizePolicy="node_width" backgroundColor="#EBEBEB" borderDistance="0.0" fontFamily="Dialog" fontSize="15" fontStyle="plain" hasLineColor="false" height="21.4609375" horizontalTextPosition="center" iconTextGap="4" modelName="internal" modelPosition="t" textColor="#000000" verticalTextPosition="bottom" visible="true" width="65.201171875" x="-7.6005859375" y="0.0">Folder 2</y:NodeLabel>
                  <y:Shape type="roundrectangle"/>
                  <y:State closed="true" closedHeight="50.0" closedWidth="50.0" innerGraphDisplayEnabled="false"/>
                  <y:Insets bottom="5" bottomF="5.0" left="5" leftF="5.0" right="5" rightF="5.0" top="5" topF="5.0"/>
                  <y:BorderInsets bottom="0" bottomF="0.0" left="0" leftF="0.0" right="0" rightF="0.0" top="0" topF="0.0"/>
                </y:GroupNode>
              </y:Realizers>
            </y:ProxyAutoBoundsNode>
          </data>
          <graph edgedefault="directed" id="n1::n1:">
            <node id="n1::n1::n0">
              <data key="d4"><![CDATA[CustomPropertyValue]]></data>
              <data key="d6"/>
              <data key="d7">
                <y:ShapeNode>
                  <y:Geometry height="30.0" width="30.0" x="0.0" y="-142.0"/>
                  <y:Fill color="#FFCC00" transparent="false"/>
                  <y:BorderStyle color="#000000" raised="false" type="line" width="1.0"/>
                  <y:NodeLabel alignment="center" autoSizePolicy="content" fontFamily="Dialog" fontSize="12" fontStyle="plain" hasBackgroundColor="false" hasLineColor="false" height="17.96875" horizontalTextPosition="center" iconTextGap="4" modelName="custom" textColor="#000000" verticalTextPosition="bottom" visible="true" width="11.634765625" x="9.1826171875" y="6.015625">5<y:LabelModel>
                      <y:SmartNodeLabelModel distance="4.0"/>
                    </y:LabelModel>
                    <y:ModelParameter>
                      <y:SmartNodeLabelModelParameter labelRatioX="0.0" labelRatioY="0.0" nodeRatioX="0.0" nodeRatioY="0.0" offsetX="0.0" offsetY="0.0" upX="0.0" upY="-1.0"/>
                    </y:ModelParameter>
                  </y:NodeLabel>
                  <y:Shape type="rectangle"/>
                </y:ShapeNode>
              </data>
            </node>
            <node id="n1::n1::n1">
              <data key="d4"><![CDATA[CustomPropertyValue]]></data>
              <data key="d6"/>
              <data key="d7">
                <y:ShapeNode>
                  <y:Geometry height="30.0" width="30.0" x="50.0" y="-142.0"/>
                  <y:Fill color="#FFCC00" transparent="false"/>
                  <y:BorderStyle color="#000000" raised="false" type="line" width="1.0"/>
                  <y:NodeLabel alignment="center" autoSizePolicy="content" fontFamily="Dialog" fontSize="12" fontStyle="plain" hasBackgroundColor="false" hasLineColor="false" height="17.96875" horizontalTextPosition="center" iconTextGap="4" modelName="custom" textColor="#000000" verticalTextPosition="bottom" visible="true" width="11.634765625" x="9.1826171875" y="6.015625">6<y:LabelModel>
                      <y:SmartNodeLabelModel distance="4.0"/>
                    </y:LabelModel>
                    <y:ModelParameter>
                      <y:SmartNodeLabelModelParameter labelRatioX="0.0" labelRatioY="0.0" nodeRatioX="0.0" nodeRatioY="0.0" offsetX="0.0" offsetY="0.0" upX="0.0" upY="-1.0"/>
                    </y:ModelParameter>
                  </y:NodeLabel>
                  <y:Shape type="rectangle"/>
                </y:ShapeNode>
              </data>
            </node>
          </graph>
        </node>
      </graph>
    </node>
    <node id="n2">
      <data key="d4"><![CDATA[CustomPropertyValue]]></data>
      <data key="d6"/>
      <data key="d7">
        <y:ShapeNode>
          <y:Geometry height="30.0" width="30.0" x="125.0" y="-142.0"/>
          <y:Fill color="#FFCC00" transparent="false"/>
          <y:BorderStyle color="#000000" raised="false" type="line" width="1.0"/>
          <y:NodeLabel alignment="center" autoSizePolicy="content" fontFamily="Dialog" fontSize="12" fontStyle="plain" hasBackgroundColor="false" hasLineColor="false" height="17.96875" horizontalTextPosition="center" iconTextGap="4" modelName="custom" textColor="#000000" verticalTextPosition="bottom" visible="true" width="11.634765625" x="9.1826171875" y="6.015625">9<y:LabelModel>
              <y:SmartNodeLabelModel distance="4.0"/>
            </y:LabelModel>
            <y:ModelParameter>
              <y:SmartNodeLabelModelParameter labelRatioX="0.0" labelRatioY="0.0" nodeRatioX="0.0" nodeRatioY="0.0" offsetX="0.0" offsetY="0.0" upX="0.0" upY="-1.0"/>
            </y:ModelParameter>
          </y:NodeLabel>
          <y:Shape type="rectangle"/>
        </y:ShapeNode>
      </data>
    </node>
    <edge id="n1::n1::e0" source="n1::n1::n0" target="n1::n1::n1">
      <data key="d10"/>
      <data key="d11">
        <y:PolyLineEdge>
          <y:Path sx="15.0" sy="-0.0" tx="-15.0" ty="-0.0"/>
          <y:LineStyle color="#000000" type="line" width="1.0"/>
          <y:Arrows source="none" target="standard"/>
          <y:BendStyle smoothed="false"/>
        </y:PolyLineEdge>
      </data>
    </edge>
    <edge id="n1::n0::e0" source="n1::n0::n1" target="n1::n0::n0">
      <data key="d10"/>
      <data key="d11">
        <y:PolyLineEdge>
          <y:Path sx="15.0" sy="-0.0" tx="-15.0" ty="-0.0"/>
          <y:LineStyle color="#000000" type="line" width="1.0"/>
          <y:Arrows source="none" target="standard"/>
          <y:BendStyle smoothed="false"/>
        </y:PolyLineEdge>
      </data>
    </edge>
    <edge id="e0" source="n1::n0::n0" target="n0">
      <data key="d10"/>
      <data key="d11">
        <y:PolyLineEdge>
          <y:Path sx="15.0" sy="-0.0" tx="-15.0" ty="-0.0"/>
          <y:LineStyle color="#000000" type="line" width="1.0"/>
          <y:Arrows source="none" target="standard"/>
          <y:BendStyle smoothed="false"/>
        </y:PolyLineEdge>
      </data>
    </edge>
    <edge id="e1" source="n1::n1::n1" target="n2">
      <data key="d10"/>
      <data key="d11">
        <y:PolyLineEdge>
          <y:Path sx="15.0" sy="-0.0" tx="-15.0" ty="-0.0"/>
          <y:LineStyle color="#000000" type="line" width="1.0"/>
          <y:Arrows source="none" target="standard"/>
          <y:BendStyle smoothed="false"/>
        </y:PolyLineEdge>
      </data>
    </edge>
  </graph>
  <data key="d8">
    <y:Resources/>
  </data>
</graphml>
"""
        # verify that nodes / attributes are correctly read when part of a group
        fh = io.BytesIO(data.encode('UTF-8'))
        G = nx.read_graphml(fh)
        data = [x for _, x in G.nodes(data=True)]
        assert_equal(len(data), 9)
        for node_data in data:
            assert_not_equal(node_data['CustomProperty'], '')


class TestWriteGraphML(BaseGraphML):
    writer = staticmethod(nx.write_graphml_lxml)

    @classmethod
    def setupClass(cls):
        try:
            import lxml.etree
        except ImportError:
            raise SkipTest('lxml.etree not available.')

    def test_write_interface(self):
        try:
            import lxml.etree
            assert_equal(nx.write_graphml, nx.write_graphml_lxml)
        except ImportError:
            assert_equal(nx.write_graphml, nx.write_graphml_xml)

    def test_write_read_simple_directed_graphml(self):
        G = self.simple_directed_graph
        G.graph['hi'] = 'there'
        fh = io.BytesIO()
        self.writer(G, fh)
        fh.seek(0)
        H = nx.read_graphml(fh)
        assert_equal(sorted(G.nodes()), sorted(H.nodes()))
        assert_equal(sorted(G.edges()), sorted(H.edges()))
        assert_equal(sorted(G.edges(data=True)), sorted(H.edges(data=True)))
        self.simple_directed_fh.seek(0)

    def test_write_read_attribute_numeric_type_graphml(self):
        from xml.etree.ElementTree import parse

        G = self.attribute_numeric_type_graph
        fh = io.BytesIO()
        self.writer(G, fh, infer_numeric_types=True)
        fh.seek(0)
        H = nx.read_graphml(fh)
        fh.seek(0)

        assert_nodes_equal(G.nodes(), H.nodes())
        assert_edges_equal(G.edges(), H.edges())
        assert_edges_equal(G.edges(data=True), H.edges(data=True))
        self.attribute_numeric_type_fh.seek(0)

        xml = parse(fh)
        # Children are the key elements, and the graph element
        children = xml.getroot().getchildren()
        assert_equal(len(children), 3)

        keys = [child.items() for child in children[:2]]

        assert_equal(len(keys), 2)
        assert_in(('attr.type', 'double'), keys[0])
        assert_in(('attr.type', 'double'), keys[1])

    def test_more_multigraph_keys(self):
        """Writing keys as edge id attributes means keys become strings.
        The original keys are stored as data, so read them back in
        if `make_str(key) == edge_id`
        This allows the adjacency to remain the same.
        """
        G = nx.MultiGraph()
        G.add_edges_from([('a', 'b', 2), ('a', 'b', 3)])
        fd, fname = tempfile.mkstemp()
        self.writer(G, fname)
        H = nx.read_graphml(fname)
        assert_true(H.is_multigraph())
        assert_edges_equal(G.edges(keys=True), H.edges(keys=True))
        assert_equal(G._adj, H._adj)
        os.close(fd)
        os.unlink(fname)

    def test_default_attribute(self):
        G = nx.Graph(name="Fred")
        G.add_node(1, label=1, color='green')
        nx.add_path(G, [0, 1, 2, 3])
        G.add_edge(1, 2, weight=3)
        G.graph['node_default'] = {'color': 'yellow'}
        G.graph['edge_default'] = {'weight': 7}
        fh = io.BytesIO()
        self.writer(G, fh)
        fh.seek(0)
        H = nx.read_graphml(fh, node_type=int)
        assert_nodes_equal(G.nodes(), H.nodes())
        assert_edges_equal(G.edges(), H.edges())
        assert_equal(G.graph, H.graph)

    def test_multigraph_to_graph(self):
        # test converting multigraph to graph if no parallel edges found
        G = nx.MultiGraph()
        G.add_edges_from([('a', 'b', 2), ('b', 'c', 3)])  # no multiedges
        fd, fname = tempfile.mkstemp()
        self.writer(G, fname)
        H = nx.read_graphml(fname)
        assert_false(H.is_multigraph())
        os.close(fd)
        os.unlink(fname)

    def test_unicode_attributes(self):
        G = nx.Graph()
        try:  # Python 3.x
            name1 = chr(2344) + chr(123) + chr(6543)
            name2 = chr(5543) + chr(1543) + chr(324)
            node_type = str
        except ValueError:  # Python 2.6+
            name1 = unichr(2344) + unichr(123) + unichr(6543)
            name2 = unichr(5543) + unichr(1543) + unichr(324)
            node_type = unicode
        G.add_edge(name1, 'Radiohead', foo=name2)
        fd, fname = tempfile.mkstemp()
        self.writer(G, fname)
        H = nx.read_graphml(fname, node_type=node_type)
        assert_equal(G._adj, H._adj)
        os.close(fd)
        os.unlink(fname)

    def test_unicode_escape(self):
        # test for handling json escaped stings in python 2 Issue #1880
        import json
        a = dict(a='{"a": "123"}')  # an object with many chars to escape
        try:  # Python 3.x
            chr(2344)
            sa = json.dumps(a)
        except ValueError:  # Python 2.6+
            sa = unicode(json.dumps(a))
        G = nx.Graph()
        G.graph['test'] = sa
        fh = io.BytesIO()
        self.writer(G, fh)
        fh.seek(0)
        H = nx.read_graphml(fh)
        assert_equal(G.graph['test'], H.graph['test'])


class TestXMLGraphML(TestWriteGraphML):
    writer = staticmethod(nx.write_graphml_xml)

    @classmethod
    def setupClass(cls):
        try:
            import xml.etree.ElementTree
        except ImportError:
            raise SkipTest('xml.etree.ElementTree not available.')
