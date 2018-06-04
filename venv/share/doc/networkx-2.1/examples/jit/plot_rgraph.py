#!/usr/bin/env python
"""
======
Rgraph
======

An example showing how to use the JavaScript InfoVis Toolkit (JIT)
JSON export

See the JIT documentation and examples at http://thejit.org

"""
__author__ = """Ollie Glass (ollieglaskovik@gmail.com)"""

import json

import matplotlib.pyplot as plt
import networkx as nx
from networkx.readwrite.json_graph import jit_data, jit_graph

# add some nodes to a graph
G = nx.Graph()

G.add_node("one", type="normal")
G.add_node("two", type="special")
G.add_node("solo")

# add edges
G.add_edge("one", "two")
G.add_edge("two", 3, type="extra special")

# convert to JIT JSON
jit_json = jit_data(G, indent=4)
print(jit_json)

X = jit_graph(json.loads(jit_json))
print("Nodes: %s" % list(X.nodes(data=True)))
print("Edges: %s" % list(X.edges(data=True)))

nx.draw(G, with_labels=True)
plt.show()
