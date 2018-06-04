#!/usr/bin/env python
"""
===========
Simple Path
===========

Draw a graph with matplotlib.
"""
import matplotlib.pyplot as plt
import networkx as nx

G = nx.path_graph(8)
nx.draw(G)
plt.show()
