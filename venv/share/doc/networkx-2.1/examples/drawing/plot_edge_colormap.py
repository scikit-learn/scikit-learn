#!/usr/bin/env python
"""
=============
Edge Colormap
=============

Draw a graph with matplotlib, color edges.
You must have matplotlib>=87.7 for this to work.
"""
# Author: Aric Hagberg (hagberg@lanl.gov)

import matplotlib.pyplot as plt
import networkx as nx

G = nx.star_graph(20)
pos = nx.spring_layout(G)
colors = range(20)
nx.draw(G, pos, node_color='#A0CBE2', edge_color=colors,
        width=4, edge_cmap=plt.cm.Blues, with_labels=False)
plt.show()
