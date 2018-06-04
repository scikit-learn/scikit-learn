#!/usr/bin/env python
"""
===============
Pygraphviz Draw
===============

An example showing how to use the interface to the pygraphviz
AGraph class to draw a graph.

Also see the pygraphviz documentation and examples at
http://pygraphviz.github.io/
"""
# Author: Aric Hagberg (hagberg@lanl.gov)

#    Copyright (C) 2006-2018 by
#    Aric Hagberg <hagberg@lanl.gov>
#    Dan Schult <dschult@colgate.edu>
#    Pieter Swart <swart@lanl.gov>
#    All rights reserved.
#    BSD license.

import networkx as nx

# plain graph

G = nx.complete_graph(5)   # start with K5 in networkx
A = nx.nx_agraph.to_agraph(G)        # convert to a graphviz graph
A.layout()            # neato layout
A.draw("k5.ps")       # write postscript in k5.ps with neato layout
