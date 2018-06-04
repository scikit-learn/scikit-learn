#!/usr/bin/env python
"""
=================
Pygraphviz Simple
=================

An example showing how to use the interface to the pygraphviz
AGraph class to convert to and from graphviz.

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
X1 = nx.nx_agraph.from_agraph(A)     # convert back to networkx (but as Graph)
X2 = nx.Graph(A)          # fancy way to do conversion
G1 = nx.Graph(X1)          # now make it a Graph

A.write('k5.dot')     # write to dot file
X3 = nx.nx_agraph.read_dot('k5.dot')  # read from dotfile
