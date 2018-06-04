#!/usr/bin/env python
"""
========================
Expected Degree Sequence
========================

Random graph from given degree sequence.
"""
# Author: Aric Hagberg (hagberg@lanl.gov)

#    Copyright (C) 2006-2018 by
#    Aric Hagberg <hagberg@lanl.gov>
#    Dan Schult <dschult@colgate.edu>
#    Pieter Swart <swart@lanl.gov>
#    All rights reserved.
#    BSD license.

import networkx as nx
from networkx.generators.degree_seq import expected_degree_graph

# make a random graph of 500 nodes with expected degrees of 50
n = 500  # n nodes
p = 0.1
w = [p * n for i in range(n)]  # w = p*n for all nodes
G = expected_degree_graph(w)  # configuration model
print("Degree histogram")
print("degree (#nodes) ****")
dh = nx.degree_histogram(G)
low = min(nx.degree(G))
for i in range(low, len(dh)):
    bar = ''.join(dh[i] * ['*'])
    print("%2s (%2s) %s" % (i, dh[i], bar))
