"""
===
Rcm
===

Cuthill-McKee ordering of matrices

The reverse Cuthill-McKee algorithm gives a sparse matrix ordering that
reduces the matrix bandwidth.
"""

# Copyright (C) 2011-2018 by
# Author:    Aric Hagberg <aric.hagberg@gmail.com>
# BSD License
import networkx as nx
from networkx.utils import reverse_cuthill_mckee_ordering
import numpy as np

# build low-bandwidth numpy matrix
G = nx.grid_2d_graph(3, 3)
rcm = list(reverse_cuthill_mckee_ordering(G))
print("ordering", rcm)

print("unordered Laplacian matrix")
A = nx.laplacian_matrix(G)
x, y = np.nonzero(A)
#print("lower bandwidth:",(y-x).max())
#print("upper bandwidth:",(x-y).max())
print("bandwidth: %d" % ((y - x).max() + (x - y).max() + 1))
print(A)

B = nx.laplacian_matrix(G, nodelist=rcm)
print("low-bandwidth Laplacian matrix")
x, y = np.nonzero(B)
#print("lower bandwidth:",(y-x).max())
#print("upper bandwidth:",(x-y).max())
print("bandwidth: %d" % ((y - x).max() + (x - y).max() + 1))
print(B)
