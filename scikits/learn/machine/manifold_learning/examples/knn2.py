"""
Demo using k-Nearest Neighbor Algorithm
"""

import numpy as np

from scikits.learn.machine.manifold_learning.regression.neighbors import \
     Neighbors, Kneighbors, Parzen

points = np.array([[0., 0., 0.],
            [1., 0., 0.],
            [0., 1., 0.],
            [0., 0., 1.],
            [1., 1., 0.],
            [1., 0., 1.],
            [0., 1., 1.],
            [1., 1., 1.]], np.float64)

n = Neighbors(points, 4, 1.2)
print n.kneighbors(points[0])
print n.parzen(points[0])

kn = Kneighbors(points, 4, 1.2)
print kn(points[0])
print kn.parzen(points[0])

pn = Parzen(points, 4, 1.2)
print pn.kneighbors(points[0])
print pn(points[0])
