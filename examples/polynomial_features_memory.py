"""
================================
Polynomial Features Memory Usage
================================

Test of small adjustment to the transform function of PolynomialFeatures.
When executing the transform function the memory usage peak during execution
only to drop after the procedure end. After investigating, the issue is with
'(X[:, None, :] ** self.powers_).prod(-1)', which first creates a huge matrix
and the reduce it after using the product function. In order to solve this, 
I simply made the transformation by example.

The results shown that for the example below, the memory consumption decreases
from 3455 MiB to 374 MiB, without increasing the execution time.

Examples
--------
>>> %timeit XP1 = p.transform_OLD(X)
1 loops, best of 3: 6.95 s per loop

>>> %memit XP1 = p.transform_OLD(X)
peak memory: 3455.32 MiB, increment: 3405.33 MiB
polynomial_features_memory
>>> %timeit XP2 = p.transform(X)
1 loops, best of 3: 6.47 s per loop

>>> %memit XP2 = p.transform(X)
peak memory: 374.59 MiB, increment: 161.95 MiB

"""
# Author: Alejandro Correa Bahnsen <al.bahnsen@gmail.com>
#
# License: BSD 3 clause

import numpy as np
from sklearn.preprocessing import PolynomialFeatures

n_samples, n_features = 2000, 20
X = np.arange(n_samples * n_features).reshape((n_samples, n_features))

p = PolynomialFeatures(degree=4)
p.fit(X)

XP1 = p.transform_OLD(X)
XP2 = p.transform(X)

print np.all(XP1 == XP2)

