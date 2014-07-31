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
1 loops, best of 3: 7.3 s per loop

>>> %memit XP1 = p.transform_OLD(X)
peak memory: 4325.37 MiB, increment: 3701.51 MiB

>>> %timeit XP2 = p.transform(X)
1 loops, best of 3: 1.75 s per loop

>>> %memit XP2 = p.transform(X)
peak memory: 1346.08 MiB, increment: 751.89 MiB

"""
# Author: Alejandro Correa Bahnsen <al.bahnsen@gmail.com>
#
# License: BSD 3 clause

import numpy as np
from sklearn.preprocessing import PolynomialFeatures

n_samples, n_features = 100000, 20
X = np.arange(-n_samples, n_samples * n_features - n_samples).reshape((n_samples, n_features))

p = PolynomialFeatures(degree=2)
p.fit(X)

XP1 = p.transform_OLD(X)
XP2 = p.transform(X)