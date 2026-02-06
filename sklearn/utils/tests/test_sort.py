# SPDX-License-Identifier: BSD-3-Clause

import numpy as np

from sklearn.utils._sorting import _py_simultaneous_sort


def test_simultaneous_sort_recursion_depth():
    dist = np.zeros(1000000, dtype=np.float64)
    ind = np.arange(dist.shape[0], dtype=np.intp)

    _py_simultaneous_sort(dist, ind, dist.shape[0])
