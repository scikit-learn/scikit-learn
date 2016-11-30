from __future__ import division, print_function, absolute_import

import numpy as np
from numpy.testing import assert_array_almost_equal
from sklearn.utils.testing import SkipTest

try:
    from scipy.sparse.csgraph import breadth_first_tree, depth_first_tree,\
        csgraph_to_dense, csgraph_from_dense
except ImportError:
    # Oldish versions of scipy don't have that
    csgraph_from_dense = None


def test_graph_breadth_first():
    if csgraph_from_dense is None:
        raise SkipTest("Old version of scipy, doesn't have csgraph.")
    csgraph = np.array([[0, 1, 2, 0, 0],
                        [1, 0, 0, 0, 3],
                        [2, 0, 0, 7, 0],
                        [0, 0, 7, 0, 1],
                        [0, 3, 0, 1, 0]])
    csgraph = csgraph_from_dense(csgraph, null_value=0)

    bfirst = np.array([[0, 1, 2, 0, 0],
                       [0, 0, 0, 0, 3],
                       [0, 0, 0, 7, 0],
                       [0, 0, 0, 0, 0],
                       [0, 0, 0, 0, 0]])

    for directed in [True, False]:
        bfirst_test = breadth_first_tree(csgraph, 0, directed)
        assert_array_almost_equal(csgraph_to_dense(bfirst_test),
                                  bfirst)


def test_graph_depth_first():
    if csgraph_from_dense is None:
        raise SkipTest("Old version of scipy, doesn't have csgraph.")
    csgraph = np.array([[0, 1, 2, 0, 0],
                        [1, 0, 0, 0, 3],
                        [2, 0, 0, 7, 0],
                        [0, 0, 7, 0, 1],
                        [0, 3, 0, 1, 0]])
    csgraph = csgraph_from_dense(csgraph, null_value=0)

    dfirst = np.array([[0, 1, 0, 0, 0],
                       [0, 0, 0, 0, 3],
                       [0, 0, 0, 0, 0],
                       [0, 0, 7, 0, 0],
                       [0, 0, 0, 1, 0]])

    for directed in [True, False]:
        dfirst_test = depth_first_tree(csgraph, 0, directed)
        assert_array_almost_equal(csgraph_to_dense(dfirst_test),
                                  dfirst)
