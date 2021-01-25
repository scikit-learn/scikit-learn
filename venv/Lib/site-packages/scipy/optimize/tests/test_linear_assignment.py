# Author: Brian M. Clapper, G. Varoquaux, Lars Buitinck
# License: BSD

from numpy.testing import assert_array_equal
from pytest import raises as assert_raises

import numpy as np

from scipy.optimize import linear_sum_assignment
from scipy.sparse.sputils import matrix


def test_linear_sum_assignment_input_validation():

    assert_raises(ValueError, linear_sum_assignment, [1, 2, 3])

    C = [[1, 2, 3], [4, 5, 6]]
    assert_array_equal(linear_sum_assignment(C),
                       linear_sum_assignment(np.asarray(C)))
    assert_array_equal(linear_sum_assignment(C),
                       linear_sum_assignment(matrix(C)))

    I = np.identity(3)
    assert_array_equal(linear_sum_assignment(I.astype(np.bool_)),
                       linear_sum_assignment(I))
    assert_raises(ValueError, linear_sum_assignment, I.astype(str))

    I[0][0] = np.nan
    assert_raises(ValueError, linear_sum_assignment, I)

    I = np.identity(3)
    I[1][1] = -np.inf
    assert_raises(ValueError, linear_sum_assignment, I)

    I = np.identity(3)
    I[:, 0] = np.inf
    assert_raises(ValueError, linear_sum_assignment, I)


def test_constant_cost_matrix():
    # Fixes #11602
    n = 8
    C = np.ones((n, n))
    row_ind, col_ind = linear_sum_assignment(C)
    assert_array_equal(row_ind, np.arange(n))
    assert_array_equal(col_ind, np.arange(n))
