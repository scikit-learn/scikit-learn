from itertools import product

from numpy.testing import assert_array_equal
import numpy as np
import pytest

from scipy.optimize import linear_sum_assignment
from scipy.sparse import csr_matrix, random
from scipy.sparse.csgraph import min_weight_full_bipartite_matching


# Tests that combine scipy.optimize.linear_sum_assignment and
# scipy.sparse.csgraph.min_weight_full_bipartite_matching
@pytest.mark.parametrize('solver_type,sign,test_case', product(
    [(linear_sum_assignment, np.array),
     (min_weight_full_bipartite_matching, csr_matrix)],
    [-1, 1],
    [
        # Square
        ([[400, 150, 400],
          [400, 450, 600],
          [300, 225, 300]],
         [150, 400, 300]),

        # Rectangular variant
        ([[400, 150, 400, 1],
          [400, 450, 600, 2],
          [300, 225, 300, 3]],
         [150, 2, 300]),

        ([[10, 10, 8],
          [9, 8, 1],
          [9, 7, 4]],
         [10, 1, 7]),

        # Square
        ([[10, 10, 8, 11],
          [9, 8, 1, 1],
          [9, 7, 4, 10]],
         [10, 1, 4]),

        # Rectangular variant
        ([[10, float("inf"), float("inf")],
          [float("inf"), float("inf"), 1],
          [float("inf"), 7, float("inf")]],
         [10, 1, 7]),
    ])
)
def test_two_methods_give_expected_result_on_small_inputs(
    solver_type, sign, test_case
):
    solver, array_type = solver_type
    cost_matrix, expected_cost = test_case
    maximize = sign == -1
    cost_matrix = sign * array_type(cost_matrix)
    expected_cost = sign * np.array(expected_cost)

    row_ind, col_ind = solver(cost_matrix, maximize=maximize)
    assert_array_equal(row_ind, np.sort(row_ind))
    assert_array_equal(expected_cost,
                       np.array(cost_matrix[row_ind, col_ind]).flatten())

    cost_matrix = cost_matrix.T
    row_ind, col_ind = solver(cost_matrix, maximize=maximize)
    assert_array_equal(row_ind, np.sort(row_ind))
    assert_array_equal(np.sort(expected_cost),
                       np.sort(np.array(
                           cost_matrix[row_ind, col_ind])).flatten())


def test_two_methods_give_same_result_on_many_sparse_inputs():
    # As opposed to the test above, here we do not spell out the expected
    # output; only assert that the two methods give the same result.
    # Concretely, the below tests 100 cases of size 100x100, out of which
    # 36 are infeasible.
    np.random.seed(1234)
    for _ in range(100):
        lsa_raises = False
        mwfbm_raises = False
        sparse = random(100, 100, density=0.06,
                        data_rvs=lambda size: np.random.randint(1, 100, size))
        # In csgraph, zeros correspond to missing edges, so we explicitly
        # replace those with infinities
        dense = np.full(sparse.shape, np.inf)
        dense[sparse.row, sparse.col] = sparse.data
        sparse = sparse.tocsr()
        try:
            row_ind, col_ind = linear_sum_assignment(dense)
            lsa_cost = dense[row_ind, col_ind].sum()
        except ValueError:
            lsa_raises = True
        try:
            row_ind, col_ind = min_weight_full_bipartite_matching(sparse)
            mwfbm_cost = sparse[row_ind, col_ind].sum()
        except ValueError:
            mwfbm_raises = True
        # Ensure that if one method raises, so does the other one.
        assert lsa_raises == mwfbm_raises
        if not lsa_raises:
            assert lsa_cost == mwfbm_cost
