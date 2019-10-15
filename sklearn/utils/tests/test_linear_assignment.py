# Author: Brian M. Clapper, G Varoquaux
# License: BSD

# TODO #0.23: Remove this test module as the methods being tested
# have been replaced by SciPy methods

import numpy as np
import pytest


@pytest.mark.filterwarnings(
  "ignore::FutureWarning")
def test_hungarian():
    from sklearn.utils.linear_assignment_ import _hungarian
    matrices = [
        # Square
        ([[400, 150, 400],
          [400, 450, 600],
          [300, 225, 300]],
         850  # expected cost
         ),

        # Rectangular variant
        ([[400, 150, 400, 1],
          [400, 450, 600, 2],
          [300, 225, 300, 3]],
         452  # expected cost
         ),

        # Square
        ([[10, 10,  8],
          [9,  8,  1],
          [9,  7,  4]],
         18
         ),

        # Rectangular variant
        ([[10, 10,  8, 11],
          [9, 8, 1, 1],
          [9, 7, 4, 10]],
         15
         ),

        # n == 2, m == 0 matrix
        ([[], []],
         0
         ),
    ]

    for cost_matrix, expected_total in matrices:
        cost_matrix = np.array(cost_matrix)
        indexes = _hungarian(cost_matrix)
        total_cost = 0
        for r, c in indexes:
            x = cost_matrix[r, c]
            total_cost += x
        assert expected_total == total_cost

        indexes = _hungarian(cost_matrix.T)
        total_cost = 0
        for c, r in indexes:
            x = cost_matrix[r, c]
            total_cost += x
        assert expected_total == total_cost
