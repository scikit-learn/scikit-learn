# Author: Brian M. Clapper, G Varoquaux
# License: BSD

import numpy as np

# XXX we should be testing the public API here
from scipy.optimize import linear_sum_assignment


def test_hungarian():
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
        r, c = linear_sum_assignment(cost_matrix)
        total_cost = cost_matrix[r, c].sum()
        assert expected_total == total_cost

        r, c = linear_sum_assignment(cost_matrix.T)
        total_cost = cost_matrix[c, r].sum()
        assert expected_total == total_cost
