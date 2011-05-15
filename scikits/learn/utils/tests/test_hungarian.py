# Author: Brian M. Clapper, G Varoquaux
# LICENSE: BSD

import numpy as np

from scikits.learn.utils.hungarian import _Hungarian

def test_hungarian():
    matrices = [
                # Square
                ([[400, 150, 400],
                  [400, 450, 600],
                  [300, 225, 300]],
                 850 # expected cost
                ),

                ## Rectangular variant
                #([[400, 150, 400, 1],
                #  [400, 450, 600, 2],
                #  [300, 225, 300, 3]],
                # 452 # expected cost
                #),

                # Square
                ([[10, 10,  8],
                  [ 9,  8,  1],
                  [ 9,  7,  4]],
                 18
                ),

                ## Rectangular variant
                #([[10, 10,  8, 11],
                #  [ 9,  8,  1, 1],
                #  [ 9,  7,  4, 10]],
                # 15
                #),
               ]

    m = _Hungarian()
    for cost_matrix, expected_total in matrices:
        print np.array(cost_matrix)
        cost_matrix = np.array(cost_matrix)
        indexes = m.compute(cost_matrix)
        total_cost = 0
        for r, c in indexes:
            x = cost_matrix[r, c]
            total_cost += x
        assert expected_total == total_cost

