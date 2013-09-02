"""
Tests for EAC clustering algorithm
"""

import pickle

import numpy as np
from collections import defaultdict

from sklearn.utils.testing import assert_equal
from sklearn.cluster.eac import _update_coassociation_matrix


def test_coassociation_matrix_building():
    """Tests that the coassociation matrix builds properly."""
    C = defaultdict(int)
    n_samples = 4
    test_labels = np.array([[0, 0, 1, 1],
                            [1, 1, 0, 0],
                            [0, 1, 1, 0],
                            [1, 1, 1, 1]])
    for labels in test_labels:
        C = _update_coassociation_matrix(C, labels)
    for i in range(n_samples):
        C[(i,i)] = len(test_labels)
    C_expected = defaultdict(int)
    # format: (i, j, C[i][j])
    data = [(0, 0, 4), (0, 1, 3), (0, 2, 1), (0, 3, 2), (1, 1, 4), (1, 2, 2),
            (1, 3, 1), (2, 2, 4), (2, 3, 3), (3, 3, 4)]
    for x, y, v in data:
        C_expected[(x, y)] = v
    assert_equal(C, C_expected)

