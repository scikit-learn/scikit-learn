"""Tests for _sketches.py."""

from __future__ import division, print_function, absolute_import
import numpy as np
from scipy.linalg import clarkson_woodruff_transform

from numpy.testing import assert_


def make_random_dense_gaussian_matrix(n_rows, n_columns, mu=0, sigma=0.01):
    """
    Make some random data with Gaussian distributed values
    """
    np.random.seed(142352345)
    res = np.random.normal(mu, sigma, n_rows*n_columns)
    return np.reshape(res, (n_rows, n_columns))


class TestClarksonWoodruffTransform(object):
    """
    Testing the Clarkson Woodruff Transform
    """
    # Big dense matrix dimensions
    n_matrix_rows = 2000
    n_matrix_columns = 100

    # Sketch matrix dimensions
    n_sketch_rows = 100

    # Error threshold
    threshold = 0.1

    dense_big_matrix = make_random_dense_gaussian_matrix(n_matrix_rows,
                                                         n_matrix_columns)

    def test_sketch_dimensions(self):
        sketch = clarkson_woodruff_transform(self.dense_big_matrix,
                                             self.n_sketch_rows)

        assert_(sketch.shape == (self.n_sketch_rows,
                                 self.dense_big_matrix.shape[1]))

    def test_sketch_rows_norm(self):
        # Given the probabilistic nature of the sketches
        # we run the 'test' multiple times and check that
        # we pass all/almost all the tries
        n_errors = 0

        seeds = [1755490010, 934377150, 1391612830, 1752708722, 2008891431,
                 1302443994, 1521083269, 1501189312, 1126232505, 1533465685]

        for seed_ in seeds:
            sketch = clarkson_woodruff_transform(self.dense_big_matrix,
                                                 self.n_sketch_rows, seed_)

            # We could use other norms (like L2)
            err = np.linalg.norm(self.dense_big_matrix) - np.linalg.norm(sketch)
            if err > self.threshold:
                n_errors += 1

        assert_(n_errors == 0)
