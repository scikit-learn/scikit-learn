import numpy as np
import scipy.sparse as sp
from numpy.testing import assert_array_almost_equal, assert_array_equal

from sklearn.datasets import make_classification
from sklearn.utils.sparsefuncs import assign_rows_csr, mean_variance_axis0


def test_mean_variance_axis0():
    X, _ = make_classification(5, 4, random_state=0)
    # Sparsify the array a little bit
    X[0, 0] = 0
    X[2, 1] = 0
    X[4, 3] = 0
    X_lil = sp.lil_matrix(X)
    X_lil[1, 0] = 0
    X[1, 0] = 0
    X_csr = sp.csr_matrix(X_lil)

    X_means, X_vars = mean_variance_axis0(X_csr)
    assert_array_almost_equal(X_means, np.mean(X, axis=0))
    assert_array_almost_equal(X_vars, np.var(X, axis=0))

    X_csc = sp.csc_matrix(X_lil)
    X_means, X_vars = mean_variance_axis0(X_csc)

    assert_array_almost_equal(X_means, np.mean(X, axis=0))
    assert_array_almost_equal(X_vars, np.var(X, axis=0))


def test_densify_rows():
    X = sp.csr_matrix([[0, 3, 0],
                       [2, 4, 0],
                       [0, 0, 0],
                       [9, 8, 7],
                       [4, 0, 5]], dtype=np.float64)
    rows = np.array([0, 2, 3], dtype=np.int)
    out = np.ones((rows.shape[0], X.shape[1]), dtype=np.float64)

    assign_rows_csr(X, rows, np.arange(out.shape[0])[::-1], out)
    assert_array_equal(out, X[rows].toarray()[::-1])
