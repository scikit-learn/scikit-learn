import numpy as np
import scipy.sparse as sp
from numpy.testing import assert_array_almost_equal

from sklearn.datasets import make_classification
from sklearn.utils.sparsefuncs import mean_variance_axis0


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
