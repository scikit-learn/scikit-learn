import numpy as np
import scipy.sparse as sp

from scipy import linalg
from numpy.testing import assert_array_almost_equal, assert_array_equal

from sklearn.datasets import make_classification
from sklearn.utils.sparsefuncs import (mean_variance_axis0,
                                       inplace_column_scale,
                                       inplace_swap_row, inplace_swap_column)
from sklearn.utils.sparsefuncs_fast import assign_rows_csr
from sklearn.utils.testing import assert_raises

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
    assert_raises(TypeError, mean_variance_axis0, X_lil)


def test_densify_rows():
    X = sp.csr_matrix([[0, 3, 0],
                       [2, 4, 0],
                       [0, 0, 0],
                       [9, 8, 7],
                       [4, 0, 5]], dtype=np.float64)
    rows = np.array([0, 2, 3], dtype=np.intp)
    out = np.ones((rows.shape[0], X.shape[1]), dtype=np.float64)

    assign_rows_csr(X, rows, 
                    np.arange(out.shape[0], dtype=np.intp)[::-1], out)
    assert_array_equal(out, X[rows].toarray()[::-1])


def test_inplace_column_scale():
    rng = np.random.RandomState(0)
    X = sp.rand(100, 200, 0.05)
    Xr = X.tocsr()
    Xc = X.tocsc()
    XA = X.toarray()
    scale = rng.rand(200)
    XA *= scale

    inplace_column_scale(Xc, scale)
    inplace_column_scale(Xr, scale)
    assert_array_almost_equal(Xr.toarray(), Xc.toarray())
    assert_array_almost_equal(XA, Xc.toarray())
    assert_array_almost_equal(XA, Xr.toarray())
    assert_raises(TypeError, inplace_column_scale, X.tolil(), scale)


def test_inplace_swap_row():
    X = np.array([[0, 3, 0],
                  [2, 4, 0],
                  [0, 0, 0],
                  [9, 8, 7],
                  [4, 0, 5]], dtype=np.float64)
    X_csr = sp.csr_matrix(X)
    X_csc = sp.csc_matrix(X)

    swap = linalg.get_blas_funcs(('swap',), (X,))
    swap = swap[0]
    X[0], X[-1] = swap(X[0], X[-1])
    inplace_swap_row(X_csr, 0, -1)
    inplace_swap_row(X_csc, 0, -1)
    assert_array_equal(X_csr.toarray(), X_csc.toarray())
    assert_array_equal(X, X_csc.toarray())
    assert_array_equal(X, X_csr.toarray())

    X[2], X[3] = swap(X[2], X[3])
    inplace_swap_row(X_csr, 2, 3)
    inplace_swap_row(X_csc, 2, 3)
    assert_array_equal(X_csr.toarray(), X_csc.toarray())
    assert_array_equal(X, X_csc.toarray())
    assert_array_equal(X, X_csr.toarray())
    assert_raises(TypeError, inplace_swap_row, X_csr.tolil())


def test_inplace_swap_column():
    X = np.array([[0, 3, 0],
                  [2, 4, 0],
                  [0, 0, 0],
                  [9, 8, 7],
                  [4, 0, 5]], dtype=np.float64)
    X_csr = sp.csr_matrix(X)
    X_csc = sp.csc_matrix(X)

    swap = linalg.get_blas_funcs(('swap',), (X,))
    swap = swap[0]
    X[:, 0], X[:, -1] = swap(X[:, 0], X[:, -1])
    inplace_swap_column(X_csr, 0, -1)
    inplace_swap_column(X_csc, 0, -1)
    assert_array_equal(X_csr.toarray(), X_csc.toarray())
    assert_array_equal(X, X_csc.toarray())
    assert_array_equal(X, X_csr.toarray())

    X[:, 0], X[:, 1] = swap(X[:, 0], X[:, 1])
    inplace_swap_column(X_csr, 0, 1)
    inplace_swap_column(X_csc, 0, 1)
    assert_array_equal(X_csr.toarray(), X_csc.toarray())
    assert_array_equal(X, X_csc.toarray())
    assert_array_equal(X, X_csr.toarray())
    assert_raises(TypeError, inplace_swap_column, X_csr.tolil())
