import numpy as np
import scipy.sparse as sp

from scipy import linalg
from numpy.testing import assert_array_almost_equal, assert_array_equal

from sklearn.datasets import make_classification
from sklearn.utils.sparsefuncs import (mean_variance_axis,
                                       inplace_column_scale,
                                       inplace_row_scale,
                                       inplace_swap_row, inplace_swap_column,
                                       min_max_axis,
                                       count_nonzero, csc_median_axis_0)
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

    X_means, X_vars = mean_variance_axis(X_csr, axis=0)
    assert_array_almost_equal(X_means, np.mean(X, axis=0))
    assert_array_almost_equal(X_vars, np.var(X, axis=0))

    X_csc = sp.csc_matrix(X_lil)
    X_means, X_vars = mean_variance_axis(X_csc, axis=0)

    assert_array_almost_equal(X_means, np.mean(X, axis=0))
    assert_array_almost_equal(X_vars, np.var(X, axis=0))
    assert_raises(TypeError, mean_variance_axis, X_lil, axis=0)

    X = X.astype(np.float32)
    X_csr = X_csr.astype(np.float32)
    X_csc = X_csr.astype(np.float32)
    X_means, X_vars = mean_variance_axis(X_csr, axis=0)
    assert_array_almost_equal(X_means, np.mean(X, axis=0))
    assert_array_almost_equal(X_vars, np.var(X, axis=0))
    X_means, X_vars = mean_variance_axis(X_csc, axis=0)
    assert_array_almost_equal(X_means, np.mean(X, axis=0))
    assert_array_almost_equal(X_vars, np.var(X, axis=0))
    assert_raises(TypeError, mean_variance_axis, X_lil, axis=0)


def test_mean_variance_illegal_axis():
    X, _ = make_classification(5, 4, random_state=0)
    # Sparsify the array a little bit
    X[0, 0] = 0
    X[2, 1] = 0
    X[4, 3] = 0
    X_csr = sp.csr_matrix(X)
    assert_raises(ValueError, mean_variance_axis, X_csr, axis=-3)
    assert_raises(ValueError, mean_variance_axis, X_csr, axis=2)
    assert_raises(ValueError, mean_variance_axis, X_csr, axis=-1)


def test_mean_variance_axis1():
    X, _ = make_classification(5, 4, random_state=0)
    # Sparsify the array a little bit
    X[0, 0] = 0
    X[2, 1] = 0
    X[4, 3] = 0
    X_lil = sp.lil_matrix(X)
    X_lil[1, 0] = 0
    X[1, 0] = 0
    X_csr = sp.csr_matrix(X_lil)

    X_means, X_vars = mean_variance_axis(X_csr, axis=1)
    assert_array_almost_equal(X_means, np.mean(X, axis=1))
    assert_array_almost_equal(X_vars, np.var(X, axis=1))

    X_csc = sp.csc_matrix(X_lil)
    X_means, X_vars = mean_variance_axis(X_csc, axis=1)

    assert_array_almost_equal(X_means, np.mean(X, axis=1))
    assert_array_almost_equal(X_vars, np.var(X, axis=1))
    assert_raises(TypeError, mean_variance_axis, X_lil, axis=1)

    X = X.astype(np.float32)
    X_csr = X_csr.astype(np.float32)
    X_csc = X_csr.astype(np.float32)
    X_means, X_vars = mean_variance_axis(X_csr, axis=1)
    assert_array_almost_equal(X_means, np.mean(X, axis=1))
    assert_array_almost_equal(X_vars, np.var(X, axis=1))
    X_means, X_vars = mean_variance_axis(X_csc, axis=1)
    assert_array_almost_equal(X_means, np.mean(X, axis=1))
    assert_array_almost_equal(X_vars, np.var(X, axis=1))
    assert_raises(TypeError, mean_variance_axis, X_lil, axis=1)


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

    X = X.astype(np.float32)
    scale = scale.astype(np.float32)
    Xr = X.tocsr()
    Xc = X.tocsc()
    XA = X.toarray()
    XA *= scale
    inplace_column_scale(Xc, scale)
    inplace_column_scale(Xr, scale)
    assert_array_almost_equal(Xr.toarray(), Xc.toarray())
    assert_array_almost_equal(XA, Xc.toarray())
    assert_array_almost_equal(XA, Xr.toarray())
    assert_raises(TypeError, inplace_column_scale, X.tolil(), scale)


def test_inplace_row_scale():
    rng = np.random.RandomState(0)
    X = sp.rand(100, 200, 0.05)
    Xr = X.tocsr()
    Xc = X.tocsc()
    XA = X.toarray()
    scale = rng.rand(100)
    XA *= scale.reshape(-1, 1)

    inplace_row_scale(Xc, scale)
    inplace_row_scale(Xr, scale)
    assert_array_almost_equal(Xr.toarray(), Xc.toarray())
    assert_array_almost_equal(XA, Xc.toarray())
    assert_array_almost_equal(XA, Xr.toarray())
    assert_raises(TypeError, inplace_column_scale, X.tolil(), scale)

    X = X.astype(np.float32)
    scale = scale.astype(np.float32)
    Xr = X.tocsr()
    Xc = X.tocsc()
    XA = X.toarray()
    XA *= scale.reshape(-1, 1)
    inplace_row_scale(Xc, scale)
    inplace_row_scale(Xr, scale)
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

    X = np.array([[0, 3, 0],
                  [2, 4, 0],
                  [0, 0, 0],
                  [9, 8, 7],
                  [4, 0, 5]], dtype=np.float32)
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

    X = np.array([[0, 3, 0],
                  [2, 4, 0],
                  [0, 0, 0],
                  [9, 8, 7],
                  [4, 0, 5]], dtype=np.float32)
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


def test_min_max_axis0():
    X = np.array([[0, 3, 0],
                  [2, -1, 0],
                  [0, 0, 0],
                  [9, 8, 7],
                  [4, 0, 5]], dtype=np.float64)
    X_csr = sp.csr_matrix(X)
    X_csc = sp.csc_matrix(X)

    mins_csr, maxs_csr = min_max_axis(X_csr, axis=0)
    assert_array_equal(mins_csr, X.min(axis=0))
    assert_array_equal(maxs_csr, X.max(axis=0))

    mins_csc, maxs_csc = min_max_axis(X_csc, axis=0)
    assert_array_equal(mins_csc, X.min(axis=0))
    assert_array_equal(maxs_csc, X.max(axis=0))

    X = X.astype(np.float32)
    X_csr = sp.csr_matrix(X)
    X_csc = sp.csc_matrix(X)
    mins_csr, maxs_csr = min_max_axis(X_csr, axis=0)
    assert_array_equal(mins_csr, X.min(axis=0))
    assert_array_equal(maxs_csr, X.max(axis=0))
    mins_csc, maxs_csc = min_max_axis(X_csc, axis=0)
    assert_array_equal(mins_csc, X.min(axis=0))
    assert_array_equal(maxs_csc, X.max(axis=0))


def test_min_max_axis1():
    X = np.array([[0, 3, 0],
                  [2, -1, 0],
                  [0, 0, 0],
                  [9, 8, 7],
                  [4, 0, 5]], dtype=np.float64)
    X_csr = sp.csr_matrix(X)
    X_csc = sp.csc_matrix(X)

    mins_csr, maxs_csr = min_max_axis(X_csr, axis=1)
    assert_array_equal(mins_csr, X.min(axis=1))
    assert_array_equal(maxs_csr, X.max(axis=1))

    mins_csc, maxs_csc = min_max_axis(X_csc, axis=1)
    assert_array_equal(mins_csc, X.min(axis=1))
    assert_array_equal(maxs_csc, X.max(axis=1))

    X = X.astype(np.float32)
    X_csr = sp.csr_matrix(X)
    X_csc = sp.csc_matrix(X)
    mins_csr, maxs_csr = min_max_axis(X_csr, axis=1)
    assert_array_equal(mins_csr, X.min(axis=1))
    assert_array_equal(maxs_csr, X.max(axis=1))
    mins_csc, maxs_csc = min_max_axis(X_csc, axis=1)
    assert_array_equal(mins_csc, X.min(axis=1))
    assert_array_equal(maxs_csc, X.max(axis=1))


def test_min_max_axis_errors():
    X = np.array([[0, 3, 0],
                  [2, -1, 0],
                  [0, 0, 0],
                  [9, 8, 7],
                  [4, 0, 5]], dtype=np.float64)
    X_csr = sp.csr_matrix(X)
    X_csc = sp.csc_matrix(X)
    assert_raises(TypeError, min_max_axis, X_csr.tolil(), axis=0)
    assert_raises(ValueError, min_max_axis, X_csr, axis=2)
    assert_raises(ValueError, min_max_axis, X_csc, axis=-3)


def test_count_nonzero():
    X = np.array([[0, 3, 0],
                  [2, -1, 0],
                  [0, 0, 0],
                  [9, 8, 7],
                  [4, 0, 5]], dtype=np.float64)
    X_csr = sp.csr_matrix(X)
    X_csc = sp.csc_matrix(X)
    X_nonzero = X != 0
    sample_weight = [.5, .2, .3, .1, .1]
    X_nonzero_weighted = X_nonzero * np.array(sample_weight)[:, None]

    for axis in [0, 1, -1, -2, None]:
        assert_array_almost_equal(count_nonzero(X_csr, axis=axis),
                                  X_nonzero.sum(axis=axis))
        assert_array_almost_equal(count_nonzero(X_csr, axis=axis,
                                                sample_weight=sample_weight),
                                  X_nonzero_weighted.sum(axis=axis))

    assert_raises(TypeError, count_nonzero, X_csc)
    assert_raises(ValueError, count_nonzero, X_csr, axis=2)


def test_csc_row_median():
    # Test csc_row_median actually calculates the median.

    # Test that it gives the same output when X is dense.
    rng = np.random.RandomState(0)
    X = rng.rand(100, 50)
    dense_median = np.median(X, axis=0)
    csc = sp.csc_matrix(X)
    sparse_median = csc_median_axis_0(csc)
    assert_array_equal(sparse_median, dense_median)

    # Test that it gives the same output when X is sparse
    X = rng.rand(51, 100)
    X[X < 0.7] = 0.0
    ind = rng.randint(0, 50, 10)
    X[ind] = -X[ind]
    csc = sp.csc_matrix(X)
    dense_median = np.median(X, axis=0)
    sparse_median = csc_median_axis_0(csc)
    assert_array_equal(sparse_median, dense_median)

    # Test for toy data.
    X = [[0, -2], [-1, -1], [1, 0], [2, 1]]
    csc = sp.csc_matrix(X)
    assert_array_equal(csc_median_axis_0(csc), np.array([0.5, -0.5]))
    X = [[0, -2], [-1, -5], [1, -3]]
    csc = sp.csc_matrix(X)
    assert_array_equal(csc_median_axis_0(csc), np.array([0., -3]))

    # Test that it raises an Error for non-csc matrices.
    assert_raises(TypeError, csc_median_axis_0, sp.csr_matrix(X))
