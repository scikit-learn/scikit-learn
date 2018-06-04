from __future__ import absolute_import

import pytest
pytest.importorskip('numpy')
pytest.importorskip('scipy')

import numpy as np
import scipy.linalg

import dask.array as da
from dask.array.linalg import tsqr, svd_compressed, qr, svd
from dask.array.utils import assert_eq, same_keys


def test_tsqr_regular_blocks():
    m, n = 20, 10
    mat = np.random.rand(m, n)
    data = da.from_array(mat, chunks=(10, n), name='A')

    q, r = tsqr(data)
    q = np.array(q)
    r = np.array(r)

    assert_eq(mat, np.dot(q, r))  # accuracy check
    assert_eq(np.eye(n, n), np.dot(q.T, q))  # q must be orthonormal
    assert_eq(r, np.triu(r))  # r must be upper triangular


def test_tsqr_irregular_blocks():
    m, n = 20, 10
    mat = np.random.rand(m, n)
    data = da.from_array(mat, chunks=(3, n), name='A')[1:]
    mat2 = mat[1:, :]

    q, r = tsqr(data)
    q = np.array(q)
    r = np.array(r)

    assert_eq(mat2, np.dot(q, r))  # accuracy check
    assert_eq(np.eye(n, n), np.dot(q.T, q))  # q must be orthonormal
    assert_eq(r, np.triu(r))  # r must be upper triangular


def test_tsqr_svd_regular_blocks():
    m, n = 20, 10
    mat = np.random.rand(m, n)
    data = da.from_array(mat, chunks=(10, n), name='A')

    u, s, vt = tsqr(data, compute_svd=True)
    u = np.array(u)
    s = np.array(s)
    vt = np.array(vt)
    usvt = np.dot(u, np.dot(np.diag(s), vt))

    s_exact = np.linalg.svd(mat)[1]

    assert_eq(mat, usvt)  # accuracy check
    assert_eq(np.eye(n, n), np.dot(u.T, u))  # u must be orthonormal
    assert_eq(np.eye(n, n), np.dot(vt, vt.T))  # v must be orthonormal
    assert_eq(s, s_exact)  # s must contain the singular values


def test_tsqr_svd_irregular_blocks():
    m, n = 20, 10
    mat = np.random.rand(m, n)
    data = da.from_array(mat, chunks=(3, n), name='A')[1:]
    mat2 = mat[1:, :]

    u, s, vt = tsqr(data, compute_svd=True)
    u = np.array(u)
    s = np.array(s)
    vt = np.array(vt)
    usvt = np.dot(u, np.dot(np.diag(s), vt))

    s_exact = np.linalg.svd(mat2)[1]

    assert_eq(mat2, usvt)  # accuracy check
    assert_eq(np.eye(n, n), np.dot(u.T, u))  # u must be orthonormal
    assert_eq(np.eye(n, n), np.dot(vt, vt.T))  # v must be orthonormal
    assert_eq(s, s_exact)  # s must contain the singular values


def test_linalg_consistent_names():
    m, n = 20, 10
    mat = np.random.rand(m, n)
    data = da.from_array(mat, chunks=(10, n), name='A')

    q1, r1 = qr(data)
    q2, r2 = qr(data)
    assert same_keys(q1, q2)
    assert same_keys(r1, r2)

    u1, s1, v1 = svd(data)
    u2, s2, v2 = svd(data)
    assert same_keys(u1, u2)
    assert same_keys(s1, s2)
    assert same_keys(v1, v2)


@pytest.mark.slow
def test_svd_compressed():
    m, n = 2000, 250
    r = 10
    np.random.seed(4321)
    mat1 = np.random.randn(m, r)
    mat2 = np.random.randn(r, n)
    mat = mat1.dot(mat2)
    data = da.from_array(mat, chunks=(500, 50))

    u, s, vt = svd_compressed(data, r, seed=4321, n_power_iter=2)
    u, s, vt = da.compute(u, s, vt)

    usvt = np.dot(u, np.dot(np.diag(s), vt))

    tol = 0.2
    assert_eq(np.linalg.norm(usvt),
              np.linalg.norm(mat),
              rtol=tol, atol=tol)  # average accuracy check

    u = u[:, :r]
    s = s[:r]
    vt = vt[:r, :]

    s_exact = np.linalg.svd(mat)[1]
    s_exact = s_exact[:r]

    assert_eq(np.eye(r, r), np.dot(u.T, u))  # u must be orthonormal
    assert_eq(np.eye(r, r), np.dot(vt, vt.T))  # v must be orthonormal
    assert_eq(s, s_exact)  # s must contain the singular values


def test_svd_compressed_deterministic():
    m, n = 30, 25
    x = da.random.RandomState(1234).random_sample(size=(m, n), chunks=(5, 5))
    u, s, vt = svd_compressed(x, 3, seed=1234)
    u2, s2, vt2 = svd_compressed(x, 3, seed=1234)

    assert all(da.compute((u == u2).all(), (s == s2).all(), (vt == vt2).all()))


def _check_lu_result(p, l, u, A):
    assert np.allclose(p.dot(l).dot(u), A)

    # check triangulars
    assert np.allclose(l, np.tril(l.compute()))
    assert np.allclose(u, np.triu(u.compute()))


def test_lu_1():
    A1 = np.array([[7, 3, -1, 2], [3, 8, 1, -4],
                  [-1, 1, 4, -1], [2, -4, -1, 6] ])

    A2 = np.array([[7,  0,  0,  0,  0,  0],
                   [0,  8,  0,  0,  0,  0],
                   [0,  0,  4,  0,  0,  0],
                   [0,  0,  0,  6,  0,  0],
                   [0,  0,  0,  0,  3,  0],
                   [0,  0,  0,  0,  0,  5]])
    # without shuffle
    for A, chunk in zip([A1, A2], [2, 2]):
        dA = da.from_array(A, chunks=(chunk, chunk))
        p, l, u = scipy.linalg.lu(A)
        dp, dl, du = da.linalg.lu(dA)
        assert_eq(p, dp)
        assert_eq(l, dl)
        assert_eq(u, du)
        _check_lu_result(dp, dl, du, A)

    A3 = np.array([[ 7,  3,  2,  1,  4,  1],
                   [ 7, 11,  5,  2,  5,  2],
                   [21, 25, 16, 10, 16,  5],
                   [21, 41, 18, 13, 16, 11],
                   [14, 46, 23, 24, 21, 22],
                   [ 0, 56, 29, 17, 14, 8]])

    # with shuffle
    for A, chunk in zip([A3], [2]):
        dA = da.from_array(A, chunks=(chunk, chunk))
        p, l, u = scipy.linalg.lu(A)
        dp, dl, du = da.linalg.lu(dA)
        _check_lu_result(dp, dl, du, A)


@pytest.mark.slow
@pytest.mark.parametrize('size', [10, 20, 30, 50])
def test_lu_2(size):
    np.random.seed(10)
    A = np.random.randint(0, 10, (size, size))

    dA = da.from_array(A, chunks=(5, 5))
    dp, dl, du = da.linalg.lu(dA)
    _check_lu_result(dp, dl, du, A)


@pytest.mark.slow
@pytest.mark.parametrize('size', [50, 100, 200])
def test_lu_3(size):
    np.random.seed(10)
    A = np.random.randint(0, 10, (size, size))

    dA = da.from_array(A, chunks=(25, 25))
    dp, dl, du = da.linalg.lu(dA)
    _check_lu_result(dp, dl, du, A)


def test_lu_errors():
    A = np.random.randint(0, 11, (10, 10, 10))
    dA = da.from_array(A, chunks=(5, 5, 5))
    pytest.raises(ValueError, lambda: da.linalg.lu(dA))

    A = np.random.randint(0, 11, (10, 8))
    dA = da.from_array(A, chunks=(5, 4))
    pytest.raises(ValueError, lambda: da.linalg.lu(dA))

    A = np.random.randint(0, 11, (20, 20))
    dA = da.from_array(A, chunks=(5, 4))
    pytest.raises(ValueError, lambda: da.linalg.lu(dA))


@pytest.mark.parametrize(('shape', 'chunk'), [(20, 10), (50, 10), (70, 20)])
def test_solve_triangular_vector(shape, chunk):
    np.random.seed(1)

    A = np.random.randint(1, 11, (shape, shape))
    b = np.random.randint(1, 11, shape)

    # upper
    Au = np.triu(A)
    dAu = da.from_array(Au, (chunk, chunk))
    db = da.from_array(b, chunk)
    res = da.linalg.solve_triangular(dAu, db)
    assert_eq(res, scipy.linalg.solve_triangular(Au, b))
    assert_eq(dAu.dot(res), b.astype(float))

    # lower
    Al = np.tril(A)
    dAl = da.from_array(Al, (chunk, chunk))
    db = da.from_array(b, chunk)
    res = da.linalg.solve_triangular(dAl, db, lower=True)
    assert_eq(res, scipy.linalg.solve_triangular(Al, b, lower=True))
    assert_eq(dAl.dot(res), b.astype(float))


@pytest.mark.parametrize(('shape', 'chunk'), [(20, 10), (50, 10), (50, 20)])
def test_solve_triangular_matrix(shape, chunk):
    np.random.seed(1)

    A = np.random.randint(1, 10, (shape, shape))
    b = np.random.randint(1, 10, (shape, 5))

    # upper
    Au = np.triu(A)
    dAu = da.from_array(Au, (chunk, chunk))
    db = da.from_array(b, (chunk, 5))
    res = da.linalg.solve_triangular(dAu, db)
    assert_eq(res, scipy.linalg.solve_triangular(Au, b))
    assert_eq(dAu.dot(res), b.astype(float))

    # lower
    Al = np.tril(A)
    dAl = da.from_array(Al, (chunk, chunk))
    db = da.from_array(b, (chunk, 5))
    res = da.linalg.solve_triangular(dAl, db, lower=True)
    assert_eq(res, scipy.linalg.solve_triangular(Al, b, lower=True))
    assert_eq(dAl.dot(res), b.astype(float))


@pytest.mark.parametrize(('shape', 'chunk'), [(20, 10), (50, 10), (50, 20)])
def test_solve_triangular_matrix2(shape, chunk):
    np.random.seed(1)

    A = np.random.randint(1, 10, (shape, shape))
    b = np.random.randint(1, 10, (shape, shape))

    # upper
    Au = np.triu(A)
    dAu = da.from_array(Au, (chunk, chunk))
    db = da.from_array(b, (chunk, chunk))
    res = da.linalg.solve_triangular(dAu, db)
    assert_eq(res, scipy.linalg.solve_triangular(Au, b))
    assert_eq(dAu.dot(res), b.astype(float))

    # lower
    Al = np.tril(A)
    dAl = da.from_array(Al, (chunk, chunk))
    db = da.from_array(b, (chunk, chunk))
    res = da.linalg.solve_triangular(dAl, db, lower=True)
    assert_eq(res, scipy.linalg.solve_triangular(Al, b, lower=True))
    assert_eq(dAl.dot(res), b.astype(float))


def test_solve_triangular_errors():
    A = np.random.randint(0, 10, (10, 10, 10))
    b = np.random.randint(1, 10, 10)
    dA = da.from_array(A, chunks=(5, 5, 5))
    db = da.from_array(b, chunks=5)
    pytest.raises(ValueError, lambda: da.linalg.solve_triangular(dA, db))

    A = np.random.randint(0, 10, (10, 10))
    b = np.random.randint(1, 10, 10)
    dA = da.from_array(A, chunks=(3, 3))
    db = da.from_array(b, chunks=5)
    pytest.raises(ValueError, lambda: da.linalg.solve_triangular(dA, db))


@pytest.mark.parametrize(('shape', 'chunk'), [(20, 10), (50, 10)])
def test_solve(shape, chunk):
    np.random.seed(1)

    A = np.random.randint(1, 10, (shape, shape))
    dA = da.from_array(A, (chunk, chunk))

    # vector
    b = np.random.randint(1, 10, shape)
    db = da.from_array(b, chunk)

    res = da.linalg.solve(dA, db)
    assert_eq(res, scipy.linalg.solve(A, b))
    assert_eq(dA.dot(res), b.astype(float))

    # tall-and-skinny matrix
    b = np.random.randint(1, 10, (shape, 5))
    db = da.from_array(b, (chunk, 5))

    res = da.linalg.solve(dA, db)
    assert_eq(res, scipy.linalg.solve(A, b))
    assert_eq(dA.dot(res), b.astype(float))

    # matrix
    b = np.random.randint(1, 10, (shape, shape))
    db = da.from_array(b, (chunk, chunk))

    res = da.linalg.solve(dA, db)
    assert_eq(res, scipy.linalg.solve(A, b))
    assert_eq(dA.dot(res), b.astype(float))


@pytest.mark.parametrize(('shape', 'chunk'), [(20, 10), (50, 10)])
def test_inv(shape, chunk):
    np.random.seed(1)

    A = np.random.randint(1, 10, (shape, shape))
    dA = da.from_array(A, (chunk, chunk))

    res = da.linalg.inv(dA)
    assert_eq(res, scipy.linalg.inv(A))
    assert_eq(dA.dot(res), np.eye(shape, dtype=float))


def _get_symmat(size):
    np.random.seed(1)
    A = np.random.randint(1, 21, (size, size))
    lA = np.tril(A)
    return lA.dot(lA.T)


@pytest.mark.parametrize(('shape', 'chunk'), [(20, 10), (30, 6)])
def test_solve_sym_pos(shape, chunk):
    np.random.seed(1)

    A = _get_symmat(shape)
    dA = da.from_array(A, (chunk, chunk))

    # vector
    b = np.random.randint(1, 10, shape)
    db = da.from_array(b, chunk)

    res = da.linalg.solve(dA, db, sym_pos=True)
    assert_eq(res, scipy.linalg.solve(A, b, sym_pos=True))
    assert_eq(dA.dot(res), b.astype(float))

    # tall-and-skinny matrix
    b = np.random.randint(1, 10, (shape, 5))
    db = da.from_array(b, (chunk, 5))

    res = da.linalg.solve(dA, db, sym_pos=True)
    assert_eq(res, scipy.linalg.solve(A, b, sym_pos=True))
    assert_eq(dA.dot(res), b.astype(float))

    # matrix
    b = np.random.randint(1, 10, (shape, shape))
    db = da.from_array(b, (chunk, chunk))

    res = da.linalg.solve(dA, db, sym_pos=True)
    assert_eq(res, scipy.linalg.solve(A, b, sym_pos=True))
    assert_eq(dA.dot(res), b.astype(float))


@pytest.mark.parametrize(('shape', 'chunk'), [(20, 10), (12, 3), (30, 3), (30, 6)])
def test_cholesky(shape, chunk):

    A = _get_symmat(shape)
    dA = da.from_array(A, (chunk, chunk))
    assert_eq(da.linalg.cholesky(dA), scipy.linalg.cholesky(A))
    assert_eq(da.linalg.cholesky(dA, lower=True), scipy.linalg.cholesky(A, lower=True))


@pytest.mark.parametrize(("nrow", "ncol", "chunk"),
                         [(20, 10, 5), (100, 10, 10)])
def test_lstsq(nrow, ncol, chunk):
    np.random.seed(1)
    A = np.random.randint(1, 20, (nrow, ncol))
    b = np.random.randint(1, 20, nrow)

    dA = da.from_array(A, (chunk, ncol))
    db = da.from_array(b, chunk)

    x, r, rank, s = np.linalg.lstsq(A, b)
    dx, dr, drank, ds = da.linalg.lstsq(dA, db)

    assert_eq(dx, x)
    assert_eq(dr, r)
    assert drank.compute() == rank
    assert_eq(ds, s)

    # reduce rank causes multicollinearity, only compare rank
    A[:, 1] = A[:, 2]
    dA = da.from_array(A, (chunk, ncol))
    db = da.from_array(b, chunk)
    x, r, rank, s = np.linalg.lstsq(A, b,
                                    rcond=np.finfo(np.double).eps * max(nrow,
                                                                        ncol))
    assert rank == ncol - 1
    dx, dr, drank, ds = da.linalg.lstsq(dA, db)
    assert drank.compute() == rank


def test_no_chunks_svd():
    x = np.random.random((100, 10))
    u, s, v = np.linalg.svd(x, full_matrices=0)

    for chunks in [((np.nan,) * 10, (10,)),
                   ((np.nan,) * 10, (np.nan,))]:
        dx = da.from_array(x, chunks=(10, 10))
        dx._chunks = chunks

        du, ds, dv = da.linalg.svd(dx)

        assert_eq(s, ds)
        assert_eq(u.dot(np.diag(s)).dot(v),
                  du.dot(da.diag(ds)).dot(dv))
        assert_eq(du.T.dot(du), np.eye(10))
        assert_eq(dv.T.dot(dv), np.eye(10))

        dx = da.from_array(x, chunks=(10, 10))
        dx._chunks = ((np.nan,) * 10, (np.nan,))
        assert_eq(abs(v), abs(dv))
        assert_eq(abs(u), abs(du))


@pytest.mark.parametrize("shape, chunks, axis", [
    [(5,), (2,), None],
    [(5,), (2,), 0],
    [(5,), (2,), (0,)],
    [(5, 6), (2, 2), None],
    [(5, 6), (2, 2), 0],
    [(5, 6), (2, 2), 1],
    [(5, 6), (2, 2), (0, 1)],
    [(5, 6), (2, 2), (1, 0)],
])
@pytest.mark.parametrize("norm", [
    None,
    1,
    -1,
    np.inf,
    -np.inf,
])
@pytest.mark.parametrize("keepdims", [
    False,
    True,
])
def test_norm_any_ndim(shape, chunks, axis, norm, keepdims):
    a = np.random.random(shape)
    d = da.from_array(a, chunks=chunks)

    a_r = np.linalg.norm(a, ord=norm, axis=axis, keepdims=keepdims)
    d_r = da.linalg.norm(d, ord=norm, axis=axis, keepdims=keepdims)

    assert_eq(a_r, d_r)


@pytest.mark.parametrize("shape, chunks, axis", [
    [(5,), (2,), None],
    [(5,), (2,), 0],
    [(5,), (2,), (0,)],
])
@pytest.mark.parametrize("norm", [
    0,
    2,
    -2,
    0.5,
])
@pytest.mark.parametrize("keepdims", [
    False,
    True,
])
def test_norm_1dim(shape, chunks, axis, norm, keepdims):
    a = np.random.random(shape)
    d = da.from_array(a, chunks=chunks)

    a_r = np.linalg.norm(a, ord=norm, axis=axis, keepdims=keepdims)
    d_r = da.linalg.norm(d, ord=norm, axis=axis, keepdims=keepdims)
    assert_eq(a_r, d_r)


@pytest.mark.parametrize("shape, chunks, axis", [
    [(5, 6), (2, 2), None],
    [(5, 6), (2, 2), (0, 1)],
    [(5, 6), (2, 2), (1, 0)],
])
@pytest.mark.parametrize("norm", [
    "fro",
    "nuc",
    2,
    -2
])
@pytest.mark.parametrize("keepdims", [
    False,
    True,
])
def test_norm_2dim(shape, chunks, axis, norm, keepdims):
    a = np.random.random(shape)
    d = da.from_array(a, chunks=chunks)

    # Need one chunk on last dimension for svd.
    if norm == "nuc" or norm == 2 or norm == -2:
        d = d.rechunk((d.chunks[0], d.shape[1]))

    a_r = np.linalg.norm(a, ord=norm, axis=axis, keepdims=keepdims)
    d_r = da.linalg.norm(d, ord=norm, axis=axis, keepdims=keepdims)

    assert_eq(a_r, d_r)
