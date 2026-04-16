#
# Created by: Pearu Peterson, April 2002
#

import math
import pytest
import numpy as np
from numpy.testing import (assert_equal, assert_almost_equal,
                           assert_array_almost_equal, assert_allclose)
from pytest import raises as assert_raises

from numpy import (arange, triu, tril, zeros, tril_indices, ones,
                   diag, append, eye, nonzero)

import scipy
from scipy.linalg import _fblas as fblas, get_blas_funcs, toeplitz, solve

try:
    from scipy.linalg import _cblas as cblas
except ImportError:
    cblas = None

REAL_DTYPES = [np.float32, np.float64]
COMPLEX_DTYPES = [np.complex64, np.complex128]
DTYPES = REAL_DTYPES + COMPLEX_DTYPES


def test_get_blas_funcs():
    # check that it returns Fortran code for arrays that are
    # fortran-ordered
    f1, f2, f3 = get_blas_funcs(
        ('axpy', 'axpy', 'axpy'),
        (np.empty((2, 2), dtype=np.complex64, order='F'),
         np.empty((2, 2), dtype=np.complex128, order='C'))
        )

    # get_blas_funcs will choose libraries depending on most generic
    # array
    assert_equal(f1.typecode, 'z')
    assert_equal(f2.typecode, 'z')
    if cblas is not None:
        assert_equal(f1.module_name, 'cblas')
        assert_equal(f2.module_name, 'cblas')

    # check defaults.
    f1 = get_blas_funcs('rotg')
    assert_equal(f1.typecode, 'd')

    # check also dtype interface
    f1 = get_blas_funcs('gemm', dtype=np.complex64)
    assert_equal(f1.typecode, 'c')
    f1 = get_blas_funcs('gemm', dtype='F')
    assert_equal(f1.typecode, 'c')

    # extended precision complex
    f1 = get_blas_funcs('gemm', dtype=np.clongdouble)
    assert_equal(f1.typecode, 'z')

    # check safe complex upcasting
    f1 = get_blas_funcs('axpy',
                        (np.empty((2, 2), dtype=np.float64),
                         np.empty((2, 2), dtype=np.complex64))
                        )
    assert_equal(f1.typecode, 'z')


def test_get_blas_funcs_alias():
    # check alias for get_blas_funcs
    f, g = get_blas_funcs(('nrm2', 'dot'), dtype=np.complex64)
    assert f.typecode == 'c'
    assert g.typecode == 'c'

    f, g, h = get_blas_funcs(('dot', 'dotc', 'dotu'), dtype=np.float64)
    assert f is g
    assert f is h


def parametrize_blas(mod, func_name, prefixes):
    if mod is None:
        return pytest.mark.skip(reason="cblas not available")
    params = []
    for prefix in prefixes:
        if 'z' in prefix:
            dtype = np.complex128
        elif 'c' in prefix:
            dtype = np.complex64
        elif 'd' in prefix:
            dtype = np.float64
        else:
            assert 's' in prefix
            dtype = np.float32

        f = getattr(mod, prefix + func_name)
        params.append(pytest.param(f, dtype, id=prefix + func_name))

    return pytest.mark.parametrize("f,dtype", params)


class TestCBLAS1Simple:
    @parametrize_blas(cblas, "axpy", "sdcz")
    def test_axpy(self, f, dtype):
        assert_array_almost_equal(f([1, 2, 3], [2, -1, 3], a=5),
                                  [7, 9, 18])
        if dtype in COMPLEX_DTYPES:
            assert_array_almost_equal(f([1, 2j, 3], [2, -1, 3], a=5),
                                      [7, 10j-1, 18])


class TestFBLAS1Simple:

    @parametrize_blas(fblas, "axpy", "sdcz")
    def test_axpy(self, f, dtype):
        assert_array_almost_equal(f([1, 2, 3], [2, -1, 3], a=5),
                                  [7, 9, 18])
        if dtype in COMPLEX_DTYPES:
            assert_array_almost_equal(f([1, 2j, 3], [2, -1, 3], a=5),
                                      [7, 10j-1, 18])

    @parametrize_blas(fblas, "copy", "sdcz")
    def test_copy(self, f, dtype):
        assert_array_almost_equal(f([3, 4, 5], [8]*3), [3, 4, 5])
        if dtype in COMPLEX_DTYPES:
            assert_array_almost_equal(f([3, 4j, 5+3j], [8]*3), [3, 4j, 5+3j])

    @parametrize_blas(fblas, "asum", ["s", "d", "sc", "dz"])
    def test_asum(self, f, dtype):
        assert_almost_equal(f([3, -4, 5]), 12)
        if dtype in COMPLEX_DTYPES:
            assert_almost_equal(f([3j, -4, 3-4j]), 14)

    @parametrize_blas(fblas, "dot", "sd")
    def test_dot(self, f, dtype):
        assert_almost_equal(f([3, -4, 5], [2, 5, 1]), -9)

    @parametrize_blas(fblas, "dotu", "cz")
    def test_dotu(self, f, dtype):
        assert_almost_equal(f([3j, -4, 3-4j], [2, 3, 1]), -9+2j)

    @parametrize_blas(fblas, "dotc", "cz")
    def test_dotc(self, f, dtype):
        assert_almost_equal(f([3j, -4, 3-4j], [2, 3j, 1]), 3-14j)

    @parametrize_blas(fblas, "nrm2", ["s", "d", "sc", "dz"])
    def test_nrm2(self, f, dtype):
        assert_almost_equal(f([3, -4, 5]), math.sqrt(50))
        if dtype in COMPLEX_DTYPES:
            assert_almost_equal(f([3j, -4, 3-4j]), math.sqrt(50))

    @parametrize_blas(fblas, "scal", ["s", "d", "cs", "zd"])
    def test_scal(self, f, dtype):
        assert_array_almost_equal(f(2, [3, -4, 5]), [6, -8, 10])
        if dtype in COMPLEX_DTYPES:
            assert_array_almost_equal(f(3, [3j, -4, 3-4j]), [9j, -12, 9-12j])

    @parametrize_blas(fblas, "swap", "sdcz")
    def test_swap(self, f, dtype):
        x, y = [2, 3, 1], [-2, 3, 7]
        x1, y1 = f(x, y)
        assert_array_almost_equal(x1, y)
        assert_array_almost_equal(y1, x)

        if dtype in COMPLEX_DTYPES:
            x, y = [2, 3j, 1], [-2, 3, 7-3j]
            x1, y1 = f(x, y)
            assert_array_almost_equal(x1, y)
            assert_array_almost_equal(y1, x)

    @parametrize_blas(fblas, "amax", ["is", "id", "ic", "iz"])
    def test_amax(self, f, dtype):
        assert_equal(f([-2, 4, 3]), 1)
        if dtype in COMPLEX_DTYPES:
            assert_equal(f([-5, 4+3j, 6]), 1)

    # XXX: need tests for rot,rotm,rotg,rotmg


class TestFBLAS2Simple:
    @parametrize_blas(fblas, "gemv", "sdcz")
    def test_gemv(self, f, dtype):
        assert_array_almost_equal(f(3, [[3]], [-4]), [-36])
        assert_array_almost_equal(f(3, [[3]], [-4], 3, [5]), [-21])
        if dtype in COMPLEX_DTYPES:
            assert_array_almost_equal(f(3j, [[3-4j]], [-4]), [-48-36j])
            assert_array_almost_equal(f(3j, [[3-4j]], [-4], 3, [5j]),
                                      [-48-21j])

    @parametrize_blas(fblas, "ger", "sd")
    def test_ger(self, f, dtype):
        assert_array_almost_equal(f(1, [1, 2], [3, 4]), [[3, 4], [6, 8]])
        assert_array_almost_equal(f(2, [1, 2, 3], [3, 4]),
                                  [[6, 8], [12, 16], [18, 24]])
        assert_array_almost_equal(f(1, [1, 2], [3, 4],
                                  a=[[1, 2], [3, 4]]), [[4, 6], [9, 12]])
        if dtype in COMPLEX_DTYPES:
            assert_array_almost_equal(f(1, [1j, 2], [3, 4]),
                                      [[3j, 4j], [6, 8]])
            assert_array_almost_equal(f(2, [1j, 2j, 3j], [3j, 4j]),
                                      [[6, 8], [12, 16], [18, 24]])

    @parametrize_blas(fblas, "geru", "cz")
    def test_geru(self, f, dtype):
        assert_array_almost_equal(f(1, [1j, 2], [3, 4]),
                                  [[3j, 4j], [6, 8]])
        assert_array_almost_equal(f(-2, [1j, 2j, 3j], [3j, 4j]),
                                  [[6, 8], [12, 16], [18, 24]])

    @parametrize_blas(fblas, "gerc", "cz")
    def test_gerc(self, f, dtype):
        assert_array_almost_equal(f(1, [1j, 2], [3, 4]),
                                  [[3j, 4j], [6, 8]])
        assert_array_almost_equal(f(2, [1j, 2j, 3j], [3j, 4j]),
                                  [[6, 8], [12, 16], [18, 24]])

    @parametrize_blas(fblas, "syr", "sdcz")
    def test_syr(self, f, dtype):
        x = np.arange(1, 5, dtype='d')
        resx = np.triu(x[:, np.newaxis] * x)
        resx_reverse = np.triu(x[::-1, np.newaxis] * x[::-1])
        y = np.linspace(0, 8.5, 17, endpoint=False)
        z = np.arange(1, 9, dtype='d').view('D')
        resz = np.triu(z[:, np.newaxis] * z)
        resz_reverse = np.triu(z[::-1, np.newaxis] * z[::-1])
        w = np.c_[np.zeros(4), z, np.zeros(4)].ravel()

        rtol = np.finfo(dtype).eps

        assert_allclose(f(1.0, x), resx, rtol=rtol)
        assert_allclose(f(1.0, x, lower=True), resx.T, rtol=rtol)
        assert_allclose(f(1.0, y, incx=2, offx=2, n=4), resx, rtol=rtol)
        # negative increments imply reversed vectors in blas
        assert_allclose(f(1.0, y, incx=-2, offx=2, n=4),
                        resx_reverse, rtol=rtol)

        if dtype in COMPLEX_DTYPES:
            assert_allclose(f(1.0, z), resz, rtol=rtol)
            assert_allclose(f(1.0, z, lower=True), resz.T, rtol=rtol)
            assert_allclose(f(1.0, w, incx=3, offx=1, n=4), resz, rtol=rtol)
            # negative increments imply reversed vectors in blas
            assert_allclose(f(1.0, w, incx=-3, offx=1, n=4),
                            resz_reverse, rtol=rtol)

            a = np.zeros((4, 4), dtype, 'F')
            b = f(1.0, z, a=a, overwrite_a=True)
            assert_allclose(a, resz, rtol=rtol)
            b = f(2.0, z, a=a)
            assert a is not b
            assert_allclose(b, 3*resz, rtol=rtol)

        else:
            a = np.zeros((4, 4), dtype, 'F')
            b = f(1.0, x, a=a, overwrite_a=True)
            assert_allclose(a, resx, rtol=rtol)
            b = f(2.0, x, a=a)
            assert a is not b
            assert_allclose(b, 3*resx, rtol=rtol)

        assert_raises(Exception, f, 1.0, x, incx=0)
        assert_raises(Exception, f, 1.0, x, offx=5)
        assert_raises(Exception, f, 1.0, x, offx=-2)
        assert_raises(Exception, f, 1.0, x, n=-2)
        assert_raises(Exception, f, 1.0, x, n=5)
        assert_raises(Exception, f, 1.0, x, lower=2)
        assert_raises(Exception, f, 1.0, x, a=np.zeros((2, 2), 'd', 'F'))

    @parametrize_blas(fblas, "her", "cz")
    def test_her(self, f, dtype):
        x = np.arange(1, 5, dtype='d')
        z = np.arange(1, 9, dtype='d').view('D')
        rehz = np.triu(z[:, np.newaxis] * z.conj())
        rehz_reverse = np.triu(z[::-1, np.newaxis] * z[::-1].conj())
        w = np.c_[np.zeros(4), z, np.zeros(4)].ravel()

        rtol = np.finfo(dtype).eps

        assert_allclose(f(1.0, z), rehz, rtol=rtol)
        assert_allclose(f(1.0, z, lower=True), rehz.T.conj(), rtol=rtol)
        assert_allclose(f(1.0, w, incx=3, offx=1, n=4), rehz, rtol=rtol)
        # negative increments imply reversed vectors in blas
        assert_allclose(f(1.0, w, incx=-3, offx=1, n=4),
                        rehz_reverse, rtol=rtol)

        a = np.zeros((4, 4), dtype, 'F')
        b = f(1.0, z, a=a, overwrite_a=True)
        assert_allclose(a, rehz, rtol=rtol)

        b = f(2.0, z, a=a)
        assert a is not b
        assert_allclose(b, 3*rehz, rtol=rtol)

        assert_raises(Exception, f, 1.0, x, incx=0)
        assert_raises(Exception, f, 1.0, x, offx=5)
        assert_raises(Exception, f, 1.0, x, offx=-2)
        assert_raises(Exception, f, 1.0, x, n=-2)
        assert_raises(Exception, f, 1.0, x, n=5)
        assert_raises(Exception, f, 1.0, x, lower=2)
        assert_raises(Exception, f, 1.0, x, a=np.zeros((2, 2), 'd', 'F'))

    @parametrize_blas(fblas, "syr2", "sd")
    def test_syr2(self, f, dtype):
        x = np.arange(1, 5, dtype='d')
        y = np.arange(5, 9, dtype='d')
        resxy = np.triu(x[:, np.newaxis] * y + y[:, np.newaxis] * x)
        resxy_reverse = np.triu(x[::-1, np.newaxis] * y[::-1]
                                + y[::-1, np.newaxis] * x[::-1])

        q = np.linspace(0, 8.5, 17, endpoint=False)
        rtol = np.finfo(dtype).eps

        assert_allclose(f(1.0, x, y), resxy, rtol=rtol)
        assert_allclose(f(1.0, x, y, n=3), resxy[:3, :3], rtol=rtol)
        assert_allclose(f(1.0, x, y, lower=True), resxy.T, rtol=rtol)

        assert_allclose(f(1.0, q, q, incx=2, offx=2, incy=2, offy=10),
                        resxy, rtol=rtol)
        assert_allclose(f(1.0, q, q, incx=2, offx=2, incy=2, offy=10, n=3),
                        resxy[:3, :3], rtol=rtol)
        # negative increments imply reversed vectors in blas
        assert_allclose(f(1.0, q, q, incx=-2, offx=2, incy=-2, offy=10),
                        resxy_reverse, rtol=rtol)

        a = np.zeros((4, 4), dtype, 'F')
        b = f(1.0, x, y, a=a, overwrite_a=True)
        assert_allclose(a, resxy, rtol=rtol)

        b = f(2.0, x, y, a=a)
        assert a is not b
        assert_allclose(b, 3*resxy, rtol=rtol)

        assert_raises(Exception, f, 1.0, x, y, incx=0)
        assert_raises(Exception, f, 1.0, x, y, offx=5)
        assert_raises(Exception, f, 1.0, x, y, offx=-2)
        assert_raises(Exception, f, 1.0, x, y, incy=0)
        assert_raises(Exception, f, 1.0, x, y, offy=5)
        assert_raises(Exception, f, 1.0, x, y, offy=-2)
        assert_raises(Exception, f, 1.0, x, y, n=-2)
        assert_raises(Exception, f, 1.0, x, y, n=5)
        assert_raises(Exception, f, 1.0, x, y, lower=2)
        assert_raises(Exception, f, 1.0, x, y, a=np.zeros((2, 2), 'd', 'F'))

    @parametrize_blas(fblas, "her2", "cz")
    def test_her2(self, f, dtype):
        x = np.arange(1, 9, dtype='d').view('D')
        y = np.arange(9, 17, dtype='d').view('D')
        resxy = x[:, np.newaxis] * y.conj() + y[:, np.newaxis] * x.conj()
        resxy = np.triu(resxy)

        resxy_reverse = x[::-1, np.newaxis] * y[::-1].conj()
        resxy_reverse += y[::-1, np.newaxis] * x[::-1].conj()
        resxy_reverse = np.triu(resxy_reverse)

        u = np.c_[np.zeros(4), x, np.zeros(4)].ravel()
        v = np.c_[np.zeros(4), y, np.zeros(4)].ravel()

        rtol = np.finfo(dtype).eps

        assert_allclose(f(1.0, x, y), resxy, rtol=rtol)
        assert_allclose(f(1.0, x, y, n=3), resxy[:3, :3], rtol=rtol)
        assert_allclose(f(1.0, x, y, lower=True), resxy.T.conj(),
                        rtol=rtol)

        assert_allclose(f(1.0, u, v, incx=3, offx=1, incy=3, offy=1),
                        resxy, rtol=rtol)
        assert_allclose(f(1.0, u, v, incx=3, offx=1, incy=3, offy=1, n=3),
                        resxy[:3, :3], rtol=rtol)
        # negative increments imply reversed vectors in blas
        assert_allclose(f(1.0, u, v, incx=-3, offx=1, incy=-3, offy=1),
                        resxy_reverse, rtol=rtol)

        a = np.zeros((4, 4), dtype, 'F')
        b = f(1.0, x, y, a=a, overwrite_a=True)
        assert_allclose(a, resxy, rtol=rtol)

        b = f(2.0, x, y, a=a)
        assert a is not b
        assert_allclose(b, 3*resxy, rtol=rtol)

        assert_raises(Exception, f, 1.0, x, y, incx=0)
        assert_raises(Exception, f, 1.0, x, y, offx=5)
        assert_raises(Exception, f, 1.0, x, y, offx=-2)
        assert_raises(Exception, f, 1.0, x, y, incy=0)
        assert_raises(Exception, f, 1.0, x, y, offy=5)
        assert_raises(Exception, f, 1.0, x, y, offy=-2)
        assert_raises(Exception, f, 1.0, x, y, n=-2)
        assert_raises(Exception, f, 1.0, x, y, n=5)
        assert_raises(Exception, f, 1.0, x, y, lower=2)
        assert_raises(Exception, f, 1.0, x, y,
                      a=np.zeros((2, 2), 'd', 'F'))

    @pytest.mark.parametrize("dtype", DTYPES)
    def test_gbmv(self, dtype):
        rng = np.random.default_rng(1234)
        n = 7
        m = 5
        kl = 1
        ku = 2
        # fake a banded matrix via toeplitz
        A = toeplitz(append(rng.random(kl+1), zeros(m-kl-1)),
                     append(rng.random(ku+1), zeros(n-ku-1)))
        A = A.astype(dtype)
        Ab = zeros((kl+ku+1, n), dtype=dtype)

        # Form the banded storage
        Ab[2, :5] = A[0, 0]  # diag
        Ab[1, 1:6] = A[0, 1]  # sup1
        Ab[0, 2:7] = A[0, 2]  # sup2
        Ab[3, :4] = A[1, 0]  # sub1

        x = rng.random(n).astype(dtype)
        y = rng.random(m).astype(dtype)
        alpha, beta = dtype(3), dtype(-5)

        func, = get_blas_funcs(('gbmv',), dtype=dtype)
        y1 = func(m=m, n=n, ku=ku, kl=kl, alpha=alpha, a=Ab,
                  x=x, y=y, beta=beta)
        y2 = alpha * A.dot(x) + beta * y
        assert_array_almost_equal(y1, y2)

        y1 = func(m=m, n=n, ku=ku, kl=kl, alpha=alpha, a=Ab,
                  x=y, y=x, beta=beta, trans=1)
        y2 = alpha * A.T.dot(y) + beta * x
        assert_array_almost_equal(y1, y2)

    @pytest.mark.parametrize("dtype", DTYPES)
    def test_sbmv_hbmv(self, dtype):
        rng = np.random.default_rng(1234)
        n = 6
        k = 2
        A = zeros((n, n), dtype=dtype)
        Ab = zeros((k+1, n), dtype=dtype)

        # Form the array and its packed banded storage
        A[arange(n), arange(n)] = rng.random(n)
        for ind2 in range(1, k+1):
            temp = rng.random(n-ind2)
            A[arange(n-ind2), arange(ind2, n)] = temp
            Ab[-1-ind2, ind2:] = temp
        A = A.astype(dtype)
        if dtype in COMPLEX_DTYPES:
            A += A.conj().T
            func, = get_blas_funcs(('hbmv',), dtype=dtype)
        else:
            A += A.T
            func, = get_blas_funcs(('sbmv',), dtype=dtype)

        Ab[-1, :] = diag(A)
        x = rng.random(n).astype(dtype)
        y = rng.random(n).astype(dtype)
        alpha, beta = dtype(1.25), dtype(3)

        y1 = func(k=k, alpha=alpha, a=Ab, x=x, y=y, beta=beta)
        y2 = alpha * A.dot(x) + beta * y
        assert_array_almost_equal(y1, y2)

    @pytest.mark.parametrize("fname,dtype", [
        *[('spmv', dtype) for dtype in REAL_DTYPES + COMPLEX_DTYPES],
        *[('hpmv', dtype) for dtype in COMPLEX_DTYPES],
    ])
    def test_spmv_hpmv(self, fname, dtype):
        rng = np.random.default_rng(1234)
        n = 3
        A = rng.random((n, n)).astype(dtype)
        if dtype in COMPLEX_DTYPES:
            A += rng.random((n, n))*1j
        A += A.T if fname == 'spmv' else A.conj().T
        c, r = tril_indices(n)
        Ap = A[r, c]
        x = rng.random(n).astype(dtype)
        y = rng.random(n).astype(dtype)
        xlong = arange(2*n).astype(dtype)
        ylong = ones(2*n).astype(dtype)
        alpha, beta = dtype(1.25), dtype(2)

        func, = get_blas_funcs((fname,), dtype=dtype)
        y1 = func(n=n, alpha=alpha, ap=Ap, x=x, y=y, beta=beta)
        y2 = alpha * A.dot(x) + beta * y
        assert_array_almost_equal(y1, y2)

        # Test inc and offsets
        y1 = func(n=n-1, alpha=alpha, beta=beta, x=xlong, y=ylong, ap=Ap,
                  incx=2, incy=2, offx=n, offy=n)
        y2 = (alpha * A[:-1, :-1]).dot(xlong[3::2]) + beta * ylong[3::2]
        assert_array_almost_equal(y1[3::2], y2)
        assert_almost_equal(y1[4], ylong[4])

    @pytest.mark.parametrize("fname,dtype", [
        *[('spr', dtype) for dtype in REAL_DTYPES + COMPLEX_DTYPES],
        *[('hpr', dtype) for dtype in COMPLEX_DTYPES],
    ])
    def test_spr_hpr(self, fname, dtype):
        rng = np.random.default_rng(1234)
        n = 3
        A = rng.random((n, n)).astype(dtype)
        if dtype in COMPLEX_DTYPES:
            A += rng.random((n, n))*1j
        A += A.T if fname == 'spr' else A.conj().T
        c, r = tril_indices(n)
        Ap = A[r, c]
        x = rng.random(n).astype(dtype)

        alpha = np.finfo(dtype).dtype.type(2.5)
        if fname == 'hpr':
            func, = get_blas_funcs(('hpr',), dtype=dtype)
            y2 = alpha * x[:, None].dot(x[None, :].conj()) + A
        else:
            func, = get_blas_funcs(('spr',), dtype=dtype)
            y2 = alpha * x[:, None].dot(x[None, :]) + A

        y1 = func(n=n, alpha=alpha, ap=Ap, x=x)
        y1f = zeros((3, 3), dtype=dtype)
        y1f[r, c] = y1
        y1f[c, r] = y1.conj() if fname == 'hpr' else y1
        assert_array_almost_equal(y1f, y2)

    @pytest.mark.parametrize("dtype", DTYPES)
    def test_spr2_hpr2(self, dtype):
        rng = np.random.default_rng(1234)
        n = 3
        A = rng.random((n, n)).astype(dtype)
        if dtype in COMPLEX_DTYPES:
            A += rng.random((n, n))*1j
            A += A.conj().T
            func, = get_blas_funcs(('hpr2',), dtype=dtype)
        else:
            A += A.T
            func, = get_blas_funcs(('spr2',), dtype=dtype)

        c, r = tril_indices(n)
        Ap = A[r, c]
        x = rng.random(n).astype(dtype)
        y = rng.random(n).astype(dtype)
        alpha = dtype(2)

        u = alpha.conj() * x[:, None].dot(y[None, :].conj())
        y2 = A + u + u.conj().T
        y1 = func(n=n, alpha=alpha, x=x, y=y, ap=Ap)
        y1f = zeros((3, 3), dtype=dtype)
        y1f[r, c] = y1
        y1f[[1, 2, 2], [0, 0, 1]] = y1[[1, 3, 4]].conj()
        assert_array_almost_equal(y1f, y2)

    @pytest.mark.parametrize("dtype", DTYPES)
    def test_tbmv(self, dtype):
        rng = np.random.default_rng(1234)
        n = 10
        k = 3
        x = rng.random(n).astype(dtype)
        A = zeros((n, n), dtype=dtype)
        # Banded upper triangular array
        for sup in range(k+1):
            A[arange(n-sup), arange(sup, n)] = rng.random(n-sup)

        # Add complex parts for c,z
        if dtype in COMPLEX_DTYPES:
            A[nonzero(A)] += 1j * rng.random((k+1)*n-(k*(k+1)//2)).astype(dtype)

        # Form the banded storage
        Ab = zeros((k+1, n), dtype=dtype)
        for row in range(k+1):
            Ab[-row-1, row:] = diag(A, k=row)
        func, = get_blas_funcs(('tbmv',), dtype=dtype)

        y1 = func(k=k, a=Ab, x=x)
        y2 = A.dot(x)
        assert_array_almost_equal(y1, y2)

        y1 = func(k=k, a=Ab, x=x, diag=1)
        A[arange(n), arange(n)] = dtype(1)
        y2 = A.dot(x)
        assert_array_almost_equal(y1, y2)

        y1 = func(k=k, a=Ab, x=x, diag=1, trans=1)
        y2 = A.T.dot(x)
        assert_array_almost_equal(y1, y2)

        y1 = func(k=k, a=Ab, x=x, diag=1, trans=2)
        y2 = A.conj().T.dot(x)
        assert_array_almost_equal(y1, y2)

    @pytest.mark.parametrize("dtype", DTYPES)
    def test_tbsv(self, dtype):
        rng = np.random.default_rng(12345)
        n = 6
        k = 3
        x = rng.random(n).astype(dtype)
        A = zeros((n, n), dtype=dtype)
        # Banded upper triangular array
        for sup in range(k+1):
            A[arange(n-sup), arange(sup, n)] = rng.random(n-sup)

        # Add complex parts for c,z
        if dtype in COMPLEX_DTYPES:
            A[nonzero(A)] += 1j * rng.random((k+1)*n-(k*(k+1)//2)).astype(dtype)

        # Form the banded storage
        Ab = zeros((k+1, n), dtype=dtype)
        for row in range(k+1):
            Ab[-row-1, row:] = diag(A, k=row)
        func, = get_blas_funcs(('tbsv',), dtype=dtype)

        y1 = func(k=k, a=Ab, x=x)
        y2 = solve(A, x)
        assert_array_almost_equal(y1, y2)

        y1 = func(k=k, a=Ab, x=x, diag=1)
        A[arange(n), arange(n)] = dtype(1)
        y2 = solve(A, x)
        assert_array_almost_equal(y1, y2)

        y1 = func(k=k, a=Ab, x=x, diag=1, trans=1)
        y2 = solve(A.T, x)
        assert_array_almost_equal(y1, y2)

        y1 = func(k=k, a=Ab, x=x, diag=1, trans=2)
        y2 = solve(A.conj().T, x)
        assert_array_almost_equal(y1, y2)

    @pytest.mark.parametrize("dtype", DTYPES)
    def test_tpmv(self, dtype):
        rng = np.random.default_rng(1234)
        n = 10
        x = rng.random(n).astype(dtype)
        # Upper triangular array
        if dtype in COMPLEX_DTYPES:
            A = triu(rng.random((n, n)) + rng.random((n, n))*1j)
        else:
            A = triu(rng.random((n, n)))

        # Form the packed storage
        c, r = tril_indices(n)
        Ap = A[r, c]
        func, = get_blas_funcs(('tpmv',), dtype=dtype)

        y1 = func(n=n, ap=Ap, x=x)
        y2 = A.dot(x)
        assert_array_almost_equal(y1, y2)

        y1 = func(n=n, ap=Ap, x=x, diag=1)
        A[arange(n), arange(n)] = dtype(1)
        y2 = A.dot(x)
        assert_array_almost_equal(y1, y2)

        y1 = func(n=n, ap=Ap, x=x, diag=1, trans=1)
        y2 = A.T.dot(x)
        assert_array_almost_equal(y1, y2)

        y1 = func(n=n, ap=Ap, x=x, diag=1, trans=2)
        y2 = A.conj().T.dot(x)
        assert_array_almost_equal(y1, y2)

    @pytest.mark.parametrize("dtype", DTYPES)
    def test_tpsv(self, dtype):
        rng = np.random.default_rng(1234)
        n = 10
        x = rng.random(n).astype(dtype)
        # Upper triangular array
        if dtype in COMPLEX_DTYPES:
            A = triu(rng.random((n, n)) + rng.random((n, n))*1j)
        else:
            A = triu(rng.random((n, n)))
        A += eye(n)
        # Form the packed storage
        c, r = tril_indices(n)
        Ap = A[r, c]
        func, = get_blas_funcs(('tpsv',), dtype=dtype)

        y1 = func(n=n, ap=Ap, x=x)
        y2 = solve(A, x)
        assert_array_almost_equal(y1, y2)

        y1 = func(n=n, ap=Ap, x=x, diag=1)
        A[arange(n), arange(n)] = dtype(1)
        y2 = solve(A, x)
        assert_array_almost_equal(y1, y2)

        y1 = func(n=n, ap=Ap, x=x, diag=1, trans=1)
        y2 = solve(A.T, x)
        assert_array_almost_equal(y1, y2)

        y1 = func(n=n, ap=Ap, x=x, diag=1, trans=2)
        y2 = solve(A.conj().T, x)
        assert_array_almost_equal(y1, y2)

    @pytest.mark.parametrize("dtype", DTYPES)
    def test_trmv(self, dtype):
        rng = np.random.default_rng(1234)
        n = 3
        A = (rng.random((n, n))+eye(n)).astype(dtype)
        x = rng.random(3).astype(dtype)
        func, = get_blas_funcs(('trmv',), dtype=dtype)

        y1 = func(a=A, x=x)
        y2 = triu(A).dot(x)
        assert_array_almost_equal(y1, y2)

        y1 = func(a=A, x=x, diag=1)
        A[arange(n), arange(n)] = dtype(1)
        y2 = triu(A).dot(x)
        assert_array_almost_equal(y1, y2)

        y1 = func(a=A, x=x, diag=1, trans=1)
        y2 = triu(A).T.dot(x)
        assert_array_almost_equal(y1, y2)

        y1 = func(a=A, x=x, diag=1, trans=2)
        y2 = triu(A).conj().T.dot(x)
        assert_array_almost_equal(y1, y2)

    @pytest.mark.parametrize("dtype", DTYPES)
    def test_trsv(self, dtype):
        rng = np.random.default_rng(1234)
        n = 15
        A = (rng.random((n, n))+eye(n)).astype(dtype)
        x = rng.random(n).astype(dtype)
        func, = get_blas_funcs(('trsv',), dtype=dtype)

        y1 = func(a=A, x=x)
        y2 = solve(triu(A), x)
        assert_array_almost_equal(y1, y2)

        y1 = func(a=A, x=x, lower=1)
        y2 = solve(tril(A), x)
        assert_array_almost_equal(y1, y2)

        y1 = func(a=A, x=x, diag=1)
        A[arange(n), arange(n)] = dtype(1)
        y2 = solve(triu(A), x)
        assert_array_almost_equal(y1, y2)

        y1 = func(a=A, x=x, diag=1, trans=1)
        y2 = solve(triu(A).T, x)
        assert_array_almost_equal(y1, y2)

        y1 = func(a=A, x=x, diag=1, trans=2)
        y2 = solve(triu(A).conj().T, x)
        assert_array_almost_equal(y1, y2)


class TestFBLAS3Simple:
    @parametrize_blas(fblas, "gemm", "sdcz")
    def test_gemm(self, f, dtype):
        assert_array_almost_equal(f(3, [3], [-4]), [[-36]])
        assert_array_almost_equal(f(3, [3], [-4], 3, [5]), [-21])
        if dtype in COMPLEX_DTYPES:
            assert_array_almost_equal(f(3j, [3-4j], [-4]), [[-48-36j]])
            assert_array_almost_equal(f(3j, [3-4j], [-4], 3, [5j]), [-48-21j])


class TestBLAS3Symm:

    def setup_method(self):
        self.a = np.array([[1., 2.],
                           [0., 1.]])
        self.b = np.array([[1., 0., 3.],
                           [0., -1., 2.]])
        self.c = np.ones((2, 3))
        self.t = np.array([[2., -1., 8.],
                           [3., 0., 9.]])

    @parametrize_blas(fblas, "symm", "sdcz")
    def test_symm(self, f, dtype):
        res = f(a=self.a, b=self.b, c=self.c, alpha=1., beta=1.)
        assert_array_almost_equal(res, self.t)

        res = f(a=self.a.T, b=self.b, lower=1, c=self.c, alpha=1., beta=1.)
        assert_array_almost_equal(res, self.t)

        res = f(a=self.a, b=self.b.T, side=1, c=self.c.T,
                alpha=1., beta=1.)
        assert_array_almost_equal(res, self.t.T)

    @parametrize_blas(fblas, "symm", "sdcz")
    def test_symm_wrong_side(self, f, dtype):
        """`side=1` means C <- B*A, hence shapes of A and B are to be
        compatible. Otherwise, f2py exception is raised.
        """
        # FIXME narrow down to _fblas.error
        with pytest.raises(Exception):
            f(a=self.a, b=self.b, alpha=1, side=1)

    @parametrize_blas(fblas, "symm", "sdcz")
    def test_symm_wrong_uplo(self, f, dtype):
        """SYMM only considers the upper/lower part of A. Hence setting
        wrong value for `lower` (default is lower=0, meaning upper triangle)
        gives a wrong result.
        """
        res = f(a=self.a, b=self.b, c=self.c, alpha=1., beta=1.)
        assert np.allclose(res, self.t)
        res = f(a=self.a, b=self.b, lower=1, c=self.c, alpha=1., beta=1.)
        assert not np.allclose(res, self.t)


class TestBLAS3Syrk:
    def setup_method(self):
        self.a = np.array([[1., 0.],
                           [0., -2.],
                           [2., 3.]])
        self.t = np.array([[1., 0., 2.],
                           [0., 4., -6.],
                           [2., -6., 13.]])
        self.tt = np.array([[5., 6.],
                            [6., 13.]])

    @parametrize_blas(fblas, "syrk", "sdcz")
    def test_syrk(self, f, dtype):
        c = f(a=self.a, alpha=1.)
        assert_array_almost_equal(np.triu(c), np.triu(self.t))

        c = f(a=self.a, alpha=1., lower=1)
        assert_array_almost_equal(np.tril(c), np.tril(self.t))

        c0 = np.ones(self.t.shape)
        c = f(a=self.a, alpha=1., beta=1., c=c0)
        assert_array_almost_equal(np.triu(c), np.triu(self.t+c0))

        c = f(a=self.a, alpha=1., trans=1)
        assert_array_almost_equal(np.triu(c), np.triu(self.tt))

    # prints '0-th dimension must be fixed to 3 but got 5',
    # FIXME: suppress?
    @parametrize_blas(fblas, "syrk", "sdcz")
    def test_syrk_wrong_c(self, f, dtype):
        # FIXME narrow down to _fblas.error
        with pytest.raises(Exception):
            f(a=self.a, alpha=1., c=np.ones((5, 8)))
        # if C is supplied, it must have compatible dimensions


class TestBLAS3Syr2k:
    def setup_method(self):
        self.a = np.array([[1., 0.],
                           [0., -2.],
                           [2., 3.]])
        self.b = np.array([[0., 1.],
                           [1., 0.],
                           [0, 1.]])
        self.t = np.array([[0., -1., 3.],
                           [-1., 0., 0.],
                           [3., 0., 6.]])
        self.tt = np.array([[0., 1.],
                            [1., 6]])

    @parametrize_blas(fblas, "syr2k", "sdcz")
    def test_syr2k(self, f, dtype):
        c = f(a=self.a, b=self.b, alpha=1.)
        assert_array_almost_equal(np.triu(c), np.triu(self.t))

        c = f(a=self.a, b=self.b, alpha=1., lower=1)
        assert_array_almost_equal(np.tril(c), np.tril(self.t))

        c0 = np.ones(self.t.shape)
        c = f(a=self.a, b=self.b, alpha=1., beta=1., c=c0)
        assert_array_almost_equal(np.triu(c), np.triu(self.t+c0))

        c = f(a=self.a, b=self.b, alpha=1., trans=1)
        assert_array_almost_equal(np.triu(c), np.triu(self.tt))

    # prints '0-th dimension must be fixed to 3 but got 5', FIXME: suppress?
    @parametrize_blas(fblas, "syr2k", "sdcz")
    def test_syr2k_wrong_c(self, f, dtype):
        with pytest.raises(Exception):
            f(a=self.a, b=self.b, alpha=1., c=np.zeros((15, 8)))
        # if C is supplied, it must have compatible dimensions


class TestSyHe:
    """Quick and simple tests for (zc)-symm, syrk, syr2k."""

    def setup_method(self):
        self.sigma_y = np.array([[0., -1.j],
                                 [1.j, 0.]])

    @parametrize_blas(fblas, "symm", "zc")
    def test_symm(self, f, dtype):
        # NB: a is symmetric w/upper diag of ONLY
        res = f(a=self.sigma_y, b=self.sigma_y, alpha=1.)
        assert_array_almost_equal(np.triu(res), np.diag([1, -1]))

    @parametrize_blas(fblas, "hemm", "zc")
    def test_hemm(self, f, dtype):
        # NB: a is hermitian w/upper diag of ONLY
        res = f(a=self.sigma_y, b=self.sigma_y, alpha=1.)
        assert_array_almost_equal(np.triu(res), np.diag([1, 1]))

    @parametrize_blas(fblas, "syrk", "zc")
    def test_syrk(self, f, dtype):
        res = f(a=self.sigma_y, alpha=1.)
        assert_array_almost_equal(np.triu(res), np.diag([-1, -1]))

    @parametrize_blas(fblas, "herk", "zc")
    def test_herk(self, f, dtype):
        res = f(a=self.sigma_y, alpha=1.)
        assert_array_almost_equal(np.triu(res), np.diag([1, 1]))

    @parametrize_blas(fblas, "syr2k", "zc")
    def test_syr2k_zr(self, f, dtype):
        res = f(a=self.sigma_y, b=self.sigma_y, alpha=1.)
        assert_array_almost_equal(np.triu(res), 2.*np.diag([-1, -1]))

    @parametrize_blas(fblas, "her2k", "zc")
    def test_her2k_zr(self, f, dtype):
        res = f(a=self.sigma_y, b=self.sigma_y, alpha=1.)
        assert_array_almost_equal(np.triu(res), 2.*np.diag([1, 1]))


class TestTRMM:
    """Quick and simple tests for *trmm."""

    def setup_method(self):
        self.a = np.array([[1., 2., ],
                           [-2., 1.]])
        self.b = np.array([[3., 4., -1.],
                           [5., 6., -2.]])

        self.a2 = np.array([[1, 1, 2, 3],
                            [0, 1, 4, 5],
                            [0, 0, 1, 6],
                            [0, 0, 0, 1]], order="f")
        self.b2 = np.array([[1, 4], [2, 5], [3, 6], [7, 8], [9, 10]],
                           order="f")

    @pytest.mark.parametrize("dtype", DTYPES)
    def test_side(self, dtype):
        trmm = get_blas_funcs("trmm", dtype=dtype)
        # Provide large A array that works for side=1 but not 0 (see gh-10841)
        assert_raises(Exception, trmm, 1.0, self.a2, self.b2)
        res = trmm(1.0, self.a2.astype(dtype), self.b2.astype(dtype),
                   side=1)
        k = self.b2.shape[1]
        assert_allclose(res, self.b2 @ self.a2[:k, :k], rtol=0.,
                        atol=100*np.finfo(dtype).eps)

    @parametrize_blas(fblas, "trmm", "sdcz")
    def test_ab(self, f, dtype):
        result = f(1., self.a, self.b)
        # default a is upper triangular
        expected = np.array([[13., 16., -5.],
                             [ 5.,  6., -2.]])
        assert_array_almost_equal(result, expected)

    @parametrize_blas(fblas, "trmm", "sdcz")
    def test_ab_lower(self, f, dtype):
        result = f(1., self.a, self.b, lower=True)
        expected = np.array([[ 3.,  4., -1.],
                             [-1., -2.,  0.]])  # now a is lower triangular
        assert_array_almost_equal(result, expected)

    @parametrize_blas(fblas, "trmm", "sdcz")
    def test_b_overwrites(self, f, dtype):
        # BLAS *trmm modifies B argument in-place.
        # Here the default is to copy, but this can be overridden
        b = self.b.astype(dtype)
        for overwr in [True, False]:
            bcopy = b.copy()
            result = f(1., self.a, bcopy, overwrite_b=overwr)
            # C-contiguous arrays are copied
            assert not bcopy.flags.f_contiguous
            assert not np.may_share_memory(bcopy, result)
            assert_equal(bcopy, b)

        bcopy = np.asfortranarray(b.copy())  # or just transpose it
        result = f(1., self.a, bcopy, overwrite_b=True)
        assert bcopy.flags.f_contiguous
        assert np.may_share_memory(bcopy, result)
        assert_array_almost_equal(bcopy, result)


@pytest.mark.parametrize("dtype", DTYPES)
def test_trsm(dtype):
    rng = np.random.default_rng(1234)
    tol = np.finfo(dtype).eps*1000
    func, = get_blas_funcs(('trsm',), dtype=dtype)

    # Test protection against size mismatches
    A = rng.random((4, 5)).astype(dtype)
    B = rng.random((4, 4)).astype(dtype)
    alpha = dtype(1)
    assert_raises(Exception, func, alpha, A, B)
    assert_raises(Exception, func, alpha, A.T, B)

    n = 8
    m = 7
    alpha = dtype(-2.5)
    if dtype in COMPLEX_DTYPES:
        A = (rng.random((m, m)) + rng.random((m, m))*1j) + eye(m)
    else:
        A = rng.random((m, m)) + eye(m)
    A = A.astype(dtype)
    Au = triu(A)
    Al = tril(A)
    B1 = rng.random((m, n)).astype(dtype)
    B2 = rng.random((n, m)).astype(dtype)

    x1 = func(alpha=alpha, a=A, b=B1)
    assert_equal(B1.shape, x1.shape)
    x2 = solve(Au, alpha*B1)
    assert_allclose(x1, x2, atol=tol)

    x1 = func(alpha=alpha, a=A, b=B1, trans_a=1)
    x2 = solve(Au.T, alpha*B1)
    assert_allclose(x1, x2, atol=tol)

    x1 = func(alpha=alpha, a=A, b=B1, trans_a=2)
    x2 = solve(Au.conj().T, alpha*B1)
    assert_allclose(x1, x2, atol=tol)

    x1 = func(alpha=alpha, a=A, b=B1, diag=1)
    Au[arange(m), arange(m)] = dtype(1)
    x2 = solve(Au, alpha*B1)
    assert_allclose(x1, x2, atol=tol)

    x1 = func(alpha=alpha, a=A, b=B2, diag=1, side=1)
    x2 = solve(Au.conj().T, alpha*B2.conj().T)
    assert_allclose(x1, x2.conj().T, atol=tol)

    x1 = func(alpha=alpha, a=A, b=B2, diag=1, side=1, lower=1)
    Al[arange(m), arange(m)] = dtype(1)
    x2 = solve(Al.conj().T, alpha*B2.conj().T)
    assert_allclose(x1, x2.conj().T, atol=tol)


@pytest.mark.xfail(run=False,
                   reason="gh-16930")
def test_gh_169309():
    x = np.repeat(10, 9)
    actual = scipy.linalg.blas.dnrm2(x, 5, 3, -1)
    expected = math.sqrt(500)
    assert_allclose(actual, expected)


def test_dnrm2_neg_incx():
    # check that dnrm2(..., incx < 0) raises
    # XXX: remove the test after the lowest supported BLAS implements
    # negative incx (new in LAPACK 3.10)
    x = np.repeat(10, 9)
    incx = -1
    with assert_raises(fblas.__fblas_error):
        scipy.linalg.blas.dnrm2(x, 5, 3, incx)
