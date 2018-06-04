#
# Created by: Pearu Peterson, September 2002
#

from __future__ import division, print_function, absolute_import

import sys
import subprocess
import time

from numpy.testing import (assert_equal, assert_array_almost_equal, assert_,
                           assert_allclose, assert_almost_equal,
                           assert_array_equal)
import pytest
from pytest import raises as assert_raises

import numpy as np
from numpy.random import rand, seed

from scipy.linalg import _flapack as flapack
from scipy.linalg import inv
from scipy.linalg import svd
from scipy.linalg.lapack import _compute_lwork

try:
    from scipy.linalg import _clapack as clapack
except ImportError:
    clapack = None
from scipy.linalg.lapack import get_lapack_funcs
from scipy.linalg.blas import get_blas_funcs

REAL_DTYPES = [np.float32, np.float64]
COMPLEX_DTYPES = [np.complex64, np.complex128]
DTYPES = REAL_DTYPES + COMPLEX_DTYPES


class TestFlapackSimple(object):

    def test_gebal(self):
        a = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
        a1 = [[1, 0, 0, 3e-4],
              [4, 0, 0, 2e-3],
              [7, 1, 0, 0],
              [0, 1, 0, 0]]
        for p in 'sdzc':
            f = getattr(flapack, p+'gebal', None)
            if f is None:
                continue
            ba, lo, hi, pivscale, info = f(a)
            assert_(not info, repr(info))
            assert_array_almost_equal(ba, a)
            assert_equal((lo, hi), (0, len(a[0])-1))
            assert_array_almost_equal(pivscale, np.ones(len(a)))

            ba, lo, hi, pivscale, info = f(a1, permute=1, scale=1)
            assert_(not info, repr(info))
            # print(a1)
            # print(ba, lo, hi, pivscale)

    def test_gehrd(self):
        a = [[-149, -50, -154],
             [537, 180, 546],
             [-27, -9, -25]]
        for p in 'd':
            f = getattr(flapack, p+'gehrd', None)
            if f is None:
                continue
            ht, tau, info = f(a)
            assert_(not info, repr(info))

    def test_trsyl(self):
        a = np.array([[1, 2], [0, 4]])
        b = np.array([[5, 6], [0, 8]])
        c = np.array([[9, 10], [11, 12]])
        trans = 'T'

        # Test single and double implementations, including most
        # of the options
        for dtype in 'fdFD':
            a1, b1, c1 = a.astype(dtype), b.astype(dtype), c.astype(dtype)
            trsyl, = get_lapack_funcs(('trsyl',), (a1,))
            if dtype.isupper():  # is complex dtype
                a1[0] += 1j
                trans = 'C'

            x, scale, info = trsyl(a1, b1, c1)
            assert_array_almost_equal(np.dot(a1, x) + np.dot(x, b1),
                                      scale * c1)

            x, scale, info = trsyl(a1, b1, c1, trana=trans, tranb=trans)
            assert_array_almost_equal(
                    np.dot(a1.conjugate().T, x) + np.dot(x, b1.conjugate().T),
                    scale * c1, decimal=4)

            x, scale, info = trsyl(a1, b1, c1, isgn=-1)
            assert_array_almost_equal(np.dot(a1, x) - np.dot(x, b1),
                                      scale * c1, decimal=4)

    def test_lange(self):
        a = np.array([
            [-149, -50, -154],
            [537, 180, 546],
            [-27, -9, -25]])

        for dtype in 'fdFD':
            for norm in 'Mm1OoIiFfEe':
                a1 = a.astype(dtype)
                if dtype.isupper():
                    # is complex dtype
                    a1[0, 0] += 1j

                lange, = get_lapack_funcs(('lange',), (a1,))
                value = lange(norm, a1)

                if norm in 'FfEe':
                    if dtype in 'Ff':
                        decimal = 3
                    else:
                        decimal = 7
                    ref = np.sqrt(np.sum(np.square(np.abs(a1))))
                    assert_almost_equal(value, ref, decimal)
                else:
                    if norm in 'Mm':
                        ref = np.max(np.abs(a1))
                    elif norm in '1Oo':
                        ref = np.max(np.sum(np.abs(a1), axis=0))
                    elif norm in 'Ii':
                        ref = np.max(np.sum(np.abs(a1), axis=1))

                    assert_equal(value, ref)


class TestLapack(object):

    def test_flapack(self):
        if hasattr(flapack, 'empty_module'):
            # flapack module is empty
            pass

    def test_clapack(self):
        if hasattr(clapack, 'empty_module'):
            # clapack module is empty
            pass


class TestLeastSquaresSolvers(object):

    def test_gels(self):
        seed(1234)
        # Test fat/tall matrix argument handling - gh-issue #8329
        for ind, dtype in enumerate(DTYPES):
            m = 10
            n = 20
            nrhs = 1
            a1 = rand(m, n).astype(dtype)
            b1 = rand(n).astype(dtype)
            gls, glslw = get_lapack_funcs(('gels', 'gels_lwork'), dtype=dtype)

            # Request of sizes
            lwork = _compute_lwork(glslw, m, n, nrhs)
            _, _, info = gls(a1, b1, lwork=lwork)
            assert_(info >= 0)
            _, _, info = gls(a1, b1, trans='TTCC'[ind], lwork=lwork)
            assert_(info >= 0)

        for dtype in REAL_DTYPES:
            a1 = np.array([[1.0, 2.0],
                          [4.0, 5.0],
                          [7.0, 8.0]], dtype=dtype)
            b1 = np.array([16.0, 17.0, 20.0], dtype=dtype)
            gels, gels_lwork, geqrf = get_lapack_funcs(
                    ('gels', 'gels_lwork', 'geqrf'), (a1, b1))

            m, n = a1.shape
            if len(b1.shape) == 2:
                nrhs = b1.shape[1]
            else:
                nrhs = 1

            # Request of sizes
            lwork = _compute_lwork(gels_lwork, m, n, nrhs)

            lqr, x, info = gels(a1, b1, lwork=lwork)
            assert_allclose(x[:-1], np.array([-14.333333333333323,
                                              14.999999999999991],
                                             dtype=dtype),
                            rtol=25*np.finfo(dtype).eps)
            lqr_truth, _, _, _ = geqrf(a1)
            assert_array_equal(lqr, lqr_truth)

        for dtype in COMPLEX_DTYPES:
            a1 = np.array([[1.0+4.0j, 2.0],
                          [4.0+0.5j, 5.0-3.0j],
                          [7.0-2.0j, 8.0+0.7j]], dtype=dtype)
            b1 = np.array([16.0, 17.0+2.0j, 20.0-4.0j], dtype=dtype)
            gels, gels_lwork, geqrf = get_lapack_funcs(
                    ('gels', 'gels_lwork', 'geqrf'), (a1, b1))

            m, n = a1.shape
            if len(b1.shape) == 2:
                nrhs = b1.shape[1]
            else:
                nrhs = 1

            # Request of sizes
            lwork = _compute_lwork(gels_lwork, m, n, nrhs)

            lqr, x, info = gels(a1, b1, lwork=lwork)
            assert_allclose(x[:-1],
                            np.array([1.161753632288328-1.901075709391912j,
                                      1.735882340522193+1.521240901196909j],
                            dtype=dtype), rtol=25*np.finfo(dtype).eps)
            lqr_truth, _, _, _ = geqrf(a1)
            assert_array_equal(lqr, lqr_truth)

    def test_gelsd(self):
        for dtype in REAL_DTYPES:
            a1 = np.array([[1.0, 2.0],
                          [4.0, 5.0],
                          [7.0, 8.0]], dtype=dtype)
            b1 = np.array([16.0, 17.0, 20.0], dtype=dtype)
            gelsd, gelsd_lwork = get_lapack_funcs(('gelsd', 'gelsd_lwork'),
                                                  (a1, b1))

            m, n = a1.shape
            if len(b1.shape) == 2:
                nrhs = b1.shape[1]
            else:
                nrhs = 1

            # Request of sizes
            work, iwork, info = gelsd_lwork(m, n, nrhs, -1)
            lwork = int(np.real(work))
            iwork_size = iwork

            x, s, rank, info = gelsd(a1, b1, lwork, iwork_size,
                                     -1, False, False)
            assert_allclose(x[:-1], np.array([-14.333333333333323,
                                             14.999999999999991], dtype=dtype),
                            rtol=25*np.finfo(dtype).eps)
            assert_allclose(s, np.array([12.596017180511966,
                                         0.583396253199685], dtype=dtype),
                            rtol=25*np.finfo(dtype).eps)

        for dtype in COMPLEX_DTYPES:
            a1 = np.array([[1.0+4.0j, 2.0],
                          [4.0+0.5j, 5.0-3.0j],
                          [7.0-2.0j, 8.0+0.7j]], dtype=dtype)
            b1 = np.array([16.0, 17.0+2.0j, 20.0-4.0j], dtype=dtype)
            gelsd, gelsd_lwork = get_lapack_funcs(('gelsd', 'gelsd_lwork'),
                                                  (a1, b1))

            m, n = a1.shape
            if len(b1.shape) == 2:
                nrhs = b1.shape[1]
            else:
                nrhs = 1

            # Request of sizes
            work, rwork, iwork, info = gelsd_lwork(m, n, nrhs, -1)
            lwork = int(np.real(work))
            rwork_size = int(rwork)
            iwork_size = iwork

            x, s, rank, info = gelsd(a1, b1, lwork, rwork_size, iwork_size,
                                     -1, False, False)
            assert_allclose(x[:-1],
                            np.array([1.161753632288328-1.901075709391912j,
                                      1.735882340522193+1.521240901196909j],
                            dtype=dtype), rtol=25*np.finfo(dtype).eps)
            assert_allclose(s,
                            np.array([13.035514762572043, 4.337666985231382],
                                     dtype=dtype), rtol=25*np.finfo(dtype).eps)

    def test_gelss(self):

        for dtype in REAL_DTYPES:
            a1 = np.array([[1.0, 2.0],
                          [4.0, 5.0],
                          [7.0, 8.0]], dtype=dtype)
            b1 = np.array([16.0, 17.0, 20.0], dtype=dtype)
            gelss, gelss_lwork = get_lapack_funcs(('gelss', 'gelss_lwork'),
                                                  (a1, b1))

            m, n = a1.shape
            if len(b1.shape) == 2:
                nrhs = b1.shape[1]
            else:
                nrhs = 1

            # Request of sizes
            work, info = gelss_lwork(m, n, nrhs, -1)
            lwork = int(np.real(work))

            v, x, s, rank, work, info = gelss(a1, b1, -1, lwork, False, False)
            assert_allclose(x[:-1], np.array([-14.333333333333323,
                                             14.999999999999991], dtype=dtype),
                            rtol=25*np.finfo(dtype).eps)
            assert_allclose(s, np.array([12.596017180511966,
                                         0.583396253199685], dtype=dtype),
                            rtol=25*np.finfo(dtype).eps)

        for dtype in COMPLEX_DTYPES:
            a1 = np.array([[1.0+4.0j, 2.0],
                          [4.0+0.5j, 5.0-3.0j],
                          [7.0-2.0j, 8.0+0.7j]], dtype=dtype)
            b1 = np.array([16.0, 17.0+2.0j, 20.0-4.0j], dtype=dtype)
            gelss, gelss_lwork = get_lapack_funcs(('gelss', 'gelss_lwork'),
                                                  (a1, b1))

            m, n = a1.shape
            if len(b1.shape) == 2:
                nrhs = b1.shape[1]
            else:
                nrhs = 1

            # Request of sizes
            work, info = gelss_lwork(m, n, nrhs, -1)
            lwork = int(np.real(work))

            v, x, s, rank, work, info = gelss(a1, b1, -1, lwork, False, False)
            assert_allclose(x[:-1],
                            np.array([1.161753632288328-1.901075709391912j,
                                      1.735882340522193+1.521240901196909j],
                                     dtype=dtype),
                            rtol=25*np.finfo(dtype).eps)
            assert_allclose(s, np.array([13.035514762572043,
                                         4.337666985231382], dtype=dtype),
                            rtol=25*np.finfo(dtype).eps)

    def test_gelsy(self):

        for dtype in REAL_DTYPES:
            a1 = np.array([[1.0, 2.0],
                          [4.0, 5.0],
                          [7.0, 8.0]], dtype=dtype)
            b1 = np.array([16.0, 17.0, 20.0], dtype=dtype)
            gelsy, gelsy_lwork = get_lapack_funcs(('gelsy', 'gelss_lwork'),
                                                  (a1, b1))

            m, n = a1.shape
            if len(b1.shape) == 2:
                nrhs = b1.shape[1]
            else:
                nrhs = 1

            # Request of sizes
            work, info = gelsy_lwork(m, n, nrhs, 10*np.finfo(dtype).eps)
            lwork = int(np.real(work))

            jptv = np.zeros((a1.shape[1], 1), dtype=np.int32)
            v, x, j, rank, info = gelsy(a1, b1, jptv, np.finfo(dtype).eps,
                                        lwork, False, False)
            assert_allclose(x[:-1], np.array([-14.333333333333323,
                                             14.999999999999991], dtype=dtype),
                            rtol=25*np.finfo(dtype).eps)

        for dtype in COMPLEX_DTYPES:
            a1 = np.array([[1.0+4.0j, 2.0],
                          [4.0+0.5j, 5.0-3.0j],
                          [7.0-2.0j, 8.0+0.7j]], dtype=dtype)
            b1 = np.array([16.0, 17.0+2.0j, 20.0-4.0j], dtype=dtype)
            gelsy, gelsy_lwork = get_lapack_funcs(('gelsy', 'gelss_lwork'),
                                                  (a1, b1))

            m, n = a1.shape
            if len(b1.shape) == 2:
                nrhs = b1.shape[1]
            else:
                nrhs = 1

            # Request of sizes
            work, info = gelsy_lwork(m, n, nrhs, 10*np.finfo(dtype).eps)
            lwork = int(np.real(work))

            jptv = np.zeros((a1.shape[1], 1), dtype=np.int32)
            v, x, j, rank, info = gelsy(a1, b1, jptv, np.finfo(dtype).eps,
                                        lwork, False, False)
            assert_allclose(x[:-1],
                            np.array([1.161753632288328-1.901075709391912j,
                                      1.735882340522193+1.521240901196909j],
                                     dtype=dtype),
                            rtol=25*np.finfo(dtype).eps)


class TestRegression(object):

    def test_ticket_1645(self):
        # Check that RQ routines have correct lwork
        for dtype in DTYPES:
            a = np.zeros((300, 2), dtype=dtype)

            gerqf, = get_lapack_funcs(['gerqf'], [a])
            assert_raises(Exception, gerqf, a, lwork=2)
            rq, tau, work, info = gerqf(a)

            if dtype in REAL_DTYPES:
                orgrq, = get_lapack_funcs(['orgrq'], [a])
                assert_raises(Exception, orgrq, rq[-2:], tau, lwork=1)
                orgrq(rq[-2:], tau, lwork=2)
            elif dtype in COMPLEX_DTYPES:
                ungrq, = get_lapack_funcs(['ungrq'], [a])
                assert_raises(Exception, ungrq, rq[-2:], tau, lwork=1)
                ungrq(rq[-2:], tau, lwork=2)


class TestDpotr(object):
    def test_gh_2691(self):
        # 'lower' argument of dportf/dpotri
        for lower in [True, False]:
            for clean in [True, False]:
                np.random.seed(42)
                x = np.random.normal(size=(3, 3))
                a = x.dot(x.T)

                dpotrf, dpotri = get_lapack_funcs(("potrf", "potri"), (a, ))

                c, info = dpotrf(a, lower, clean=clean)
                dpt = dpotri(c, lower)[0]

                if lower:
                    assert_allclose(np.tril(dpt), np.tril(inv(a)))
                else:
                    assert_allclose(np.triu(dpt), np.triu(inv(a)))


class TestDlasd4(object):
    def test_sing_val_update(self):

        sigmas = np.array([4., 3., 2., 0])
        m_vec = np.array([3.12, 5.7, -4.8, -2.2])

        M = np.hstack((np.vstack((np.diag(sigmas[0:-1]),
                       np.zeros((1, len(m_vec) - 1)))), m_vec[:, np.newaxis]))
        SM = svd(M, full_matrices=False, compute_uv=False, overwrite_a=False,
                 check_finite=False)

        it_len = len(sigmas)
        sgm = np.concatenate((sigmas[::-1], (sigmas[0] +
                              it_len*np.sqrt(np.sum(np.power(m_vec, 2))),)))
        mvc = np.concatenate((m_vec[::-1], (0,)))

        lasd4 = get_lapack_funcs('lasd4', (sigmas,))

        roots = []
        for i in range(0, it_len):
            res = lasd4(i, sgm, mvc)
            roots.append(res[1])

            assert_((res[3] <= 0), "LAPACK root finding dlasd4 failed to find \
                                    the singular value %i" % i)
        roots = np.array(roots)[::-1]

        assert_((not np.any(np.isnan(roots)), "There are NaN roots"))
        assert_allclose(SM, roots, atol=100*np.finfo(np.float64).eps,
                        rtol=100*np.finfo(np.float64).eps)


def test_lartg():
    for dtype in 'fdFD':
        lartg = get_lapack_funcs('lartg', dtype=dtype)

        f = np.array(3, dtype)
        g = np.array(4, dtype)

        if np.iscomplexobj(g):
            g *= 1j

        cs, sn, r = lartg(f, g)

        assert_allclose(cs, 3.0/5.0)
        assert_allclose(r, 5.0)

        if np.iscomplexobj(g):
            assert_allclose(sn, -4.0j/5.0)
            assert_(type(r) == complex)
            assert_(type(cs) == float)
        else:
            assert_allclose(sn, 4.0/5.0)


def test_rot():
    # srot, drot from blas and crot and zrot from lapack.

    for dtype in 'fdFD':
        c = 0.6
        s = 0.8

        u = np.ones(4, dtype) * 3
        v = np.ones(4, dtype) * 4
        atol = 10**-(np.finfo(dtype).precision-1)

        if dtype in 'fd':
            rot = get_blas_funcs('rot', dtype=dtype)
            f = 4
        else:
            rot = get_lapack_funcs('rot', dtype=dtype)
            s *= -1j
            v *= 1j
            f = 4j

        assert_allclose(rot(u, v, c, s), [[5, 5, 5, 5],
                                          [0, 0, 0, 0]], atol=atol)
        assert_allclose(rot(u, v, c, s, n=2), [[5, 5, 3, 3],
                                               [0, 0, f, f]], atol=atol)
        assert_allclose(rot(u, v, c, s, offx=2, offy=2),
                        [[3, 3, 5, 5], [f, f, 0, 0]], atol=atol)
        assert_allclose(rot(u, v, c, s, incx=2, offy=2, n=2),
                        [[5, 3, 5, 3], [f, f, 0, 0]], atol=atol)
        assert_allclose(rot(u, v, c, s, offx=2, incy=2, n=2),
                        [[3, 3, 5, 5], [0, f, 0, f]], atol=atol)
        assert_allclose(rot(u, v, c, s, offx=2, incx=2, offy=2, incy=2, n=1),
                        [[3, 3, 5, 3], [f, f, 0, f]], atol=atol)
        assert_allclose(rot(u, v, c, s, incx=-2, incy=-2, n=2),
                        [[5, 3, 5, 3], [0, f, 0, f]], atol=atol)

        a, b = rot(u, v, c, s, overwrite_x=1, overwrite_y=1)
        assert_(a is u)
        assert_(b is v)
        assert_allclose(a, [5, 5, 5, 5], atol=atol)
        assert_allclose(b, [0, 0, 0, 0], atol=atol)


def test_larfg_larf():
    np.random.seed(1234)
    a0 = np.random.random((4, 4))
    a0 = a0.T.dot(a0)

    a0j = np.random.random((4, 4)) + 1j*np.random.random((4, 4))
    a0j = a0j.T.conj().dot(a0j)

    # our test here will be to do one step of reducing a hermetian matrix to
    # tridiagonal form using householder transforms.

    for dtype in 'fdFD':
        larfg, larf = get_lapack_funcs(['larfg', 'larf'], dtype=dtype)

        if dtype in 'FD':
            a = a0j.copy()
        else:
            a = a0.copy()

        # generate a householder transform to clear a[2:,0]
        alpha, x, tau = larfg(a.shape[0]-1, a[1, 0], a[2:, 0])

        # create expected output
        expected = np.zeros_like(a[:, 0])
        expected[0] = a[0, 0]
        expected[1] = alpha

        # assemble householder vector
        v = np.zeros_like(a[1:, 0])
        v[0] = 1.0
        v[1:] = x

        # apply transform from the left
        a[1:, :] = larf(v, tau.conjugate(), a[1:, :], np.zeros(a.shape[1]))

        # apply transform from the right
        a[:, 1:] = larf(v, tau, a[:,1:], np.zeros(a.shape[0]), side='R')

        assert_allclose(a[:, 0], expected, atol=1e-5)
        assert_allclose(a[0, :], expected, atol=1e-5)


@pytest.mark.xslow
def test_sgesdd_lwork_bug_workaround():
    # Test that SGESDD lwork is sufficiently large for LAPACK.
    #
    # This checks that workaround around an apparent LAPACK bug
    # actually works. cf. gh-5401
    #
    # xslow: requires 1GB+ of memory

    p = subprocess.Popen([sys.executable, '-c',
                          'import numpy as np; '
                          'from scipy.linalg import svd; '
                          'a = np.zeros([9537, 9537], dtype=np.float32); '
                          'svd(a)'],
                         stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT)

    # Check if it an error occurred within 5 sec; the computation can
    # take substantially longer, and we will not wait for it to finish
    for j in range(50):
        time.sleep(0.1)
        if p.poll() is not None:
            returncode = p.returncode
            break
    else:
        # Didn't exit in time -- probably entered computation.  The
        # error is raised before entering computation, so things are
        # probably OK.
        returncode = 0
        p.terminate()

    assert_equal(returncode, 0,
                 "Code apparently failed: " + p.stdout.read())


class TestSytrd(object):
    def test_sytrd(self):
        for dtype in REAL_DTYPES:
            # Assert that a 0x0 matrix raises an error
            A = np.zeros((0, 0), dtype=dtype)
            sytrd, sytrd_lwork = \
                get_lapack_funcs(('sytrd', 'sytrd_lwork'), (A,))
            assert_raises(ValueError, sytrd, A)

            # Tests for n = 1 currently fail with
            # ```
            # ValueError: failed to create intent(cache|hide)|optional array--
            # must have defined dimensions but got (0,)
            # ```
            # This is a NumPy issue
            # <https://github.com/numpy/numpy/issues/9617>.
            # TODO once the issue has been resolved, test for n=1

            # some upper triangular array
            n = 3
            A = np.zeros((n, n), dtype=dtype)
            A[np.triu_indices_from(A)] = \
                np.arange(1, n*(n+1)//2+1, dtype=dtype)

            # query lwork
            lwork, info = sytrd_lwork(n)
            assert_equal(info, 0)

            # check lower=1 behavior (shouldn't do much since the matrix is
            # upper triangular)
            data, d, e, tau, info = sytrd(A, lower=1, lwork=lwork)
            assert_equal(info, 0)

            assert_allclose(data, A, atol=5*np.finfo(dtype).eps, rtol=1.0)
            assert_allclose(d, np.diag(A))
            assert_allclose(e, 0.0)
            assert_allclose(tau, 0.0)

            # and now for the proper test (lower=0 is the default)
            data, d, e, tau, info = sytrd(A, lwork=lwork)
            assert_equal(info, 0)

            # assert Q^T*A*Q = tridiag(e, d, e)

            # build tridiagonal matrix
            T = np.zeros_like(A, dtype=dtype)
            k = np.arange(A.shape[0])
            T[k, k] = d
            k2 = np.arange(A.shape[0]-1)
            T[k2+1, k2] = e
            T[k2, k2+1] = e

            # build Q
            Q = np.eye(n, n, dtype=dtype)
            for i in range(n-1):
                v = np.zeros(n, dtype=dtype)
                v[:i] = data[:i, i+1]
                v[i] = 1.0
                H = np.eye(n, n, dtype=dtype) - tau[i] * np.outer(v, v)
                Q = np.dot(H, Q)

            # Make matrix fully symmetric
            i_lower = np.tril_indices(n, -1)
            A[i_lower] = A.T[i_lower]

            QTAQ = np.dot(Q.T, np.dot(A, Q))

            # disable rtol here since some values in QTAQ and T are very close
            # to 0.
            assert_allclose(QTAQ, T, atol=5*np.finfo(dtype).eps, rtol=1.0)


class TestHetrd(object):
    def test_hetrd(self):
        for real_dtype, complex_dtype in zip(REAL_DTYPES, COMPLEX_DTYPES):
            # Assert that a 0x0 matrix raises an error
            A = np.zeros((0, 0), dtype=complex_dtype)
            hetrd, hetrd_lwork = \
                get_lapack_funcs(('hetrd', 'hetrd_lwork'), (A,))
            assert_raises(ValueError, hetrd, A)

            # Tests for n = 1 currently fail with
            # ```
            # ValueError: failed to create intent(cache|hide)|optional array--
            # must have defined dimensions but got (0,)
            # ```
            # This is a NumPy issue
            # <https://github.com/numpy/numpy/issues/9617>.
            # TODO once the issue has been resolved, test for n=1

            # some upper triangular array
            n = 3
            A = np.zeros((n, n), dtype=complex_dtype)
            A[np.triu_indices_from(A)] = (
                np.arange(1, n*(n+1)//2+1, dtype=real_dtype)
                + 1j * np.arange(1, n*(n+1)//2+1, dtype=real_dtype)
                )
            np.fill_diagonal(A, np.real(np.diag(A)))

            # query lwork
            lwork, info = hetrd_lwork(n)
            assert_equal(info, 0)

            # check lower=1 behavior (shouldn't do much since the matrix is
            # upper triangular)
            data, d, e, tau, info = hetrd(A, lower=1, lwork=lwork)
            assert_equal(info, 0)

            assert_allclose(data, A, atol=5*np.finfo(real_dtype).eps, rtol=1.0)

            assert_allclose(d, np.real(np.diag(A)))
            assert_allclose(e, 0.0)
            assert_allclose(tau, 0.0)

            # and now for the proper test (lower=0 is the default)
            data, d, e, tau, info = hetrd(A, lwork=lwork)
            assert_equal(info, 0)

            # assert Q^T*A*Q = tridiag(e, d, e)

            # build tridiagonal matrix
            T = np.zeros_like(A, dtype=real_dtype)
            k = np.arange(A.shape[0], dtype=int)
            T[k, k] = d
            k2 = np.arange(A.shape[0]-1, dtype=int)
            T[k2+1, k2] = e
            T[k2, k2+1] = e

            # build Q
            Q = np.eye(n, n, dtype=complex_dtype)
            for i in range(n-1):
                v = np.zeros(n, dtype=complex_dtype)
                v[:i] = data[:i, i+1]
                v[i] = 1.0
                H = np.eye(n, n, dtype=complex_dtype) \
                    - tau[i] * np.outer(v, np.conj(v))
                Q = np.dot(H, Q)

            # Make matrix fully Hermetian
            i_lower = np.tril_indices(n, -1)
            A[i_lower] = np.conj(A.T[i_lower])

            QHAQ = np.dot(np.conj(Q.T), np.dot(A, Q))

            # disable rtol here since some values in QTAQ and T are very close
            # to 0.
            assert_allclose(
                QHAQ, T, atol=10*np.finfo(real_dtype).eps, rtol=1.0
                )


def test_gglse():
    # Example data taken from NAG manual
    for ind, dtype in enumerate(DTYPES):
        # DTYPES = <s,d,c,z> gglse
        func, func_lwork = get_lapack_funcs(('gglse', 'gglse_lwork'),
                                            dtype=dtype)
        lwork = _compute_lwork(func_lwork, m=6, n=4, p=2)
        # For <s,d>gglse
        if ind < 2:
            a = np.array([[-0.57, -1.28, -0.39, 0.25],
                          [-1.93, 1.08, -0.31, -2.14],
                          [2.30, 0.24, 0.40, -0.35],
                          [-1.93, 0.64, -0.66, 0.08],
                          [0.15, 0.30, 0.15, -2.13],
                          [-0.02, 1.03, -1.43, 0.50]], dtype=dtype)
            c = np.array([-1.50, -2.14, 1.23, -0.54, -1.68, 0.82], dtype=dtype)
            d = np.array([0., 0.], dtype=dtype)
        # For <s,d>gglse
        else:
            a = np.array([[0.96-0.81j, -0.03+0.96j, -0.91+2.06j, -0.05+0.41j],
                          [-0.98+1.98j, -1.20+0.19j, -0.66+0.42j, -0.81+0.56j],
                          [0.62-0.46j, 1.01+0.02j, 0.63-0.17j, -1.11+0.60j],
                          [0.37+0.38j, 0.19-0.54j, -0.98-0.36j, 0.22-0.20j],
                          [0.83+0.51j, 0.20+0.01j, -0.17-0.46j, 1.47+1.59j],
                          [1.08-0.28j, 0.20-0.12j, -0.07+1.23j, 0.26+0.26j]])
            c = np.array([[-2.54+0.09j],
                          [1.65-2.26j],
                          [-2.11-3.96j],
                          [1.82+3.30j],
                          [-6.41+3.77j],
                          [2.07+0.66j]])
            d = np.zeros(2, dtype=dtype)

        b = np.array([[1., 0., -1., 0.], [0., 1., 0., -1.]], dtype=dtype)

        _, _, _, result, _ = func(a, b, c, d, lwork=lwork)
        if ind < 2:
            expected = np.array([0.48904455,
                                 0.99754786,
                                 0.48904455,
                                 0.99754786])
        else:
            expected = np.array([1.08742917-1.96205783j,
                                 -0.74093902+3.72973919j,
                                 1.08742917-1.96205759j,
                                 -0.74093896+3.72973895j])
        assert_array_almost_equal(result, expected, decimal=4)


def test_sycon_hecon():
    seed(1234)
    for ind, dtype in enumerate(DTYPES+COMPLEX_DTYPES):
        # DTYPES + COMPLEX DTYPES = <s,d,c,z> sycon + <c,z>hecon
        n = 10
        # For <s,d,c,z>sycon
        if ind < 4:
            func_lwork = get_lapack_funcs('sytrf_lwork', dtype=dtype)
            funcon, functrf = get_lapack_funcs(('sycon', 'sytrf'), dtype=dtype)
            A = (rand(n, n)).astype(dtype)
        # For <c,z>hecon
        else:
            func_lwork = get_lapack_funcs('hetrf_lwork', dtype=dtype)
            funcon, functrf = get_lapack_funcs(('hecon', 'hetrf'), dtype=dtype)
            A = (rand(n, n) + rand(n, n)*1j).astype(dtype)

        # Since sycon only refers to upper/lower part, conj() is safe here.
        A = (A + A.conj().T)/2 + 2*np.eye(n, dtype=dtype)

        anorm = np.linalg.norm(A, 1)
        lwork = _compute_lwork(func_lwork, n)
        ldu, ipiv, _ = functrf(A, lwork=lwork, lower=1)
        rcond, _ = funcon(a=ldu, ipiv=ipiv, anorm=anorm, lower=1)
        # The error is at most 1-fold
        assert_(abs(1/rcond - np.linalg.cond(A, p=1))*rcond < 1)


def test_sygst():
    seed(1234)
    for ind, dtype in enumerate(REAL_DTYPES):
        # DTYPES = <s,d> sygst
        n = 10

        potrf, sygst, syevd, sygvd = get_lapack_funcs(('potrf', 'sygst', 'syevd', 'sygvd'), dtype=dtype)

        A = rand(n, n).astype(dtype)
        A = (A + A.T)/2
        # B must be positive definite
        B = rand(n, n).astype(dtype)
        B = (B + B.T)/2 + 2 * np.eye(n, dtype=dtype)

        # Perform eig (sygvd)
        _, eig_gvd, info = sygvd(A, B)
        assert_(info == 0)

        # Convert to std problem potrf
        b, info = potrf(B)
        assert_(info == 0)
        a, info = sygst(A, b)
        assert_(info == 0)

        eig, _, info = syevd(a)
        assert_(info == 0)
        assert_allclose(eig, eig_gvd, rtol=1e-4)


def test_hegst():
    seed(1234)
    for ind, dtype in enumerate(COMPLEX_DTYPES):
        # DTYPES = <c,z> hegst
        n = 10

        potrf, hegst, heevd, hegvd = get_lapack_funcs(('potrf', 'hegst', 'heevd', 'hegvd'), dtype=dtype)

        A = rand(n, n).astype(dtype) + 1j * rand(n, n).astype(dtype)
        A = (A + A.conj().T)/2
        # B must be positive definite
        B = rand(n, n).astype(dtype) + 1j * rand(n, n).astype(dtype)
        B = (B + B.conj().T)/2 + 2 * np.eye(n, dtype=dtype)

        # Perform eig (hegvd)
        _, eig_gvd, info = hegvd(A, B)
        assert_(info == 0)

        # Convert to std problem potrf
        b, info = potrf(B)
        assert_(info == 0)
        a, info = hegst(A, b)
        assert_(info == 0)

        eig, _, info = heevd(a)
        assert_(info == 0)
        assert_allclose(eig, eig_gvd, rtol=1e-4)
