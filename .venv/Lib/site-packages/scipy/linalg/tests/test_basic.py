import os
import platform
import itertools
import warnings

import numpy as np
from numpy import (arange, array, dot, zeros, identity, conjugate, transpose,
                   float32)

from numpy.testing import (assert_equal, assert_almost_equal, assert_,
                           assert_array_almost_equal, assert_allclose,
                           assert_array_equal)
import pytest
from pytest import raises as assert_raises

from scipy.linalg import (solve, inv, det, lstsq, pinv, pinvh, norm,
                          solve_banded, solveh_banded, solve_triangular,
                          solve_circulant, circulant, LinAlgError, block_diag,
                          matrix_balance, qr, LinAlgWarning)

from scipy.linalg._testutils import assert_no_overwrite
from scipy._lib._testutils import check_free_memory, IS_MUSL
from scipy.linalg.blas import HAS_ILP64
from scipy.conftest import skip_xp_invalid_arg

REAL_DTYPES = (np.float32, np.float64, np.longdouble)
COMPLEX_DTYPES = (np.complex64, np.complex128, np.clongdouble)
DTYPES = REAL_DTYPES + COMPLEX_DTYPES


parametrize_overwrite_arg = pytest.mark.parametrize(
    "overwrite_kw", [{"overwrite_a": True}, {"overwrite_a": False}, {}],
    ids=["True", "False", "None"]
)


parametrize_overwrite_b_arg = pytest.mark.parametrize(
    "overwrite_b_kw", [{"overwrite_b": True}, {"overwrite_b": False}, {}],
    ids=["True", "False", "None"]
)


def _eps_cast(dtyp):
    """Get the epsilon for dtype, possibly downcast to BLAS types."""
    dt = dtyp
    if dt == np.longdouble:
        dt = np.float64
    elif dt == np.clongdouble:
        dt = np.complex128
    return np.finfo(dt).eps


class TestSolveBanded:

    def test_real(self):
        a = array([[1.0, 20, 0, 0],
                   [-30, 4, 6, 0],
                   [2, 1, 20, 2],
                   [0, -1, 7, 14]])
        ab = array([[0.0, 20, 6, 2],
                    [1, 4, 20, 14],
                    [-30, 1, 7, 0],
                    [2, -1, 0, 0]])
        l, u = 2, 1
        b4 = array([10.0, 0.0, 2.0, 14.0])
        b4by1 = b4.reshape(-1, 1)
        b4by2 = array([[2, 1],
                       [-30, 4],
                       [2, 3],
                       [1, 3]])
        b4by4 = array([[1, 0, 0, 0],
                       [0, 0, 0, 1],
                       [0, 1, 0, 0],
                       [0, 1, 0, 0]])
        for b in [b4, b4by1, b4by2, b4by4]:
            x = solve_banded((l, u), ab, b)
            assert_array_almost_equal(dot(a, x), b)

    def test_complex(self):
        a = array([[1.0, 20, 0, 0],
                   [-30, 4, 6, 0],
                   [2j, 1, 20, 2j],
                   [0, -1, 7, 14]])
        ab = array([[0.0, 20, 6, 2j],
                    [1, 4, 20, 14],
                    [-30, 1, 7, 0],
                    [2j, -1, 0, 0]])
        l, u = 2, 1
        b4 = array([10.0, 0.0, 2.0, 14.0j])
        b4by1 = b4.reshape(-1, 1)
        b4by2 = array([[2, 1],
                       [-30, 4],
                       [2, 3],
                       [1, 3]])
        b4by4 = array([[1, 0, 0, 0],
                       [0, 0, 0, 1j],
                       [0, 1, 0, 0],
                       [0, 1, 0, 0]])
        for b in [b4, b4by1, b4by2, b4by4]:
            x = solve_banded((l, u), ab, b)
            assert_array_almost_equal(dot(a, x), b)

    def test_tridiag_real(self):
        ab = array([[0.0, 20, 6, 2],
                   [1, 4, 20, 14],
                   [-30, 1, 7, 0]])
        a = np.diag(ab[0, 1:], 1) + np.diag(ab[1, :], 0) + np.diag(
                                                                ab[2, :-1], -1)
        b4 = array([10.0, 0.0, 2.0, 14.0])
        b4by1 = b4.reshape(-1, 1)
        b4by2 = array([[2, 1],
                       [-30, 4],
                       [2, 3],
                       [1, 3]])
        b4by4 = array([[1, 0, 0, 0],
                       [0, 0, 0, 1],
                       [0, 1, 0, 0],
                       [0, 1, 0, 0]])
        for b in [b4, b4by1, b4by2, b4by4]:
            x = solve_banded((1, 1), ab, b)
            assert_array_almost_equal(dot(a, x), b)

    def test_tridiag_complex(self):
        ab = array([[0.0, 20, 6, 2j],
                   [1, 4, 20, 14],
                   [-30, 1, 7, 0]])
        a = np.diag(ab[0, 1:], 1) + np.diag(ab[1, :], 0) + np.diag(
                                                               ab[2, :-1], -1)
        b4 = array([10.0, 0.0, 2.0, 14.0j])
        b4by1 = b4.reshape(-1, 1)
        b4by2 = array([[2, 1],
                       [-30, 4],
                       [2, 3],
                       [1, 3]])
        b4by4 = array([[1, 0, 0, 0],
                       [0, 0, 0, 1],
                       [0, 1, 0, 0],
                       [0, 1, 0, 0]])
        for b in [b4, b4by1, b4by2, b4by4]:
            x = solve_banded((1, 1), ab, b)
            assert_array_almost_equal(dot(a, x), b)

    def test_check_finite(self):
        a = array([[1.0, 20, 0, 0],
                   [-30, 4, 6, 0],
                   [2, 1, 20, 2],
                   [0, -1, 7, 14]])
        ab = array([[0.0, 20, 6, 2],
                    [1, 4, 20, 14],
                    [-30, 1, 7, 0],
                    [2, -1, 0, 0]])
        l, u = 2, 1
        b4 = array([10.0, 0.0, 2.0, 14.0])
        x = solve_banded((l, u), ab, b4, check_finite=False)
        assert_array_almost_equal(dot(a, x), b4)

    def test_bad_shape(self):
        ab = array([[0.0, 20, 6, 2],
                    [1, 4, 20, 14],
                    [-30, 1, 7, 0],
                    [2, -1, 0, 0]])
        l, u = 2, 1
        bad = array([1.0, 2.0, 3.0, 4.0]).reshape(-1, 4)
        assert_raises(ValueError, solve_banded, (l, u), ab, bad)
        assert_raises(ValueError, solve_banded, (l, u), ab, [1.0, 2.0])

        # Values of (l,u) are not compatible with ab.
        assert_raises(ValueError, solve_banded, (1, 1), ab, [1.0, 2.0])

    def test_1x1(self):
        # gh-8906 noted that the case of A@x = b with 1x1 A was handled
        # incorrectly; check that this is resolved. Typical case:
        # nupper == nlower == 0
        # A = [[2]]
        b = array([[1., 2., 3.]])
        ref = array([[0.5, 1.0, 1.5]])
        x = solve_banded((0, 0), [[2]], b)
        assert_allclose(x, ref, rtol=1e-15)

        # However, the user *can* represent the same system with garbage rows
        # in `ab`. Test the case with `nupper == 1, nlower == 1`.
        x = solve_banded((1, 1), [[0], [2], [0]], b)
        assert_allclose(x, ref, rtol=1e-15)
        assert_equal(x.dtype, np.dtype('f8'))
        assert_array_equal(b, [[1.0, 2.0, 3.0]])

    def test_native_list_arguments(self):
        a = [[1.0, 20, 0, 0],
             [-30, 4, 6, 0],
             [2, 1, 20, 2],
             [0, -1, 7, 14]]
        ab = [[0.0, 20, 6, 2],
              [1, 4, 20, 14],
              [-30, 1, 7, 0],
              [2, -1, 0, 0]]
        l, u = 2, 1
        b = [10.0, 0.0, 2.0, 14.0]
        x = solve_banded((l, u), ab, b)
        assert_array_almost_equal(dot(a, x), b)

    @pytest.mark.parametrize('dt_ab', [int, float, np.float32, complex, np.complex64])
    @pytest.mark.parametrize('dt_b', [int, float, np.float32, complex, np.complex64])
    def test_empty(self, dt_ab, dt_b):
        # ab contains one empty row corresponding to the diagonal
        ab = np.array([[]], dtype=dt_ab)
        b = np.array([], dtype=dt_b)
        x = solve_banded((0, 0), ab, b)

        assert x.shape == (0,)
        assert x.dtype == solve(np.eye(1, dtype=dt_ab), np.ones(1, dtype=dt_b)).dtype

        b = np.empty((0, 0), dtype=dt_b)
        x = solve_banded((0, 0), ab, b)

        assert x.shape == (0, 0)
        assert x.dtype == solve(np.eye(1, dtype=dt_ab), np.ones(1, dtype=dt_b)).dtype


class TestSolveHBanded:

    def test_01_upper(self):
        # Solve
        # [ 4 1 2 0]     [1]
        # [ 1 4 1 2] X = [4]
        # [ 2 1 4 1]     [1]
        # [ 0 2 1 4]     [2]
        # with the RHS as a 1D array.
        ab = array([[0.0, 0.0, 2.0, 2.0],
                    [-99, 1.0, 1.0, 1.0],
                    [4.0, 4.0, 4.0, 4.0]])
        b = array([1.0, 4.0, 1.0, 2.0])
        x = solveh_banded(ab, b)
        assert_array_almost_equal(x, [0.0, 1.0, 0.0, 0.0])

    def test_02_upper(self):
        # Solve
        # [ 4 1 2 0]     [1 6]
        # [ 1 4 1 2] X = [4 2]
        # [ 2 1 4 1]     [1 6]
        # [ 0 2 1 4]     [2 1]
        #
        ab = array([[0.0, 0.0, 2.0, 2.0],
                    [-99, 1.0, 1.0, 1.0],
                    [4.0, 4.0, 4.0, 4.0]])
        b = array([[1.0, 6.0],
                   [4.0, 2.0],
                   [1.0, 6.0],
                   [2.0, 1.0]])
        x = solveh_banded(ab, b)
        expected = array([[0.0, 1.0],
                          [1.0, 0.0],
                          [0.0, 1.0],
                          [0.0, 0.0]])
        assert_array_almost_equal(x, expected)

    def test_03_upper(self):
        # Solve
        # [ 4 1 2 0]     [1]
        # [ 1 4 1 2] X = [4]
        # [ 2 1 4 1]     [1]
        # [ 0 2 1 4]     [2]
        # with the RHS as a 2D array with shape (3,1).
        ab = array([[0.0, 0.0, 2.0, 2.0],
                    [-99, 1.0, 1.0, 1.0],
                    [4.0, 4.0, 4.0, 4.0]])
        b = array([1.0, 4.0, 1.0, 2.0]).reshape(-1, 1)
        x = solveh_banded(ab, b)
        assert_array_almost_equal(x, array([0., 1., 0., 0.]).reshape(-1, 1))

    def test_01_lower(self):
        # Solve
        # [ 4 1 2 0]     [1]
        # [ 1 4 1 2] X = [4]
        # [ 2 1 4 1]     [1]
        # [ 0 2 1 4]     [2]
        #
        ab = array([[4.0, 4.0, 4.0, 4.0],
                    [1.0, 1.0, 1.0, -99],
                    [2.0, 2.0, 0.0, 0.0]])
        b = array([1.0, 4.0, 1.0, 2.0])
        x = solveh_banded(ab, b, lower=True)
        assert_array_almost_equal(x, [0.0, 1.0, 0.0, 0.0])

    def test_02_lower(self):
        # Solve
        # [ 4 1 2 0]     [1 6]
        # [ 1 4 1 2] X = [4 2]
        # [ 2 1 4 1]     [1 6]
        # [ 0 2 1 4]     [2 1]
        #
        ab = array([[4.0, 4.0, 4.0, 4.0],
                    [1.0, 1.0, 1.0, -99],
                    [2.0, 2.0, 0.0, 0.0]])
        b = array([[1.0, 6.0],
                   [4.0, 2.0],
                   [1.0, 6.0],
                   [2.0, 1.0]])
        x = solveh_banded(ab, b, lower=True)
        expected = array([[0.0, 1.0],
                          [1.0, 0.0],
                          [0.0, 1.0],
                          [0.0, 0.0]])
        assert_array_almost_equal(x, expected)

    def test_01_float32(self):
        # Solve
        # [ 4 1 2 0]     [1]
        # [ 1 4 1 2] X = [4]
        # [ 2 1 4 1]     [1]
        # [ 0 2 1 4]     [2]
        #
        ab = array([[0.0, 0.0, 2.0, 2.0],
                    [-99, 1.0, 1.0, 1.0],
                    [4.0, 4.0, 4.0, 4.0]], dtype=float32)
        b = array([1.0, 4.0, 1.0, 2.0], dtype=float32)
        x = solveh_banded(ab, b)
        assert_array_almost_equal(x, [0.0, 1.0, 0.0, 0.0])

    def test_02_float32(self):
        # Solve
        # [ 4 1 2 0]     [1 6]
        # [ 1 4 1 2] X = [4 2]
        # [ 2 1 4 1]     [1 6]
        # [ 0 2 1 4]     [2 1]
        #
        ab = array([[0.0, 0.0, 2.0, 2.0],
                    [-99, 1.0, 1.0, 1.0],
                    [4.0, 4.0, 4.0, 4.0]], dtype=float32)
        b = array([[1.0, 6.0],
                   [4.0, 2.0],
                   [1.0, 6.0],
                   [2.0, 1.0]], dtype=float32)
        x = solveh_banded(ab, b)
        expected = array([[0.0, 1.0],
                          [1.0, 0.0],
                          [0.0, 1.0],
                          [0.0, 0.0]])
        assert_array_almost_equal(x, expected)

    def test_01_complex(self):
        # Solve
        # [ 4 -j  2  0]     [2-j]
        # [ j  4 -j  2] X = [4-j]
        # [ 2  j  4 -j]     [4+j]
        # [ 0  2  j  4]     [2+j]
        #
        ab = array([[0.0, 0.0, 2.0, 2.0],
                    [-99, -1.0j, -1.0j, -1.0j],
                    [4.0, 4.0, 4.0, 4.0]])
        b = array([2-1.0j, 4.0-1j, 4+1j, 2+1j])
        x = solveh_banded(ab, b)
        assert_array_almost_equal(x, [0.0, 1.0, 1.0, 0.0])

    def test_02_complex(self):
        # Solve
        # [ 4 -j  2  0]     [2-j 2+4j]
        # [ j  4 -j  2] X = [4-j -1-j]
        # [ 2  j  4 -j]     [4+j 4+2j]
        # [ 0  2  j  4]     [2+j j]
        #
        ab = array([[0.0, 0.0, 2.0, 2.0],
                    [-99, -1.0j, -1.0j, -1.0j],
                    [4.0, 4.0, 4.0, 4.0]])
        b = array([[2-1j, 2+4j],
                   [4.0-1j, -1-1j],
                   [4.0+1j, 4+2j],
                   [2+1j, 1j]])
        x = solveh_banded(ab, b)
        expected = array([[0.0, 1.0j],
                          [1.0, 0.0],
                          [1.0, 1.0],
                          [0.0, 0.0]])
        assert_array_almost_equal(x, expected)

    def test_tridiag_01_upper(self):
        # Solve
        # [ 4 1 0]     [1]
        # [ 1 4 1] X = [4]
        # [ 0 1 4]     [1]
        # with the RHS as a 1D array.
        ab = array([[-99, 1.0, 1.0], [4.0, 4.0, 4.0]])
        b = array([1.0, 4.0, 1.0])
        x = solveh_banded(ab, b)
        assert_array_almost_equal(x, [0.0, 1.0, 0.0])

    def test_tridiag_02_upper(self):
        # Solve
        # [ 4 1 0]     [1 4]
        # [ 1 4 1] X = [4 2]
        # [ 0 1 4]     [1 4]
        #
        ab = array([[-99, 1.0, 1.0],
                    [4.0, 4.0, 4.0]])
        b = array([[1.0, 4.0],
                   [4.0, 2.0],
                   [1.0, 4.0]])
        x = solveh_banded(ab, b)
        expected = array([[0.0, 1.0],
                          [1.0, 0.0],
                          [0.0, 1.0]])
        assert_array_almost_equal(x, expected)

    def test_tridiag_03_upper(self):
        # Solve
        # [ 4 1 0]     [1]
        # [ 1 4 1] X = [4]
        # [ 0 1 4]     [1]
        # with the RHS as a 2D array with shape (3,1).
        ab = array([[-99, 1.0, 1.0], [4.0, 4.0, 4.0]])
        b = array([1.0, 4.0, 1.0]).reshape(-1, 1)
        x = solveh_banded(ab, b)
        assert_array_almost_equal(x, array([0.0, 1.0, 0.0]).reshape(-1, 1))

    def test_tridiag_01_lower(self):
        # Solve
        # [ 4 1 0]     [1]
        # [ 1 4 1] X = [4]
        # [ 0 1 4]     [1]
        #
        ab = array([[4.0, 4.0, 4.0],
                    [1.0, 1.0, -99]])
        b = array([1.0, 4.0, 1.0])
        x = solveh_banded(ab, b, lower=True)
        assert_array_almost_equal(x, [0.0, 1.0, 0.0])

    def test_tridiag_02_lower(self):
        # Solve
        # [ 4 1 0]     [1 4]
        # [ 1 4 1] X = [4 2]
        # [ 0 1 4]     [1 4]
        #
        ab = array([[4.0, 4.0, 4.0],
                    [1.0, 1.0, -99]])
        b = array([[1.0, 4.0],
                   [4.0, 2.0],
                   [1.0, 4.0]])
        x = solveh_banded(ab, b, lower=True)
        expected = array([[0.0, 1.0],
                          [1.0, 0.0],
                          [0.0, 1.0]])
        assert_array_almost_equal(x, expected)

    def test_tridiag_01_float32(self):
        # Solve
        # [ 4 1 0]     [1]
        # [ 1 4 1] X = [4]
        # [ 0 1 4]     [1]
        #
        ab = array([[-99, 1.0, 1.0], [4.0, 4.0, 4.0]], dtype=float32)
        b = array([1.0, 4.0, 1.0], dtype=float32)
        x = solveh_banded(ab, b)
        assert_array_almost_equal(x, [0.0, 1.0, 0.0])

    def test_tridiag_02_float32(self):
        # Solve
        # [ 4 1 0]     [1 4]
        # [ 1 4 1] X = [4 2]
        # [ 0 1 4]     [1 4]
        #
        ab = array([[-99, 1.0, 1.0],
                    [4.0, 4.0, 4.0]], dtype=float32)
        b = array([[1.0, 4.0],
                   [4.0, 2.0],
                   [1.0, 4.0]], dtype=float32)
        x = solveh_banded(ab, b)
        expected = array([[0.0, 1.0],
                          [1.0, 0.0],
                          [0.0, 1.0]])
        assert_array_almost_equal(x, expected)

    def test_tridiag_01_complex(self):
        # Solve
        # [ 4 -j 0]     [ -j]
        # [ j 4 -j] X = [4-j]
        # [ 0 j  4]     [4+j]
        #
        ab = array([[-99, -1.0j, -1.0j], [4.0, 4.0, 4.0]])
        b = array([-1.0j, 4.0-1j, 4+1j])
        x = solveh_banded(ab, b)
        assert_array_almost_equal(x, [0.0, 1.0, 1.0])

    def test_tridiag_02_complex(self):
        # Solve
        # [ 4 -j 0]     [ -j    4j]
        # [ j 4 -j] X = [4-j  -1-j]
        # [ 0 j  4]     [4+j   4  ]
        #
        ab = array([[-99, -1.0j, -1.0j],
                    [4.0, 4.0, 4.0]])
        b = array([[-1j, 4.0j],
                   [4.0-1j, -1.0-1j],
                   [4.0+1j, 4.0]])
        x = solveh_banded(ab, b)
        expected = array([[0.0, 1.0j],
                          [1.0, 0.0],
                          [1.0, 1.0]])
        assert_array_almost_equal(x, expected)

    def test_check_finite(self):
        # Solve
        # [ 4 1 0]     [1]
        # [ 1 4 1] X = [4]
        # [ 0 1 4]     [1]
        # with the RHS as a 1D array.
        ab = array([[-99, 1.0, 1.0], [4.0, 4.0, 4.0]])
        b = array([1.0, 4.0, 1.0])
        x = solveh_banded(ab, b, check_finite=False)
        assert_array_almost_equal(x, [0.0, 1.0, 0.0])

    def test_bad_shapes(self):
        ab = array([[-99, 1.0, 1.0],
                    [4.0, 4.0, 4.0]])
        b = array([[1.0, 4.0],
                   [4.0, 2.0]])
        assert_raises(ValueError, solveh_banded, ab, b)
        assert_raises(ValueError, solveh_banded, ab, [1.0, 2.0])
        assert_raises(ValueError, solveh_banded, ab, [1.0])

    def test_1x1(self):
        x = solveh_banded([[1]], [[1, 2, 3]])
        assert_array_equal(x, [[1.0, 2.0, 3.0]])
        assert_equal(x.dtype, np.dtype('f8'))

    def test_native_list_arguments(self):
        # Same as test_01_upper, using python's native list.
        ab = [[0.0, 0.0, 2.0, 2.0],
              [-99, 1.0, 1.0, 1.0],
              [4.0, 4.0, 4.0, 4.0]]
        b = [1.0, 4.0, 1.0, 2.0]
        x = solveh_banded(ab, b)
        assert_array_almost_equal(x, [0.0, 1.0, 0.0, 0.0])

    @pytest.mark.parametrize('dt_ab', [int, float, np.float32, complex, np.complex64])
    @pytest.mark.parametrize('dt_b', [int, float, np.float32, complex, np.complex64])
    def test_empty(self, dt_ab, dt_b):
        # ab contains one empty row corresponding to the diagonal
        ab = np.array([[]], dtype=dt_ab)
        b = np.array([], dtype=dt_b)
        x = solveh_banded(ab, b)

        assert x.shape == (0,)
        assert x.dtype == solve(np.eye(1, dtype=dt_ab), np.ones(1, dtype=dt_b)).dtype

        b = np.empty((0, 0), dtype=dt_b)
        x = solveh_banded(ab, b)

        assert x.shape == (0, 0)
        assert x.dtype == solve(np.eye(1, dtype=dt_ab), np.ones(1, dtype=dt_b)).dtype


class TestSolve:
    def test_20Feb04_bug(self):
        a = [[1, 1], [1.0, 0]]  # ok
        x0 = solve(a, [1, 0j])
        assert_array_almost_equal(dot(a, x0), [1, 0])

        # gives failure with clapack.zgesv(..,rowmajor=0)
        a = [[1, 1], [1.2, 0]]
        b = [1, 0j]
        x0 = solve(a, b)
        assert_array_almost_equal(dot(a, x0), [1, 0])

    def test_simple(self):
        a = [[1, 20], [-30, 4]]
        for b in ([[1, 0], [0, 1]],
                  [1, 0],
                  [[2, 1], [-30, 4]]
                  ):
            x = solve(a, b)
            assert_array_almost_equal(dot(a, x), b)

    def test_simple_complex(self):
        a = array([[5, 2], [2j, 4]], 'D')
        for b in ([1j, 0],
                  [[1j, 1j], [0, 2]],
                  [1, 0j],
                  array([1, 0], 'D'),
                  ):
            x = solve(a, b)
            assert_array_almost_equal(dot(a, x), b)

    def test_simple_pos(self):
        a = [[2, 3], [3, 5]]
        for lower in [0, 1]:
            for b in ([[1, 0], [0, 1]],
                      [1, 0]
                      ):
                x = solve(a, b, assume_a='pos', lower=lower)
                assert_array_almost_equal(dot(a, x), b)

    def test_simple_pos_complexb(self):
        a = [[5, 2], [2, 4]]
        for b in ([1j, 0],
                  [[1j, 1j], [0, 2]],
                  ):
            x = solve(a, b, assume_a='pos')
            assert_array_almost_equal(dot(a, x), b)

    def test_simple_sym(self):
        a = [[2, 3], [3, -5]]
        for lower in [0, 1]:
            for b in ([[1, 0], [0, 1]],
                      [1, 0]
                      ):
                x = solve(a, b, assume_a='sym', lower=lower)
                assert_array_almost_equal(dot(a, x), b)

    def test_simple_sym_complexb(self):
        a = [[5, 2], [2, -4]]
        for b in ([1j, 0],
                  [[1j, 1j], [0, 2]]
                  ):
            x = solve(a, b, assume_a='sym')
            assert_array_almost_equal(dot(a, x), b)

    def test_simple_sym_complex(self):
        a = [[5, 2+1j], [2+1j, -4]]
        for b in ([1j, 0],
                  [1, 0],
                  [[1j, 1j], [0, 2]]
                  ):
            x = solve(a, b, assume_a='sym')
            assert_array_almost_equal(dot(a, x), b)

    def test_simple_her_actuallysym(self):
        a = [[2, 3], [3, -5]]
        for lower in [0, 1]:
            for b in ([[1, 0], [0, 1]],
                      [1, 0],
                      [1j, 0],
                      ):
                x = solve(a, b, assume_a='her', lower=lower)
                assert_array_almost_equal(dot(a, x), b)

    def test_simple_her(self):
        a = [[5, 2+1j], [2-1j, -4]]
        for b in ([1j, 0],
                  [1, 0],
                  [[1j, 1j], [0, 2]]
                  ):
            x = solve(a, b, assume_a='her')
            assert_array_almost_equal(dot(a, x), b)

    def test_nils_20Feb04(self):
        rng = np.random.default_rng(1234)
        n = 2
        A = rng.random([n, n])+rng.random([n, n])*1j
        X = zeros((n, n), 'D')
        Ainv = inv(A)
        R = identity(n)+identity(n)*0j
        for i in arange(0, n):
            r = R[:, i]
            X[:, i] = solve(A, r)
        assert_array_almost_equal(X, Ainv)

    def test_random(self):
        rng = np.random.default_rng(1234)
        n = 20
        a = rng.random([n, n])
        for i in range(n):
            a[i, i] = 20*(.1+a[i, i])
        for i in range(4):
            b = rng.random([n, 3])
            x = solve(a, b)
            assert_array_almost_equal(dot(a, x), b)

    def test_random_complex(self):
        rng = np.random.default_rng(1234)
        n = 20
        a = rng.random([n, n]) + 1j * rng.random([n, n])
        for i in range(n):
            a[i, i] = 20*(.1+a[i, i])
        for i in range(2):
            b = rng.random([n, 3])
            x = solve(a, b)
            assert_array_almost_equal(dot(a, x), b)

    def test_random_sym(self):
        rng = np.random.default_rng(1234)
        n = 20
        a = rng.random([n, n])
        for i in range(n):
            a[i, i] = abs(20*(.1+a[i, i]))
            for j in range(i):
                a[i, j] = a[j, i]
        for i in range(4):
            b = rng.random([n])
            x = solve(a, b, assume_a="pos")
            assert_array_almost_equal(dot(a, x), b)

    def test_random_sym_complex(self):
        rng = np.random.default_rng(1234)
        n = 20
        a = rng.random([n, n])
        a = a + 1j*rng.random([n, n])
        for i in range(n):
            a[i, i] = abs(20*(.1+a[i, i]))
            for j in range(i):
                a[i, j] = conjugate(a[j, i])
        b = rng.random([n])+2j*rng.random([n])
        for i in range(2):
            x = solve(a, b, assume_a="pos")
            assert_array_almost_equal(dot(a, x), b)

    def test_check_finite(self):
        a = [[1, 20], [-30, 4]]
        for b in ([[1, 0], [0, 1]], [1, 0],
                  [[2, 1], [-30, 4]]):
            x = solve(a, b, check_finite=False)
            assert_array_almost_equal(dot(a, x), b)

    def test_scalar_a_and_1D_b(self):
        a = 1
        b = [1, 2, 3]
        x = solve(a, b)
        assert_array_almost_equal(x.ravel(), b)
        assert_(x.shape == (3,), 'Scalar_a_1D_b test returned wrong shape')

    def test_simple2(self):
        a = np.array([[1.80, 2.88, 2.05, -0.89],
                      [525.00, -295.00, -95.00, -380.00],
                      [1.58, -2.69, -2.90, -1.04],
                      [-1.11, -0.66, -0.59, 0.80]])

        b = np.array([[9.52, 18.47],
                      [2435.00, 225.00],
                      [0.77, -13.28],
                      [-6.22, -6.21]])

        x = solve(a, b)
        assert_array_almost_equal(x, np.array([[1., -1, 3, -5],
                                               [3, 2, 4, 1]]).T)

    def test_simple_complex2(self):
        a = np.array([[-1.34+2.55j, 0.28+3.17j, -6.39-2.20j, 0.72-0.92j],
                      [-1.70-14.10j, 33.10-1.50j, -1.50+13.40j, 12.90+13.80j],
                      [-3.29-2.39j, -1.91+4.42j, -0.14-1.35j, 1.72+1.35j],
                      [2.41+0.39j, -0.56+1.47j, -0.83-0.69j, -1.96+0.67j]])

        b = np.array([[26.26+51.78j, 31.32-6.70j],
                      [64.30-86.80j, 158.60-14.20j],
                      [-5.75+25.31j, -2.15+30.19j],
                      [1.16+2.57j, -2.56+7.55j]])

        x = solve(a, b)
        assert_array_almost_equal(x, np. array([[1+1.j, -1-2.j],
                                                [2-3.j, 5+1.j],
                                                [-4-5.j, -3+4.j],
                                                [6.j, 2-3.j]]))

    @pytest.mark.parametrize("assume_a", ['her', 'sym'])
    def test_symmetric_hermitian(self, assume_a):
        # An upper triangular matrix will be used for symmetric/hermitian matrix a
        a = np.array([[-1.84, 0.11-0.11j, -1.78-1.18j, 3.91-1.50j],
                      [0, -4.63, -1.84+0.03j, 2.21+0.21j],
                      [0, 0, -8.87, 1.58-0.90j],
                      [0, 0, 0, -1.36]])
        b = np.array([[2.98-10.18j, 28.68-39.89j],
                      [-9.58+3.88j, -24.79-8.40j],
                      [-0.77-16.05j, 4.23-70.02j],
                      [7.79+5.48j, -35.39+18.01j]])

        a2 = a.T if assume_a == 'sym' else a.conj().T  # for testing `lower`
        a3 = a + a2                                    # for reference solution
        a3[np.arange(4), np.arange(4)] = np.diag(a)
        ref = solve(a3, b, assume_a='general')

        x = solve(a, b, assume_a=assume_a)
        assert_array_almost_equal(x, ref)
        # Also transpose(/conjugate) `a` and test for lower triangular data
        # This also tests gh-22265 resolution; otherwise, a warning would be emitted
        x = solve(a2, b, assume_a=assume_a, lower=True)
        assert_array_almost_equal(x, ref)

    def test_pos_and_sym(self):
        A = np.arange(1, 10).reshape(3, 3)
        x = solve(np.tril(A)/9, np.ones(3), assume_a='pos')
        assert_array_almost_equal(x, [9., 1.8, 1.])
        x = solve(np.tril(A)/9, np.ones(3), assume_a='sym')
        assert_array_almost_equal(x, [9., 1.8, 1.])

    def test_singularity(self):
        a = np.array([[1, 0, 0, 0, 0, 0, 1, 0, 1],
                      [1, 1, 1, 0, 0, 0, 1, 0, 1],
                      [0, 1, 1, 0, 0, 0, 1, 0, 1],
                      [1, 0, 1, 1, 1, 1, 0, 0, 0],
                      [1, 0, 1, 1, 1, 1, 0, 0, 0],
                      [1, 0, 1, 1, 1, 1, 0, 0, 0],
                      [1, 0, 1, 1, 1, 1, 0, 0, 0],
                      [1, 1, 1, 1, 1, 1, 1, 1, 1],
                      [1, 1, 1, 1, 1, 1, 1, 1, 1]])
        b = np.arange(9)[:, None]
        assert_raises(LinAlgError, solve, a, b)

    @pytest.mark.parametrize('structure',
                             ('diagonal', 'tridiagonal', 'lower triangular',
                              'upper triangular', 'symmetric', 'hermitian',
                              'positive definite', 'general', 'banded', None))
    def test_ill_condition_warning(self, structure):
        rng = np.random.default_rng(234859349452)
        n = 10
        d = np.logspace(0, 50, n)
        A = np.diag(d)
        b = rng.random(size=n)
        message = "(Ill-conditioned matrix|An ill-conditioned matrix)"
        with pytest.warns(LinAlgWarning, match=message):
            solve(A, b, assume_a=structure)

    @pytest.mark.parametrize('structure',
                             ('diagonal', 'tridiagonal', 'lower triangular',
                              'upper triangular', 'symmetric', 'hermitian',
                              'positive definite', 'general', None))
    def test_exactly_singular_gh22263(self, structure):
        n = 10
        A = np.zeros((n, n))
        b = np.ones(n)
        with (pytest.raises(LinAlgError, match="singular"), np.errstate(all='ignore')):
            solve(A, b, assume_a=structure)

    @pytest.mark.parametrize('b', [0, 1, [0, 1]])
    def test_singular_scalar(self, b):
        # regression test for gh-24355: scalar a=0 is singular
        # thus should raise the same error 

        with pytest.raises(LinAlgError):
            a = np.zeros((1, 1))
            solve(a, b)

        with pytest.raises(LinAlgError):
            solve(0, b)

        with pytest.raises(LinAlgError):
            solve([[0]], b)

    def test_multiple_rhs(self):
        a = np.eye(2)
        rng = np.random.default_rng(1234)
        b = rng.random((2, 12))
        x = solve(a, b)
        assert_array_almost_equal(x, b)

    def test_transposed_keyword(self):
        A = np.arange(9).reshape(3, 3) + 1
        x = solve(np.tril(A)/9, np.ones(3), transposed=True)
        assert_array_almost_equal(x, [1.2, 0.2, 1])
        x = solve(np.tril(A)/9, np.ones(3), transposed=False)
        assert_array_almost_equal(x, [9, -5.4, -1.2])

    @pytest.mark.skip(reason="1. why? 2. deprecate the kwarg altogether?")
    def test_transposed_notimplemented(self):
        a = np.eye(3).astype(complex)
        with assert_raises(NotImplementedError):
            solve(a, a, transposed=True)

    def test_nonsquare_a(self):
        assert_raises(ValueError, solve, [1, 2], 1)

    def test_size_mismatch_with_1D_b(self):
        assert_array_almost_equal(solve(np.eye(3), np.ones(3)), np.ones(3))
        assert_raises(ValueError, solve, np.eye(3), np.ones(4))

    def test_assume_a_keyword(self):
        assert_raises(ValueError, solve, 1, 1, assume_a='zxcv')

    @pytest.mark.parametrize("size", [10, 100])
    @pytest.mark.parametrize("assume_a", ['gen', 'sym', 'pos', 'her', 'tridiagonal'])
    @pytest.mark.parametrize(
        "dtype", [np.float32, np.float64, np.complex64, np.complex128]
    )
    def test_all_type_size_routine_combinations(self, size, dtype, assume_a):
        rng = np.random.default_rng(1234)
        is_complex = dtype in (np.complex64, np.complex128)

        a = rng.standard_normal((size, size)).astype(dtype)
        b = rng.standard_normal(size).astype(dtype)
        if is_complex:
            a += (1j*rng.standard_normal((size, size))).astype(dtype)

        if assume_a == 'sym':  # Can still be complex but only symmetric
            a = a + a.T
        elif assume_a == 'her':  # Handle hermitian matrices here instead
            a = a + a.T.conj()
        elif assume_a == 'pos':
            a = a.T.conj() @ a + 0.1*np.eye(size)
        elif assume_a == 'tridiagonal':
            a = (np.diag(np.diag(a)) +
                 np.diag(np.diag(a, 1), 1) +
                 np.diag(np.diag(a, -1), -1)
            )

        tol = 1e-12 if dtype in (np.float64, np.complex128) else 1e-6

        if assume_a in ['gen', 'sym', 'her']:
            # We revert the tolerance from before
            #   4b4a6e7c34fa4060533db38f9a819b98fa81476c
            if dtype in (np.float32, np.complex64):
                tol *= 10

        x = solve(a, b, assume_a=assume_a)
        assert_allclose(a @ x, b, atol=tol * size, rtol=tol * size)

        if assume_a == 'sym' and not is_complex:
            x = solve(a, b, assume_a=assume_a, transposed=True)
            assert_allclose(a @ x, b, atol=tol * size, rtol=tol * size)

    @pytest.mark.parametrize('dt_a', [int, float, np.float32, complex, np.complex64])
    @pytest.mark.parametrize('dt_b', [int, float, np.float32, complex, np.complex64])
    def test_empty(self, dt_a, dt_b):
        a = np.empty((0, 0), dtype=dt_a)
        b = np.empty(0, dtype=dt_b)
        x = solve(a, b)

        assert x.size == 0
        dt_nonempty = solve(np.eye(2, dtype=dt_a), np.ones(2, dtype=dt_b)).dtype
        assert x.dtype == dt_nonempty
        assert x.shape == np.linalg.solve(a, b).shape

        a = np.ones((3, 0, 2, 2), dtype=dt_a)
        b = np.ones((2, 4), dtype=dt_b)
        x = solve(a, b)
        assert x.shape == (3, 0, 2, 4)
        assert x.dtype == dt_nonempty

    def test_empty_rhs(self):
        a = np.eye(2)
        b = [[], []]
        x = solve(a, b)
        assert_(x.size == 0, 'Returned array is not empty')
        assert_(x.shape == (2, 0), 'Returned empty array shape is wrong')

    @pytest.mark.parametrize('dtype', [np.float64, np.complex128])
    @pytest.mark.parametrize('assume_a', ['diagonal', 'tridiagonal', 'banded',
                                          'lower triangular', 'upper triangular',
                                          'pos', 'positive definite',
                                          'symmetric', 'hermitian', 'banded',
                                          'general', 'sym', 'her', 'gen'])
    @pytest.mark.parametrize('nrhs', [(), (5,)])
    @pytest.mark.parametrize('transposed', [True, False])
    @pytest.mark.parametrize('overwrite', [True, False])
    @pytest.mark.parametrize('fortran', [True, False])
    def test_structure_detection(self, dtype, assume_a, nrhs, transposed,
                                 overwrite, fortran):
        rng = np.random.default_rng(982345982439826)
        n = 5 if not assume_a == 'banded' else 20
        b = rng.random(size=(n,) + nrhs)
        A = rng.random(size=(n, n))

        if np.issubdtype(dtype, np.complexfloating):
            b = b + rng.random(size=(n,) + nrhs) * 1j
            A = A + rng.random(size=(n, n)) * 1j

        if assume_a == 'diagonal':
            A = np.diag(np.diag(A))
        elif assume_a == 'lower triangular':
            A = np.tril(A)
        elif assume_a == 'upper triangular':
            A = np.triu(A)
        elif assume_a == 'tridiagonal':
            A = (np.diag(np.diag(A))
                 + np.diag(np.diag(A, -1), -1)
                 + np.diag(np.diag(A, 1), 1))
        elif assume_a == 'banded':
            A = np.triu(np.tril(A, 2), -1)
        elif assume_a in {'symmetric', 'sym'}:
            A = A + A.T
        elif assume_a in {'hermitian', 'her'}:
            A = A + A.conj().T
        elif assume_a in {'positive definite', 'pos'}:
            A = A @ A.T.conj()

        if fortran:
            A = np.asfortranarray(A)

        A_copy = A.copy(order='A')
        b_copy = b.copy()

        if np.issubdtype(dtype, np.complexfloating) and transposed:
            message = "scipy.linalg.solve can currently..."
            with pytest.raises(NotImplementedError, match=message):
                solve(A, b, overwrite_a=overwrite, overwrite_b=overwrite,
                      transposed=transposed)
            return

        res = solve(A, b, overwrite_a=overwrite, overwrite_b=overwrite,
                    transposed=transposed, assume_a=assume_a)

        # Check that solution this solution is *correct*
        ref = np.linalg.solve(A_copy.T if transposed else A_copy, b_copy)
        assert_allclose(res, ref)

        # Check that `solve` correctly identifies the structure and returns
        # *exactly* the same solution whether `assume_a` is specified or not
        if assume_a != 'banded':  # structure detection removed for banded
            assert_allclose(
                solve(A_copy, b_copy, transposed=transposed), res, atol=1e-15
            )

        # Check that overwrite was respected
        if not overwrite:
            assert_equal(A, A_copy)
            assert_equal(b, b_copy)

    @pytest.mark.skipif(
        np.__version__ < '2', reason="solve chokes on b.ndim == 1 in numpy < 2"
    )
    @pytest.mark.parametrize(
        "assume_a",
        [
            None, "diagonal", "general", "upper triangular", "lower triangular", "pos",
        ]
    )
    def test_vs_np_solve(self, assume_a):
        e = np.eye(2)
        a = np.arange(1, 4*3*2 + 1).reshape((4, 3, 2, 1, 1)) * e

        b = np.ones(2)
        assert_allclose(solve(a, b, assume_a=assume_a), np.linalg.solve(a, b))

        b = np.ones((2, 1))
        assert_allclose(solve(a, b, assume_a=assume_a), np.linalg.solve(a, b))

        b = np.ones((2, 2)) * [1, 2]
        assert_allclose(solve(a, b, assume_a=assume_a), np.linalg.solve(a, b))

    def test_pos_lower(self):
        # regression test for
        # https://github.com/scipy/scipy/pull/23071#issuecomment-3085826112
        rng = np.random.default_rng(0)
        a = rng.normal(size=(4, 4))
        a = np.tril(np.matmul(a, np.conj(a.T)))  # lower triangle of hermitian array
        b = rng.normal(size=(4, 2))
        out = solve(a, b, assume_a='pos', lower=True)

        aa = a + a.T - np.diag(np.diag(a))   # the full hermitian array
        result_np = np.linalg.solve(aa, b)
        assert_allclose(out, result_np, atol=1e-15)

        # repeat with uplo='U'
        out = solve(a.T, b, assume_a='pos', lower=False)
        assert_allclose(out, result_np, atol=1e-15)

    def test_pos_fails_sym_complex(self):
        # regression test for the `solve` analog of gh-24359
        # the matrix is 1) symmetric not hermitian, and 2) not positive definite:
        a = np.asarray([[ 182.56985285-64.28859483j, -177.24879835+11.0780499j ],
                        [-177.24879835+11.0780499j ,  177.24879835-11.0780499j ]])
        b = np.eye(2)

        ainv = solve(a, b)
        assert_allclose(ainv @ a, np.eye(2), atol=1e-14)

        ainv_sym = solve(a, b, assume_a="sym")
        assert_allclose(ainv_sym, ainv, atol=1e-14)

        # Specifying assume_a="pos" disables the structure detection, and directly
        # calls LAPACK routines zportf and zpotri.
        # Since zportf(a) does not error out, neither does solve.
        ainv_chol = solve(a, b, assume_a="pos")
        assert not np.allclose(ainv, ainv_chol, atol=1e-14)

        # Setting assume_a="pos" with a non-pos def matrix returned nonsense.
        # This is at least consistent with inv.
        ainv_inv = inv(a, assume_a="pos")
        assert_allclose(ainv_chol, ainv_inv, atol=1e-14)

    def test_readonly(self):
        a = np.eye(3)
        a.flags.writeable = False
        b = np.ones(3)
        x = solve(a, b)
        assert_allclose(x, b, atol=1e-14)

    @parametrize_overwrite_arg
    def test_batch_negative_stride(self, overwrite_kw):
        a = np.arange(3*8).reshape(2, 3, 2, 2)
        a = a[:, ::-1, :, :]
        b = np.ones(2)
        x = solve(a, b, **overwrite_kw)
        assert x.shape == a.shape[:-1]
        assert_allclose(a @ x[..., None] - b, 0, atol=1e-14)

        # use b with a negative stride now
        b = np.ones((2, 4))[:, ::-1]
        x = solve(a, b, **overwrite_kw)
        assert x.shape == a.shape[:-1] + (b.shape[-1],)
        assert_allclose(a @ x - b, 0, atol=1e-14)

    @parametrize_overwrite_arg
    def test_core_negative_stride(self, overwrite_kw):
        a = np.arange(3*8).reshape(2, 3, 2, 2)
        a = a[:, :, ::-1, :]
        b = np.ones(2)
        x = solve(a, b, **overwrite_kw)

        assert x.shape == a.shape[:-1]
        assert_allclose(a @ x[..., None] - b, 0, atol=1e-14)

        # use b with a negative stride now
        b = np.ones((2, 4))[::-1, :]
        x = solve(a, b, **overwrite_kw)
        assert x.shape == a.shape[:-1] + (b.shape[-1],)
        assert_allclose(a @ x - b, 0, atol=1e-14)

    @parametrize_overwrite_arg
    def test_core_non_contiguous(self, overwrite_kw):
        a = np.arange(3*8*2).reshape(2, 3, 2, 4)
        a = a[..., ::2]
        b = np.ones(2)
        x = solve(a, b, **overwrite_kw)
        assert x.shape == a.shape[:-1]
        assert_allclose(a @ x[..., None] - b, 0, atol=1e-14)

        # use strided b now
        b = np.ones(4)[::2]
        x = solve(a, b, **overwrite_kw)
        assert x.shape == a.shape[:-1]
        assert_allclose(a @ x[..., None] - b, 0, atol=1e-14)

    @parametrize_overwrite_arg
    def test_batch_non_contiguous(self, overwrite_kw):
        a = np.arange(3*8*2).reshape(2, 6, 2, 2)
        a = a[:, ::2, ...]
        b = np.ones(2)
        x = solve(a, b, **overwrite_kw)
        assert x.shape == a.shape[:-1]
        assert_allclose(a @ x[..., None] - b, 0, atol=1e-14)

        # use strided b now
        b = np.ones((2, 6))[:, ::2]
        x = solve(a, b, **overwrite_kw)
        assert x.shape == a.shape[:-1] + (b.shape[-1],)
        assert_allclose(a @ x - b, 0, atol=1e-14)

    @parametrize_overwrite_arg
    def test_batch_weird_strides(self, overwrite_kw):
        a = np.arange(3*8*2).reshape(2, 3, 2, 2, 2)
        a = a.transpose(1, 3, 4, 0, 2)

        b = np.ones(2)
        x = solve(a, b, **overwrite_kw)
        assert x.shape == a.shape[:-1]
        assert_allclose(a @ x[..., None] - b, 0, atol=1e-14)

    @parametrize_overwrite_arg
    @parametrize_overwrite_b_arg
    @pytest.mark.parametrize('a_dtype', [int, float])
    @pytest.mark.parametrize('a_order', ['C', 'F'])
    @pytest.mark.parametrize('b_dtype', [int, float])
    @pytest.mark.parametrize('b_order', ['C', 'F'])
    @pytest.mark.parametrize('b_ndim', [1, 2])    # XXX ndim > 2
    @pytest.mark.parametrize('transposed', [True, False])
    def test_overwrite_args(
        self, overwrite_kw, overwrite_b_kw, a_dtype, a_order,
        b_dtype, b_order, b_ndim, transposed
    ):
        n = 3
        a = np.arange(1, n**2 + 1).reshape(n, n) + np.eye(n)
        a = a.astype(a_dtype, order=a_order)

        b = np.arange(n)
        if b_ndim > 1:
            b = np.stack([b*j for j in range(b_ndim)]).T
        b = b.astype(b_dtype, order=b_order)

        a_ref = a.copy()
        b_ref = b.copy()

        # solve and check that the solution is correct for all parameters
        x = solve(a, b, **overwrite_kw, **overwrite_b_kw, transposed=transposed)
        a_or_aT = a_ref.T if transposed else a_ref
        assert_allclose(a_or_aT @ x, b_ref, atol=1e-14)

        # now check that it worked in-place where expected
        overwrite_a = overwrite_kw.get('overwrite_a', False)
        a_inplace = overwrite_a and (a.dtype != int) and a.flags['F_CONTIGUOUS']

        overwrite_b = overwrite_b_kw.get('overwrite_b', False)
        b_inplace = overwrite_b and (b.dtype != int) and b.flags['F_CONTIGUOUS']

        assert np.shares_memory(x, b) == b_inplace

        assert (b == b_ref).all() != b_inplace
        assert (a == a_ref).all() != a_inplace

    def test_posdef_not_posdef(self):
        # the `b` matrix is invertible but not positive definite
        a = np.arange(9).reshape(3, 3)
        A = a + a.T + np.eye(3)
        b = np.ones(3)

        # cholesky solver fails, and the routine falls back to the general inverse
        x0 = solve(A, b)
        assert_allclose(A @ x0, b, atol=1e-14)

        # but it does not fall back if `assume_a` is given
        with assert_raises(LinAlgError):
            solve(A, b, assume_a='pos')

    def test_diagonal(self):
        a = np.stack([np.triu(np.ones((3, 3))), np.diag(np.arange(1, 4))])
        b = np.ones(3)
        x = solve(a, b)

        # basic diagonal solve
        assert_allclose(x[1, ...], 1 / np.arange(1, 4), atol=1e-14)

        # ill-conditioned inputs warn
        a = np.asarray([[1e30, 0], [0, 1]])
        b = np.ones(2)
        with pytest.warns(LinAlgWarning):
            solve(a, b, assume_a="diagonal")

        # singular input raises
        a = np.asarray([[0, 0], [0, 1]])
        b = np.ones(2)
        with pytest.raises(LinAlgError):
            solve(a, b, assume_a="diagonal")

    def test_tridiagonal(self):
        n = 4
        a = -2*np.diag(np.ones(n)) + np.diag(np.ones(3), 1) + np.diag(np.ones(3), -1)
        a = np.stack([np.triu(np.ones((n, n))), a])
        b = np.ones(4)
        x = solve(a, b)

        # basic tridiag solve
        assert_allclose(x[1, ...], np.asarray([-2., -3., -3., -2.]), atol=1e-15)

        # ill-conditioned inputs warn
        a[1, 0, 0] = 1e20
        with pytest.warns(LinAlgWarning):
            solve(a, b, assume_a="tridiagonal")

        # singular inputss raise
        a[1, 0, 0] = a[1, 0, 1] = 0
        with pytest.raises(LinAlgError):
            solve(a, b, assume_a="tridiagonal")


class TestSolveTriangular:

    def test_simple(self):
        """
        solve_triangular on a simple 2x2 matrix.
        """
        A = array([[1, 0], [1, 2]])
        b = [1, 1]
        sol = solve_triangular(A, b, lower=True)
        assert_array_almost_equal(sol, [1, 0])

        # check that it works also for non-contiguous matrices
        sol = solve_triangular(A.T, b, lower=False)
        assert_array_almost_equal(sol, [.5, .5])

        # and that it gives the same result as trans=1
        sol = solve_triangular(A, b, lower=True, trans=1)
        assert_array_almost_equal(sol, [.5, .5])

        b = identity(2)
        sol = solve_triangular(A, b, lower=True, trans=1)
        assert_array_almost_equal(sol, [[1., -.5], [0, 0.5]])

    def test_simple_complex(self):
        """
        solve_triangular on a simple 2x2 complex matrix
        """
        A = array([[1+1j, 0], [1j, 2]])
        b = identity(2)
        sol = solve_triangular(A, b, lower=True, trans=1)
        assert_array_almost_equal(sol, [[.5-.5j, -.25-.25j], [0, 0.5]])

        # check other option combinations with complex rhs
        b = np.diag([1+1j, 1+2j])
        sol = solve_triangular(A, b, lower=True, trans=0)
        assert_array_almost_equal(sol, [[1, 0], [-0.5j, 0.5+1j]])

        sol = solve_triangular(A, b, lower=True, trans=1)
        assert_array_almost_equal(sol, [[1, 0.25-0.75j], [0, 0.5+1j]])

        sol = solve_triangular(A, b, lower=True, trans=2)
        assert_array_almost_equal(sol, [[1j, -0.75-0.25j], [0, 0.5+1j]])

        sol = solve_triangular(A.T, b, lower=False, trans=0)
        assert_array_almost_equal(sol, [[1, 0.25-0.75j], [0, 0.5+1j]])

        sol = solve_triangular(A.T, b, lower=False, trans=1)
        assert_array_almost_equal(sol, [[1, 0], [-0.5j, 0.5+1j]])

        sol = solve_triangular(A.T, b, lower=False, trans=2)
        assert_array_almost_equal(sol, [[1j, 0], [-0.5, 0.5+1j]])

    def test_check_finite(self):
        """
        solve_triangular on a simple 2x2 matrix.
        """
        A = array([[1, 0], [1, 2]])
        b = [1, 1]
        sol = solve_triangular(A, b, lower=True, check_finite=False)
        assert_array_almost_equal(sol, [1, 0])

    @pytest.mark.parametrize('dt_a', [int, float, np.float32, complex, np.complex64])
    @pytest.mark.parametrize('dt_b', [int, float, np.float32, complex, np.complex64])
    def test_empty(self, dt_a, dt_b):
        a = np.empty((0, 0), dtype=dt_a)
        b = np.empty(0, dtype=dt_b)
        x = solve_triangular(a, b)

        assert x.size == 0
        dt_nonempty = solve_triangular(
            np.eye(2, dtype=dt_a), np.ones(2, dtype=dt_b)
        ).dtype
        assert x.dtype == dt_nonempty

    def test_empty_rhs(self):
        a = np.eye(2)
        b = [[], []]
        x = solve_triangular(a, b)
        assert_(x.size == 0, 'Returned array is not empty')
        assert_(x.shape == (2, 0), 'Returned empty array shape is wrong')


class TestInv:
    def test_simple(self):
        a = [[1, 2], [3, 4]]
        a_inv = inv(a)
        assert_array_almost_equal(dot(a, a_inv), np.eye(2))
        a = [[1, 2, 3], [4, 5, 6], [7, 8, 10]]
        a_inv = inv(a)
        assert_array_almost_equal(dot(a, a_inv), np.eye(3))

    def test_random(self):
        rng = np.random.default_rng(1234)
        n = 20
        for i in range(4):
            a = rng.random([n, n])
            for i in range(n):
                a[i, i] = 20*(.1+a[i, i])
            a_inv = inv(a)
            assert_array_almost_equal(dot(a, a_inv),
                                      identity(n))

    def test_simple_complex(self):
        a = [[1, 2], [3, 4j]]
        a_inv = inv(a)
        assert_array_almost_equal(dot(a, a_inv), [[1, 0], [0, 1]])

    def test_random_complex(self):
        rng = np.random.default_rng(1234)
        n = 20
        for i in range(4):
            a = rng.random([n, n])+2j*rng.random([n, n])
            for i in range(n):
                a[i, i] = 20*(.1+a[i, i])
            a_inv = inv(a)
            assert_array_almost_equal(dot(a, a_inv),
                                      identity(n))

    def test_check_finite(self):
        a = [[1, 2], [3, 4]]
        a_inv = inv(a, check_finite=False)
        assert_array_almost_equal(dot(a, a_inv), [[1, 0], [0, 1]])

    @pytest.mark.parametrize('dt', [int, float, np.float32, complex, np.complex64])
    def test_empty(self, dt):
        a = np.empty((0, 0), dtype=dt)
        a_inv = inv(a)
        assert a_inv.size == 0
        assert a_inv.dtype == inv(np.eye(2, dtype=dt)).dtype

        a = np.ones((3, 0, 2, 2), dtype=dt)
        a_inv = inv(a)
        assert a_inv.shape == (3, 0, 2, 2)

        a = np.ones((3, 1, 0, 0), dtype=dt)
        a_inv = inv(a)
        assert a_inv.shape == (3, 1, 0, 0)

    @parametrize_overwrite_arg
    def test_overwrite_a(self, overwrite_kw):
        n = 3
        a0 = np.arange(1, n**2 + 1).reshape(n, n) + np.eye(n)

        # int arrays are copied internally
        a = a0.copy()
        a_inv = inv(a, **overwrite_kw)
        assert_allclose(a_inv @ a, np.eye(n), atol=1e-14)
        assert_equal(a, a0)
        assert not np.shares_memory(a, a_inv)

        # float C ordered arrays are copied, too
        a = a0.copy().astype(float)
        a_inv = inv(a, **overwrite_kw)
        assert_allclose(a_inv @ a0, np.eye(n), atol=1e-14)
        assert_equal(a, a0)
        assert not np.shares_memory(a, a_inv)

        # 2D F-ordered arrays of LAPACK-compatible dtypes: inv works inplace.
        # IOW, the output is always the inverse, and the original input may be
        # destroyed, depending on the `overwrite_a` kwarg value
        a = a0.astype(float).copy(order='F')
        a_inv = inv(a, **overwrite_kw)
        assert_allclose(a_inv @ a0, np.eye(n), atol=1e-14)

        overwrite_a = overwrite_kw.get("overwrite_a", False)
        assert (a == a0).all() != overwrite_a
        assert np.shares_memory(a, a_inv) == overwrite_a

    @pytest.mark.parametrize(
        "dtyp", [np.float16, np.float32, np.longdouble, np.clongdouble]
    )
    def test_dtypes(self, dtyp):
        # backwards compat: inv(float16)->float32 ; inv(clongdouble)->complex128 etc
        a = np.arange(4).reshape(2, 2).astype(dtyp)

        a_inv = inv(a)
        assert_allclose(a @ a_inv, np.eye(a.shape[0]), atol=100*np.finfo(a.dtype).eps)

        dt_map = {
            'e': 'f',  # float16 -> float32
            'f': 'f',
            'g': 'd',  # longdouble -> float64
            'G': 'D'   # clongdouble -> complex128
        }
        assert a_inv.dtype.char == dt_map[a.dtype.char]

    def test_readonly(self):
        a = np.eye(3)
        a.flags.writeable = False

        a_inv = inv(a)
        assert_allclose(a_inv, a, atol=1e-14)

    @pytest.mark.parametrize('dt', [int, float, np.float32, complex, np.complex64])
    def test_batch_core_1x1(self, dt):
        a = np.arange(3*2, dtype=dt).reshape(3, 2, 1, 1) + 1
        a_inv = inv(a)
        assert a_inv.shape == a.shape
        assert_allclose(a @ a_inv, 1.)

    @parametrize_overwrite_arg
    def test_batch_zero_stride(self, overwrite_kw):
        a = np.arange(3*2*2, dtype=float).reshape(3, 2, 2)
        aa = a[None, ...]
        a_inv = inv(aa, **overwrite_kw)
        assert a_inv.shape == aa.shape
        assert_allclose(aa @ a_inv, np.broadcast_to(np.eye(2), aa.shape), atol=2e-14)

        aa = a[:, None, ...]
        a_inv = inv(aa, **overwrite_kw)
        assert a_inv.shape == aa.shape
        assert_allclose(aa @ a_inv, np.broadcast_to(np.eye(2), aa.shape), atol=2e-14)

    @parametrize_overwrite_arg
    def test_batch_negative_stride(self, overwrite_kw):
        a = np.arange(3*8).reshape(2, 3, 2, 2)
        a = a[:, ::-1, :, :]
        a_inv = inv(a, **overwrite_kw)
        assert a_inv.shape == a.shape
        assert_allclose(a @ a_inv, np.broadcast_to(np.eye(2), a.shape), atol=5e-14)

    @parametrize_overwrite_arg
    def test_core_negative_stride(self, overwrite_kw):
        a = np.arange(3*8).reshape(2, 3, 2, 2)
        a = a[:, :, ::-1, :]
        a_inv = inv(a, **overwrite_kw)
        assert a_inv.shape == a.shape
        assert_allclose(a @ a_inv, np.broadcast_to(np.eye(2), a.shape), atol=5e-14)

    @parametrize_overwrite_arg
    def test_core_non_contiguous(self, overwrite_kw):
        a = np.arange(3*8*2).reshape(2, 3, 2, 4)
        a = a[..., ::2]
        a_inv = inv(a, **overwrite_kw)
        assert a_inv.shape == (2, 3, 2, 2)
        assert_allclose(a @ a_inv, np.broadcast_to(np.eye(2), a.shape), atol=5e-14)

    @parametrize_overwrite_arg
    def test_batch_non_contiguous(self, overwrite_kw):
        a = np.arange(3*8*2).reshape(2, 6, 2, 2)
        a = a[:, ::2, ...]
        a_inv = inv(a, **overwrite_kw)
        assert a_inv.shape == (2, 3, 2, 2)
        assert_allclose(a @ a_inv, np.broadcast_to(np.eye(2), a.shape), atol=2e-13)

    @parametrize_overwrite_arg
    def test_singular(self, overwrite_kw):
        # 2D case: A singular matrix: raise

        with assert_raises(LinAlgError):
            inv(np.ones((2, 2)))

        # batched case: If all slices are singlar, raise
        with assert_raises(LinAlgError):
            inv(np.ones((3, 2, 2)))

        # XXX: shall we make this behavior configurable somehow?
        # A "keep-going" option would be this:
        # if some of the slices are singular and some are not,
        # - singular slices are filled with nans
        # - non-singular slices are inverted
        # - there is no error
        a = np.stack((np.ones((2, 2), dtype=complex), np.arange(4).reshape(2, 2)))
        with assert_raises(LinAlgError):
            inv(a)

        # this would be true for a "keep-going" option
        # assert np.isnan(a_inv[0, ...]).all()
        # assert_allclose(a_inv[1, ...] @ a[1, ...],  np.eye(2), atol=1e-14)

    def test_ill_cond(self):
        a = np.diag([1., 1e-20])
        with pytest.warns(LinAlgWarning):
            inv(a)

        a2 = np.stack([np.diag([1., 1e-20]), np.diag([1, 1]), np.diag([1, 1e-20])])
        with pytest.warns(LinAlgWarning):
            inv(a2)

    def test_wrong_assume_a(self):
        with assert_raises(KeyError):
            inv(np.eye(2), assume_a="kaboom")

    def test_posdef(self):
        x = np.arange(25, dtype=float).reshape(5, 5)
        y = x + x.T
        y += 21*np.eye(5)

        y_inv0 = inv(y)
        y_inv1 = inv(y, assume_a="pos")

        assert_allclose(y_inv1, y_inv0, atol=1e-15)

        # check that the lower triangle is not referenced for `lower=False`
        mask = np.where(1 - np.tri(*y.shape, -1) == 0, np.nan, 1)
        y_inv2 = inv(y*mask, check_finite=False, assume_a="pos", lower=False)
        assert_allclose(y_inv2, y_inv0, atol=1e-15)

        # repeat with the upper triangle
        y_inv3 = inv(y*mask.T, check_finite=False, assume_a="pos", lower=True)
        assert_allclose(y_inv3, y_inv0, atol=1e-15)

    @pytest.mark.parametrize('complex_', [False, True])
    def test_posdef_not_posdef(self, complex_):
        # the `b` matrix is invertible but not pos definite: test the "sym" fallback
        a = np.arange(9).reshape(3, 3)
        b = a + a.T + np.eye(3)
        if complex_:
            b = b + 1j*b

        # cholesky solver fails, and the routine falls back to the symmetric inverse
        b_inv0 = inv(b)
        assert_allclose(b_inv0 @ b, np.eye(3), atol=3e-15)

        # but it does not fall back if `assume_a` is given
        with assert_raises(LinAlgError):
            inv(b, assume_a='pos')

        # test posdef fallback to the hermitian solver, too
        if complex_:
            a = np.arange(9).reshape(3, 3)
            a = a + 1j*a
            b = a + a.T.conj() + np.eye(3)
            assert_allclose(inv(b) @ b, np.eye(3), atol=3e-15)

    def test_pos_fails_sym_complex(self):
        # regression test for gh-24359
        # the matrix is 1) symmetric not hermitian, and 2) not positive definite:
        a = np.asarray([[ 182.56985285-64.28859483j, -177.24879835+11.0780499j ],
                        [-177.24879835+11.0780499j ,  177.24879835-11.0780499j ]])

        ainv = inv(a)
        assert_allclose(ainv @ a, np.eye(2), atol=1e-14)

        ainv_sym = inv(a, assume_a="sym")
        assert_allclose(ainv_sym, ainv, atol=1e-14)

        # Specifying assume_a="pos" disables the structure detection, and directly
        # calls LAPACK routines zportf and zpotri.
        # Since zportf(a) does not error out, neither does inv
        ainv_chol = inv(a, assume_a="pos")
        assert not np.allclose(ainv, ainv_chol, atol=1e-14)

        # Setting assume_a="pos" with a non-pos def matrix returned nonsense.
        # This is at least consistent with solve.
        ainv_slv = solve(a, np.eye(2), assume_a="pos")
        assert_allclose(ainv_chol, ainv_slv, atol=1e-14)

        # Repeat it for bunch of simple cases to cover more branches
        # Real symmetric, positive definite
        a = np.eye(4) + np.ones(4)
        res = inv(a)
        assert_allclose(res @ a, np.eye(4), atol=1e-14)

        # Real symmetric, NOT positive definite
        a = -np.eye(4) + np.ones(4)
        res = inv(a)
        assert_allclose(res @ a, np.eye(4), atol=1e-14)

        # Real, not symmetric
        a = -np.eye(4) + np.ones(4)
        a[0, -1] = 2.
        res = inv(a)
        assert_allclose(res @ a, np.eye(4), atol=1e-14)

        # | Test                                  | is_symm | is_herm | pos def |
        # |---------------------------------------|---------|---------|---------|
        # | Complex, both sym+herm, pos def       |    1    |    1    |   yes   |
        # | Complex, symmetric only               |    1    |    0    |    -    |
        # | Complex, both sym+herm, NOT pos def   |    1    |    1    |   no    |
        # | Complex, neither                      |    0    |    0    |    -    |
        # | Complex, hermitian only, pos def      |    0    |    1    |   yes   |
        # | Complex, hermitian only, NOT pos def  |    0    |    1    |   no    |

        # Complex, both symmetric and hermitian, positive definite
        a = (np.eye(4) + np.ones(4)).astype(np.complex128)
        res = inv(a)
        assert_allclose(res @ a, np.eye(4), atol=1e-14)

        # Complex, symmetric only (not hermitian)
        a = (np.eye(4)*1.0j + np.ones(4)).astype(np.complex128)
        res = inv(a)
        assert_allclose(res @ a, np.eye(4), atol=1e-14)

        # Complex, both symmetric and hermitian, NOT positive definite
        a = (-np.eye(4) + np.ones(4)).astype(np.complex128)
        res = inv(a)
        assert_allclose(res @ a, np.eye(4), atol=1e-14)

        # Complex, neither symmetric nor hermitian
        a = (-np.eye(4) + np.ones(4)).astype(np.complex128)
        a[0, -1] = 2.
        res = inv(a)
        assert_allclose(res @ a, np.eye(4), atol=1e-14)

        # Complex, hermitian only, positive definite
        a = np.array([[2, 1+1j], [1-1j, 2]], dtype=np.complex128)
        res = inv(a)
        assert_allclose(res @ a, np.eye(2), atol=1e-14)

        # Complex, hermitian only, NOT positive definite
        a = np.array([[-1, 1+1j], [1-1j, -1]], dtype=np.complex128)
        res = inv(a)
        assert_allclose(res @ a, np.eye(2), atol=1e-14)

    @pytest.mark.parametrize('complex_', [False, True])
    @pytest.mark.parametrize('sym_herm', ['sym', 'her'])
    def test_sym_her(self, complex_, sym_herm):
        # test "sym" and "her" modes
        a = np.arange(9).reshape(3, 3)
        if complex_:
            a = a + 1j*a

        if sym_herm == "sym":
            b = a + a.T
        else:   # sym_herm == "herm":
            b = a + a.T.conj()

        b = b + np.eye(3)

        b_inv0 = np.linalg.inv(b)
        assert_allclose(b_inv0 @ b, np.eye(3), atol=1e-14)

        b_inv1 = inv(b, assume_a=sym_herm)
        assert_allclose(b_inv0, b_inv1, atol=1e-15)

        # check that the "other" triangle is not referenced
        mask = np.where(1 - np.tri(*a.shape, -1) == 0, np.nan, 1)
        b_inv2 = inv(b*mask, check_finite=False, assume_a=sym_herm, lower=False)
        assert_allclose(b_inv2, b_inv0, atol=1e-15)

        # repeat with the upper triangle
        b_inv3 = inv(b*mask.T, check_finite=False, assume_a=sym_herm, lower=True)
        assert_allclose(b_inv3, b_inv0, atol=1e-15)

    def test_triangular_1(self):
        x = np.arange(25, dtype=float).reshape(5, 5)
        y = x + x.T
        y += 21*np.eye(5)
        y_inv0 = inv(y, assume_a='upper triangular')

        # check that upper triangular differs from posdef
        y_inv_posdef = inv(y, assume_a='pos')
        assert not np.allclose(y_inv0, y_inv_posdef)

    def test_triangular_2(self):
        y = np.ones(25, dtype=float).reshape(5, 5)

        y_inv_0_u = inv(np.triu(y))
        assert_allclose(y_inv_0_u @ np.triu(y), np.eye(5), atol=1e-15)

        y_inv_1_u = inv(y, assume_a='upper triangular')
        assert_allclose(y_inv_1_u @ np.triu(y), np.eye(5), atol=1e-15)

        # check that the lower triangle is not referenced for "upper triangular"
        mask = np.where(1 - np.tri(*y.shape, -1) == 0, np.nan, 1)
        y_inv_2_u = inv(y*mask, check_finite=False, assume_a='upper triangular')
        assert_allclose(y_inv_2_u @ np.triu(y), np.eye(5), atol=1e-15)

        # repeat for the lower traingular matrix
        y_inv_0_l = inv(np.tril(y))
        assert_allclose(y_inv_0_l @ np.tril(y), np.eye(5), atol=1e-15)

        y_inv_1_l = inv(y, assume_a='lower triangular')
        assert_allclose(y_inv_1_l @ np.tril(y), np.eye(5), atol=1e-15)

        # check that the lower triangle is not referenced for "lower triangular"
        mask = np.where(1 - np.tri(*y.shape, -1) == 0, np.nan, 1)
        y_inv_2_l = inv(y*mask.T, check_finite=False, assume_a='lower triangular')
        assert_allclose(y_inv_2_l @ np.tril(y), np.eye(5), atol=1e-15)

    def test_diagonal(self):
        a = np.stack([np.triu(np.ones((3, 3))), np.diag(np.arange(1, 4))])
        inv_a = inv(a)

        # basic diagonal invert
        assert_allclose(inv_a[1], np.diag(1 / np.arange(1, 4)), atol=1e-14)

        # ill-conditioned inputs warn
        a = np.asarray([[1e30, 0], [0, 1]])
        with pytest.warns(LinAlgWarning):
            inv(a, assume_a="diagonal")

        # singular input raises
        a = np.asarray([[0, 0], [0, 1]])
        with pytest.raises(LinAlgError):
            inv(a, assume_a="diagonal")


class TestDet:
    def test_1x1_all_singleton_dims(self):
        a = np.array([[1]])
        deta = det(a)
        assert deta.dtype.char == 'd'
        assert np.isscalar(deta)
        assert deta == 1.
        a = np.array([[[[1]]]], dtype='f')
        deta = det(a)
        assert deta.dtype.char == 'd'
        assert deta.shape == (1, 1)
        assert_equal(deta, [[1.0]])
        a = np.array([[[1 + 3.j]]], dtype=np.complex64)
        deta = det(a)
        assert deta.dtype.char == 'D'
        assert deta.shape == (1,)
        assert_equal(deta, [1.+3.j])

    def test_1by1_stacked_input_output(self):
        rng = np.random.default_rng(1680305949878959)
        a = rng.random([4, 5, 1, 1], dtype=np.float32)
        deta = det(a)
        assert deta.dtype.char == 'd'
        assert deta.shape == (4, 5)
        assert_allclose(deta, np.squeeze(a))

        a = rng.random([4, 5, 1, 1], dtype=np.float32)*np.complex64(1.j)
        deta = det(a)
        assert deta.dtype.char == 'D'
        assert deta.shape == (4, 5)
        assert_allclose(deta, np.squeeze(a))

    @pytest.mark.parametrize('shape', [[2, 2], [20, 20], [3, 2, 20, 20]])
    def test_simple_det_shapes_real_complex(self, shape):
        rng = np.random.default_rng(1680305949878959)
        a = rng.uniform(-1., 1., size=shape)
        d1, d2 = det(a), np.linalg.det(a)
        assert_allclose(d1, d2)

        b = rng.uniform(-1., 1., size=shape)*1j
        b += rng.uniform(-0.5, 0.5, size=shape)
        d3, d4 = det(b), np.linalg.det(b)
        assert_allclose(d3, d4)

    def test_for_known_det_values(self):
        # Hadamard8
        a = np.array([[1, 1, 1, 1, 1, 1, 1, 1],
                      [1, -1, 1, -1, 1, -1, 1, -1],
                      [1, 1, -1, -1, 1, 1, -1, -1],
                      [1, -1, -1, 1, 1, -1, -1, 1],
                      [1, 1, 1, 1, -1, -1, -1, -1],
                      [1, -1, 1, -1, -1, 1, -1, 1],
                      [1, 1, -1, -1, -1, -1, 1, 1],
                      [1, -1, -1, 1, -1, 1, 1, -1]])
        assert_allclose(det(a), 4096.)

        # consecutive number array always singular
        assert_allclose(det(np.arange(25).reshape(5, 5)), 0.)

        # simple anti-diagonal block array
        # Upper right has det (-2+1j) and lower right has (-2-1j)
        # det(a) = - (-2+1j) (-2-1j) = 5.
        a = np.array([[0.+0.j, 0.+0.j, 0.-1.j, 1.-1.j],
                      [0.+0.j, 0.+0.j, 1.+0.j, 0.-1.j],
                      [0.+1.j, 1.+1.j, 0.+0.j, 0.+0.j],
                      [1.+0.j, 0.+1.j, 0.+0.j, 0.+0.j]], dtype=np.complex64)
        assert_allclose(det(a), 5.+0.j)

        # Fiedler companion complexified
        # >>> a = scipy.linalg.fiedler_companion(np.arange(1, 10))
        a = np.array([[-2., -3., 1., 0., 0., 0., 0., 0.],
                      [1., 0., 0., 0., 0., 0., 0., 0.],
                      [0., -4., 0., -5., 1., 0., 0., 0.],
                      [0., 1., 0., 0., 0., 0., 0., 0.],
                      [0., 0., 0., -6., 0., -7., 1., 0.],
                      [0., 0., 0., 1., 0., 0., 0., 0.],
                      [0., 0., 0., 0., 0., -8., 0., -9.],
                      [0., 0., 0., 0., 0., 1., 0., 0.]])*1.j
        assert_allclose(det(a), 9.)

    # g and G dtypes are handled differently in windows and other platforms
    @pytest.mark.parametrize('typ', [x for x in np.typecodes['All'][:20]
                                     if x not in 'gG'])
    def test_sample_compatible_dtype_input(self, typ):
        rng = np.random.default_rng(1680305949878959)
        n = 4
        a = rng.random([n, n]).astype(typ)  # value is not important
        assert isinstance(det(a), (np.float64 | np.complex128))

    def test_incompatible_dtype_input(self):
        # Double backslashes needed for escaping pytest regex.
        msg = 'cannot be cast to float\\(32, 64\\)'

        for c, t in zip('SUO', ['bytes8', 'str32', 'object']):
            with assert_raises(TypeError, match=msg):
                det(np.array([['a', 'b']]*2, dtype=c))
        with assert_raises(TypeError, match=msg):
            det(np.array([[b'a', b'b']]*2, dtype='V'))
        with assert_raises(TypeError, match=msg):
            det(np.array([[100, 200]]*2, dtype='datetime64[s]'))
        with assert_raises(TypeError, match=msg):
            det(np.array([[100, 200]]*2, dtype='timedelta64[s]'))

    def test_empty_edge_cases(self):
        assert_allclose(det(np.empty([0, 0])), 1.)
        assert_allclose(det(np.empty([0, 0, 0])), np.array([]))
        assert_allclose(det(np.empty([3, 0, 0])), np.array([1., 1., 1.]))
        with assert_raises(ValueError, match='Last 2 dimensions'):
            det(np.empty([0, 0, 3]))
        with assert_raises(ValueError, match='at least two-dimensional'):
            det(np.array([]))
        with assert_raises(ValueError, match='Last 2 dimensions'):
            det(np.array([[]]))
        with assert_raises(ValueError, match='Last 2 dimensions'):
            det(np.array([[[]]]))

    @pytest.mark.parametrize('dt', [int, float, np.float32, complex, np.complex64])
    def test_empty_dtype(self, dt):
        a = np.empty((0, 0), dtype=dt)
        d = det(a)
        assert d.shape == ()
        assert d.dtype == det(np.eye(2, dtype=dt)).dtype

        a = np.empty((3, 0, 0), dtype=dt)
        d = det(a)
        assert d.shape == (3,)
        assert d.dtype == det(np.zeros((3, 1, 1), dtype=dt)).dtype

    def test_overwrite_a(self):
        # If all conditions are met then input should be overwritten;
        #   - dtype is one of 'fdFD'
        #   - C-contiguous
        #   - writeable
        a = np.arange(9).reshape(3, 3).astype(np.float32)
        ac = a.copy()
        deta = det(ac, overwrite_a=True)
        assert_allclose(deta, 0.)
        assert not (a == ac).all()

    def test_readonly_array(self):
        a = np.array([[2., 0., 1.], [5., 3., -1.], [1., 1., 1.]])
        a.setflags(write=False)
        # overwrite_a will be overridden
        assert_allclose(det(a, overwrite_a=True), 10.)

    def test_simple_check_finite(self):
        a = [[1, 2], [3, np.inf]]
        with assert_raises(ValueError, match='array must not contain'):
            det(a)


def direct_lstsq(a, b, cmplx=0):
    at = transpose(a)
    if cmplx:
        at = conjugate(at)
    a1 = dot(at, a)
    b1 = dot(at, b)
    return solve(a1, b1)


class TestLstsq:
    lapack_drivers = ('gelsd', 'gelss', 'gelsy', None)

    def test_simple_exact(self):
        for dtype in REAL_DTYPES:
            a = np.array([[1, 20], [-30, 4]], dtype=dtype)
            for lapack_driver in TestLstsq.lapack_drivers:
                for overwrite in (True, False):
                    for bt in (((1, 0), (0, 1)), (1, 0),
                               ((2, 1), (-30, 4))):
                        # Store values in case they are overwritten
                        # later
                        a1 = a.copy()
                        b = np.array(bt, dtype=dtype)
                        b1 = b.copy()
                        out = lstsq(a1, b1,
                                    lapack_driver=lapack_driver,
                                    overwrite_a=overwrite,
                                    overwrite_b=overwrite)
                        x = out[0]
                        r = out[2]
                        assert_(r == 2,
                                f'expected efficient rank 2, got {r}')
                        assert_allclose(dot(a, x), b,
                                        atol=25 * _eps_cast(a1.dtype),
                                        rtol=25 * _eps_cast(a1.dtype),
                                        err_msg=f"driver: {lapack_driver}")

    def test_simple_overdet(self):
        for dtype in REAL_DTYPES:
            a = np.array([[1, 2], [4, 5], [3, 4]], dtype=dtype)
            b = np.array([1, 2, 3], dtype=dtype)
            for lapack_driver in TestLstsq.lapack_drivers:
                for overwrite in (True, False):
                    # Store values in case they are overwritten later
                    a1 = a.copy()
                    b1 = b.copy()
                    out = lstsq(a1, b1, lapack_driver=lapack_driver,
                                overwrite_a=overwrite,
                                overwrite_b=overwrite)
                    x = out[0]
                    if lapack_driver == 'gelsy':
                        residuals = np.sum((b - a.dot(x))**2)
                    else:
                        residuals = out[1]
                    r = out[2]
                    assert_(r == 2, f'expected efficient rank 2, got {r}')
                    assert_allclose(abs((dot(a, x) - b)**2).sum(axis=0),
                                    residuals,
                                    rtol=25 * _eps_cast(a1.dtype),
                                    atol=25 * _eps_cast(a1.dtype),
                                    err_msg=f"driver: {lapack_driver}")
                    assert_allclose(x, (-0.428571428571429, 0.85714285714285),
                                    rtol=25 * _eps_cast(a1.dtype),
                                    atol=25 * _eps_cast(a1.dtype),
                                    err_msg=f"driver: {lapack_driver}")

    def test_simple_overdet_complex(self):
        for dtype in COMPLEX_DTYPES:
            a = np.array([[1+2j, 2], [4, 5], [3, 4]], dtype=dtype)
            b = np.array([1, 2+4j, 3], dtype=dtype)
            for lapack_driver in TestLstsq.lapack_drivers:
                for overwrite in (True, False):
                    # Store values in case they are overwritten later
                    a1 = a.copy()
                    b1 = b.copy()
                    out = lstsq(a1, b1, lapack_driver=lapack_driver,
                                overwrite_a=overwrite,
                                overwrite_b=overwrite)

                    x = out[0]
                    if lapack_driver == 'gelsy':
                        res = b - a.dot(x)
                        residuals = np.sum(res * res.conj())
                    else:
                        residuals = out[1]
                    r = out[2]
                    assert_(r == 2, f'expected efficient rank 2, got {r}')
                    assert_allclose(abs((dot(a, x) - b)**2).sum(axis=0),
                                    residuals,
                                    rtol=25 * _eps_cast(a1.dtype),
                                    atol=25 * _eps_cast(a1.dtype),
                                    err_msg=f"driver: {lapack_driver}")
                    assert_allclose(
                                x, (-0.4831460674157303 + 0.258426966292135j,
                                    0.921348314606741 + 0.292134831460674j),
                                rtol=25 * _eps_cast(a1.dtype),
                                atol=25 * _eps_cast(a1.dtype),
                                err_msg=f"driver: {lapack_driver}")

    def test_simple_underdet(self):
        for dtype in REAL_DTYPES:
            a = np.array([[1, 2, 3], [4, 5, 6]], dtype=dtype)
            b = np.array([1, 2], dtype=dtype)
            for lapack_driver in TestLstsq.lapack_drivers:
                for overwrite in (True, False):
                    # Store values in case they are overwritten later
                    a1 = a.copy()
                    b1 = b.copy()
                    out = lstsq(a1, b1, lapack_driver=lapack_driver,
                                overwrite_a=overwrite,
                                overwrite_b=overwrite)

                    x = out[0]
                    r = out[2]
                    assert_(r == 2, f'expected efficient rank 2, got {r}')
                    assert_allclose(x, (-0.055555555555555, 0.111111111111111,
                                        0.277777777777777),
                                    rtol=25 * _eps_cast(a1.dtype),
                                    atol=25 * _eps_cast(a1.dtype),
                                    err_msg=f"driver: {lapack_driver}")

    @pytest.mark.parametrize("dtype", REAL_DTYPES)
    @pytest.mark.parametrize("n", (20, 200))
    @pytest.mark.parametrize("lapack_driver", lapack_drivers)
    @pytest.mark.parametrize("overwrite", (True, False))
    def test_random_exact(self, dtype, n, lapack_driver, overwrite):
        rng = np.random.RandomState(1234)

        a = np.asarray(rng.random([n, n]), dtype=dtype)
        for i in range(n):
            a[i, i] = 20 * (0.1 + a[i, i])
        for i in range(4):
            b = np.asarray(rng.random([n, 3]), dtype=dtype)
            # Store values in case they are overwritten later
            a1 = a.copy()
            b1 = b.copy()
            out = lstsq(a1, b1,
                        lapack_driver=lapack_driver,
                        overwrite_a=overwrite,
                        overwrite_b=overwrite)
            x = out[0]
            r = out[2]
            assert_(r == n, f'expected efficient rank {n}, '
                    f'got {r}')
            if dtype is np.float32:
                assert_allclose(
                          dot(a, x), b,
                          rtol=500 * _eps_cast(a1.dtype),
                          atol=500 * _eps_cast(a1.dtype),
                          err_msg=f"driver: {lapack_driver}")
            else:
                assert_allclose(
                          dot(a, x), b,
                          rtol=1000 * _eps_cast(a1.dtype),
                          atol=1000 * _eps_cast(a1.dtype),
                          err_msg=f"driver: {lapack_driver}")

    @pytest.mark.skipif(IS_MUSL, reason="may segfault on Alpine, see gh-17630")
    @pytest.mark.parametrize("dtype", COMPLEX_DTYPES)
    @pytest.mark.parametrize("n", (20, 200))
    @pytest.mark.parametrize("lapack_driver", lapack_drivers)
    @pytest.mark.parametrize("overwrite", (True, False))
    def test_random_complex_exact(self, dtype, n, lapack_driver, overwrite):
        rng = np.random.RandomState(1234)

        a = np.asarray(rng.random([n, n]) + 1j*rng.random([n, n]),
                       dtype=dtype)
        for i in range(n):
            a[i, i] = 20 * (0.1 + a[i, i])
        for i in range(2):
            b = np.asarray(rng.random([n, 3]), dtype=dtype)
            # Store values in case they are overwritten later
            a1 = a.copy()
            b1 = b.copy()
            out = lstsq(a1, b1, lapack_driver=lapack_driver,
                        overwrite_a=overwrite,
                        overwrite_b=overwrite)
            x = out[0]
            r = out[2]
            assert_(r == n, f'expected efficient rank {n}, '
                    f'got {r}')
            if dtype is np.complex64:
                assert_allclose(
                          dot(a, x), b,
                          rtol=400 * _eps_cast(a1.dtype),
                          atol=400 * _eps_cast(a1.dtype),
                          err_msg=f"driver: {lapack_driver}")
            else:
                assert_allclose(
                          dot(a, x), b,
                          rtol=1000 * _eps_cast(a1.dtype),
                          atol=1000 * _eps_cast(a1.dtype),
                          err_msg=f"driver: {lapack_driver}")

    def test_random_overdet(self):
        rng = np.random.RandomState(1234)
        for dtype in REAL_DTYPES:
            for (n, m) in ((20, 15), (200, 2)):
                for lapack_driver in TestLstsq.lapack_drivers:
                    for overwrite in (True, False):
                        a = np.asarray(rng.random([n, m]), dtype=dtype)
                        for i in range(m):
                            a[i, i] = 20 * (0.1 + a[i, i])
                        for i in range(4):
                            b = np.asarray(rng.random([n, 3]), dtype=dtype)
                            # Store values in case they are overwritten later
                            a1 = a.copy()
                            b1 = b.copy()
                            out = lstsq(a1, b1,
                                        lapack_driver=lapack_driver,
                                        overwrite_a=overwrite,
                                        overwrite_b=overwrite)
                            x = out[0]
                            r = out[2]
                            assert_(r == m, f'expected efficient rank {m}, '
                                    f'got {r}')
                            assert_allclose(
                                          x, direct_lstsq(a, b, cmplx=0),
                                          rtol=25 * _eps_cast(a1.dtype),
                                          atol=25 * _eps_cast(a1.dtype),
                                          err_msg=f"driver: {lapack_driver}")

    def test_random_complex_overdet(self):
        rng = np.random.RandomState(1234)
        for dtype in COMPLEX_DTYPES:
            for (n, m) in ((20, 15), (200, 2)):
                for lapack_driver in TestLstsq.lapack_drivers:
                    for overwrite in (True, False):
                        a = np.asarray(rng.random([n, m]) + 1j*rng.random([n, m]),
                                       dtype=dtype)
                        for i in range(m):
                            a[i, i] = 20 * (0.1 + a[i, i])
                        for i in range(2):
                            b = np.asarray(rng.random([n, 3]), dtype=dtype)
                            # Store values in case they are overwritten
                            # later
                            a1 = a.copy()
                            b1 = b.copy()
                            out = lstsq(a1, b1,
                                        lapack_driver=lapack_driver,
                                        overwrite_a=overwrite,
                                        overwrite_b=overwrite)
                            x = out[0]
                            r = out[2]
                            assert_(r == m, f'expected efficient rank {m}, '
                                    f'got {r}')
                            assert_allclose(
                                      x, direct_lstsq(a, b, cmplx=1),
                                      rtol=25 * _eps_cast(a1.dtype),
                                      atol=25 * _eps_cast(a1.dtype),
                                      err_msg=f"driver: {lapack_driver}")

    def test_check_finite(self):
        with warnings.catch_warnings():
            # On (some) OSX this tests triggers a warning (gh-7538)
            warnings.filterwarnings("ignore",
                                    "internal gelsd driver lwork query error,.*"
                                    "Falling back to 'gelss' driver.", RuntimeWarning)

        at = np.array(((1, 20), (-30, 4)))
        for dtype, bt, lapack_driver, overwrite, check_finite in \
            itertools.product(REAL_DTYPES,
                              (((1, 0), (0, 1)), (1, 0), ((2, 1), (-30, 4))),
                              TestLstsq.lapack_drivers,
                              (True, False),
                              (True, False)):

            a = at.astype(dtype)
            b = np.array(bt, dtype=dtype)
            # Store values in case they are overwritten
            # later
            a1 = a.copy()
            b1 = b.copy()
            out = lstsq(a1, b1, lapack_driver=lapack_driver,
                        check_finite=check_finite, overwrite_a=overwrite,
                        overwrite_b=overwrite)
            x = out[0]
            r = out[2]
            assert_(r == 2, f'expected efficient rank 2, got {r}')
            assert_allclose(dot(a, x), b,
                            rtol=25 * _eps_cast(a.dtype),
                            atol=25 * _eps_cast(a.dtype),
                            err_msg=f"driver: {lapack_driver}")

    def test_empty(self):
        for a_shape, b_shape in (((0, 2), (0,)),
                                 ((0, 4), (0, 2)),
                                 ((4, 0), (4,)),
                                 ((4, 0), (4, 2))):
            b = np.ones(b_shape)
            x, residues, rank, s = lstsq(np.zeros(a_shape), b)
            assert_equal(x, np.zeros((a_shape[1],) + b_shape[1:]))
            residues_should_be = (np.empty((0,)) if a_shape[1]
                                  else np.linalg.norm(b, axis=0)**2)
            assert_equal(residues, residues_should_be)
            assert_(rank == 0, 'expected rank 0')
            assert_equal(s, np.empty((0,)))

    @pytest.mark.parametrize('dt_a', [int, float, np.float32, complex, np.complex64])
    @pytest.mark.parametrize('dt_b', [int, float, np.float32, complex, np.complex64])
    def test_empty_dtype(self, dt_a, dt_b):
        a = np.empty((0, 0), dtype=dt_a)
        b = np.empty(0, dtype=dt_b)
        x, residues, rank, s = lstsq(a, b)

        assert x.size == 0
        dt_nonempty = lstsq(np.eye(2, dtype=dt_a), np.ones(2, dtype=dt_b))[0].dtype
        assert x.dtype == dt_nonempty


class TestPinv:
    def test_simple_real(self):
        a = array([[1, 2, 3], [4, 5, 6], [7, 8, 10]], dtype=float)
        a_pinv = pinv(a)
        assert_array_almost_equal(dot(a, a_pinv), np.eye(3))

    def test_simple_complex(self):
        a = (array([[1, 2, 3], [4, 5, 6], [7, 8, 10]],
             dtype=float) + 1j * array([[10, 8, 7], [6, 5, 4], [3, 2, 1]],
                                       dtype=float))
        a_pinv = pinv(a)
        assert_array_almost_equal(dot(a, a_pinv), np.eye(3))

    def test_simple_singular(self):
        a = array([[1, 2, 3], [4, 5, 6], [7, 8, 9]], dtype=float)
        a_pinv = pinv(a)
        expected = array([[-6.38888889e-01, -1.66666667e-01, 3.05555556e-01],
                          [-5.55555556e-02, 1.30136518e-16, 5.55555556e-02],
                          [5.27777778e-01, 1.66666667e-01, -1.94444444e-01]])
        assert_array_almost_equal(a_pinv, expected)

    def test_simple_cols(self):
        a = array([[1, 2, 3], [4, 5, 6]], dtype=float)
        a_pinv = pinv(a)
        expected = array([[-0.94444444, 0.44444444],
                          [-0.11111111, 0.11111111],
                          [0.72222222, -0.22222222]])
        assert_array_almost_equal(a_pinv, expected)

    def test_simple_rows(self):
        a = array([[1, 2], [3, 4], [5, 6]], dtype=float)
        a_pinv = pinv(a)
        expected = array([[-1.33333333, -0.33333333, 0.66666667],
                          [1.08333333, 0.33333333, -0.41666667]])
        assert_array_almost_equal(a_pinv, expected)

    def test_check_finite(self):
        a = array([[1, 2, 3], [4, 5, 6.], [7, 8, 10]])
        a_pinv = pinv(a, check_finite=False)
        assert_array_almost_equal(dot(a, a_pinv), np.eye(3))

    def test_native_list_argument(self):
        a = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
        a_pinv = pinv(a)
        expected = array([[-6.38888889e-01, -1.66666667e-01, 3.05555556e-01],
                          [-5.55555556e-02, 1.30136518e-16, 5.55555556e-02],
                          [5.27777778e-01, 1.66666667e-01, -1.94444444e-01]])
        assert_array_almost_equal(a_pinv, expected)

    def test_atol_rtol(self):
        rng = np.random.default_rng(1234)
        n = 12
        # get a random ortho matrix for shuffling
        q, _ = qr(rng.random((n, n)))
        a_m = np.arange(35.0).reshape(7, 5)
        a = a_m.copy()
        a[0, 0] = 0.001
        atol = 1e-5
        rtol = 0.05
        # svds of a_m is ~ [116.906, 4.234, tiny, tiny, tiny]
        # svds of a is ~ [116.906, 4.234, 4.62959e-04, tiny, tiny]
        # Just abs cutoff such that we arrive at a_modified
        a_p = pinv(a_m, atol=atol, rtol=0.)
        adiff1 = a @ a_p @ a - a
        adiff2 = a_m @ a_p @ a_m - a_m
        # Now adiff1 should be around atol value while adiff2 should be
        # relatively tiny
        assert_allclose(np.linalg.norm(adiff1), 5e-4, atol=5.e-4)
        assert_allclose(np.linalg.norm(adiff2), 5e-14, atol=5.e-14)

        # Now do the same but remove another sv ~4.234 via rtol
        a_p = pinv(a_m, atol=atol, rtol=rtol)
        adiff1 = a @ a_p @ a - a
        adiff2 = a_m @ a_p @ a_m - a_m
        assert_allclose(np.linalg.norm(adiff1), 4.233, rtol=0.01)
        assert_allclose(np.linalg.norm(adiff2), 4.233, rtol=0.01)

    @pytest.mark.parametrize('dt', [float, np.float32, complex, np.complex64])
    def test_empty(self, dt):
        a = np.empty((0, 0), dtype=dt)
        a_pinv = pinv(a)
        assert a_pinv.size == 0
        assert a_pinv.dtype == pinv(np.eye(2, dtype=dt)).dtype


class TestPinvSymmetric:
    def test_simple_real(self):
        a = array([[1, 2, 3], [4, 5, 6], [7, 8, 10]], dtype=float)
        a = np.dot(a, a.T)
        a_pinv = pinvh(a)
        assert_array_almost_equal(np.dot(a, a_pinv), np.eye(3))

    def test_nonpositive(self):
        a = array([[1, 2, 3], [4, 5, 6], [7, 8, 9]], dtype=float)
        a = np.dot(a, a.T)
        u, s, vt = np.linalg.svd(a)
        s[0] *= -1
        a = np.dot(u * s, vt)  # a is now symmetric non-positive and singular
        a_pinv = pinv(a)
        a_pinvh = pinvh(a)
        assert_array_almost_equal(a_pinv, a_pinvh)

    def test_simple_complex(self):
        a = (array([[1, 2, 3], [4, 5, 6], [7, 8, 10]],
             dtype=float) + 1j * array([[10, 8, 7], [6, 5, 4], [3, 2, 1]],
                                       dtype=float))
        a = np.dot(a, a.conj().T)
        a_pinv = pinvh(a)
        assert_array_almost_equal(np.dot(a, a_pinv), np.eye(3))

    def test_native_list_argument(self):
        a = array([[1, 2, 3], [4, 5, 6], [7, 8, 10]], dtype=float)
        a = np.dot(a, a.T)
        a_pinv = pinvh(a.tolist())
        assert_array_almost_equal(np.dot(a, a_pinv), np.eye(3))

    def test_zero_eigenvalue(self):
        # https://github.com/scipy/scipy/issues/12515
        # the SYEVR eigh driver may give the zero eigenvalue > eps
        a = np.array([[1, -1, 0], [-1, 2, -1], [0, -1, 1]])
        p = pinvh(a)
        assert_allclose(p @ a @ p, p, atol=1e-15)
        assert_allclose(a @ p @ a, a, atol=1e-15)

    def test_atol_rtol(self):
        rng = np.random.default_rng(1234)
        n = 12
        # get a random ortho matrix for shuffling
        q, _ = qr(rng.random((n, n)))
        a = np.diag([4, 3, 2, 1, 0.99e-4, 0.99e-5] + [0.99e-6]*(n-6))
        a = q.T @ a @ q
        a_m = np.diag([4, 3, 2, 1, 0.99e-4, 0.] + [0.]*(n-6))
        a_m = q.T @ a_m @ q
        atol = 1e-5
        rtol = (4.01e-4 - 4e-5)/4
        # Just abs cutoff such that we arrive at a_modified
        a_p = pinvh(a, atol=atol, rtol=0.)
        adiff1 = a @ a_p @ a - a
        adiff2 = a_m @ a_p @ a_m - a_m
        # Now adiff1 should dance around atol value since truncation
        # while adiff2 should be relatively tiny
        assert_allclose(norm(adiff1), atol, rtol=0.1)
        assert_allclose(norm(adiff2), 1e-12, atol=1e-11)

        # Now do the same but through rtol cancelling atol value
        a_p = pinvh(a, atol=atol, rtol=rtol)
        adiff1 = a @ a_p @ a - a
        adiff2 = a_m @ a_p @ a_m - a_m
        # adiff1 and adiff2 should be elevated to ~1e-4 due to mismatch
        assert_allclose(norm(adiff1), 1e-4, rtol=0.1)
        assert_allclose(norm(adiff2), 1e-4, rtol=0.1)

    @pytest.mark.parametrize('dt', [float, np.float32, complex, np.complex64])
    def test_empty(self, dt):
        a = np.empty((0, 0), dtype=dt)
        a_pinv = pinvh(a)
        assert a_pinv.size == 0
        assert a_pinv.dtype == pinv(np.eye(2, dtype=dt)).dtype


@pytest.mark.parametrize('scale', (1e-20, 1., 1e20))
@pytest.mark.parametrize('pinv_', (pinv, pinvh))
def test_auto_rcond(scale, pinv_):
    x = np.array([[1, 0], [0, 1e-10]]) * scale
    expected = np.diag(1. / np.diag(x))
    x_inv = pinv_(x)
    assert_allclose(x_inv, expected)


class TestVectorNorms:

    def test_types(self):
        for dtype in np.typecodes['AllFloat']:
            x = np.array([1, 2, 3], dtype=dtype)
            tol = max(1e-15, np.finfo(dtype).eps.real * 20)
            assert_allclose(norm(x), np.sqrt(14), rtol=tol)
            assert_allclose(norm(x, 2), np.sqrt(14), rtol=tol)

        for dtype in np.typecodes['Complex']:
            x = np.array([1j, 2j, 3j], dtype=dtype)
            tol = max(1e-15, np.finfo(dtype).eps.real * 20)
            assert_allclose(norm(x), np.sqrt(14), rtol=tol)
            assert_allclose(norm(x, 2), np.sqrt(14), rtol=tol)

    def test_overflow(self):
        # unlike numpy's norm, this one is
        # safer on overflow
        a = array([1e20], dtype=float32)
        assert_almost_equal(norm(a), a)

    def test_stable(self):
        # more stable than numpy's norm
        a = array([1e4] + [1]*10000, dtype=float32)
        try:
            # snrm in double precision; we obtain the same as for float64
            # -- large atol needed due to varying blas implementations
            assert_allclose(norm(a) - 1e4, 0.5, atol=1e-2)
        except AssertionError:
            # snrm implemented in single precision, == np.linalg.norm result
            msg = ": Result should equal either 0.0 or 0.5 (depending on " \
                  "implementation of snrm2)."
            assert_almost_equal(norm(a) - 1e4, 0.0, err_msg=msg)

    def test_zero_norm(self):
        assert_equal(norm([1, 0, 3], 0), 2)
        assert_equal(norm([1, 2, 3], 0), 3)

    def test_axis_kwd(self):
        a = np.array([[[2, 1], [3, 4]]] * 2, 'd')
        assert_allclose(norm(a, axis=1), [[3.60555128, 4.12310563]] * 2)
        assert_allclose(norm(a, 1, axis=1), [[5.] * 2] * 2)

    def test_keepdims_kwd(self):
        a = np.array([[[2, 1], [3, 4]]] * 2, 'd')
        b = norm(a, axis=1, keepdims=True)
        assert_allclose(b, [[[3.60555128, 4.12310563]]] * 2)
        assert_(b.shape == (2, 1, 2))
        assert_allclose(norm(a, 1, axis=2, keepdims=True), [[[3.], [7.]]] * 2)

    @pytest.mark.skipif(not HAS_ILP64, reason="64-bit BLAS required")
    def test_large_vector(self):
        check_free_memory(free_mb=17000)
        x = np.zeros([2**31], dtype=np.float64)
        x[-1] = 1
        res = norm(x)
        del x
        assert_allclose(res, 1.0)


class TestMatrixNorms:

    def test_matrix_norms(self):
        # Not all of these are matrix norms in the most technical sense.
        rng = np.random.default_rng(1234)
        for n, m in (1, 1), (1, 3), (3, 1), (4, 4), (4, 5), (5, 4):
            for t in np.float32, np.float64, np.complex64, np.complex128, np.int64:
                A = 10 * rng.standard_normal((n, m)).astype(t)
                if np.issubdtype(A.dtype, np.complexfloating):
                    A += 10j * rng.standard_normal((n, m))
                    t_high = np.complex128
                else:
                    t_high = np.float64
                for order in (None, 'fro', 1, -1, 2, -2, np.inf, -np.inf):
                    actual = norm(A, ord=order)
                    desired = np.linalg.norm(A, ord=order)
                    # SciPy may return higher precision matrix norms.
                    # This is a consequence of using LAPACK.
                    if not np.allclose(actual, desired):
                        desired = np.linalg.norm(A.astype(t_high), ord=order)
                        assert_allclose(actual, desired)

    def test_axis_kwd(self):
        a = np.array([[[2, 1], [3, 4]]] * 2, 'd')
        b = norm(a, ord=np.inf, axis=(1, 0))
        c = norm(np.swapaxes(a, 0, 1), ord=np.inf, axis=(0, 1))
        d = norm(a, ord=1, axis=(0, 1))
        assert_allclose(b, c)
        assert_allclose(c, d)
        assert_allclose(b, d)
        assert_(b.shape == c.shape == d.shape)
        b = norm(a, ord=1, axis=(1, 0))
        c = norm(np.swapaxes(a, 0, 1), ord=1, axis=(0, 1))
        d = norm(a, ord=np.inf, axis=(0, 1))
        assert_allclose(b, c)
        assert_allclose(c, d)
        assert_allclose(b, d)
        assert_(b.shape == c.shape == d.shape)

    def test_keepdims_kwd(self):
        a = np.arange(120, dtype='d').reshape(2, 3, 4, 5)
        b = norm(a, ord=np.inf, axis=(1, 0), keepdims=True)
        c = norm(a, ord=1, axis=(0, 1), keepdims=True)
        assert_allclose(b, c)
        assert_(b.shape == c.shape)

    def test_empty(self):
        a = np.empty((0, 0))
        assert_allclose(norm(a), 0.)
        assert_allclose(norm(a, axis=0), np.zeros((0,)))
        assert_allclose(norm(a, keepdims=True), np.zeros((1, 1)))

        a = np.empty((0, 3))
        assert_allclose(norm(a), 0.)
        assert_allclose(norm(a, axis=0), np.zeros((3,)))
        assert_allclose(norm(a, keepdims=True), np.zeros((1, 1)))


class TestOverwrite:
    def test_solve(self):
        assert_no_overwrite(solve, [(3, 3), (3,)])

    def test_solve_triangular(self):
        assert_no_overwrite(solve_triangular, [(3, 3), (3,)])

    def test_solve_banded(self):
        assert_no_overwrite(lambda ab, b: solve_banded((2, 1), ab, b),
                            [(4, 6), (6,)])

    def test_solveh_banded(self):
        assert_no_overwrite(solveh_banded, [(2, 6), (6,)])

    def test_inv(self):
        assert_no_overwrite(inv, [(3, 3)])

    def test_det(self):
        assert_no_overwrite(det, [(3, 3)])

    def test_lstsq(self):
        assert_no_overwrite(lstsq, [(3, 2), (3,)])

    def test_pinv(self):
        assert_no_overwrite(pinv, [(3, 3)])

    def test_pinvh(self):
        assert_no_overwrite(pinvh, [(3, 3)])


class TestSolveCirculant:

    def test_basic1(self):
        c = np.array([1, 2, 3, 5])
        b = np.array([1, -1, 1, 0])
        x = solve_circulant(c, b)
        y = solve(circulant(c), b)
        assert_allclose(x, y)

    def test_basic2(self):
        # b is a 2-d matrix.
        c = np.array([1, 2, -3, -5])
        b = np.arange(12).reshape(4, 3)
        x = solve_circulant(c, b)
        y = solve(circulant(c), b)
        assert_allclose(x, y)

    def test_basic3(self):
        # b is a 3-d matrix.
        c = np.array([1, 2, -3, -5])
        b = np.arange(24).reshape(4, 3, 2)
        x = solve_circulant(c, b)
        y = solve(circulant(c), b.reshape(4, -1)).reshape(b.shape)
        assert_allclose(x, y)

    def test_complex(self):
        # Complex b and c
        c = np.array([1+2j, -3, 4j, 5])
        b = np.arange(8).reshape(4, 2) + 0.5j
        x = solve_circulant(c, b)
        y = solve(circulant(c), b)
        assert_allclose(x, y)

    def test_random_b_and_c(self):
        # Random b and c
        rng = np.random.RandomState(54321)
        c = rng.standard_normal(50)
        b = rng.standard_normal(50)
        x = solve_circulant(c, b)
        y = solve(circulant(c), b)
        assert_allclose(x, y)

    def test_singular(self):
        # c gives a singular circulant matrix.
        c = np.array([1, 1, 0, 0])
        b = np.array([1, 2, 3, 4])
        x = solve_circulant(c, b, singular='lstsq')
        y, res, rnk, s = lstsq(circulant(c), b)
        assert_allclose(x, y)
        assert_raises(LinAlgError, solve_circulant, x, y)

    def test_axis_args(self):
        # Test use of caxis, baxis and outaxis.

        # c has shape (2, 1, 4)
        c = np.array([[[-1, 2.5, 3, 3.5]], [[1, 6, 6, 6.5]]])

        # b has shape (3, 4)
        b = np.array([[0, 0, 1, 1], [1, 1, 0, 0], [1, -1, 0, 0]])

        x = solve_circulant(c, b, baxis=1)
        assert_equal(x.shape, (4, 2, 3))
        expected = np.empty_like(x)
        expected[:, 0, :] = solve(circulant(c[0].ravel()), b.T)
        expected[:, 1, :] = solve(circulant(c[1].ravel()), b.T)
        assert_allclose(x, expected)

        x = solve_circulant(c, b, baxis=1, outaxis=-1)
        assert_equal(x.shape, (2, 3, 4))
        assert_allclose(np.moveaxis(x, -1, 0), expected)

        # np.swapaxes(c, 1, 2) has shape (2, 4, 1); b.T has shape (4, 3).
        x = solve_circulant(np.swapaxes(c, 1, 2), b.T, caxis=1)
        assert_equal(x.shape, (4, 2, 3))
        assert_allclose(x, expected)

    def test_native_list_arguments(self):
        # Same as test_basic1 using python's native list.
        c = [1, 2, 3, 5]
        b = [1, -1, 1, 0]
        x = solve_circulant(c, b)
        y = solve(circulant(c), b)
        assert_allclose(x, y)

    @pytest.mark.parametrize('dt_c', [int, float, np.float32, complex, np.complex64])
    @pytest.mark.parametrize('dt_b', [int, float, np.float32, complex, np.complex64])
    def test_empty(self, dt_c, dt_b):
        c = np.array([], dtype=dt_c)
        b = np.array([], dtype=dt_b)
        x = solve_circulant(c, b)
        assert x.shape == (0,)
        assert x.dtype == solve_circulant(np.arange(3, dtype=dt_c),
                                          np.ones(3, dtype=dt_b)).dtype

        b = np.empty((0, 0), dtype=dt_b)
        x1 = solve_circulant(c, b)
        assert x1.shape == (0, 0)
        assert x1.dtype == x.dtype


class TestMatrix_Balance:
    @skip_xp_invalid_arg
    def test_string_arg(self):
        assert_raises(ValueError, matrix_balance, 'Some string for fail')

    def test_infnan_arg(self):
        assert_raises(ValueError, matrix_balance,
                      np.array([[1, 2], [3, np.inf]]))
        assert_raises(ValueError, matrix_balance,
                      np.array([[1, 2], [3, np.nan]]))

    def test_scaling(self):
        _, y = matrix_balance(np.array([[1000, 1], [1000, 0]]))
        # Pre/post LAPACK 3.5.0 gives the same result up to an offset
        # since in each case col norm is x1000 greater and
        # 1000 / 32 ~= 1 * 32 hence balanced with 2 ** 5.
        assert_allclose(np.diff(np.log2(np.diag(y))), [5])

    def test_scaling_order(self):
        A = np.array([[1, 0, 1e-4], [1, 1, 1e-2], [1e4, 1e2, 1]])
        x, y = matrix_balance(A)
        assert_allclose(solve(y, A).dot(y), x)

    def test_separate(self):
        _, (y, z) = matrix_balance(np.array([[1000, 1], [1000, 0]]),
                                   separate=1)
        assert_equal(np.diff(np.log2(y)), [5])
        assert_allclose(z, np.arange(2))

    def test_permutation(self):
        A = block_diag(np.ones((2, 2)), np.tril(np.ones((2, 2))),
                       np.ones((3, 3)))
        x, (y, z) = matrix_balance(A, separate=1)
        assert_allclose(y, np.ones_like(y))
        assert_allclose(z, np.array([0, 1, 6, 5, 4, 3, 2]))

    def test_perm_and_scaling(self):
        # Matrix with its diagonal removed
        cases = (  # Case 0
                 np.array([[0., 0., 0., 0., 0.000002],
                           [0., 0., 0., 0., 0.],
                           [2., 2., 0., 0., 0.],
                           [2., 2., 0., 0., 0.],
                           [0., 0., 0.000002, 0., 0.]]),
                 #  Case 1 user reported GH-7258
                 np.array([[-0.5, 0., 0., 0.],
                           [0., -1., 0., 0.],
                           [1., 0., -0.5, 0.],
                           [0., 1., 0., -1.]]),
                 #  Case 2 user reported GH-7258
                 np.array([[-3., 0., 1., 0.],
                           [-1., -1., -0., 1.],
                           [-3., -0., -0., 0.],
                           [-1., -0., 1., -1.]])
                 )

        for A in cases:
            x, y = matrix_balance(A)
            x, (s, p) = matrix_balance(A, separate=1)
            ip = np.empty_like(p)
            ip[p] = np.arange(A.shape[0])
            assert_allclose(y, np.diag(s)[ip, :])
            assert_allclose(solve(y, A).dot(y), x)

    @pytest.mark.parametrize('dt', [int, float, np.float32, complex, np.complex64])
    def test_empty(self, dt):
        a = np.empty((0, 0), dtype=dt)
        b, t = matrix_balance(a)

        assert b.size == 0
        assert t.size == 0

        b_n, t_n = matrix_balance(np.eye(2, dtype=dt))
        assert b.dtype == b_n.dtype
        assert t.dtype == t_n.dtype

        b, (scale, perm) = matrix_balance(a, separate=True)
        assert b.size == 0
        assert scale.size == 0
        assert perm.size == 0

        b_n, (scale_n, perm_n) = matrix_balance(a, separate=True)
        assert b.dtype == b_n.dtype
        assert scale.dtype == scale_n.dtype
        assert perm.dtype == perm_n.dtype


class TestDTypes:
    """Check backwards compatibility for dtypes vs scipy 1.16."""

    def get_arr2D(self, tcode):
        # return a valid 2D array for the typecode
        if tcode == 'M':
            return np.eye(2, dtype='datetime64[ms]')
        elif tcode == 'V':
            return np.asarray([[b'a', b'b'], [b'c', b'd']], dtype='V')
        else:
            return np.eye(2, dtype=tcode)

    def get_arr1D(self, tcode):
        # return a valid 1D array for the typecode
        if tcode == 'M':
            return np.ones(2, dtype='datetime64[ms]')
        elif tcode == 'V':
            return np.asarray([b'a', b'b'], dtype='V')
        else:
            return np.ones(2, dtype=tcode)

    @pytest.mark.parametrize("tcode", np.typecodes['All'])
    def test_inv(self, tcode):
        # check backwards compat vs scipy 1.16
        a = self.get_arr2D(tcode)
        if tcode in 'SUVO':
            # raises
            with pytest.raises(ValueError):
                inv(a)
        else:
            # passes
            inv(a)

    @pytest.mark.parametrize("tcode", np.typecodes['All'])
    def test_det(self, tcode):
        a = self.get_arr2D(tcode)

        is_arm = platform.machine() == 'arm64'
        is_windows = os.name == 'nt'

        failing_tcodes = 'SUVOmM'
        if not (is_arm or is_windows):
            failing_tcodes += 'gG'

        if tcode in failing_tcodes:
            # raises
            with pytest.raises(TypeError):
                det(a)
        else:
            # passes
            det(a)

    @pytest.mark.filterwarnings("ignore:Casting complex values")
    @pytest.mark.parametrize("tcode_a", np.typecodes['All'])
    @pytest.mark.parametrize("tcode_b", np.typecodes['All'])
    def test_solve(self, tcode_a, tcode_b):
        a = self.get_arr2D(tcode_a)
        b = self.get_arr1D(tcode_b)

        can_combine = True
        try:
            np.result_type(tcode_a, tcode_b)
        except TypeError:
            can_combine = False

        if not can_combine:
            # np.exceptions.DTypePromotionError subclasses TypeError
            with pytest.raises(TypeError):
                solve(a, b)
        elif tcode_a in 'SUVO' or tcode_b in 'VO':
            with pytest.raises(ValueError):
                solve(a, b)
        else:
            solve(a, b)
