"""Test functions for the sparse.linalg._krylov_funm module."""
from functools import partial

import numpy as np
import pytest
from numpy.testing import (assert_allclose)
import scipy.sparse
import scipy.linalg
from scipy.sparse.linalg import aslinearoperator
from scipy.linalg import (expm, cosm, coshm, sinm, sinhm)

from scipy.sparse.linalg._funm_multiply_krylov import funm_multiply_krylov

REAL_DTYPES = (np.float32, np.float64)
COMPLEX_DTYPES = (np.complex64, np.complex128)
DTYPES = REAL_DTYPES + COMPLEX_DTYPES

# This the phi_1 function from exponential integrators: phi_1 (z) = (e^z - 1) / z
def custom(X):
    return scipy.linalg.solve(X, expm(X) - np.eye(X.shape[0]))

FUNCS = (expm, cosm, coshm, sinm, sinhm, custom)

class TestKrylovFunmv:

    def test_krylov_funm_zero_vector(self):
        n = 20
        A = np.zeros((n, n))
        b = np.zeros(n)
        observed = funm_multiply_krylov(expm, A, b)
        expected = np.zeros(n)
        assert_allclose(observed, expected)

    @pytest.mark.parametrize("f", FUNCS)
    def test_funm_multiply_krylov_nonhermitian_dense(self, f):
        rng = np.random.default_rng(1738151906092735)
        n = 60
        nsamples = 10

        for i in range(nsamples):
            A = rng.standard_normal((n, n))
            b = rng.standard_normal(n)

            fA = f(A)
            expected = fA @ b
            observed = funm_multiply_krylov(f, A, b)
            assert_allclose(observed, expected, rtol = 1E-6, atol = 1E-8)
            observed = funm_multiply_krylov(f, aslinearoperator(A), b)
            assert_allclose(observed, expected, rtol = 1E-6, atol = 1E-8)

    @pytest.mark.parametrize("f", FUNCS)
    def test_funm_multiply_krylov_nonhermitian_sparse(self, f, num_parallel_threads):

        rng = np.random.default_rng(1738151906092735)
        n = 100
        nsamples = 1 + 9 // num_parallel_threads  # Very slow otherwise

        for i in range(nsamples):
            D = scipy.sparse.diags(rng.standard_normal(n))
            A = scipy.sparse.random_array((n, n), density = 0.01, rng = rng) + D
            denseA = A.todense()
            b = rng.standard_normal(n)

            fA = f(denseA)
            expected = fA @ b
            observed = funm_multiply_krylov(f, A, b)
            assert_allclose(observed, expected, rtol = 1E-6, atol = 1E-8)
            observed = funm_multiply_krylov(f, aslinearoperator(A), b)
            assert_allclose(observed, expected, rtol = 1E-6, atol = 1E-8)

    @pytest.mark.parametrize("f", FUNCS)
    def test_funm_multiply_krylov_hermitian_dense(self, f):

        rng = np.random.default_rng(1738151906092735)
        n = 60
        nsamples = 10

        for i in range(nsamples):
            R = np.triu(rng.standard_normal((n, n)))
            A = R.T + R
            b = rng.standard_normal(n)

            fA = f(A)
            expected = fA @ b
            observed = funm_multiply_krylov(f, A, b, assume_a = 'her')
            assert_allclose(observed, expected, rtol = 1E-6, atol = 1E-8)
            observed = funm_multiply_krylov(f, aslinearoperator(A), b, assume_a = 'her')
            assert_allclose(observed, expected, rtol = 1E-6, atol = 1E-8)

    @pytest.mark.parametrize("f", FUNCS)
    def test_funm_multiply_krylov_hermitian_sparse(self, f, num_parallel_threads):

        rng = np.random.default_rng(1738151906092735)
        n = 100
        nsamples = 1 + 9 // num_parallel_threads  # Very slow otherwise

        for i in range(nsamples):
            D = scipy.sparse.diags(rng.standard_normal(n))
            A = scipy.sparse.random_array((n, n), density = 0.01, rng = rng)
            R = scipy.sparse.triu(A)
            A = R + R.T + D
            denseA = A.todense()
            b = rng.standard_normal(n)

            fA = f(denseA)
            expected = fA @ b
            observed = funm_multiply_krylov(f, A, b, assume_a = 'her')
            assert_allclose(observed, expected, rtol = 1E-6, atol = 1E-8)

            observed = funm_multiply_krylov(f, aslinearoperator(A), b, assume_a = 'her')
            assert_allclose(observed, expected, rtol = 1E-6, atol = 1E-8)

    def test_funm_multiply_krylov_breakdown(self):
        rng = np.random.default_rng(1738151906092735)

        # From test_iterative
        A = np.array([[0, 0, 0, 0, 0, 1, -1, -0, -0, -0, -0],
                      [0, 0, 0, 0, 0, 2, -0, -1, -0, -0, -0],
                      [0, 0, 0, 0, 0, 2, -0, -0, -1, -0, -0],
                      [0, 0, 0, 0, 0, 2, -0, -0, -0, -1, -0],
                      [0, 0, 0, 0, 0, 1, -0, -0, -0, -0, -1],
                      [1, 2, 2, 2, 1, 0, -0, -0, -0, -0, -0],
                      [-1, 0, 0, 0, 0, 0, -1, -0, -0, -0, -0],
                      [0, -1, 0, 0, 0, 0, -0, -1, -0, -0, -0],
                      [0, 0, -1, 0, 0, 0, -0, -0, -1, -0, -0],
                      [0, 0, 0, -1, 0, 0, -0, -0, -0, -1, -0],
                      [0, 0, 0, 0, -1, 0, -0, -0, -0, -0, -1]], dtype = float)
        b = rng.standard_normal(A.shape[0])

        fA = expm(A)
        expected = fA @ b
        observed = funm_multiply_krylov(expm, A, b, restart_every_m = 40)
        assert_allclose(observed, expected)

    def test_funm_multiply_krylov_invalid_input(self):
            A = np.array([[1, 2], [3, 4]])  # Non-hermitian matrix
            b = np.array([1.0, 2.0])  # Ensure 'b' is a 1D array of floats

            # Test for invalid 'b' (not 1D)
            b_invalid = np.array([[1.0], [2.0]])  # 2D array
            with pytest.raises(ValueError, 
                    match="argument 'b' must be a 1D array."): 
                funm_multiply_krylov(np.exp, A, b_invalid)

            # Test for invalid restart parameter
            with pytest.raises(ValueError, 
                    match="argument 'restart_every_m' must be positive."): 
                funm_multiply_krylov(np.exp, A, b, restart_every_m=0)

            # Test for invalid max_restarts
            with pytest.raises(ValueError, 
                    match="argument 'max_restarts' must be positive."): 
                funm_multiply_krylov(np.exp, A, b, max_restarts=0)

            # Test for invalid 'assume_a' string
            with pytest.raises(ValueError, 
                    match="is not a recognized matrix structure"): 
                funm_multiply_krylov(np.exp, A, b, assume_a='invalid')


@pytest.mark.parametrize("dtype_a", DTYPES)
@pytest.mark.parametrize("dtype_b", DTYPES)
def test_funm_multiply_krylov_types(dtype_a, dtype_b):
    assert_allclose_ = (partial(assert_allclose, rtol = 1.8e-3, atol = 1e-5)
                        if {dtype_a, dtype_b} else assert_allclose)

    rng = np.random.default_rng(1738151906092735)
    n = 50

    if dtype_a in REAL_DTYPES:
        A = rng.random([n, n]).astype(dtype_a)
    else:
        A = (rng.random([n, n]) + 1j * rng.random([n, n])).astype(dtype_a)

    if dtype_b in REAL_DTYPES:
        b = (2 * rng.random(n)).astype(dtype_b)
    else:
        b = (rng.random(n) + 1j * rng.random(n)).astype(dtype_b)

        expA = expm(A)
        expected = expA @ b
        observed = funm_multiply_krylov(expm, A, b)
        assert_allclose_(observed, expected)
        observed = funm_multiply_krylov(expm, aslinearoperator(A), b)
        assert_allclose_(observed, expected)
