from itertools import product, permutations

import numpy as np
import pytest
from numpy.testing import assert_allclose
from pytest import raises as assert_raises

from scipy.linalg import orthogonal_procrustes
from scipy.sparse._sputils import matrix
from scipy._lib._array_api import make_xp_test_case, xp_assert_close
from scipy.conftest import skip_xp_invalid_arg

def _centered(A, xp):
    mu = xp.mean(A, axis=0)
    return A - mu, mu

@make_xp_test_case(orthogonal_procrustes)
class TestOrthogonalProcrustes:
    def test_orthogonal_procrustes_ndim_too_small(self, xp):
        rng = np.random.RandomState(1234)
        A = xp.asarray(rng.randn(3))
        B = xp.asarray(rng.randn(3))
        assert_raises(ValueError, orthogonal_procrustes, A, B)

    def test_orthogonal_procrustes_shape_mismatch(self, xp):
        rng = np.random.RandomState(1234)
        shapes = ((3, 3), (3, 4), (4, 3), (4, 4))
        for a, b in permutations(shapes, 2):
            A = xp.asarray(rng.randn(*a))
            B = xp.asarray(rng.randn(*b))
            assert_raises(ValueError, orthogonal_procrustes, A, B)

    def test_orthogonal_procrustes_checkfinite_exception(self, xp):
        rng = np.random.RandomState(1234)
        m, n = 2, 3
        A_good = rng.randn(m, n)
        B_good = rng.randn(m, n)
        for bad_value in np.inf, -np.inf, np.nan:
            A_bad = A_good.copy()
            A_bad[1, 2] = bad_value
            B_bad = B_good.copy()
            B_bad[1, 2] = bad_value
            for A, B in ((A_good, B_bad), (A_bad, B_good), (A_bad, B_bad)):
                assert_raises(ValueError, orthogonal_procrustes, xp.asarray(A),
                              xp.asarray(B))

    def test_orthogonal_procrustes_scale_invariance(self, xp):
        rng = np.random.RandomState(1234)
        m, n = 4, 3
        for i in range(3):
            A_orig = xp.asarray(rng.randn(m, n))
            B_orig = xp.asarray(rng.randn(m, n))
            R_orig, s = orthogonal_procrustes(A_orig, B_orig)
            for A_scale in np.square(rng.randn(3)):
                for B_scale in np.square(rng.randn(3)):
                    R, s = orthogonal_procrustes(A_orig * xp.asarray(A_scale),
                                                 B_orig * xp.asarray(B_scale))
                    xp_assert_close(R, R_orig)

    @skip_xp_invalid_arg()
    def test_orthogonal_procrustes_array_conversion(self):
        rng = np.random.RandomState(1234)
        for m, n in ((6, 4), (4, 4), (4, 6)):
            A_arr = rng.randn(m, n)
            B_arr = rng.randn(m, n)
            As = (A_arr, A_arr.tolist(), matrix(A_arr))
            Bs = (B_arr, B_arr.tolist(), matrix(B_arr))
            R_arr, s = orthogonal_procrustes(A_arr, B_arr)
            AR_arr = A_arr.dot(R_arr)
            for A, B in product(As, Bs):
                R, s = orthogonal_procrustes(A, B)
                AR = A_arr.dot(R)
                assert_allclose(AR, AR_arr)

    def test_orthogonal_procrustes(self, xp):
        rng = np.random.RandomState(1234)
        for m, n in ((6, 4), (4, 4), (4, 6)):
            # Sample a random target matrix.
            B = xp.asarray(rng.randn(m, n))
            # Sample a random orthogonal matrix
            # by computing eigh of a sampled symmetric matrix.
            X = xp.asarray(rng.randn(n, n))
            w, V = xp.linalg.eigh(X.T + X)
            xp_assert_close(xp.linalg.inv(V), V.T)
            # Compute a matrix with a known orthogonal transformation that gives B.
            A = B @ V.T
            # Check that an orthogonal transformation from A to B can be recovered.
            R, s = orthogonal_procrustes(A, B)
            xp_assert_close(xp.linalg.inv(R), R.T)
            xp_assert_close(A @ R, B)
            # Create a perturbed input matrix.
            A_perturbed = A + 1e-2 * xp.asarray(rng.randn(m, n))
            # Check that the orthogonal procrustes function can find an orthogonal
            # transformation that is better than the orthogonal transformation
            # computed from the original input matrix.
            R_prime, s = orthogonal_procrustes(A_perturbed, B)
            xp_assert_close(xp.linalg.inv(R_prime), R_prime.T)
            # Compute the naive and optimal transformations of the perturbed input.
            naive_approx = A_perturbed @ R
            optim_approx = A_perturbed @ R_prime
            # Compute the Frobenius norm errors of the matrix approximations.
            naive_approx_error = xp.linalg.matrix_norm(naive_approx - B, ord='fro')
            optim_approx_error = xp.linalg.matrix_norm(optim_approx - B, ord='fro')
            # Check that the orthogonal Procrustes approximation is better.
            assert xp.all(optim_approx_error < naive_approx_error)

    def test_orthogonal_procrustes_exact_example(self, xp):
        # Check a small application.
        # It uses translation, scaling, reflection, and rotation.
        #
        #         |
        #   a  b  |
        #         |
        #   d  c  |        w
        #         |
        # --------+--- x ----- z ---
        #         |
        #         |        y
        #         |
        #
        A_orig = xp.asarray([[-3, 3], [-2, 3], [-2, 2], [-3, 2]], dtype=xp.float64)
        B_orig = xp.asarray([[3, 2], [1, 0], [3, -2], [5, 0]], dtype=xp.float64)
        A, A_mu = _centered(A_orig, xp)
        B, B_mu = _centered(B_orig, xp)
        R, s = orthogonal_procrustes(A, B)
        scale = s / xp.linalg.matrix_norm(A)**2
        B_approx = scale * A @ R + B_mu
        xp_assert_close(B_approx, B_orig, atol=1e-8)

    def test_orthogonal_procrustes_stretched_example(self, xp):
        # Try again with a target with a stretched y axis.
        A_orig = xp.asarray([[-3, 3], [-2, 3], [-2, 2], [-3, 2]], dtype=xp.float64)
        B_orig = xp.asarray([[3, 40], [1, 0], [3, -40], [5, 0]], dtype=xp.float64)
        A, A_mu = _centered(A_orig, xp)
        B, B_mu = _centered(B_orig, xp)
        R, s = orthogonal_procrustes(A, B)
        scale = s / xp.linalg.matrix_norm(A)**2
        B_approx = scale * A @ R + B_mu
        expected = xp.asarray([[3, 21], [-18, 0], [3, -21], [24, 0]], dtype=xp.float64)
        xp_assert_close(B_approx, expected, atol=1e-8)
        # Check disparity symmetry.
        expected_disparity = xp.asarray(0.4501246882793018, dtype=xp.float64)[()]
        AB_disparity = (xp.linalg.matrix_norm(B_approx - B_orig)
                        / xp.linalg.matrix_norm(B))**2
        xp_assert_close(AB_disparity, expected_disparity)
        R, s = orthogonal_procrustes(B, A)
        scale = s / xp.linalg.matrix_norm(B)**2
        A_approx = scale * B @ R + A_mu
        BA_disparity = (xp.linalg.matrix_norm(A_approx - A_orig)
                        / xp.linalg.matrix_norm(A))**2
        xp_assert_close(BA_disparity, expected_disparity)

    def test_orthogonal_procrustes_skbio_example(self, xp):
        # This transformation is also exact.
        # It uses translation, scaling, and reflection.
        #
        #   |
        #   | a
        #   | b
        #   | c d
        # --+---------
        #   |
        #   |       w
        #   |
        #   |       x
        #   |
        #   |   z   y
        #   |
        #
        A_orig = xp.asarray([[4, -2], [4, -4], [4, -6], [2, -6]], dtype=xp.float64)
        B_orig = xp.asarray([[1, 3], [1, 2], [1, 1], [2, 1]], dtype=xp.float64)
        B_standardized = xp.asarray([[-0.13363062, 0.6681531],
                                     [-0.13363062, 0.13363062],
                                     [-0.13363062, -0.40089186],
                                     [0.40089186, -0.40089186]], dtype=xp.float64)
        A, A_mu = _centered(A_orig, xp)
        B, B_mu = _centered(B_orig, xp)
        R, s = orthogonal_procrustes(A, B)
        scale = s / xp.linalg.matrix_norm(A)**2
        B_approx = scale * A @ R + B_mu
        xp_assert_close(B_approx, B_orig)
        xp_assert_close(B / xp.linalg.matrix_norm(B), B_standardized)

    def test_empty(self, xp):
        a = xp.empty((0, 0))
        r, s = orthogonal_procrustes(a, a)
        xp_assert_close(r, xp.empty((0, 0)))

        a = xp.empty((0, 3))
        r, s = orthogonal_procrustes(a, a)
        xp_assert_close(r, xp.eye(3))

    @pytest.mark.parametrize('shape', [(4, 5), (5, 5), (5, 4)])
    def test_unitary(self, shape, xp):
        # gh-12071 added support for unitary matrices; check that it
        # works as intended.
        m, n = shape
        rng = np.random.default_rng(589234981235)
        A = xp.asarray(rng.random(shape) + rng.random(shape) * 1j)
        Q = xp.asarray(rng.random((n, n)) + rng.random((n, n)) * 1j)
        Q, _ = xp.linalg.qr(Q)
        B = A @ Q
        R, scale = orthogonal_procrustes(A, B)
        xp_assert_close(R @ xp.conj(R).T, xp.eye(n, dtype=xp.complex128), atol=1e-14)
        xp_assert_close(A @ Q, B)
        if shape != (4, 5):  # solution is unique
            xp_assert_close(R, Q)
        _, s, _ = xp.linalg.svd(xp.conj(A).T @ B)
        xp_assert_close(scale, xp.sum(s))
