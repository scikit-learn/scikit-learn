#******************************************************************************
#   Copyright (C) 2013 Kenneth L. Ho
#   Redistribution and use in source and binary forms, with or without
#   modification, are permitted provided that the following conditions are met:
#
#   Redistributions of source code must retain the above copyright notice, this
#   list of conditions and the following disclaimer. Redistributions in binary
#   form must reproduce the above copyright notice, this list of conditions and
#   the following disclaimer in the documentation and/or other materials
#   provided with the distribution.
#
#   None of the names of the copyright holders may be used to endorse or
#   promote products derived from this software without specific prior written
#   permission.
#
#   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
#   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
#   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
#   ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
#   LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
#   CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
#   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
#   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
#   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
#   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
#   POSSIBILITY OF SUCH DAMAGE.
#******************************************************************************

import scipy.linalg.interpolative as pymatrixid
import numpy as np
from scipy.linalg import hilbert, svdvals, norm
from scipy.sparse.linalg import aslinearoperator
from scipy.linalg.interpolative import interp_decomp
import time
import itertools

from numpy.testing import assert_, assert_allclose
from pytest import raises as assert_raises


def _debug_print(s):
    if 0:
        print(s)


class TestInterpolativeDecomposition(object):
    def test_id(self):
        for dtype in [np.float64, np.complex128]:
            self.check_id(dtype)

    def check_id(self, dtype):
        # Test ID routines on a Hilbert matrix.

        # set parameters
        n = 300
        eps = 1e-12

        # construct Hilbert matrix
        A = hilbert(n).astype(dtype)
        if np.issubdtype(dtype, np.complexfloating):
            A = A * (1 + 1j)
        L = aslinearoperator(A)

        # find rank
        S = np.linalg.svd(A, compute_uv=False)
        try:
            rank = np.nonzero(S < eps)[0][0]
        except IndexError:
            rank = n

        # print input summary
        _debug_print("Hilbert matrix dimension:        %8i" % n)
        _debug_print("Working precision:               %8.2e" % eps)
        _debug_print("Rank to working precision:       %8i" % rank)

        # set print format
        fmt = "%8.2e (s) / %5s"

        # test real ID routines
        _debug_print("-----------------------------------------")
        _debug_print("Real ID routines")
        _debug_print("-----------------------------------------")

        # fixed precision
        _debug_print("Calling iddp_id / idzp_id  ...",)
        t0 = time.time()
        k, idx, proj = pymatrixid.interp_decomp(A, eps, rand=False)
        t = time.time() - t0
        B = pymatrixid.reconstruct_matrix_from_id(A[:, idx[:k]], idx, proj)
        _debug_print(fmt % (t, np.allclose(A, B, eps)))
        assert_(np.allclose(A, B, eps))

        _debug_print("Calling iddp_aid / idzp_aid ...",)
        t0 = time.time()
        k, idx, proj = pymatrixid.interp_decomp(A, eps)
        t = time.time() - t0
        B = pymatrixid.reconstruct_matrix_from_id(A[:, idx[:k]], idx, proj)
        _debug_print(fmt % (t, np.allclose(A, B, eps)))
        assert_(np.allclose(A, B, eps))

        _debug_print("Calling iddp_rid / idzp_rid ...",)
        t0 = time.time()
        k, idx, proj = pymatrixid.interp_decomp(L, eps)
        t = time.time() - t0
        B = pymatrixid.reconstruct_matrix_from_id(A[:, idx[:k]], idx, proj)
        _debug_print(fmt % (t, np.allclose(A, B, eps)))
        assert_(np.allclose(A, B, eps))

        # fixed rank
        k = rank

        _debug_print("Calling iddr_id / idzr_id  ...",)
        t0 = time.time()
        idx, proj = pymatrixid.interp_decomp(A, k, rand=False)
        t = time.time() - t0
        B = pymatrixid.reconstruct_matrix_from_id(A[:, idx[:k]], idx, proj)
        _debug_print(fmt % (t, np.allclose(A, B, eps)))
        assert_(np.allclose(A, B, eps))

        _debug_print("Calling iddr_aid / idzr_aid ...",)
        t0 = time.time()
        idx, proj = pymatrixid.interp_decomp(A, k)
        t = time.time() - t0
        B = pymatrixid.reconstruct_matrix_from_id(A[:, idx[:k]], idx, proj)
        _debug_print(fmt % (t, np.allclose(A, B, eps)))
        assert_(np.allclose(A, B, eps))

        _debug_print("Calling iddr_rid / idzr_rid ...",)
        t0 = time.time()
        idx, proj = pymatrixid.interp_decomp(L, k)
        t = time.time() - t0
        B = pymatrixid.reconstruct_matrix_from_id(A[:, idx[:k]], idx, proj)
        _debug_print(fmt % (t, np.allclose(A, B, eps)))
        assert_(np.allclose(A, B, eps))

        # check skeleton and interpolation matrices
        idx, proj = pymatrixid.interp_decomp(A, k, rand=False)
        P = pymatrixid.reconstruct_interp_matrix(idx, proj)
        B = pymatrixid.reconstruct_skel_matrix(A, k, idx)
        assert_(np.allclose(B, A[:,idx[:k]], eps))
        assert_(np.allclose(B.dot(P), A, eps))

        # test SVD routines
        _debug_print("-----------------------------------------")
        _debug_print("SVD routines")
        _debug_print("-----------------------------------------")

        # fixed precision
        _debug_print("Calling iddp_svd / idzp_svd ...",)
        t0 = time.time()
        U, S, V = pymatrixid.svd(A, eps, rand=False)
        t = time.time() - t0
        B = np.dot(U, np.dot(np.diag(S), V.T.conj()))
        _debug_print(fmt % (t, np.allclose(A, B, eps)))
        assert_(np.allclose(A, B, eps))

        _debug_print("Calling iddp_asvd / idzp_asvd...",)
        t0 = time.time()
        U, S, V = pymatrixid.svd(A, eps)
        t = time.time() - t0
        B = np.dot(U, np.dot(np.diag(S), V.T.conj()))
        _debug_print(fmt % (t, np.allclose(A, B, eps)))
        assert_(np.allclose(A, B, eps))

        _debug_print("Calling iddp_rsvd / idzp_rsvd...",)
        t0 = time.time()
        U, S, V = pymatrixid.svd(L, eps)
        t = time.time() - t0
        B = np.dot(U, np.dot(np.diag(S), V.T.conj()))
        _debug_print(fmt % (t, np.allclose(A, B, eps)))
        assert_(np.allclose(A, B, eps))

        # fixed rank
        k = rank

        _debug_print("Calling iddr_svd / idzr_svd ...",)
        t0 = time.time()
        U, S, V = pymatrixid.svd(A, k, rand=False)
        t = time.time() - t0
        B = np.dot(U, np.dot(np.diag(S), V.T.conj()))
        _debug_print(fmt % (t, np.allclose(A, B, eps)))
        assert_(np.allclose(A, B, eps))

        _debug_print("Calling iddr_asvd / idzr_asvd ...",)
        t0 = time.time()
        U, S, V = pymatrixid.svd(A, k)
        t = time.time() - t0
        B = np.dot(U, np.dot(np.diag(S), V.T.conj()))
        _debug_print(fmt % (t, np.allclose(A, B, eps)))
        assert_(np.allclose(A, B, eps))

        _debug_print("Calling iddr_rsvd / idzr_rsvd ...",)
        t0 = time.time()
        U, S, V = pymatrixid.svd(L, k)
        t = time.time() - t0
        B = np.dot(U, np.dot(np.diag(S), V.T.conj()))
        _debug_print(fmt % (t, np.allclose(A, B, eps)))
        assert_(np.allclose(A, B, eps))

        # ID to SVD
        idx, proj = pymatrixid.interp_decomp(A, k, rand=False)
        Up, Sp, Vp = pymatrixid.id_to_svd(A[:, idx[:k]], idx, proj)
        B = U.dot(np.diag(S).dot(V.T.conj()))
        assert_(np.allclose(A, B, eps))

        # Norm estimates
        s = svdvals(A)
        norm_2_est = pymatrixid.estimate_spectral_norm(A)
        assert_(np.allclose(norm_2_est, s[0], 1e-6))

        B = A.copy()
        B[:,0] *= 1.2
        s = svdvals(A - B)
        norm_2_est = pymatrixid.estimate_spectral_norm_diff(A, B)
        assert_(np.allclose(norm_2_est, s[0], 1e-6))

        # Rank estimates
        B = np.array([[1, 1, 0], [0, 0, 1], [0, 0, 1]], dtype=dtype)
        for M in [A, B]:
            ML = aslinearoperator(M)

            rank_tol = 1e-9
            rank_np = np.linalg.matrix_rank(M, norm(M, 2)*rank_tol)
            rank_est = pymatrixid.estimate_rank(M, rank_tol)
            rank_est_2 = pymatrixid.estimate_rank(ML, rank_tol)

            assert_(rank_est >= rank_np)
            assert_(rank_est <= rank_np + 10)

            assert_(rank_est_2 >= rank_np - 4)
            assert_(rank_est_2 <= rank_np + 4)

    def test_rand(self):
        pymatrixid.seed('default')
        assert_(np.allclose(pymatrixid.rand(2), [0.8932059, 0.64500803], 1e-4))

        pymatrixid.seed(1234)
        x1 = pymatrixid.rand(2)
        assert_(np.allclose(x1, [0.7513823, 0.06861718], 1e-4))

        np.random.seed(1234)
        pymatrixid.seed()
        x2 = pymatrixid.rand(2)

        np.random.seed(1234)
        pymatrixid.seed(np.random.rand(55))
        x3 = pymatrixid.rand(2)

        assert_allclose(x1, x2)
        assert_allclose(x1, x3)

    def test_badcall(self):
        A = hilbert(5).astype(np.float32)
        assert_raises(ValueError, pymatrixid.interp_decomp, A, 1e-6, rand=False)

    def test_rank_too_large(self):
        # svd(array, k) should not segfault
        a = np.ones((4, 3))
        with assert_raises(ValueError):
            pymatrixid.svd(a, 4)

    def test_full_rank(self):
        eps = 1.0e-12

        # fixed precision
        A = np.random.rand(16, 8)
        k, idx, proj = pymatrixid.interp_decomp(A, eps)
        assert_(k == A.shape[1])

        P = pymatrixid.reconstruct_interp_matrix(idx, proj)
        B = pymatrixid.reconstruct_skel_matrix(A, k, idx)
        assert_allclose(A, B.dot(P))

        # fixed rank
        idx, proj = pymatrixid.interp_decomp(A, k)

        P = pymatrixid.reconstruct_interp_matrix(idx, proj)
        B = pymatrixid.reconstruct_skel_matrix(A, k, idx)
        assert_allclose(A, B.dot(P))

    def test_bug_9793(self):
        dtypes = [np.float_, np.complex_]
        rands = [True, False]
        epss = [1, 0.1]

        for dtype, eps, rand in itertools.product(dtypes, epss, rands):
            A = np.array([[-1, -1, -1, 0, 0, 0],
                          [0, 0, 0, 1, 1, 1],
                          [1, 0, 0, 1, 0, 0],
                          [0, 1, 0, 0, 1, 0],
                          [0, 0, 1, 0, 0, 1]],
                         dtype=dtype, order="C")
            B = A.copy()
            interp_decomp(A.T, eps, rand=rand)
            assert_(np.array_equal(A, B))
