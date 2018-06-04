""" Test functions for the sparse.linalg.eigen.lobpcg module
"""
from __future__ import division, print_function, absolute_import

import itertools

import numpy as np
from numpy.testing import (assert_almost_equal, assert_equal,
        assert_allclose, assert_array_less, assert_)

from scipy import ones, rand, r_, diag, linalg, eye
from scipy.linalg import eig, eigh, toeplitz
import scipy.sparse
from scipy.sparse.linalg.eigen.lobpcg import lobpcg


def ElasticRod(n):
    # Fixed-free elastic rod
    L = 1.0
    le = L/n
    rho = 7.85e3
    S = 1.e-4
    E = 2.1e11
    mass = rho*S*le/6.
    k = E*S/le
    A = k*(diag(r_[2.*ones(n-1),1])-diag(ones(n-1),1)-diag(ones(n-1),-1))
    B = mass*(diag(r_[4.*ones(n-1),2])+diag(ones(n-1),1)+diag(ones(n-1),-1))
    return A,B


def MikotaPair(n):
    # Mikota pair acts as a nice test since the eigenvalues
    # are the squares of the integers n, n=1,2,...
    x = np.arange(1,n+1)
    B = diag(1./x)
    y = np.arange(n-1,0,-1)
    z = np.arange(2*n-1,0,-2)
    A = diag(z)-diag(y,-1)-diag(y,1)
    return A,B


def compare_solutions(A,B,m):
    n = A.shape[0]

    np.random.seed(0)

    V = rand(n,m)
    X = linalg.orth(V)

    eigs,vecs = lobpcg(A, X, B=B, tol=1e-5, maxiter=30)
    eigs.sort()

    w,v = eig(A,b=B)
    w.sort()

    assert_almost_equal(w[:int(m/2)],eigs[:int(m/2)],decimal=2)


def test_Small():
    A,B = ElasticRod(10)
    compare_solutions(A,B,10)
    A,B = MikotaPair(10)
    compare_solutions(A,B,10)


def test_ElasticRod():
    A,B = ElasticRod(100)
    compare_solutions(A,B,20)


def test_MikotaPair():
    A,B = MikotaPair(100)
    compare_solutions(A,B,20)


def test_trivial():
    n = 5
    X = ones((n, 1))
    A = eye(n)
    compare_solutions(A, None, n)


def test_regression():
    # https://mail.python.org/pipermail/scipy-user/2010-October/026944.html
    n = 10
    X = np.ones((n, 1))
    A = np.identity(n)
    w, V = lobpcg(A, X)
    assert_allclose(w, [1])


def test_diagonal():
    # This test was moved from '__main__' in lobpcg.py.
    # Coincidentally or not, this is the same eigensystem
    # required to reproduce arpack bug
    # http://forge.scilab.org/index.php/p/arpack-ng/issues/1397/
    # even using the same n=100.

    np.random.seed(1234)

    # The system of interest is of size n x n.
    n = 100

    # We care about only m eigenpairs.
    m = 4

    # Define the generalized eigenvalue problem Av = cBv
    # where (c, v) is a generalized eigenpair,
    # and where we choose A to be the diagonal matrix whose entries are 1..n
    # and where B is chosen to be the identity matrix.
    vals = np.arange(1, n+1, dtype=float)
    A = scipy.sparse.diags([vals], [0], (n, n))
    B = scipy.sparse.eye(n)

    # Let the preconditioner M be the inverse of A.
    M = scipy.sparse.diags([np.reciprocal(vals)], [0], (n, n))

    # Pick random initial vectors.
    X = np.random.rand(n, m)

    # Require that the returned eigenvectors be in the orthogonal complement
    # of the first few standard basis vectors.
    m_excluded = 3
    Y = np.eye(n, m_excluded)

    eigs, vecs = lobpcg(A, X, B, M=M, Y=Y, tol=1e-4, maxiter=40, largest=False)

    assert_allclose(eigs, np.arange(1+m_excluded, 1+m_excluded+m))
    _check_eigen(A, eigs, vecs, rtol=1e-3, atol=1e-3)


def _check_eigen(M, w, V, rtol=1e-8, atol=1e-14):
    mult_wV = np.multiply(w, V)
    dot_MV = M.dot(V)
    assert_allclose(mult_wV, dot_MV, rtol=rtol, atol=atol)


def _check_fiedler(n, p):
    # This is not necessarily the recommended way to find the Fiedler vector.
    np.random.seed(1234)
    col = np.zeros(n)
    col[1] = 1
    A = toeplitz(col)
    D = np.diag(A.sum(axis=1))
    L = D - A
    # Compute the full eigendecomposition using tricks, e.g.
    # http://www.cs.yale.edu/homes/spielman/561/2009/lect02-09.pdf
    tmp = np.pi * np.arange(n) / n
    analytic_w = 2 * (1 - np.cos(tmp))
    analytic_V = np.cos(np.outer(np.arange(n) + 1/2, tmp))
    _check_eigen(L, analytic_w, analytic_V)
    # Compute the full eigendecomposition using eigh.
    eigh_w, eigh_V = eigh(L)
    _check_eigen(L, eigh_w, eigh_V)
    # Check that the first eigenvalue is near zero and that the rest agree.
    assert_array_less(np.abs([eigh_w[0], analytic_w[0]]), 1e-14)
    assert_allclose(eigh_w[1:], analytic_w[1:])

    # Check small lobpcg eigenvalues.
    X = analytic_V[:, :p]
    lobpcg_w, lobpcg_V = lobpcg(L, X, largest=False)
    assert_equal(lobpcg_w.shape, (p,))
    assert_equal(lobpcg_V.shape, (n, p))
    _check_eigen(L, lobpcg_w, lobpcg_V)
    assert_array_less(np.abs(np.min(lobpcg_w)), 1e-14)
    assert_allclose(np.sort(lobpcg_w)[1:], analytic_w[1:p])

    # Check large lobpcg eigenvalues.
    X = analytic_V[:, -p:]
    lobpcg_w, lobpcg_V = lobpcg(L, X, largest=True)
    assert_equal(lobpcg_w.shape, (p,))
    assert_equal(lobpcg_V.shape, (n, p))
    _check_eigen(L, lobpcg_w, lobpcg_V)
    assert_allclose(np.sort(lobpcg_w), analytic_w[-p:])

    # Look for the Fiedler vector using good but not exactly correct guesses.
    fiedler_guess = np.concatenate((np.ones(n//2), -np.ones(n-n//2)))
    X = np.vstack((np.ones(n), fiedler_guess)).T
    lobpcg_w, lobpcg_V = lobpcg(L, X, largest=False)
    # Mathematically, the smaller eigenvalue should be zero
    # and the larger should be the algebraic connectivity.
    lobpcg_w = np.sort(lobpcg_w)
    assert_allclose(lobpcg_w, analytic_w[:2], atol=1e-14)


def test_fiedler_small_8():
    # This triggers the dense path because 8 < 2*5.
    _check_fiedler(8, 2)


def test_fiedler_large_12():
    # This does not trigger the dense path, because 2*5 <= 12.
    _check_fiedler(12, 2)


def test_hermitian():
    np.random.seed(1234)

    sizes = [3, 10, 50]
    ks = [1, 3, 10, 50]
    gens = [True, False]

    for size, k, gen in itertools.product(sizes, ks, gens):
        if k > size:
            continue

        H = np.random.rand(size, size) + 1.j * np.random.rand(size, size)
        H = 10 * np.eye(size) + H + H.T.conj()

        X = np.random.rand(size, k)

        if not gen:
            B = np.eye(size)
            w, v = lobpcg(H, X, maxiter=5000)
            w0, v0 = eigh(H)
        else:
            B = np.random.rand(size, size) + 1.j * np.random.rand(size, size)
            B = 10 * np.eye(size) + B.dot(B.T.conj())
            w, v = lobpcg(H, X, B, maxiter=5000)
            w0, v0 = eigh(H, B)

        for wx, vx in zip(w, v.T):
            # Check eigenvector
            assert_allclose(np.linalg.norm(H.dot(vx) - B.dot(vx) * wx) / np.linalg.norm(H.dot(vx)),
                            0, atol=5e-4, rtol=0)

            # Compare eigenvalues
            j = np.argmin(abs(w0 - wx))
            assert_allclose(wx, w0[j], rtol=1e-4)
