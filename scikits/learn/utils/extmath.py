"""
Extended math utilities.
"""
# Authors: G. Varoquaux, A. Gramfort, A. Passos, O. Grisel
# License: BSD

import sys
import math

import numpy as np

#XXX: We should have a function with numpy's slogdet API
def _fast_logdet(A):
    """
    Compute log(det(A)) for A symmetric
    Equivalent to : np.log(np.linalg.det(A))
    but more robust
    It returns -Inf if det(A) is non positive or is not defined.
    """
    # XXX: Should be implemented as in numpy, using ATLAS
    # http://projects.scipy.org/numpy/browser/trunk/numpy/linalg/linalg.py#L1559
    from scipy import linalg
    ld = np.sum(np.log(np.diag(A)))
    a = np.exp(ld/A.shape[0])
    d = np.linalg.det(A/a)
    ld += np.log(d)
    if not np.isfinite(ld):
        return -np.inf
    return ld

def _fast_logdet_numpy(A):
    """
    Compute log(det(A)) for A symmetric
    Equivalent to : np.log(nl.det(A))
    but more robust
    It returns -Inf if det(A) is non positive or is not defined.
    """
    from scipy import linalg
    sign, ld = np.linalg.slogdet(A)
    if not sign > 0:
        return -np.inf
    return ld


# Numpy >= 1.5 provides a fast logdet
if hasattr(np.linalg, 'slogdet'):
    fast_logdet = _fast_logdet_numpy
else:
    fast_logdet = _fast_logdet

if sys.version_info[1] < 6:
    # math.factorial is only available in 2.6
    def factorial(x) :
        # simple recursive implementation
        if x == 0: return 1
        return x * factorial(x-1)
else:
    factorial = math.factorial


if sys.version_info[1] < 6:
    def combinations(seq, r=None):
        """Generator returning combinations of items from sequence <seq>
        taken <r> at a time. Order is not significant. If <r> is not given,
        the entire sequence is returned.
        """
        if r == None:
            r = len(seq)
        if r <= 0:
            yield []
        else:
            for i in xrange(len(seq)):
                for cc in combinations(seq[i+1:], r-1):
                    yield [seq[i]]+cc

else:
    import itertools
    combinations = itertools.combinations


def density(w, **kwargs):
    """Compute density of a sparse vector

    Return a value between 0 and 1
    """
    d = 0 if w is None else float((w != 0).sum()) / w.size
    return d


def safe_sparse_dot(a, b, dense_output=False):
    """Dot product that handle the sparse matrix case correctly"""
    from scipy import sparse
    if sparse.issparse(a) or sparse.issparse(b):
        ret = a * b
        if dense_output and hasattr(ret, "toarray"):
            ret = ret.toarray()
        return ret
    else:
        return np.dot(a,b)


def fast_svd(M, k, p=None, q=0, transpose='auto', rng=0):
    """Computes the k-truncated randomized SVD

    Parameters
    ===========
    M: ndarray or sparse matrix
        Matrix to decompose

    k: int
        Number of singular values and vectors to extract.

    p: int (default is k)
        Additional number of samples of the range of M to ensure proper
        conditioning. See the notes below.

    q: int (default is 0)
        Number of power iterations (can be used to deal with very noisy
        problems).

    transpose: True, False or 'auto' (default)
        Whether the algorithm should be applied to M.T instead of M. The
        result should approximately be the same. The 'auto' mode will
        trigger the transposition if M.shape[1] > M.shape[0] since this
        implementation of randomized SVD tend to be a little faster in that
        case).

    rng: RandomState or an int seed (0 by default)
        A random number generator instance to make behavior

    Notes
    =====
    This algorithm finds the exact truncated singular values decomposition
    using randomization to speed up the computations. It is particularly
    fast on large matrices on which you whish to extract only a small
    number of components.

    (k + p) should be strictly higher than the rank of M. This can be
    checked by ensuring that the lowest extracted singular value is on
    the order of the machine precision of floating points.

    References
    ==========
    Finding structure with randomness: Stochastic algorithms for constructing
    approximate matrix decompositions
    Halko, et al., 2009 (arXiv:909)

    A randomized algorithm for the decomposition of matrices
    Per-Gunnar Martinsson, Vladimir Rokhlin and Mark Tygert
    """
    # lazy import of scipy sparse, because it is very slow.
    from scipy import sparse

    if p == None:
        p = k

    if rng is None:
        rng = np.random.RandomState()
    elif isinstance(rng, int):
        rng = np.random.RandomState(rng)

    n_samples, n_features = M.shape

    if transpose == 'auto' and n_samples > n_features:
        transpose = True
    if transpose:
        # this implementation is a bit faster with smaller shape[1]
        M = M.T

   # generating random gaussian vectors r with shape: (M.shape[1], k + p)
    r = rng.normal(size=(M.shape[1], k + p))

    # sampling the range of M using by linear projection of r
    Y = safe_sparse_dot(M, r)
    del r

    # apply q power iterations on Y to make to further 'imprint' the top
    # singular values of M in Y
    for i in xrange(q):
        Y = safe_sparse_dot(M, safe_sparse_dot(M.T, Y))

    # extracting an orthonormal basis of the M range samples
    from .fixes import qr_economic
    Q, R = qr_economic(Y)
    del R

    # project M to the (k + p) dimensional space using the basis vectors
    B = safe_sparse_dot(Q.T, M)

    # compute the SVD on the thin matrix: (k + p) wide
    from scipy import linalg
    Uhat, s, V = linalg.svd(B, full_matrices=False)
    del B
    U = np.dot(Q, Uhat)

    if transpose:
        # transpose back the results according to the input convention
        return V[:k, :].T, s[:k], U[:, :k].T
    else:
        return U[:, :k], s[:k], V[:k, :]

