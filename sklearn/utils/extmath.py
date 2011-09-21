"""
Extended math utilities.
"""
# Authors: G. Varoquaux, A. Gramfort, A. Passos, O. Grisel
# License: BSD

import math

from . import check_random_state
import numpy as np

from scipy import linalg

def norm(v):
    v = np.asarray(v)
    __nrm2, = linalg.get_blas_funcs(['nrm2'], [v])
    return __nrm2(v)

def _fast_logdet(A):
    """
    Compute log(det(A)) for A symmetric
    Equivalent to : np.log(np.linalg.det(A))
    but more robust
    It returns -Inf if det(A) is non positive or is not defined.
    """
    # XXX: Should be implemented as in numpy, using ATLAS
    # http://projects.scipy.org/numpy/browser/trunk/numpy/linalg/linalg.py#L1559
    ld = np.sum(np.log(np.diag(A)))
    a = np.exp(ld / A.shape[0])
    d = np.linalg.det(A / a)
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
    sign, ld = np.linalg.slogdet(A)
    if not sign > 0:
        return -np.inf
    return ld


# Numpy >= 1.5 provides a fast logdet
if hasattr(np.linalg, 'slogdet'):
    fast_logdet = _fast_logdet_numpy
else:
    fast_logdet = _fast_logdet

try:
    factorial = math.factorial
except AttributeError:
     # math.factorial is only available in Python >= 2.6
     import operator
     def factorial(n):
         # don't use reduce operator or 2to3 will fail.
         # ripped from http://www.joelbdalley.com/page.pl?38
         # Ensure that n is a Natural number
         n = abs(int(n))
         if n < 1: n = 1

         # Store n! in variable x
         x = 1

         # Compute n!
         for i in range(1, n + 1):
             x = i * x

         # Return n!
         return x


try:
    import itertools
    combinations = itertools.combinations
except AttributeError:
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


def density(w, **kwargs):
    """Compute density of a sparse vector

    Return a value between 0 and 1
    """
    if hasattr(w, "tocsr"):
        d = float(w.data.size) / w.size
    else:
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


def fast_svd(M, k, p=None, q=0, transpose='auto', random_state=0):
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

    random_state: RandomState or an int seed (0 by default)
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
    if p == None:
        p = k

    random_state = check_random_state(random_state)
    n_samples, n_features = M.shape

    if transpose == 'auto' and n_samples > n_features:
        transpose = True
    if transpose:
        # this implementation is a bit faster with smaller shape[1]
        M = M.T

   # generating random gaussian vectors r with shape: (M.shape[1], k + p)
    r = random_state.normal(size=(M.shape[1], k + p))

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


def logsum(arr, axis=0):
    """ Computes the sum of arr assuming arr is in the log domain.

    Returns log(sum(exp(arr))) while minimizing the possibility of
    over/underflow.

    Examples
    ========

    >>> import numpy as np
    >>> from sklearn.utils.extmath import logsum
    >>> a = np.arange(10)
    >>> np.log(np.sum(np.exp(a)))
    9.4586297444267107
    >>> logsum(a)
    9.4586297444267107
    """
    arr = np.rollaxis(arr, axis)
    # Use the max to normalize, as with the log this is what accumulates
    # the less errors
    vmax = arr.max(axis=0)
    out = np.log(np.sum(np.exp(arr - vmax), axis=0))
    out += vmax
    return out


def weighted_mode(a, w, axis=0):
    """Returns an array of the weighted modal (most common) value in the
    passed array.

    If there is more than one such value, only the first is returned.
    The bin-count for the modal bins is also returned.

    This is an extension of the algorithm in scipy.stats.mode.

    Parameters
    ----------
    a : array_like
        n-dimensional array of which to find mode(s).
    w : array_like
        n-dimensional array of weights for each value
    axis : int, optional
        Axis along which to operate. Default is 0, i.e. the first axis.

    Returns
    -------
    vals : ndarray
        Array of modal values.
    score : ndarray
        Array of weighted counts for each mode.

    Examples
    --------
    >>> from sklearn.utils.extmath import weighted_mode
    >>> x = [4, 1, 4, 2, 4, 2]
    >>> weights = [1, 1, 1, 1, 1, 1]
    >>> weighted_mode(x, weights)
    (array([ 4.]), array([ 3.]))

    The value 4 appears three times: with uniform weights, the result is
    simply the mode of the distribution.

    >>> weights = [1, 3, 0.5, 1.5, 1, 2] # deweight the 4's
    >>> weighted_mode(x, weights)
    (array([ 2.]), array([ 3.5]))

    The value 2 has the highest score: it appears twice with weights of
    1.5 and 2: the sum of these is 3.

    See Also
    --------
    scipy.stats.mode
    """
    if axis is None:
        a = np.ravel(a)
        w = np.ravel(w)
        axis = 0
    else:
        a = np.asarray(a)
        w = np.asarray(w)
        axis = axis

    if a.shape != w.shape:
        w = np.zeros(a.shape, dtype=w.dtype) + w

    scores = np.unique(np.ravel(a))       # get ALL unique values
    testshape = list(a.shape)
    testshape[axis] = 1
    oldmostfreq = np.zeros(testshape)
    oldcounts = np.zeros(testshape)
    for score in scores:
        template = np.zeros(a.shape)
        ind = (a == score)
        template[ind] = w[ind]
        counts = np.expand_dims(np.sum(template, axis), axis)
        mostfrequent = np.where(counts > oldcounts, score, oldmostfreq)
        oldcounts = np.maximum(counts, oldcounts)
        oldmostfreq = mostfrequent
    return mostfrequent, oldcounts
