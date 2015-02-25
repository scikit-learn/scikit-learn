"""
Extended math utilities.
"""
# Authors: Gael Varoquaux
#          Alexandre Gramfort
#          Alexandre T. Passos
#          Olivier Grisel
#          Lars Buitinck
#          Stefan van der Walt
#          Kyle Kastner
# License: BSD 3 clause

from __future__ import division
from functools import partial
import warnings

import numpy as np
from scipy import linalg
from scipy.sparse import issparse

from . import check_random_state, deprecated
from .fixes import np_version
from ._logistic_sigmoid import _log_logistic_sigmoid
from ..externals.six.moves import xrange
from .sparsefuncs_fast import csr_row_norms
from .validation import check_array, NonBLASDotWarning


def norm(x):
    """Compute the Euclidean or Frobenius norm of x.

    Returns the Euclidean norm when x is a vector, the Frobenius norm when x
    is a matrix (2-d array). More precise than sqrt(squared_norm(x)).
    """
    x = np.asarray(x)
    nrm2, = linalg.get_blas_funcs(['nrm2'], [x])
    return nrm2(x)


# Newer NumPy has a ravel that needs less copying.
if np_version < (1, 7, 1):
    _ravel = np.ravel
else:
    _ravel = partial(np.ravel, order='K')


def squared_norm(x):
    """Squared Euclidean or Frobenius norm of x.

    Returns the Euclidean norm when x is a vector, the Frobenius norm when x
    is a matrix (2-d array). Faster than norm(x) ** 2.
    """
    x = _ravel(x)
    return np.dot(x, x)


def row_norms(X, squared=False):
    """Row-wise (squared) Euclidean norm of X.

    Equivalent to np.sqrt((X * X).sum(axis=1)), but also supports CSR sparse
    matrices and does not create an X.shape-sized temporary.

    Performs no input validation.
    """
    if issparse(X):
        norms = csr_row_norms(X)
    else:
        norms = np.einsum('ij,ij->i', X, X)

    if not squared:
        np.sqrt(norms, norms)
    return norms


def fast_logdet(A):
    """Compute log(det(A)) for A symmetric

    Equivalent to : np.log(nl.det(A)) but more robust.
    It returns -Inf if det(A) is non positive or is not defined.
    """
    sign, ld = np.linalg.slogdet(A)
    if not sign > 0:
        return -np.inf
    return ld


def _impose_f_order(X):
    """Helper Function"""
    # important to access flags instead of calling np.isfortran,
    # this catches corner cases.
    if X.flags.c_contiguous:
        return check_array(X.T, copy=False, order='F'), True
    else:
        return check_array(X, copy=False, order='F'), False


def _fast_dot(A, B):
    if B.shape[0] != A.shape[A.ndim - 1]:  # check adopted from '_dotblas.c'
        raise ValueError

    if A.dtype != B.dtype or any(x.dtype not in (np.float32, np.float64)
                                 for x in [A, B]):
        warnings.warn('Data must be of same type. Supported types '
                      'are 32 and 64 bit float. '
                      'Falling back to np.dot.', NonBLASDotWarning)
        raise ValueError

    if min(A.shape) == 1 or min(B.shape) == 1 or A.ndim != 2 or B.ndim != 2:
        raise ValueError

    # scipy 0.9 compliant API
    dot = linalg.get_blas_funcs(['gemm'], (A, B))[0]
    A, trans_a = _impose_f_order(A)
    B, trans_b = _impose_f_order(B)
    return dot(alpha=1.0, a=A, b=B, trans_a=trans_a, trans_b=trans_b)


def _have_blas_gemm():
    try:
        linalg.get_blas_funcs(['gemm'])
        return True
    except (AttributeError, ValueError):
        warnings.warn('Could not import BLAS, falling back to np.dot')
        return False


# Only use fast_dot for older NumPy; newer ones have tackled the speed issue.
if np_version < (1, 7, 2) and _have_blas_gemm():
    def fast_dot(A, B):
        """Compute fast dot products directly calling BLAS.

        This function calls BLAS directly while warranting Fortran contiguity.
        This helps avoiding extra copies `np.dot` would have created.
        For details see section `Linear Algebra on large Arrays`:
        http://wiki.scipy.org/PerformanceTips

        Parameters
        ----------
        A, B: instance of np.ndarray
            Input arrays. Arrays are supposed to be of the same dtype and to
            have exactly 2 dimensions. Currently only floats are supported.
            In case these requirements aren't met np.dot(A, B) is returned
            instead. To activate the related warning issued in this case
            execute the following lines of code:

            >> import warnings
            >> from sklearn.utils.validation import NonBLASDotWarning
            >> warnings.simplefilter('always', NonBLASDotWarning)
        """
        try:
            return _fast_dot(A, B)
        except ValueError:
            # Maltyped or malformed data.
            return np.dot(A, B)
else:
    fast_dot = np.dot


def density(w, **kwargs):
    """Compute density of a sparse vector

    Return a value between 0 and 1
    """
    if hasattr(w, "toarray"):
        d = float(w.nnz) / (w.shape[0] * w.shape[1])
    else:
        d = 0 if w is None else float((w != 0).sum()) / w.size
    return d


def safe_sparse_dot(a, b, dense_output=False):
    """Dot product that handle the sparse matrix case correctly

    Uses BLAS GEMM as replacement for numpy.dot where possible
    to avoid unnecessary copies.
    """
    if issparse(a) or issparse(b):
        ret = a * b
        if dense_output and hasattr(ret, "toarray"):
            ret = ret.toarray()
        return ret
    else:
        return fast_dot(a, b)


def randomized_range_finder(A, size, n_iter, random_state=None):
    """Computes an orthonormal matrix whose range approximates the range of A.

    Parameters
    ----------
    A: 2D array
        The input data matrix
    size: integer
        Size of the return array
    n_iter: integer
        Number of power iterations used to stabilize the result
    random_state: RandomState or an int seed (0 by default)
        A random number generator instance

    Returns
    -------
    Q: 2D array
        A (size x size) projection matrix, the range of which
        approximates well the range of the input matrix A.

    Notes
    -----

    Follows Algorithm 4.3 of
    Finding structure with randomness: Stochastic algorithms for constructing
    approximate matrix decompositions
    Halko, et al., 2009 (arXiv:909) http://arxiv.org/pdf/0909.4061
    """
    random_state = check_random_state(random_state)

    # generating random gaussian vectors r with shape: (A.shape[1], size)
    R = random_state.normal(size=(A.shape[1], size))

    # sampling the range of A using by linear projection of r
    Y = safe_sparse_dot(A, R)
    del R

    # perform power iterations with Y to further 'imprint' the top
    # singular vectors of A in Y
    for i in xrange(n_iter):
        Y = safe_sparse_dot(A, safe_sparse_dot(A.T, Y))

    # extracting an orthonormal basis of the A range samples
    Q, R = linalg.qr(Y, mode='economic')
    return Q


def randomized_svd(M, n_components, n_oversamples=10, n_iter=0,
                   transpose='auto', flip_sign=True, random_state=0):
    """Computes a truncated randomized SVD

    Parameters
    ----------
    M: ndarray or sparse matrix
        Matrix to decompose

    n_components: int
        Number of singular values and vectors to extract.

    n_oversamples: int (default is 10)
        Additional number of random vectors to sample the range of M so as
        to ensure proper conditioning. The total number of random vectors
        used to find the range of M is n_components + n_oversamples.

    n_iter: int (default is 0)
        Number of power iterations (can be used to deal with very noisy
        problems).

    transpose: True, False or 'auto' (default)
        Whether the algorithm should be applied to M.T instead of M. The
        result should approximately be the same. The 'auto' mode will
        trigger the transposition if M.shape[1] > M.shape[0] since this
        implementation of randomized SVD tend to be a little faster in that
        case).

    flip_sign: boolean, (True by default)
        The output of a singular value decomposition is only unique up to a
        permutation of the signs of the singular vectors. If `flip_sign` is
        set to `True`, the sign ambiguity is resolved by making the largest
        loadings for each component in the left singular vectors positive.

    random_state: RandomState or an int seed (0 by default)
        A random number generator instance to make behavior

    Notes
    -----
    This algorithm finds a (usually very good) approximate truncated
    singular value decomposition using randomization to speed up the
    computations. It is particularly fast on large matrices on which
    you wish to extract only a small number of components.

    References
    ----------
    * Finding structure with randomness: Stochastic algorithms for constructing
      approximate matrix decompositions
      Halko, et al., 2009 http://arxiv.org/abs/arXiv:0909.4061

    * A randomized algorithm for the decomposition of matrices
      Per-Gunnar Martinsson, Vladimir Rokhlin and Mark Tygert
    """
    random_state = check_random_state(random_state)
    n_random = n_components + n_oversamples
    n_samples, n_features = M.shape

    if transpose == 'auto' and n_samples > n_features:
        transpose = True
    if transpose:
        # this implementation is a bit faster with smaller shape[1]
        M = M.T

    Q = randomized_range_finder(M, n_random, n_iter, random_state)

    # project M to the (k + p) dimensional space using the basis vectors
    B = safe_sparse_dot(Q.T, M)

    # compute the SVD on the thin matrix: (k + p) wide
    Uhat, s, V = linalg.svd(B, full_matrices=False)
    del B
    U = np.dot(Q, Uhat)

    if flip_sign:
        U, V = svd_flip(U, V)

    if transpose:
        # transpose back the results according to the input convention
        return V[:n_components, :].T, s[:n_components], U[:, :n_components].T
    else:
        return U[:, :n_components], s[:n_components], V[:n_components, :]


def logsumexp(arr, axis=0):
    """Computes the sum of arr assuming arr is in the log domain.

    Returns log(sum(exp(arr))) while minimizing the possibility of
    over/underflow.

    Examples
    --------

    >>> import numpy as np
    >>> from sklearn.utils.extmath import logsumexp
    >>> a = np.arange(10)
    >>> np.log(np.sum(np.exp(a)))
    9.4586297444267107
    >>> logsumexp(a)
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
    """Returns an array of the weighted modal (most common) value in a

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


def pinvh(a, cond=None, rcond=None, lower=True):
    """Compute the (Moore-Penrose) pseudo-inverse of a hermetian matrix.

    Calculate a generalized inverse of a symmetric matrix using its
    eigenvalue decomposition and including all 'large' eigenvalues.

    Parameters
    ----------
    a : array, shape (N, N)
        Real symmetric or complex hermetian matrix to be pseudo-inverted

    cond : float or None, default None
        Cutoff for 'small' eigenvalues.
        Singular values smaller than rcond * largest_eigenvalue are considered
        zero.

        If None or -1, suitable machine precision is used.

    rcond : float or None, default None (deprecated)
        Cutoff for 'small' eigenvalues.
        Singular values smaller than rcond * largest_eigenvalue are considered
        zero.

        If None or -1, suitable machine precision is used.
 
    lower : boolean
        Whether the pertinent array data is taken from the lower or upper
        triangle of a. (Default: lower)

    Returns
    -------
    B : array, shape (N, N)

    Raises
    ------
    LinAlgError
        If eigenvalue does not converge

    Examples
    --------
    >>> import numpy as np
    >>> a = np.random.randn(9, 6)
    >>> a = np.dot(a, a.T)
    >>> B = pinvh(a)
    >>> np.allclose(a, np.dot(a, np.dot(B, a)))
    True
    >>> np.allclose(B, np.dot(B, np.dot(a, B)))
    True

    """
    a = np.asarray_chkfinite(a)
    s, u = linalg.eigh(a, lower=lower)

    if rcond is not None:
        cond = rcond
    if cond in [None, -1]:
        t = u.dtype.char.lower()
        factor = {'f': 1E3, 'd': 1E6}
        cond = factor[t] * np.finfo(t).eps

    # unlike svd case, eigh can lead to negative eigenvalues
    above_cutoff = (abs(s) > cond * np.max(abs(s)))
    psigma_diag = np.zeros_like(s)
    psigma_diag[above_cutoff] = 1.0 / s[above_cutoff]

    return np.dot(u * psigma_diag, np.conjugate(u).T)


def cartesian(arrays, out=None):
    """Generate a cartesian product of input arrays.

    Parameters
    ----------
    arrays : list of array-like
        1-D arrays to form the cartesian product of.
    out : ndarray
        Array to place the cartesian product in.

    Returns
    -------
    out : ndarray
        2-D array of shape (M, len(arrays)) containing cartesian products
        formed of input arrays.

    Examples
    --------
    >>> cartesian(([1, 2, 3], [4, 5], [6, 7]))
    array([[1, 4, 6],
           [1, 4, 7],
           [1, 5, 6],
           [1, 5, 7],
           [2, 4, 6],
           [2, 4, 7],
           [2, 5, 6],
           [2, 5, 7],
           [3, 4, 6],
           [3, 4, 7],
           [3, 5, 6],
           [3, 5, 7]])

    """
    arrays = [np.asarray(x) for x in arrays]
    shape = (len(x) for x in arrays)
    dtype = arrays[0].dtype

    ix = np.indices(shape)
    ix = ix.reshape(len(arrays), -1).T

    if out is None:
        out = np.empty_like(ix, dtype=dtype)

    for n, arr in enumerate(arrays):
        out[:, n] = arrays[n][ix[:, n]]

    return out


def svd_flip(u, v, u_based_decision=True):
    """Sign correction to ensure deterministic output from SVD.

    Adjusts the columns of u and the rows of v such that the loadings in the
    columns in u that are largest in absolute value are always positive.

    Parameters
    ----------
    u, v : ndarray
        u and v are the output of `linalg.svd` or 
        `sklearn.utils.extmath.randomized_svd`, with matching inner dimensions
        so one can compute `np.dot(u * s, v)`.

    u_based_decision : boolean, (default=True)
        If True, use the columns of u as the basis for sign flipping. Otherwise,
        use the rows of v. The choice of which variable to base the decision on
        is generally algorithm dependent.


    Returns
    -------
    u_adjusted, v_adjusted : arrays with the same dimensions as the input.

    """
    if u_based_decision:
        # columns of u, rows of v
        max_abs_cols = np.argmax(np.abs(u), axis=0)
        signs = np.sign(u[max_abs_cols, xrange(u.shape[1])])
        u *= signs
        v *= signs[:, np.newaxis]
    else:
        # rows of v, columns of u
        max_abs_rows = np.argmax(np.abs(v), axis=1)
        signs = np.sign(v[xrange(v.shape[0]), max_abs_rows])
        u *= signs
        v *= signs[:, np.newaxis]
    return u, v


@deprecated('to be removed in 0.17; use scipy.special.expit or log_logistic')
def logistic_sigmoid(X, log=False, out=None):
    """Logistic function, ``1 / (1 + e ** (-x))``, or its log."""
    from .fixes import expit
    fn = log_logistic if log else expit
    return fn(X, out)


def log_logistic(X, out=None):
    """Compute the log of the logistic function, ``log(1 / (1 + e ** -x))``.

    This implementation is numerically stable because it splits positive and
    negative values::

        -log(1 + exp(-x_i))     if x_i > 0
        x_i - log(1 + exp(x_i)) if x_i <= 0

    For the ordinary logistic function, use ``sklearn.utils.fixes.expit``.

    Parameters
    ----------
    X: array-like, shape (M, N)
        Argument to the logistic function

    out: array-like, shape: (M, N), optional:
        Preallocated output array.

    Returns
    -------
    out: array, shape (M, N)
        Log of the logistic function evaluated at every point in x

    Notes
    -----
    See the blog post describing this implementation:
    http://fa.bianp.net/blog/2013/numerical-optimizers-for-logistic-regression/
    """
    is_1d = X.ndim == 1
    X = check_array(X, dtype=np.float)

    n_samples, n_features = X.shape

    if out is None:
        out = np.empty_like(X)

    _log_logistic_sigmoid(n_samples, n_features, X, out)

    if is_1d:
        return np.squeeze(out)
    return out


def safe_min(X):
    """Returns the minimum value of a dense or a CSR/CSC matrix.

    Adapated from http://stackoverflow.com/q/13426580

    """
    if issparse(X):
        if len(X.data) == 0:
            return 0
        m = X.data.min()
        return m if X.getnnz() == X.size else min(m, 0)
    else:
        return X.min()


def make_nonnegative(X, min_value=0):
    """Ensure `X.min()` >= `min_value`."""
    min_ = safe_min(X)
    if min_ < min_value:
        if issparse(X):
            raise ValueError("Cannot make the data matrix"
                             " nonnegative because it is sparse."
                             " Adding a value to every entry would"
                             " make it no longer sparse.")
        X = X + (min_value - min_)
    return X


def _batch_mean_variance_update(X, old_mean, old_variance, old_sample_count):
    """Calculate an average mean update and a Youngs and Cramer variance update.

    From the paper "Algorithms for computing the sample variance: analysis and
    recommendations", by Chan, Golub, and LeVeque.

    Parameters
    ----------
    X : array-like, shape (n_samples, n_features)
        Data to use for variance update

    old_mean : array-like, shape: (n_features,)

    old_variance : array-like, shape: (n_features,)

    old_sample_count : int

    Returns
    -------
    updated_mean : array, shape (n_features,)

    updated_variance : array, shape (n_features,)

    updated_sample_count : int

    References
    ----------
    T. Chan, G. Golub, R. LeVeque. Algorithms for computing the sample variance:
        recommendations, The American Statistician, Vol. 37, No. 3, pp. 242-247

    """
    new_sum = X.sum(axis=0)
    new_variance = X.var(axis=0) * X.shape[0]
    old_sum = old_mean * old_sample_count
    n_samples = X.shape[0]
    updated_sample_count = old_sample_count + n_samples
    partial_variance = old_sample_count / (n_samples * updated_sample_count) * (
        n_samples / old_sample_count * old_sum - new_sum) ** 2
    unnormalized_variance = old_variance * old_sample_count + new_variance + \
        partial_variance
    return ((old_sum + new_sum) / updated_sample_count,
            unnormalized_variance / updated_sample_count,
            updated_sample_count)


def _deterministic_vector_sign_flip(u):
    """Modify the sign of vectors for reproducibility
    
    Flips the sign of elements of all the vectors (rows of u) such that
    the absolute maximum element of each vector is positive.

    Parameters
    ----------
    u : ndarray
        Array with vectors as its rows.

    Returns
    -------
    u_flipped : ndarray with same shape as u
        Array with the sign flipped vectors as its rows.
    """
    max_abs_rows = np.argmax(np.abs(u), axis=1)
    signs = np.sign(u[range(u.shape[0]), max_abs_rows])
    u *= signs[:, np.newaxis]
    return u
