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
#          Giorgio Patrini
# License: BSD 3 clause

import warnings

import numpy as np
from scipy import linalg, sparse

from . import check_random_state
from ._logistic_sigmoid import _log_logistic_sigmoid
from .sparsefuncs_fast import csr_row_norms
from .validation import check_array
from .validation import _deprecate_positional_args


def squared_norm(x):
    """Squared Euclidean or Frobenius norm of x.

    Faster than norm(x) ** 2.

    Parameters
    ----------
    x : array-like

    Returns
    -------
    float
        The Euclidean norm when x is a vector, the Frobenius norm when x
        is a matrix (2-d array).
    """
    x = np.ravel(x, order='K')
    if np.issubdtype(x.dtype, np.integer):
        warnings.warn('Array type is integer, np.dot may overflow. '
                      'Data should be float type to avoid this issue',
                      UserWarning)
    return np.dot(x, x)


def row_norms(X, squared=False):
    """Row-wise (squared) Euclidean norm of X.

    Equivalent to np.sqrt((X * X).sum(axis=1)), but also supports sparse
    matrices and does not create an X.shape-sized temporary.

    Performs no input validation.

    Parameters
    ----------
    X : array-like
        The input array.
    squared : bool, default=False
        If True, return squared norms.

    Returns
    -------
    array-like
        The row-wise (squared) Euclidean norm of X.
    """
    if sparse.issparse(X):
        if not isinstance(X, sparse.csr_matrix):
            X = sparse.csr_matrix(X)
        norms = csr_row_norms(X)
    else:
        norms = np.einsum('ij,ij->i', X, X)

    if not squared:
        np.sqrt(norms, norms)
    return norms


def fast_logdet(A):
    """Compute log(det(A)) for A symmetric.

    Equivalent to : np.log(nl.det(A)) but more robust.
    It returns -Inf if det(A) is non positive or is not defined.

    Parameters
    ----------
    A : array-like
        The matrix.
    """
    sign, ld = np.linalg.slogdet(A)
    if not sign > 0:
        return -np.inf
    return ld


def density(w, **kwargs):
    """Compute density of a sparse vector.

    Parameters
    ----------
    w : array-like
        The sparse vector.

    Returns
    -------
    float
        The density of w, between 0 and 1.
    """
    if hasattr(w, "toarray"):
        d = float(w.nnz) / (w.shape[0] * w.shape[1])
    else:
        d = 0 if w is None else float((w != 0).sum()) / w.size
    return d


@_deprecate_positional_args
def safe_sparse_dot(a, b, *, dense_output=False):
    """Dot product that handle the sparse matrix case correctly.

    Parameters
    ----------
    a : {ndarray, sparse matrix}
    b : {ndarray, sparse matrix}
    dense_output : bool, default=False
        When False, ``a`` and ``b`` both being sparse will yield sparse output.
        When True, output will always be a dense array.

    Returns
    -------
    dot_product : {ndarray, sparse matrix}
        Sparse if ``a`` and ``b`` are sparse and ``dense_output=False``.
    """
    if a.ndim > 2 or b.ndim > 2:
        if sparse.issparse(a):
            # sparse is always 2D. Implies b is 3D+
            # [i, j] @ [k, ..., l, m, n] -> [i, k, ..., l, n]
            b_ = np.rollaxis(b, -2)
            b_2d = b_.reshape((b.shape[-2], -1))
            ret = a @ b_2d
            ret = ret.reshape(a.shape[0], *b_.shape[1:])
        elif sparse.issparse(b):
            # sparse is always 2D. Implies a is 3D+
            # [k, ..., l, m] @ [i, j] -> [k, ..., l, j]
            a_2d = a.reshape(-1, a.shape[-1])
            ret = a_2d @ b
            ret = ret.reshape(*a.shape[:-1], b.shape[1])
        else:
            ret = np.dot(a, b)
    else:
        ret = a @ b

    if (sparse.issparse(a) and sparse.issparse(b)
            and dense_output and hasattr(ret, "toarray")):
        return ret.toarray()
    return ret


@_deprecate_positional_args
def randomized_range_finder(A, *, size, n_iter,
                            power_iteration_normalizer='auto',
                            random_state=None):
    """Computes an orthonormal matrix whose range approximates the range of A.

    Parameters
    ----------
    A : 2D array
        The input data matrix.

    size : int
        Size of the return array.

    n_iter : int
        Number of power iterations used to stabilize the result.

    power_iteration_normalizer : {'auto', 'QR', 'LU', 'none'}, default='auto'
        Whether the power iterations are normalized with step-by-step
        QR factorization (the slowest but most accurate), 'none'
        (the fastest but numerically unstable when `n_iter` is large, e.g.
        typically 5 or larger), or 'LU' factorization (numerically stable
        but can lose slightly in accuracy). The 'auto' mode applies no
        normalization if `n_iter` <= 2 and switches to LU otherwise.

        .. versionadded:: 0.18

    random_state : int, RandomState instance or None, default=None
        The seed of the pseudo random number generator to use when shuffling
        the data, i.e. getting the random vectors to initialize the algorithm.
        Pass an int for reproducible results across multiple function calls.
        See :term:`Glossary <random_state>`.

    Returns
    -------
    Q : ndarray
        A (size x size) projection matrix, the range of which
        approximates well the range of the input matrix A.

    Notes
    -----

    Follows Algorithm 4.3 of
    Finding structure with randomness: Stochastic algorithms for constructing
    approximate matrix decompositions
    Halko, et al., 2009 (arXiv:909) https://arxiv.org/pdf/0909.4061.pdf

    An implementation of a randomized algorithm for principal component
    analysis
    A. Szlam et al. 2014
    """
    random_state = check_random_state(random_state)

    # Generating normal random vectors with shape: (A.shape[1], size)
    Q = random_state.normal(size=(A.shape[1], size))
    if A.dtype.kind == 'f':
        # Ensure f32 is preserved as f32
        Q = Q.astype(A.dtype, copy=False)

    # Deal with "auto" mode
    if power_iteration_normalizer == 'auto':
        if n_iter <= 2:
            power_iteration_normalizer = 'none'
        else:
            power_iteration_normalizer = 'LU'

    # Perform power iterations with Q to further 'imprint' the top
    # singular vectors of A in Q
    for i in range(n_iter):
        if power_iteration_normalizer == 'none':
            Q = safe_sparse_dot(A, Q)
            Q = safe_sparse_dot(A.T, Q)
        elif power_iteration_normalizer == 'LU':
            Q, _ = linalg.lu(safe_sparse_dot(A, Q), permute_l=True)
            Q, _ = linalg.lu(safe_sparse_dot(A.T, Q), permute_l=True)
        elif power_iteration_normalizer == 'QR':
            Q, _ = linalg.qr(safe_sparse_dot(A, Q), mode='economic')
            Q, _ = linalg.qr(safe_sparse_dot(A.T, Q), mode='economic')

    # Sample the range of A using by linear projection of Q
    # Extract an orthonormal basis
    Q, _ = linalg.qr(safe_sparse_dot(A, Q), mode='economic')
    return Q


@_deprecate_positional_args
def randomized_svd(M, n_components, *, n_oversamples=10, n_iter='auto',
                   power_iteration_normalizer='auto', transpose='auto',
                   flip_sign=True, random_state=0):
    """Computes a truncated randomized SVD.

    Parameters
    ----------
    M : {ndarray, sparse matrix}
        Matrix to decompose.

    n_components : int
        Number of singular values and vectors to extract.

    n_oversamples : int, default=10
        Additional number of random vectors to sample the range of M so as
        to ensure proper conditioning. The total number of random vectors
        used to find the range of M is n_components + n_oversamples. Smaller
        number can improve speed but can negatively impact the quality of
        approximation of singular vectors and singular values.

    n_iter : int or 'auto', default='auto'
        Number of power iterations. It can be used to deal with very noisy
        problems. When 'auto', it is set to 4, unless `n_components` is small
        (< .1 * min(X.shape)) `n_iter` in which case is set to 7.
        This improves precision with few components.

        .. versionchanged:: 0.18

    power_iteration_normalizer : {'auto', 'QR', 'LU', 'none'}, default='auto'
        Whether the power iterations are normalized with step-by-step
        QR factorization (the slowest but most accurate), 'none'
        (the fastest but numerically unstable when `n_iter` is large, e.g.
        typically 5 or larger), or 'LU' factorization (numerically stable
        but can lose slightly in accuracy). The 'auto' mode applies no
        normalization if `n_iter` <= 2 and switches to LU otherwise.

        .. versionadded:: 0.18

    transpose : bool or 'auto', default='auto'
        Whether the algorithm should be applied to M.T instead of M. The
        result should approximately be the same. The 'auto' mode will
        trigger the transposition if M.shape[1] > M.shape[0] since this
        implementation of randomized SVD tend to be a little faster in that
        case.

        .. versionchanged:: 0.18

    flip_sign : bool, default=True
        The output of a singular value decomposition is only unique up to a
        permutation of the signs of the singular vectors. If `flip_sign` is
        set to `True`, the sign ambiguity is resolved by making the largest
        loadings for each component in the left singular vectors positive.

    random_state : int, RandomState instance or None, default=0
        The seed of the pseudo random number generator to use when shuffling
        the data, i.e. getting the random vectors to initialize the algorithm.
        Pass an int for reproducible results across multiple function calls.
        See :term:`Glossary <random_state>`.

    Notes
    -----
    This algorithm finds a (usually very good) approximate truncated
    singular value decomposition using randomization to speed up the
    computations. It is particularly fast on large matrices on which
    you wish to extract only a small number of components. In order to
    obtain further speed up, `n_iter` can be set <=2 (at the cost of
    loss of precision).

    References
    ----------
    * Finding structure with randomness: Stochastic algorithms for constructing
      approximate matrix decompositions
      Halko, et al., 2009 https://arxiv.org/abs/0909.4061

    * A randomized algorithm for the decomposition of matrices
      Per-Gunnar Martinsson, Vladimir Rokhlin and Mark Tygert

    * An implementation of a randomized algorithm for principal component
      analysis
      A. Szlam et al. 2014
    """
    if isinstance(M, (sparse.lil_matrix, sparse.dok_matrix)):
        warnings.warn("Calculating SVD of a {} is expensive. "
                      "csr_matrix is more efficient.".format(
                          type(M).__name__),
                      sparse.SparseEfficiencyWarning)

    random_state = check_random_state(random_state)
    n_random = n_components + n_oversamples
    n_samples, n_features = M.shape

    if n_iter == 'auto':
        # Checks if the number of iterations is explicitly specified
        # Adjust n_iter. 7 was found a good compromise for PCA. See #5299
        n_iter = 7 if n_components < .1 * min(M.shape) else 4

    if transpose == 'auto':
        transpose = n_samples < n_features
    if transpose:
        # this implementation is a bit faster with smaller shape[1]
        M = M.T

    Q = randomized_range_finder(
        M, size=n_random, n_iter=n_iter,
        power_iteration_normalizer=power_iteration_normalizer,
        random_state=random_state)

    # project M to the (k + p) dimensional space using the basis vectors
    B = safe_sparse_dot(Q.T, M)

    # compute the SVD on the thin matrix: (k + p) wide
    Uhat, s, Vt = linalg.svd(B, full_matrices=False)

    del B
    U = np.dot(Q, Uhat)

    if flip_sign:
        if not transpose:
            U, Vt = svd_flip(U, Vt)
        else:
            # In case of transpose u_based_decision=false
            # to actually flip based on u and not v.
            U, Vt = svd_flip(U, Vt, u_based_decision=False)

    if transpose:
        # transpose back the results according to the input convention
        return Vt[:n_components, :].T, s[:n_components], U[:, :n_components].T
    else:
        return U[:, :n_components], s[:n_components], Vt[:n_components, :]


@_deprecate_positional_args
def weighted_mode(a, w, *, axis=0):
    """Returns an array of the weighted modal (most common) value in a.

    If there is more than one such value, only the first is returned.
    The bin-count for the modal bins is also returned.

    This is an extension of the algorithm in scipy.stats.mode.

    Parameters
    ----------
    a : array-like
        n-dimensional array of which to find mode(s).
    w : array-like
        n-dimensional array of weights for each value.
    axis : int, default=0
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
    (array([4.]), array([3.]))

    The value 4 appears three times: with uniform weights, the result is
    simply the mode of the distribution.

    >>> weights = [1, 3, 0.5, 1.5, 1, 2]  # deweight the 4's
    >>> weighted_mode(x, weights)
    (array([2.]), array([3.5]))

    The value 2 has the highest score: it appears twice with weights of
    1.5 and 2: the sum of these is 3.5.

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

    if a.shape != w.shape:
        w = np.full(a.shape, w, dtype=w.dtype)

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


def cartesian(arrays, out=None):
    """Generate a cartesian product of input arrays.

    Parameters
    ----------
    arrays : list of array-like
        1-D arrays to form the cartesian product of.
    out : ndarray, default=None
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

    Notes
    -----
    This function may not be used on more than 32 arrays
    because the underlying numpy functions do not support it.
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
    u : ndarray
        u and v are the output of `linalg.svd` or
        :func:`~sklearn.utils.extmath.randomized_svd`, with matching inner
        dimensions so one can compute `np.dot(u * s, v)`.

    v : ndarray
        u and v are the output of `linalg.svd` or
        :func:`~sklearn.utils.extmath.randomized_svd`, with matching inner
        dimensions so one can compute `np.dot(u * s, v)`.
        The input v should really be called vt to be consistent with scipy's
        ouput.

    u_based_decision : bool, default=True
        If True, use the columns of u as the basis for sign flipping.
        Otherwise, use the rows of v. The choice of which variable to base the
        decision on is generally algorithm dependent.


    Returns
    -------
    u_adjusted, v_adjusted : arrays with the same dimensions as the input.

    """
    if u_based_decision:
        # columns of u, rows of v
        max_abs_cols = np.argmax(np.abs(u), axis=0)
        signs = np.sign(u[max_abs_cols, range(u.shape[1])])
        u *= signs
        v *= signs[:, np.newaxis]
    else:
        # rows of v, columns of u
        max_abs_rows = np.argmax(np.abs(v), axis=1)
        signs = np.sign(v[range(v.shape[0]), max_abs_rows])
        u *= signs
        v *= signs[:, np.newaxis]
    return u, v


def log_logistic(X, out=None):
    """Compute the log of the logistic function, ``log(1 / (1 + e ** -x))``.

    This implementation is numerically stable because it splits positive and
    negative values::

        -log(1 + exp(-x_i))     if x_i > 0
        x_i - log(1 + exp(x_i)) if x_i <= 0

    For the ordinary logistic function, use ``scipy.special.expit``.

    Parameters
    ----------
    X : array-like of shape (M, N) or (M,)
        Argument to the logistic function.

    out : array-like of shape (M, N) or (M,), default=None
        Preallocated output array.

    Returns
    -------
    out : ndarray of shape (M, N) or (M,)
        Log of the logistic function evaluated at every point in x.

    Notes
    -----
    See the blog post describing this implementation:
    http://fa.bianp.net/blog/2013/numerical-optimizers-for-logistic-regression/
    """
    is_1d = X.ndim == 1
    X = np.atleast_2d(X)
    X = check_array(X, dtype=np.float64)

    n_samples, n_features = X.shape

    if out is None:
        out = np.empty_like(X)

    _log_logistic_sigmoid(n_samples, n_features, X, out)

    if is_1d:
        return np.squeeze(out)
    return out


def softmax(X, copy=True):
    """
    Calculate the softmax function.

    The softmax function is calculated by
    np.exp(X) / np.sum(np.exp(X), axis=1)

    This will cause overflow when large values are exponentiated.
    Hence the largest value in each row is subtracted from each data
    point to prevent this.

    Parameters
    ----------
    X : array-like of float of shape (M, N)
        Argument to the logistic function.

    copy : bool, default=True
        Copy X or not.

    Returns
    -------
    out : ndarray of shape (M, N)
        Softmax function evaluated at every point in x.
    """
    if copy:
        X = np.copy(X)
    max_prob = np.max(X, axis=1).reshape((-1, 1))
    X -= max_prob
    np.exp(X, X)
    sum_prob = np.sum(X, axis=1).reshape((-1, 1))
    X /= sum_prob
    return X


def make_nonnegative(X, min_value=0):
    """Ensure `X.min()` >= `min_value`.

    Parameters
    ----------
    X : array-like
        The matrix to make non-negative.
    min_value : float, default=0
        The threshold value.

    Returns
    -------
    array-like
        The thresholded array.

    Raises
    ------
    ValueError
        When X is sparse.
    """
    min_ = X.min()
    if min_ < min_value:
        if sparse.issparse(X):
            raise ValueError("Cannot make the data matrix"
                             " nonnegative because it is sparse."
                             " Adding a value to every entry would"
                             " make it no longer sparse.")
        X = X + (min_value - min_)
    return X


# Use at least float64 for the accumulating functions to avoid precision issue
# see https://github.com/numpy/numpy/issues/9393. The float64 is also retained
# as it is in case the float overflows
def _safe_accumulator_op(op, x, *args, **kwargs):
    """
    This function provides numpy accumulator functions with a float64 dtype
    when used on a floating point input. This prevents accumulator overflow on
    smaller floating point dtypes.

    Parameters
    ----------
    op : function
        A numpy accumulator function such as np.mean or np.sum.
    x : ndarray
        A numpy array to apply the accumulator function.
    *args : positional arguments
        Positional arguments passed to the accumulator function after the
        input x.
    **kwargs : keyword arguments
        Keyword arguments passed to the accumulator function.

    Returns
    -------
    result
        The output of the accumulator function passed to this function.
    """
    if np.issubdtype(x.dtype, np.floating) and x.dtype.itemsize < 8:
        result = op(x, *args, **kwargs, dtype=np.float64)
    else:
        result = op(x, *args, **kwargs)
    return result


def _incremental_weighted_mean_and_var(X, sample_weight,
                                       last_mean,
                                       last_variance,
                                       last_weight_sum):
    """Calculate weighted mean and weighted variance incremental update.

    .. versionadded:: 0.24

    Parameters
    ----------
    X : array-like of shape (n_samples, n_features)
        Data to use for mean and variance update.

    sample_weight : array-like of shape (n_samples,) or None
        Sample weights. If None, then samples are equally weighted.

    last_mean : array-like of shape (n_features,)
        Mean before the incremental update.

    last_variance : array-like of shape (n_features,) or None
        Variance before the incremental update.
        If None, variance update is not computed (in case scaling is not
        required).

    last_weight_sum : array-like of shape (n_features,)
        Sum of weights before the incremental update.

    Returns
    -------
    updated_mean : array of shape (n_features,)

    updated_variance : array of shape (n_features,) or None
        If None, only mean is computed.

    updated_weight_sum : array of shape (n_features,)

    Notes
    -----
    NaNs in `X` are ignored.

    `last_mean` and `last_variance` are statistics computed at the last step
    by the function. Both must be initialized to 0.0.
    The mean is always required (`last_mean`) and returned (`updated_mean`),
    whereas the variance can be None (`last_variance` and `updated_variance`).

    For further details on the algorithm to perform the computation in a
    numerically stable way, see [Finch2009]_, Sections 4 and 5.

    References
    ----------
    .. [Finch2009] `Tony Finch,
       "Incremental calculation of weighted mean and variance",
       University of Cambridge Computing Service, February 2009.
       <https://fanf2.user.srcf.net/hermes/doc/antiforgery/stats.pdf>`_

    """
    # last = stats before the increment
    # new = the current increment
    # updated = the aggregated stats
    if sample_weight is None:
        return _incremental_mean_and_var(X, last_mean, last_variance,
                                         last_weight_sum)
    nan_mask = np.isnan(X)
    sample_weight_T = np.reshape(sample_weight, (1, -1))
    # new_weight_sum with shape (n_features,)
    new_weight_sum = np.dot(sample_weight_T,
                            ~nan_mask).ravel().astype(np.float64)
    total_weight_sum = _safe_accumulator_op(np.sum, sample_weight, axis=0)

    X_0 = np.where(nan_mask, 0, X)
    new_mean = np.average(X_0,
                          weights=sample_weight, axis=0).astype(np.float64)
    new_mean *= total_weight_sum / new_weight_sum
    updated_weight_sum = last_weight_sum + new_weight_sum
    updated_mean = (
            (last_weight_sum * last_mean + new_weight_sum * new_mean)
            / updated_weight_sum)

    if last_variance is None:
        updated_variance = None
    else:
        X_0 = np.where(nan_mask, 0, (X-new_mean)**2)
        new_variance =\
            _safe_accumulator_op(
                np.average, X_0, weights=sample_weight, axis=0)
        new_variance *= total_weight_sum / new_weight_sum
        new_term = (
                new_weight_sum *
                (new_variance +
                 (new_mean - updated_mean) ** 2))
        last_term = (
                last_weight_sum *
                (last_variance +
                 (last_mean - updated_mean) ** 2))
        updated_variance = (new_term + last_term) / updated_weight_sum

    return updated_mean, updated_variance, updated_weight_sum


def _incremental_mean_and_var(X, last_mean, last_variance, last_sample_count):
    """Calculate mean update and a Youngs and Cramer variance update.

    last_mean and last_variance are statistics computed at the last step by the
    function. Both must be initialized to 0.0. In case no scaling is required
    last_variance can be None. The mean is always required and returned because
    necessary for the calculation of the variance. last_n_samples_seen is the
    number of samples encountered until now.

    From the paper "Algorithms for computing the sample variance: analysis and
    recommendations", by Chan, Golub, and LeVeque.

    Parameters
    ----------
    X : array-like of shape (n_samples, n_features)
        Data to use for variance update.

    last_mean : array-like of shape (n_features,)

    last_variance : array-like of shape (n_features,)

    last_sample_count : array-like of shape (n_features,)

    Returns
    -------
    updated_mean : ndarray of shape (n_features,)

    updated_variance : ndarray of shape (n_features,)
        If None, only mean is computed.

    updated_sample_count : ndarray of shape (n_features,)

    Notes
    -----
    NaNs are ignored during the algorithm.

    References
    ----------
    T. Chan, G. Golub, R. LeVeque. Algorithms for computing the sample
        variance: recommendations, The American Statistician, Vol. 37, No. 3,
        pp. 242-247

    Also, see the sparse implementation of this in
    `utils.sparsefuncs.incr_mean_variance_axis` and
    `utils.sparsefuncs_fast.incr_mean_variance_axis0`
    """
    # old = stats until now
    # new = the current increment
    # updated = the aggregated stats
    last_sum = last_mean * last_sample_count
    new_sum = _safe_accumulator_op(np.nansum, X, axis=0)

    new_sample_count = np.sum(~np.isnan(X), axis=0)
    updated_sample_count = last_sample_count + new_sample_count

    updated_mean = (last_sum + new_sum) / updated_sample_count

    if last_variance is None:
        updated_variance = None
    else:
        new_unnormalized_variance = (
            _safe_accumulator_op(np.nanvar, X, axis=0) * new_sample_count)
        last_unnormalized_variance = last_variance * last_sample_count

        with np.errstate(divide='ignore', invalid='ignore'):
            last_over_new_count = last_sample_count / new_sample_count
            updated_unnormalized_variance = (
                last_unnormalized_variance + new_unnormalized_variance +
                last_over_new_count / updated_sample_count *
                (last_sum / last_over_new_count - new_sum) ** 2)

        zeros = last_sample_count == 0
        updated_unnormalized_variance[zeros] = new_unnormalized_variance[zeros]
        updated_variance = updated_unnormalized_variance / updated_sample_count

    return updated_mean, updated_variance, updated_sample_count


def _deterministic_vector_sign_flip(u):
    """Modify the sign of vectors for reproducibility.

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


def stable_cumsum(arr, axis=None, rtol=1e-05, atol=1e-08):
    """Use high precision for cumsum and check that final value matches sum.

    Parameters
    ----------
    arr : array-like
        To be cumulatively summed as flat.
    axis : int, default=None
        Axis along which the cumulative sum is computed.
        The default (None) is to compute the cumsum over the flattened array.
    rtol : float, default=1e-05
        Relative tolerance, see ``np.allclose``.
    atol : float, default=1e-08
        Absolute tolerance, see ``np.allclose``.
    """
    out = np.cumsum(arr, axis=axis, dtype=np.float64)
    expected = np.sum(arr, axis=axis, dtype=np.float64)
    if not np.all(np.isclose(out.take(-1, axis=axis), expected, rtol=rtol,
                             atol=atol, equal_nan=True)):
        warnings.warn('cumsum was found to be unstable: '
                      'its last element does not correspond to sum',
                      RuntimeWarning)
    return out
