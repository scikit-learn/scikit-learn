"""Utilities for input validation"""
# Authors: Olivier Grisel and Gael Varoquaux and others (please update me)
# License: BSD 3

import numpy as np
from scipy import sparse
import warnings

from .fixes import safe_copy

def assert_all_finite(X):
    """Throw a ValueError if X contains NaN or infinity.

    Input MUST be an np.ndarray instance or a scipy.sparse matrix."""

    # First try an O(n) time, O(1) space solution for the common case that
    # there everything is finite; fall back to O(n) space np.isfinite to
    # prevent false positives from overflow in sum method.
    if X.dtype.char in np.typecodes['AllFloat'] and not np.isfinite(X.sum()) \
      and not np.isfinite(X.data if sparse.issparse(X) else X).all():
            raise ValueError("array contains NaN or infinity")


def safe_asarray(X, dtype=None, order=None):
    """Convert X to an array or sparse matrix.

    Prevents copying X when possible; sparse matrices are passed through."""
    if not sparse.issparse(X):
        X = np.asarray(X, dtype, order)
    assert_all_finite(X)
    return X


def as_float_array(X, copy=True):
    """Converts an array-like to an array of floats

    The new dtype will be np.float32 or np.float64, depending on the original
    type. The function can create a copy or modify the argument depending
    on the argument copy.

    Parameters
    ----------
    X : {array-like, sparse matrix}

    copy : bool, optional
        If True, a copy of X will be created. If False, a copy may still be
        returned if X's dtype is not a floating point type.

    Returns
    -------
    XT : {array, sparse matrix}
        An array of type np.float
    """
    if isinstance(X, np.matrix) or (not isinstance(X, np.ndarray)
                                    and not sparse.issparse(X)):
        return safe_asarray(X, dtype=np.float64)
    elif sparse.issparse(X) and X.dtype in [np.float32, np.float64]:
        return X.copy() if copy else X
    elif X.dtype in [np.float32, np.float64]:  # is numpy array
        return X.copy('F' if X.flags['F_CONTIGUOUS'] else 'C') if copy else X
    else:
        return X.astype(np.float32 if X.dtype == np.int32 else np.float64)


def array2d(X, dtype=None, order=None, copy=False):
    """Returns at least 2-d array with data from X"""
    if sparse.issparse(X):
        raise TypeError('A sparse matrix was passed, but dense data '
                        'is required. Use X.todense() to convert to dense.')
    X_2d = np.asarray(np.atleast_2d(X), dtype=dtype, order=order)
    if X is X_2d and copy:
        X_2d = safe_copy(X_2d)
    return X_2d


def atleast2d_or_csc(X, dtype=None, order=None, copy=False):
    """Like numpy.atleast_2d, but converts sparse matrices to CSC format

    Also, converts np.matrix to np.ndarray.
    """
    if sparse.issparse(X):
        # Note: order is ignored because CSR matrices hold data in 1-d arrays
        if dtype is None or X.dtype == dtype:
            X = X.tocsc()
        else:
            X = sparse.csc_matrix(X, dtype=dtype)
    else:
        X = array2d(X, dtype=dtype, order=order, copy=copy)
    assert_all_finite(X)
    return X


def atleast2d_or_csr(X, dtype=None, order=None, copy=False):
    """Like numpy.atleast_2d, but converts sparse matrices to CSR format

    Also, converts np.matrix to np.ndarray.
    """
    if sparse.issparse(X):
        # Note: order is ignored because CSR matrices hold data in 1-d arrays
        if dtype is None or X.dtype == dtype:
            X = X.tocsr()
        else:
            X = sparse.csr_matrix(X, dtype=dtype)
    else:
        X = array2d(X, dtype=dtype, order=order, copy=copy)
    assert_all_finite(X)
    return X


def _num_samples(x):
    """Return number of samples in array-like x."""
    if not hasattr(x, '__len__') and not hasattr(x, 'shape'):
        raise TypeError("Expected sequence or array-like, got %r" % x)
    return x.shape[0] if hasattr(x, 'shape') else len(x)


def check_arrays(*arrays, **options):
    """Checked that all arrays have consistent first dimensions

    Parameters
    ----------
    *arrays : sequence of arrays or scipy.sparse matrices with same shape[0]
        Python lists or tuples occurring in arrays are converted to 1D numpy
        arrays.

    sparse_format : 'csr', 'csc' or 'dense', None by default
        If not None, any scipy.sparse matrix is converted to
        Compressed Sparse Rows or Compressed Sparse Columns representations.
        If 'dense', an error is raised when a sparse array is
        passed.

    copy : boolean, False by default
        If copy is True, ensure that returned arrays are copies of the original
        (if not already converted to another format earlier in the process).

    check_ccontiguous : boolean, False by default
        Check that the arrays are C contiguous

    dtype : a numpy dtype instance, None by default
        Enforce a specific dtype.
    """
    sparse_format = options.pop('sparse_format', None)
    if sparse_format not in (None, 'csr', 'csc', 'dense'):
        raise ValueError('Unexpected sparse format: %r' % sparse_format)
    copy = options.pop('copy', False)
    check_ccontiguous = options.pop('check_ccontiguous', False)
    dtype = options.pop('dtype', None)
    if options:
        raise TypeError("Unexpected keyword arguments: %r" % options.keys())

    if len(arrays) == 0:
        return None

    n_samples = _num_samples(arrays[0])

    checked_arrays = []
    for array in arrays:
        array_orig = array
        if array is None:
            # special case: ignore optional y=None kwarg pattern
            checked_arrays.append(array)
            continue

        size = _num_samples(array)

        if size != n_samples:
            raise ValueError("Found array with dim %d. Expected %d" % (
                size, n_samples))

        if sparse.issparse(array):
            if sparse_format == 'csr':
                array = array.tocsr()
            elif sparse_format == 'csc':
                array = array.tocsc()
            elif sparse_format == 'dense':
                raise TypeError('A sparse matrix was passed, but dense data '
                    'is required. Use X.todense() to convert to dense.')
            if check_ccontiguous:
                array.data = np.ascontiguousarray(array.data, dtype=dtype)
            else:
                array.data = np.asarray(array.data, dtype=dtype)
        else:
            if check_ccontiguous:
                array = np.ascontiguousarray(array, dtype=dtype)
            else:
                array = np.asarray(array, dtype=dtype)

        if copy and array is array_orig:
            array = array.copy()
        checked_arrays.append(array)

    return checked_arrays


def warn_if_not_float(X, estimator='This algorithm'):
    """Warning utility function to check that data type is floating point"""
    if not isinstance(estimator, basestring):
        estimator = estimator.__class__.__name__
    if X.dtype.kind != 'f':
        warnings.warn("%s assumes floating point values as input, "
                      "got %s" % (estimator, X.dtype))


def check_random_state(seed):
    """Turn seed into a np.random.RandomState instance

    If seed is None, return the RandomState singleton used by np.random.
    If seed is an int, return a new RandomState instance seeded with seed.
    If seed is already a RandomState instance, return it.
    Otherwise raise ValueError.
    """
    if seed is None or seed is np.random:
        return np.random.mtrand._rand
    if isinstance(seed, (int, np.integer)):
        return np.random.RandomState(seed)
    if isinstance(seed, np.random.RandomState):
        return seed
    raise ValueError('%r cannot be used to seed a numpy.random.RandomState'
                     ' instance' % seed)
