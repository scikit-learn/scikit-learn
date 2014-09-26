"""
The :mod:`sklearn.utils` module includes various utilities.
"""
from collections import Sequence

import numpy as np
from scipy.sparse import issparse
import warnings

from .murmurhash import murmurhash3_32
from .validation import (as_float_array,
                         assert_all_finite, warn_if_not_float,
                         check_random_state, column_or_1d, check_array,
                         check_consistent_length, check_X_y, indexable)
from .class_weight import compute_class_weight


__all__ = ["murmurhash3_32", "as_float_array",
           "assert_all_finite", "check_array",
           "warn_if_not_float",
           "check_random_state",
           "compute_class_weight",
           "column_or_1d", "safe_indexing",
           "check_consistent_length", "check_X_y", 'indexable']


class deprecated(object):
    """Decorator to mark a function or class as deprecated.

    Issue a warning when the function is called/the class is instantiated and
    adds a warning to the docstring.

    The optional extra argument will be appended to the deprecation message
    and the docstring. Note: to use this with the default value for extra, put
    in an empty of parentheses:

    >>> from sklearn.utils import deprecated
    >>> deprecated() # doctest: +ELLIPSIS
    <sklearn.utils.deprecated object at ...>

    >>> @deprecated()
    ... def some_function(): pass
    """

    # Adapted from http://wiki.python.org/moin/PythonDecoratorLibrary,
    # but with many changes.

    def __init__(self, extra=''):
        """
        Parameters
        ----------
        extra: string
          to be added to the deprecation messages

        """
        self.extra = extra

    def __call__(self, obj):
        if isinstance(obj, type):
            return self._decorate_class(obj)
        else:
            return self._decorate_fun(obj)

    def _decorate_class(self, cls):
        msg = "Class %s is deprecated" % cls.__name__
        if self.extra:
            msg += "; %s" % self.extra

        # FIXME: we should probably reset __new__ for full generality
        init = cls.__init__

        def wrapped(*args, **kwargs):
            warnings.warn(msg, category=DeprecationWarning)
            return init(*args, **kwargs)
        cls.__init__ = wrapped

        wrapped.__name__ = '__init__'
        wrapped.__doc__ = self._update_doc(init.__doc__)
        wrapped.deprecated_original = init

        return cls

    def _decorate_fun(self, fun):
        """Decorate function fun"""

        msg = "Function %s is deprecated" % fun.__name__
        if self.extra:
            msg += "; %s" % self.extra

        def wrapped(*args, **kwargs):
            warnings.warn(msg, category=DeprecationWarning)
            return fun(*args, **kwargs)

        wrapped.__name__ = fun.__name__
        wrapped.__dict__ = fun.__dict__
        wrapped.__doc__ = self._update_doc(fun.__doc__)

        return wrapped

    def _update_doc(self, olddoc):
        newdoc = "DEPRECATED"
        if self.extra:
            newdoc = "%s: %s" % (newdoc, self.extra)
        if olddoc:
            newdoc = "%s\n\n%s" % (newdoc, olddoc)
        return newdoc


def safe_mask(X, mask):
    """Return a mask which is safe to use on X.

    Parameters
    ----------
    X : {array-like, sparse matrix}
        Data on which to apply mask.

    mask: array
        Mask to be used on X.

    Returns
    -------
        mask
    """
    mask = np.asarray(mask)
    if np.issubdtype(mask.dtype, np.int):
        return mask

    if hasattr(X, "toarray"):
        ind = np.arange(mask.shape[0])
        mask = ind[mask]
    return mask


def safe_indexing(X, indices):
    """Return items or rows from X using indices.

    Allows simple indexing of lists or arrays.

    Parameters
    ----------
    X : array-like, sparse-matrix, list.
        Data from which to sample rows or items.

    indices : array-like, list
        Indices according to which X will be subsampled.
    """
    if hasattr(X, "iloc"):
        # Pandas Dataframes and Series
        return X.iloc[indices]
    elif hasattr(X, "shape"):
        if hasattr(X, 'take') and (hasattr(indices, 'dtype') and
                                   indices.dtype.kind == 'i'):
            # This is often substantially faster than X[indices]
            return X.take(indices, axis=0)
        else:
            return X[indices]
    else:
        return [X[idx] for idx in indices]


def resample(*arrays, **options):
    """Resample arrays or sparse matrices in a consistent way

    The default strategy implements one step of the bootstrapping
    procedure.

    Parameters
    ----------
    `*arrays` : sequence of arrays or scipy.sparse matrices with same shape[0]

    replace : boolean, True by default
        Implements resampling with replacement. If False, this will implement
        (sliced) random permutations.

    n_samples : int, None by default
        Number of samples to generate. If left to None this is
        automatically set to the first dimension of the arrays.

    random_state : int or RandomState instance
        Control the shuffling for reproducible behavior.

    Returns
    -------
    Sequence of resampled views of the collections. The original arrays are
    not impacted.

    Examples
    --------
    It is possible to mix sparse and dense arrays in the same run::

      >>> X = [[1., 0.], [2., 1.], [0., 0.]]
      >>> y = np.array([0, 1, 2])

      >>> from scipy.sparse import coo_matrix
      >>> X_sparse = coo_matrix(X)

      >>> from sklearn.utils import resample
      >>> X, X_sparse, y = resample(X, X_sparse, y, random_state=0)
      >>> X
      array([[ 1.,  0.],
             [ 2.,  1.],
             [ 1.,  0.]])

      >>> X_sparse                   # doctest: +ELLIPSIS +NORMALIZE_WHITESPACE
      <3x2 sparse matrix of type '<... 'numpy.float64'>'
          with 4 stored elements in Compressed Sparse Row format>

      >>> X_sparse.toarray()
      array([[ 1.,  0.],
             [ 2.,  1.],
             [ 1.,  0.]])

      >>> y
      array([0, 1, 0])

      >>> resample(y, n_samples=2, random_state=0)
      array([0, 1])


    See also
    --------
    :class:`sklearn.cross_validation.Bootstrap`
    :func:`sklearn.utils.shuffle`
    """
    random_state = check_random_state(options.pop('random_state', None))
    replace = options.pop('replace', True)
    max_n_samples = options.pop('n_samples', None)
    if options:
        raise ValueError("Unexpected kw arguments: %r" % options.keys())

    if len(arrays) == 0:
        return None

    first = arrays[0]
    n_samples = first.shape[0] if hasattr(first, 'shape') else len(first)

    if max_n_samples is None:
        max_n_samples = n_samples

    if max_n_samples > n_samples:
        raise ValueError("Cannot sample %d out of arrays with dim %d" % (
            max_n_samples, n_samples))

    check_consistent_length(*arrays)
    arrays = [check_array(x, accept_sparse='csr', ensure_2d=False,
                          allow_nd=True) for x in arrays]

    if replace:
        indices = random_state.randint(0, n_samples, size=(max_n_samples,))
    else:
        indices = np.arange(n_samples)
        random_state.shuffle(indices)
        indices = indices[:max_n_samples]

    resampled_arrays = []

    for array in arrays:
        array = array[indices]
        resampled_arrays.append(array)

    if len(resampled_arrays) == 1:
        # syntactic sugar for the unit argument case
        return resampled_arrays[0]
    else:
        return resampled_arrays


def shuffle(*arrays, **options):
    """Shuffle arrays or sparse matrices in a consistent way

    This is a convenience alias to ``resample(*arrays, replace=False)`` to do
    random permutations of the collections.

    Parameters
    ----------
    `*arrays` : sequence of arrays or scipy.sparse matrices with same shape[0]

    random_state : int or RandomState instance
        Control the shuffling for reproducible behavior.

    n_samples : int, None by default
        Number of samples to generate. If left to None this is
        automatically set to the first dimension of the arrays.

    Returns
    -------
    Sequence of shuffled views of the collections. The original arrays are
    not impacted.

    Examples
    --------
    It is possible to mix sparse and dense arrays in the same run::

      >>> X = [[1., 0.], [2., 1.], [0., 0.]]
      >>> y = np.array([0, 1, 2])

      >>> from scipy.sparse import coo_matrix
      >>> X_sparse = coo_matrix(X)

      >>> from sklearn.utils import shuffle
      >>> X, X_sparse, y = shuffle(X, X_sparse, y, random_state=0)
      >>> X
      array([[ 0.,  0.],
             [ 2.,  1.],
             [ 1.,  0.]])

      >>> X_sparse                   # doctest: +ELLIPSIS +NORMALIZE_WHITESPACE
      <3x2 sparse matrix of type '<... 'numpy.float64'>'
          with 3 stored elements in Compressed Sparse Row format>

      >>> X_sparse.toarray()
      array([[ 0.,  0.],
             [ 2.,  1.],
             [ 1.,  0.]])

      >>> y
      array([2, 1, 0])

      >>> shuffle(y, n_samples=2, random_state=0)
      array([0, 1])

    See also
    --------
    :func:`sklearn.utils.resample`
    """
    options['replace'] = False
    return resample(*arrays, **options)


def safe_sqr(X, copy=True):
    """Element wise squaring of array-likes and sparse matrices.

    Parameters
    ----------
    X : array like, matrix, sparse matrix

    Returns
    -------
    X ** 2 : element wise square
    """
    X = check_array(X, accept_sparse=['csr', 'csc', 'coo'])
    if issparse(X):
        if copy:
            X = X.copy()
        X.data **= 2
    else:
        if copy:
            X = X ** 2
        else:
            X **= 2
    return X


def gen_batches(n, batch_size):
    """Generator to create slices containing batch_size elements, from 0 to n.

    The last slice may contain less than batch_size elements, when batch_size
    does not divide n.

    Examples
    --------
    >>> from sklearn.utils import gen_batches
    >>> list(gen_batches(7, 3))
    [slice(0, 3, None), slice(3, 6, None), slice(6, 7, None)]
    >>> list(gen_batches(6, 3))
    [slice(0, 3, None), slice(3, 6, None)]
    >>> list(gen_batches(2, 3))
    [slice(0, 2, None)]
    """
    start = 0
    for _ in range(int(n // batch_size)):
        end = start + batch_size
        yield slice(start, end)
        start = end
    if start < n:
        yield slice(start, n)


def gen_even_slices(n, n_packs, n_samples=None):
    """Generator to create n_packs slices going up to n.

    Pass n_samples when the slices are to be used for sparse matrix indexing;
    slicing off-the-end raises an exception, while it works for NumPy arrays.

    Examples
    --------
    >>> from sklearn.utils import gen_even_slices
    >>> list(gen_even_slices(10, 1))
    [slice(0, 10, None)]
    >>> list(gen_even_slices(10, 10))                     #doctest: +ELLIPSIS
    [slice(0, 1, None), slice(1, 2, None), ..., slice(9, 10, None)]
    >>> list(gen_even_slices(10, 5))                      #doctest: +ELLIPSIS
    [slice(0, 2, None), slice(2, 4, None), ..., slice(8, 10, None)]
    >>> list(gen_even_slices(10, 3))
    [slice(0, 4, None), slice(4, 7, None), slice(7, 10, None)]
    """
    start = 0
    for pack_num in range(n_packs):
        this_n = n // n_packs
        if pack_num < n % n_packs:
            this_n += 1
        if this_n > 0:
            end = start + this_n
            if n_samples is not None:
                end = min(n_samples, end)
            yield slice(start, end, None)
            start = end


def tosequence(x):
    """Cast iterable x to a Sequence, avoiding a copy if possible."""
    if isinstance(x, np.ndarray):
        return np.asarray(x)
    elif isinstance(x, Sequence):
        return x
    else:
        return list(x)


class ConvergenceWarning(Warning):
    "Custom warning to capture convergence problems"
