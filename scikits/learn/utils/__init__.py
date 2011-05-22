
import numpy as np
import scipy.sparse as sp


def safe_asanyarray(X, dtype=None, order=None):
    if sp.issparse(X):
        return X
        #return type(X)(X, dtype)
    else:
        return np.asanyarray(X, dtype, order)


def check_random_state(seed):
    """Turn seed into a np.random.RandomState instance

    If seed is None, return the RandomState singleton used by np.random.
    If seed is an int, return a new RandomState instance seeded with seed.
    If seed is already a RandomState instance, return it.
    Otherwise raise ValueError.
    """
    if seed is None or seed is np.random:
        return np.random.mtrand._rand
    if isinstance(seed, int):
        return np.random.RandomState(seed)
    if isinstance(seed, np.random.RandomState):
        return seed
    raise ValueError('%r cannot be used to seed a numpy.random.RandomState'
                     ' instance' % seed)


def check_arrays(*arrays, **options):
    """Checked that all arrays have consistent first dimensions

    Parameters
    ----------
    *arrays : sequence of arrays or scipy.sparse matrices with same shape[0]
        Python lists or tuples occurring in arrays are converted to 1D numpy
        arrays.

    force_csr : boolean, False by default
        If force_csr is True, any scipy.sparse matrix is converted to
        Compressed Sparse Row representation.
    """
    force_csr = options.pop('force_csr', False)
    if options:
        raise ValueError("Unexpected kw arguments: %r" % options.keys())

    if len(arrays) == 0:
        return None

    first = arrays[0]
    if not hasattr(first, '__len__') and not hasattr(first, 'shape'):
        raise ValueError("Expected python sequence or array, got %r" % first)
    n_samples = first.shape[0] if hasattr(first, 'shape') else len(first)

    checked_arrays = []
    for array in arrays:
        if array is None:
            # special case: ignore optional y=None kwarg pattern
            checked_arrays.append(array)
            continue

        if not hasattr(array, '__len__') and not hasattr(array, 'shape'):
            raise ValueError("Expected python sequence or array, got %r"
                             % array)
        size = array.shape[0] if hasattr(array, 'shape') else len(array)

        if size != n_samples:
            raise ValueError("Found array with dim %d. Expected %d" % (
                size, n_samples))

        if hasattr(array, 'tocsr'):
            if force_csr:
                array = array.tocsr()
        else:
            array = np.asanyarray(array)

        checked_arrays.append(array)

    return checked_arrays


def resample(*arrays, **options):
    """Resample arrays or sparse matrices in a consistent way

    The default strategy implements one step of the bootstrapping
    procedure.

    Parameters
    ----------
    *arrays : sequence of arrays or scipy.sparse matrices with same shape[0]

    replace : boolean, True by default
        Implements resampling with replacement. If False, this will implement
        (sliced) random permutations.

    n_samples : int, None by default
        Number of samples to generate. If left to None this is
        automatically set to the first dimension of the arrays.

    random_state : int or RandomState instance
        Control the shuffling for reproducible behavior.

    Return
    ------
    Sequence of resampled views of the collections. The original arrays are
    not impacted.

    Example
    -------
    It is possible to mix sparse and dense arrays in the same run::

      >>> X = [[1., 0.], [2., 1.], [0., 0.]]
      >>> y = np.array([0, 1, 2])

      >>> from scipy.sparse import coo_matrix
      >>> X_sparse = coo_matrix(X)

      >>> X, X_sparse, y = resample(X, X_sparse, y, random_state=0)
      >>> X
      array([[ 1.,  0.],
             [ 2.,  1.],
             [ 1.,  0.]])

      >>> X_sparse                            # doctest: +NORMALIZE_WHITESPACE
      <3x2 sparse matrix of type '<type 'numpy.float64'>'
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
    :class:`scikits.learn.cross_val.Bootstrap`
    :func:`scikits.learn.utils.shuffle`
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

    arrays = check_arrays(*arrays, force_csr=True)

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

    This is a convenience alias to resample(*arrays, replace=False) to do
    random permutations of the collections.

    Parameters
    ----------
    *arrays : sequence of arrays or scipy.sparse matrices with same shape[0]

    random_state : int or RandomState instance
        Control the shuffling for reproducible behavior.

    n_samples : int, None by default
        Number of samples to generate. If left to None this is
        automatically set to the first dimension of the arrays.

    Return
    ------
    Sequence of shuffled views of the collections. The original arrays are
    not impacted.

    Example
    -------
    It is possible to mix sparse and dense arrays in the same run::

      >>> X = [[1., 0.], [2., 1.], [0., 0.]]
      >>> y = np.array([0, 1, 2])

      >>> from scipy.sparse import coo_matrix
      >>> X_sparse = coo_matrix(X)

      >>> X, X_sparse, y = shuffle(X, X_sparse, y, random_state=0)
      >>> X
      array([[ 0.,  0.],
             [ 2.,  1.],
             [ 1.,  0.]])

      >>> X_sparse                            # doctest: +NORMALIZE_WHITESPACE
      <3x2 sparse matrix of type '<type 'numpy.float64'>'
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
    :func:`scikits.learn.utils.resample`
    """
    options['replace'] = False
    return resample(*arrays, **options)
