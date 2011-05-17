
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


def shuffle(*args, **kwargs):
    """Shuffle arrays or sparse matrices in a consistent way

    Parameters
    ----------
    *args : sequence of arrays or scipy.sparse matrices with same shape[0]

    random_state : int or RandomState instance
        Control the shuffling for reproducible behavior.

    copy : boolean (True by default)
        If False, perform inplace shuffling unless the collection is a python
        list or scipy.sparse matrix in a format other than CSR.

    Return
    ------
    Sequence of shuffled collections.

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

    """
    copy = kwargs.get('copy', True)
    random_state = check_random_state(kwargs.get('random_state'))

    if len(args) == 0:
        return args

    first = args[0]
    n_samples = first.shape[0] if hasattr(first, 'shape') else len(first)
    indices = np.arange(n_samples)
    random_state.shuffle(indices)

    shuffled = []

    for arg in args:
        if hasattr(arg, 'tocsr'):
            arg = arg.tocsr()
        else:
            arg = np.asanyarray(arg)
        if copy:
            arg = arg.copy()
        arg = arg[indices]
        shuffled.append(arg)

    return shuffled
