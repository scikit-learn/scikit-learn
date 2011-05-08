
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


def inplace_row_normalize(X, norm=1):
    """Inplace row normalization of the 2D array or scipy.sparse matrix

    Parameters
    ----------
      X: 2D array-like or scipy sparse matrix

      norm: 1 or 2 (default is 1)
        the norm to use for row normalization
    """
    if norm not in (1, 2):
        raise ValueError("%r is not a supported norm" % norm)

    if hasattr(X, 'tocsr'):
        # inplace row normalization for sparse matrices

        # TODO: compare the speed with the existing cython implementation
        # available in the preprocessing package

        # lazy import of scipy.sparse for performance
        from scipy.sparse.sparsetools import csr_scale_rows

        # ensure the sparse matrix is in Compressed Sparse Rows format
        X = X.tocsr()

        if norm == 1:
            norms = np.abs(X).sum(axis=1).A.flatten()
        elif norm == 2:
            norms = np.sqrt(X.multiply(X).sum(axis=1)).A.flatten()
        norms[norms == 0.0] = 1.0

        # inplace rescaling of the CSR matrix
        csr_scale_rows(X.shape[0], X.shape[1], X.indptr, X.indices, X.data,
                       1.0 / norms)
    else:
        # in-place row normalization for ndarray
        if norm == 1:
            norms = np.abs(X).sum(axis=1)
        elif norm == 2:
            norms = np.sqrt((X ** 2).sum(axis=1))
        norms[norms == 0.0] = 1.0
        X /= norms[:, np.newaxis]

    return X
