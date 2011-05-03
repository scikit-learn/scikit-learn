
import numpy as np
import scipy.sparse as sp

def safe_asanyarray(X, dtype=None, order=None):
    if sp.issparse(X):
        return X
        #return type(X)(X, dtype)
    else:
        return np.asanyarray(X, dtype, order)

def make_rng(seed):
    """Turn seed into a np.random.RandomState instance

    If seed is None, return the np.random singleton.
    If seed is an int, return a new RandomState instance seeded with seed.
    If seed is already a RandomState instance, return it.
    Otherwise raise ValueError.
    """
    if seed is None or seed is np.random:
        return np.random
    if isinstance(seed, int):
        return np.random.RandomState(seed)
    if isinstance(seed, np.random.RandomState):
        return seed
    raise ValueError('%r cannot be used to seed a numpy.random.RandomState'
                     ' instance' % seed)
