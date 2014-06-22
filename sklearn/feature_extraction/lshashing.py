"""
Locality Sensitive Hashing Algorithms
-------------------------------------
"""
# Author: Maheshakya Wijewardena <maheshakya.10@cse.mrt.ac.lk>

import numpy as np
from abc import ABCMeta, abstractmethod
from ..externals.six import with_metaclass
from ..utils import check_random_state

__all__ = ["RandomProjections"]


class BaseHash(with_metaclass(ABCMeta)):
    """
    Base class for LSH algorithms.

    Warning: This class should not be used directly. Use derived classes
    instead.

    Parameters
    ----------
    n_dim: int
        Number of dimensions in the data set which is being hashed.

    hash_size: int
        Length of the hash

    seed: float, optional (defualt=1)
        Seed to initialize pseudo random hash functions.

    """
    @abstractmethod
    def __init__(self, n_dim=None, hash_size=None, random_state=1):
        if n_dim is None or hash_size is None:
            raise ValueError("n_dim or hash_size cannot be None.")

        self.n_dim = n_dim
        self.random_state = random_state
        self.hash_size = hash_size

    def generate_hash_function(self):
        """Generates a hash function"""

    def do_hash(self, input_point=None, hash_function=None):
        """Performs hashing on the input_point with hash_function"""


class RandomProjections(BaseHash):
    """
    Performs random projections [1] as an LSH algorithm

    References
    ----------

    .. [1] Wikipedia, "Random Projection", Avaliable[online]:
           http://en.wikipedia.org/wiki/Locality-sensitive_hashing#Random_projection
    """
    def __init__(self, n_dim, hash_size, random_state=1):
        super(RandomProjections, self).__init__(n_dim=n_dim,
                                                hash_size=hash_size,
                                                random_state=random_state)

    def generate_hash_function(self):
        """
        Generates hyperplanes of shape (hash_size, n_dim) from standard
        normal distribution.
        """
        random_state = check_random_state(self.random_state)
        return random_state.randn(self.hash_size, self.n_dim)

    def do_hash(self, input_point=None, hash_function=None):
        """
        Does hashing on the data point with the provided hash_function.
        """
        if input_point is None or hash_function is None:
            raise ValueError("input_point or hash_function cannot be None.")

        projections = np.dot(hash_function, input_point)
        return "".join(['1' if i > 0 else '0' for i in projections])
