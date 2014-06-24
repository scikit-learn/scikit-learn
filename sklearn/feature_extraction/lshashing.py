"""
Locality Sensitive Hashing Algorithms
-------------------------------------
"""
# Author: Maheshakya Wijewardena <maheshakya.10@cse.mrt.ac.lk>

import numpy as np
from abc import ABCMeta, abstractmethod
from ..externals.six import with_metaclass
from ..utils import check_random_state

from ..random_projection import GaussianRandomProjection

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
    Performs random projections [1] as an LSH algorithm. This uses Gaussian
    Random Projections from `sklearn.random_projections` module.

    References
    ----------

    .. [1] Wikipedia, "Random Projection", Avaliable[online]:
           http://en.wikipedia.org/wiki/Locality-sensitive_hashing#Random_projection
    """
    def __init__(self, n_dim, hash_size, random_state=1):
        super(RandomProjections, self).__init__(n_dim=n_dim,
                                                hash_size=hash_size,
                                                random_state=random_state)

    def _generate_hash_function(self):
        """
        Fits a `GaussianRandomProjections` with `n_components=hash_size
        and n_features=n_dim.
        """
        random_state = check_random_state(self.random_state)
        grp = GaussianRandomProjection(n_components=self.hash_size,
                                       random_state=random_state.randint(0,
                                                                         10))
        X = np.zeros((2, self.n_dim), dtype=float)
        grp.fit(X)
        return grp

    def do_hash(self, input_array=None):
        """
        Does hashing on an array of data points.
        This creates a binary hash by getting the dot product of
        input_point and hash_function then transforming the projection
        into a binary string array based on the sign(positive/negative)
        of the projection.

        Parameters
        ----------

        input_array: array_like, shape (n_samples, n_features)
            A matrix of dimensions (n_samples, n_features), which is being
            hashed.
        """
        if input_array is None:
            raise ValueError("input_array cannot be None.")

        grp = self._generate_hash_function()
        res = grp.transform(input_array)

        bin_hashes = np.empty(res.shape[0], dtype='|S'+str(self.hash_size))
        for i in range(res.shape[0]):
            bin_hashes[i] = "".join(map(str, np.array(res[i] > 0, dtype=int)))

        return bin_hashes, grp.components_
