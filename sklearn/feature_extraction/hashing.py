# Author: Lars Buitinck <L.J.Buitinck@uva.nl>
# License: 3-clause BSD.

import itertools
import numbers

import numpy as np
import scipy.sparse as sp

from . import _hashing
from ..base import BaseEstimator, TransformerMixin


def _iteritems(d):
    """Like d.iteritems, but accepts any collections.Mapping."""
    return d.iteritems() if hasattr(d, "iteritems") else d.items()


class FeatureHasher(BaseEstimator, TransformerMixin):
    """Implements feature hashing, aka the hashing trick.

    This class turns sequences of symbolic feature names (strings) into
    scipy.sparse matrices, using a hash function to compute the matrix column
    corresponding to a name. The hash function employed is the signed 32-bit
    version of Murmurhash3.

    Feature names of type byte string are used as-is. Unicode strings are
    converted to UTF-8 first, but no Unicode normalization is done.

    This class is a low-memory alternative to DictVectorizer and
    CountVectorizer, intended for large-scale (online) learning and situations
    where memory is tight, e.g. when running prediction code on embedded
    devices.

    Parameters
    ----------
    n_features : integer, optional
        The number of features (columns) in the output matrices. Small numbers
        of features are likely to cause hash collisions, but large numbers
        will cause larger coefficient dimensions in linear learners.
    dtype : NumPy type, optional
        The type of feature values. Passed to scipy.sparse matrix constructors
        as the dtype argument. Do not set this to bool, np.boolean or any
        unsigned integer type.
    input_type : string, optional
        Either "dict" (the default) to accept dictionaries over
        (feature_name, value); "pair" to accept pairs of (feature_name, value);
        or "string" to accept single strings.
        feature_name should be a string, while value should be a number.
        In the case of "string", a value of 1 is implied.
        The feature_name is hashed to find the appropriate column for the
        feature. The value's sign might be flipped in the output (but see
        non_negative, below).
    non_negative : boolean, optional
        Whether output matrices should contain non-negative values only;
        effectively calls abs on the matrix prior to returning it.
        When True, output values will be multinomially distributed.
        When False, output values will be normally distributed (Gaussian) with
        mean 0, assuming a good hash function.

    """

    def __init__(self, n_features=(2 ** 20), input_type="dict",
                 dtype=np.float64, non_negative=False):
        self._validate_params(n_features, input_type)

        self.dtype = dtype
        self.input_type = input_type
        self.n_features = n_features
        self.non_negative = non_negative

    @staticmethod
    def _validate_params(n_features, input_type):
        # strangely, np.int16 instances are not instances of Integral,
        # while np.int64 instances are...
        if not isinstance(n_features, (numbers.Integral, np.integer)):
            raise TypeError("n_features must be integral, got %r (%s)."
                            % (n_features, type(n_features)))
        elif n_features < 1 or n_features >= 2**31:
            raise ValueError("Invalid number of features (%d)." % n_features)

        if input_type not in ("dict", "pair", "string"):
            raise ValueError("input_type must be 'dict', 'pair' or 'string',"
                             " got %r." % input_type)

    def fit(self, X=None, y=None):
        """No-op.

        This method doesn't do anything. It exists purely for compatibility
        with the scikit-learn transformer API.

        Returns
        -------
        self : FeatureHasher

        """
        # repeat input validation for grid search (which calls set_params)
        self._validate_params(self.n_features, self.input_type)
        return self

    def transform(self, raw_X, y=None):
        """Transform a sequence of instances to a scipy.sparse matrix.

        Parameters
        ----------
        raw_X : iterable over iterable over raw features, length = n_samples
            Samples. Each sample must be iterable an (e.g., a list or tuple)
            containing/generating feature names (and optionally values, see
            the input_type constructor argument) which will be hashed.
            raw_X need not support the len function, so it can be the result
            of a generator; n_samples is determined on the fly.
        y : (ignored)

        Returns
        -------
        X : scipy.sparse matrix, shape = (n_samples, self.n_features)
            Feature matrix, for use with estimators or further transformers.

        """
        raw_X = iter(raw_X)
        if self.input_type == "dict":
            raw_X = (_iteritems(d) for d in raw_X)
        elif self.input_type == "string":
            raw_X = (((f, 1) for f in x) for x in raw_X)
        indices, indptr, values = \
            _hashing.transform(raw_X, self.n_features, self.dtype)
        n_samples = indptr.shape[0] - 1

        if n_samples == 0:
            raise ValueError("Cannot vectorize empty sequence.")

        X = sp.csr_matrix((values, indices, indptr), dtype=self.dtype,
                          shape=(n_samples, self.n_features))
        X.sum_duplicates()  # also sorts the indices
        if self.non_negative:
            np.abs(X.data, X.data)
        return X
