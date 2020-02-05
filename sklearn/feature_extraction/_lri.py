# Author: Lars Buitinck
# License: BSD 3 clause

import numbers

import numpy as np
import scipy.sparse as sp

from ..utils import IS_PYPY
from ..base import BaseEstimator, TransformerMixin

if not IS_PYPY:
    from ._lri_fast import transform as _lri_transform
else:
    def _lri_transform(*args, **kwargs):
        raise NotImplementedError(
                'FeatureLightweightRandomIndexing is not compatible with '
                'PyPy (see '
                'https://github.com/scikit-learn/scikit-learn/issues/11540 '
                'for the status updates).')


def _iteritems(d):
    """Like d.iteritems, but accepts any collections.Mapping."""
    return d.iteritems() if hasattr(d, "iteritems") else d.items()


class FeatureLightweightRandomIndexing(TransformerMixin, BaseEstimator):
    """Applies lightweight random indexing [1] to a collection of texts.

    Derived from :ref:`FeatureHasher`, it computes two hashes for every
    feature (avoiding collisions). The resulting indexing vectors for features
    have good properties of orthogonality and information preservation in
    reduced feature spaces.

    Parameters
    ----------
    n_features : int, default=2**20
        The number of features (columns) in the output matrices. Small numbers
        of features are likely to cause hash collisions, but large numbers
        will cause larger coefficient dimensions in linear learners.
    input_type : {"dict", "pair"}, default="dict"
        Either "dict" (the default) to accept dictionaries over
        (feature_name, value); "pair" to accept pairs of (feature_name, value);
        or "string" to accept single strings.
        feature_name should be a string, while value should be a number.
        In the case of "string", a value of 1 is implied.
        The feature_name is hashed to find the appropriate column for the
        feature. The value's sign might be flipped in the output (but see
        non_negative, below).
    dtype : numpy dtype, default=np.float64
        The type of feature values. Passed to scipy.sparse matrix constructors
        as the dtype argument. Do not set this to bool, np.boolean or any
        unsigned integer type.

    Examples
    --------
    >>> from sklearn.feature_extraction import FeatureLightweightRandomIndexing
    >>> h = FeatureLightweightRandomIndexing(n_features=10)
    >>> D = [{'dog': 1, 'cat':2, 'elephant':4},{'dog': 2, 'run': 5}]
    >>> f = h.transform(D)
    >>> f.toarray()
    array([[ 0.,  2., -4., -1., -1.,  0.,  0.,  0., -4.,  2.],
           [ 5.,  0.,  0., -2., -7.,  0.,  0.,  0.,  0.,  0.]])

    See also
    --------
    DictVectorizer : vectorizes string-valued features using a hash table.
    FeatureHasher : implements the hashing trick
    sklearn.preprocessing.OneHotEncoder : handles nominal/categorical features.

    References
    ----------
    .. [1] Lightweight random indexing for polylingual text classification,
    A. Moreo, A. Esuli & F. Sebastiani, JAIR, 2016.
    """

    def __init__(self, n_features=(2 ** 20), input_type="dict",
                 dtype=np.float64):
        self._validate_params(n_features, input_type)

        self.dtype = dtype
        self.input_type = input_type
        self.n_features = n_features

    @staticmethod
    def _validate_params(n_features, input_type):
        # strangely, np.int16 instances are not instances of Integral,
        # while np.int64 instances are...
        if not isinstance(n_features, numbers.Integral):
            raise TypeError("n_features must be integral, got %r (%s)."
                            % (n_features, type(n_features)))
        elif n_features < 1 or n_features >= np.iinfo(np.int32).max + 1:
            raise ValueError("Invalid number of features (%d)." % n_features)

        if input_type not in ("dict", "pair", "string"):
            raise ValueError("input_type must be 'dict', 'pair' or 'string',"
                             " got %r." % input_type)

    def fit(self, X=None, y=None):
        """No-op.

        This method doesn't do anything. It exists purely for compatibility
        with the scikit-learn transformer API.

        Parameters
        ----------
        X : ndarray

        Returns
        -------
        self : FeatureLightweightRandomIndexing

        """
        # repeat input validation for grid search (which calls set_params)
        self._validate_params(self.n_features, self.input_type)
        return self

    def transform(self, raw_X):
        """Transform a sequence of instances to a scipy.sparse matrix.

        Parameters
        ----------
        raw_X : iterable over iterable over raw features, length = n_samples
            Samples. Each sample must be iterable an (e.g., a list or tuple)
            containing/generating feature names (and optionally values, see
            the input_type constructor argument) which will be hashed.
            raw_X need not support the len function, so it can be the result
            of a generator; n_samples is determined on the fly.

        Returns
        -------
        X : sparse matrix of shape (n_samples, n_features)
            Feature matrix, for use with estimators or further transformers.

        """
        raw_X = iter(raw_X)
        if self.input_type == "dict":
            raw_X = (_iteritems(d) for d in raw_X)
        elif self.input_type == "string":
            raw_X = (((f, 1) for f in x) for x in raw_X)
        indices, indptr, values = \
            _lri_transform(raw_X, self.n_features, self.dtype, seed=0)
        n_samples = indptr.shape[0] - 1

        if n_samples == 0:
            raise ValueError("Cannot vectorize empty sequence.")

        X = sp.csr_matrix((values, indices, indptr), dtype=self.dtype,
                          shape=(n_samples, self.n_features))
        X.sum_duplicates()  # also sorts the indices

        return X

    def _more_tags(self):
        return {'X_types': [self.input_type]}
