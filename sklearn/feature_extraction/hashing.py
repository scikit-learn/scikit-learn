# Author: Lars Buitinck
# License: BSD 3 clause

import numbers
import warnings

import numpy as np
import scipy.sparse as sp

from ..base import BaseEstimator, TransformerMixin

from ..externals import six

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
    Feature values must be (finite) numbers.

    This class is a low-memory alternative to DictVectorizer and
    CountVectorizer, intended for large-scale (online) learning and situations
    where memory is tight, e.g. when running prediction code on embedded
    devices.

    Read more in the :ref:`User Guide <feature_hashing>`.

    Parameters
    ----------
    n_features : integer, optional
        The number of features (columns) in the output matrices. Small numbers
        of features are likely to cause hash collisions, but large numbers
        will cause larger coefficient dimensions in linear learners.
    input_type : string, optional, default "dict"
        Either "dict" (the default) to accept dictionaries over
        (feature_name, value); "pair" to accept pairs of (feature_name, value);
        or "string" to accept single strings.
        feature_name should be a string, while value should be a number.
        In the case of "string", a value of 1 is implied.
        The feature_name is hashed to find the appropriate column for the
        feature. The value's sign might be flipped in the output (but see
        non_negative, below).
    dtype : numpy type, optional, default np.float64
        The type of feature values. Passed to scipy.sparse matrix constructors
        as the dtype argument. Do not set this to bool, np.boolean or any
        unsigned integer type.
    alternate_sign : boolean, optional, default True
        When True, an alternating sign is added to the features as to
        approximately conserve the inner product in the hashed space even for
        small n_features. This approach is similar to sparse random projection.

    non_negative : boolean, optional, default False
        When True, an absolute value is applied to the features matrix prior to
        returning it. When used in conjunction with alternate_sign=True, this
        significantly reduces the inner product preservation property.

        .. deprecated:: 0.19
            This option will be removed in 0.21.


    Examples
    --------
    >>> from sklearn.feature_extraction import FeatureHasher
    >>> h = FeatureHasher(n_features=10)
    >>> D = [{'dog': 1, 'cat':2, 'elephant':4},{'dog': 2, 'run': 5}]
    >>> f = h.transform(D)
    >>> f.toarray()
    array([[ 0.,  0., -4., -1.,  0.,  0.,  0.,  0.,  0.,  2.],
           [ 0.,  0.,  0., -2., -5.,  0.,  0.,  0.,  0.,  0.]])

    See also
    --------
    DictVectorizer : vectorizes string-valued features using a hash table.
    sklearn.preprocessing.OneHotEncoder : handles nominal/categorical features
        encoded as columns of integers.
    """

    def __init__(self, n_features=(2 ** 20), input_type="dict",
                 dtype=np.float64, alternate_sign=True, non_negative=False):
        self._validate_params(n_features, input_type)
        if non_negative:
            warnings.warn("the option non_negative=True has been deprecated"
                          " in 0.19 and will be removed"
                          " in version 0.21.", DeprecationWarning)

        self.dtype = dtype
        self.input_type = input_type
        self.n_features = n_features
        self.alternate_sign = alternate_sign
        self.non_negative = non_negative

    @staticmethod
    def _validate_params(n_features, input_type):
        # strangely, np.int16 instances are not instances of Integral,
        # while np.int64 instances are...
        if not isinstance(n_features, (numbers.Integral, np.integer)):
            raise TypeError("n_features must be integral, got %r (%s)."
                            % (n_features, type(n_features)))
        elif n_features < 1 or n_features >= 2 ** 31:
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
        X : array-like

        Returns
        -------
        self : FeatureHasher

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
        X : scipy.sparse matrix, shape = (n_samples, self.n_features)
            Feature matrix, for use with estimators or further transformers.

        """
        raw_X = iter(raw_X)
        if self.input_type == "dict":
            raw_X = (_iteritems(d) for d in raw_X)
        elif self.input_type == "string":
            raw_X = (((f, 1) for f in x) for x in raw_X)
        indices, indptr, values = \
            _transform(raw_X, self.n_features, self.dtype,
                               self.alternate_sign)
        n_samples = indptr.shape[0] - 1

        if n_samples == 0:
            raise ValueError("Cannot vectorize empty sequence.")

        X = sp.csr_matrix((values, indices, indptr), dtype=self.dtype,
                          shape=(n_samples, self.n_features))
        X.sum_duplicates()  # also sorts the indices

        if self.non_negative:
            np.abs(X.data, X.data)
        return X


import sys
import array
import numpy as np

from sklearn.utils.murmurhash import murmurhash3_bytes_s32
from sklearn.utils.fixes import sp_version

def _transform(raw_X, n_features, dtype, alternate_sign=1):
    """Guts of FeatureHasher.transform.

    Returns
    -------
    n_samples : integer
    indices, indptr, values : lists
        For constructing a scipy.sparse.csr_matrix.

    """
    assert n_features > 0

    indices = []
    if sys.version_info >= (3, 3):
        indices_array_dtype = "q"
        indices_np_dtype = np.longlong
    else:
        # On Windows with PY2.7 long int would still correspond to 32 bit.
        indices_array_dtype = "l"
        indices_np_dtype = np.int_

    indptr = [0]

    # Since Python array does not understand Numpy dtypes, we grow the indices
    # and values arrays ourselves. Use a Py_ssize_t capacity for safety.
    capacity = 8192     # arbitrary
    size = 0
    values = np.empty(capacity, dtype=dtype)

    for x in raw_X:
        for f, v in x:
            if isinstance(v, six.string_types):
                f = "%s%s%s" % (f, '=', v)
                value = 1
            else:
                value = v

            if value == 0:
                continue

            if isinstance(f, six.text_type):
                f = f.encode("utf-8")
            # Need explicit type check because Murmurhash does not propagate
            # all exceptions. Add "except *" there?
            elif not isinstance(f, six.binary_type):
                raise TypeError("feature names must be strings")

            h = murmurhash3_bytes_s32(f, 0)

            indices.append(abs(h) % n_features)
            # improve inner product preservation in the hashed space
            if alternate_sign:
                value *= (h >= 0) * 2 - 1
            values[size] = value
            size += 1

            if size == capacity:
                capacity *= 2
                # can't use resize member because there might be multiple
                # references to the arrays due to Cython's error checking
                values = np.resize(values, capacity)

        indptr.append(size)

    indices = array.array("i", indices)
    indptr = array.array(indices_array_dtype, indptr)
    indices_a = np.frombuffer(indices, dtype=np.int32)
    indptr_a = np.frombuffer(indptr, dtype=indices_np_dtype)

    if indptr[-1] > 2147483648:  # = 2**31
        if sp_version < (0, 14):
            raise ValueError(('sparse CSR array has {} non-zero '
                              'elements and requires 64 bit indexing, '
                              ' which is unsupported with scipy {}. '
                              'Please upgrade to scipy >=0.14')
                             .format(indptr[-1], '.'.join(sp_version)))
        # both indices and indptr have the same dtype in CSR arrays
        indices_a = indices_a.astype(np.int64)
    else:
        indptr_a = indptr_a.astype(np.int32)

    return (indices_a, indptr_a, values[:size])
