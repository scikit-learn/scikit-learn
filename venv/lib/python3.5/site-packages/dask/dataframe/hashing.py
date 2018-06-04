from __future__ import absolute_import, division, print_function

import numpy as np
import pandas as pd

from .utils import PANDAS_VERSION

# In Pandas 0.19.2, a function to hash pandas objects was introduced. Object
# arrays are assumed to be strings, and are hashed with a cython implementation
# of siphash. However, the version in 0.19.2 hashes categoricals based on their
# integer codes, instead of taking into account the values they represent. This
# is fixed in pandas > 0.19.2. To support versions 0.19.0 and up, we do do the
# following:
#
# - For versions > 0.19.2, we use the provided `hash_pandas_object` function.
# - For 0.19.0 through 0.19.2, we copy the definition of `hash_pandas_object`
#   from pandas master (will be released as 0.20.0).
# - For 0.19.0 and 0.19.1, we use python's `hash` builtin to hash strings.
# - For 0.19.2, we use the `hash_object_array` method provided in pandas
#   (an implementation of siphash)
#
# When dask drops support for pandas <= 0.19.2, all this can be removed.

# XXX: Pandas uses release branches > 0.19.0, which doesn't play well with
# versioneer, since the release tags aren't ancestors of master. As such, we
# need to use this hacky awfulness to check if we're > 0.19.2.
if PANDAS_VERSION >= '0.20.0':
    from pandas.util import hash_pandas_object
elif PANDAS_VERSION not in ('0.19.1', '0.19.2') and PANDAS_VERSION > '0.19.0+460':
    from pandas.tools.hashing import hash_pandas_object
else:
    from pandas.types.common import (is_categorical_dtype, is_numeric_dtype,
                                     is_datetime64_dtype, is_timedelta64_dtype)
    from pandas.lib import is_bool_array

    if PANDAS_VERSION == '0.19.2':
        from pandas._hash import hash_object_array
    else:  # 0.19.0 and 0.19.1
        def hash_object_array(x, hash_key, encoding):
            return np.array([hash(i) for i in x], dtype=np.uint64)

    # 16 byte long hashing key
    _default_hash_key = '0123456789123456'

    def hash_pandas_object(obj, index=True, encoding='utf8', hash_key=None,
                           categorize=True):
        if hash_key is None:
            hash_key = _default_hash_key

        def adder(h, hashed_to_add):
            h = np.multiply(h, np.uint(3), h)
            return np.add(h, hashed_to_add, h)

        if isinstance(obj, pd.Index):
            h = hash_array(obj.values, encoding, hash_key,
                           categorize).astype('uint64')
            h = pd.Series(h, index=obj, dtype='uint64')
        elif isinstance(obj, pd.Series):
            h = hash_array(obj.values, encoding, hash_key,
                           categorize).astype('uint64')
            if index:
                h = adder(h, hash_pandas_object(obj.index,
                                                index=False,
                                                encoding=encoding,
                                                hash_key=hash_key,
                                                categorize=categorize).values)
            h = pd.Series(h, index=obj.index, dtype='uint64')
        elif isinstance(obj, pd.DataFrame):
            cols = obj.iteritems()
            first_series = next(cols)[1]
            h = hash_array(first_series.values, encoding,
                           hash_key, categorize).astype('uint64')
            for _, col in cols:
                h = adder(h, hash_array(col.values, encoding, hash_key,
                                        categorize))
            if index:
                h = adder(h, hash_pandas_object(obj.index,
                                                index=False,
                                                encoding=encoding,
                                                hash_key=hash_key,
                                                categorize=categorize).values)

            h = pd.Series(h, index=obj.index, dtype='uint64')
        else:
            raise TypeError("Unexpected type for hashing %s" % type(obj))
        return h

    def _hash_categorical(c, encoding, hash_key):
        hashed = hash_array(c.categories.values, encoding, hash_key,
                            categorize=False)
        mask = c.isnull()
        if len(hashed):
            result = hashed.take(c.codes)
        else:
            result = np.zeros(len(mask), dtype='uint64')

        if mask.any():
            result[mask] = np.iinfo(np.uint64).max

        return result

    def hash_array(vals, encoding='utf8', hash_key=None, categorize=True):
        if hash_key is None:
            hash_key = _default_hash_key

        # For categoricals, we hash the categories, then remap the codes to the
        # hash values. (This check is above the complex check so that we don't
        # ask numpy if categorical is a subdtype of complex, as it will choke.
        if is_categorical_dtype(vals.dtype):
            return _hash_categorical(vals, encoding, hash_key)

        # we'll be working with everything as 64-bit values, so handle this
        # 128-bit value early
        if np.issubdtype(vals.dtype, np.complex128):
            return hash_array(vals.real) + 23 * hash_array(vals.imag)

        # First, turn whatever array this is into unsigned 64-bit ints, if we
        # can manage it.
        if is_bool_array(vals):
            vals = vals.astype('u8')
        elif ((is_datetime64_dtype(vals) or
               is_timedelta64_dtype(vals) or
               is_numeric_dtype(vals)) and vals.dtype.itemsize <= 8):
            vals = vals.view('u{}'.format(vals.dtype.itemsize)).astype('u8')
        else:
            # With repeated values, its MUCH faster to categorize object
            # dtypes, then hash and rename categories. We allow skipping the
            # categorization when the values are known/likely to be unique.
            if categorize:
                codes, categories = pd.factorize(vals, sort=False)
                cat = pd.Categorical(codes, pd.Index(categories),
                                     ordered=False, fastpath=True)
                return _hash_categorical(cat, encoding, hash_key)

            vals = hash_object_array(vals, hash_key, encoding)

        # Then, redistribute these 64-bit ints within the space of 64-bit ints
        vals ^= vals >> 30
        vals *= np.uint64(0xbf58476d1ce4e5b9)
        vals ^= vals >> 27
        vals *= np.uint64(0x94d049bb133111eb)
        vals ^= vals >> 31
        return vals
