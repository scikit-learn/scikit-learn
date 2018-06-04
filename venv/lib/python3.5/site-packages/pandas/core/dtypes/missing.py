"""
missing types & inference
"""
import numpy as np
from pandas._libs import lib, missing as libmissing
from pandas._libs.tslib import NaT, iNaT
from .generic import (ABCMultiIndex, ABCSeries,
                      ABCIndexClass, ABCGeneric,
                      ABCExtensionArray)
from .common import (is_string_dtype, is_datetimelike,
                     is_datetimelike_v_numeric, is_float_dtype,
                     is_datetime64_dtype, is_datetime64tz_dtype,
                     is_timedelta64_dtype, is_interval_dtype,
                     is_period_dtype,
                     is_complex_dtype,
                     is_string_like_dtype, is_bool_dtype,
                     is_integer_dtype, is_dtype_equal,
                     is_extension_array_dtype,
                     needs_i8_conversion, _ensure_object,
                     pandas_dtype,
                     is_scalar,
                     is_object_dtype,
                     is_integer,
                     _TD_DTYPE,
                     _NS_DTYPE)
from .inference import is_list_like

isposinf_scalar = libmissing.isposinf_scalar
isneginf_scalar = libmissing.isneginf_scalar


def isna(obj):
    """
    Detect missing values for an array-like object.

    This function takes a scalar or array-like object and indictates
    whether values are missing (``NaN`` in numeric arrays, ``None`` or ``NaN``
    in object arrays, ``NaT`` in datetimelike).

    Parameters
    ----------
    obj : scalar or array-like
        Object to check for null or missing values.

    Returns
    -------
    bool or array-like of bool
        For scalar input, returns a scalar boolean.
        For array input, returns an array of boolean indicating whether each
        corresponding element is missing.

    See Also
    --------
    notna : boolean inverse of pandas.isna.
    Series.isna : Detetct missing values in a Series.
    DataFrame.isna : Detect missing values in a DataFrame.
    Index.isna : Detect missing values in an Index.

    Examples
    --------
    Scalar arguments (including strings) result in a scalar boolean.

    >>> pd.isna('dog')
    False

    >>> pd.isna(np.nan)
    True

    ndarrays result in an ndarray of booleans.

    >>> array = np.array([[1, np.nan, 3], [4, 5, np.nan]])
    >>> array
    array([[ 1., nan,  3.],
           [ 4.,  5., nan]])
    >>> pd.isna(array)
    array([[False,  True, False],
           [False, False,  True]])

    For indexes, an ndarray of booleans is returned.

    >>> index = pd.DatetimeIndex(["2017-07-05", "2017-07-06", None,
    ...                           "2017-07-08"])
    >>> index
    DatetimeIndex(['2017-07-05', '2017-07-06', 'NaT', '2017-07-08'],
                  dtype='datetime64[ns]', freq=None)
    >>> pd.isna(index)
    array([False, False,  True, False])

    For Series and DataFrame, the same type is returned, containing booleans.

    >>> df = pd.DataFrame([['ant', 'bee', 'cat'], ['dog', None, 'fly']])
    >>> df
         0     1    2
    0  ant   bee  cat
    1  dog  None  fly
    >>> pd.isna(df)
           0      1      2
    0  False  False  False
    1  False   True  False

    >>> pd.isna(df[1])
    0    False
    1     True
    Name: 1, dtype: bool
    """
    return _isna(obj)


isnull = isna


def _isna_new(obj):
    if is_scalar(obj):
        return libmissing.checknull(obj)
    # hack (for now) because MI registers as ndarray
    elif isinstance(obj, ABCMultiIndex):
        raise NotImplementedError("isna is not defined for MultiIndex")
    elif isinstance(obj, (ABCSeries, np.ndarray, ABCIndexClass,
                          ABCExtensionArray)):
        return _isna_ndarraylike(obj)
    elif isinstance(obj, ABCGeneric):
        return obj._constructor(obj._data.isna(func=isna))
    elif isinstance(obj, list):
        return _isna_ndarraylike(np.asarray(obj, dtype=object))
    elif hasattr(obj, '__array__'):
        return _isna_ndarraylike(np.asarray(obj))
    else:
        return obj is None


def _isna_old(obj):
    """Detect missing values. Treat None, NaN, INF, -INF as null.

    Parameters
    ----------
    arr: ndarray or object value

    Returns
    -------
    boolean ndarray or boolean
    """
    if is_scalar(obj):
        return libmissing.checknull_old(obj)
    # hack (for now) because MI registers as ndarray
    elif isinstance(obj, ABCMultiIndex):
        raise NotImplementedError("isna is not defined for MultiIndex")
    elif isinstance(obj, (ABCSeries, np.ndarray, ABCIndexClass)):
        return _isna_ndarraylike_old(obj)
    elif isinstance(obj, ABCGeneric):
        return obj._constructor(obj._data.isna(func=_isna_old))
    elif isinstance(obj, list):
        return _isna_ndarraylike_old(np.asarray(obj, dtype=object))
    elif hasattr(obj, '__array__'):
        return _isna_ndarraylike_old(np.asarray(obj))
    else:
        return obj is None


_isna = _isna_new


def _use_inf_as_na(key):
    """Option change callback for na/inf behaviour
    Choose which replacement for numpy.isnan / -numpy.isfinite is used.

    Parameters
    ----------
    flag: bool
        True means treat None, NaN, INF, -INF as null (old way),
        False means None and NaN are null, but INF, -INF are not null
        (new way).

    Notes
    -----
    This approach to setting global module values is discussed and
    approved here:

    * http://stackoverflow.com/questions/4859217/
      programmatically-creating-variables-in-python/4859312#4859312
    """
    from pandas.core.config import get_option
    flag = get_option(key)
    if flag:
        globals()['_isna'] = _isna_old
    else:
        globals()['_isna'] = _isna_new


def _isna_ndarraylike(obj):
    values = getattr(obj, 'values', obj)
    dtype = values.dtype

    if is_extension_array_dtype(obj):
        if isinstance(obj, (ABCIndexClass, ABCSeries)):
            values = obj._values
        else:
            values = obj
        result = values.isna()
    elif is_interval_dtype(values):
        # TODO(IntervalArray): remove this if block
        from pandas import IntervalIndex
        result = IntervalIndex(obj).isna()
    elif is_string_dtype(dtype):
        # Working around NumPy ticket 1542
        shape = values.shape

        if is_string_like_dtype(dtype):
            # object array of strings
            result = np.zeros(values.shape, dtype=bool)
        else:
            # object array of non-strings
            result = np.empty(shape, dtype=bool)
            vec = libmissing.isnaobj(values.ravel())
            result[...] = vec.reshape(shape)

    elif needs_i8_conversion(obj):
        # this is the NaT pattern
        result = values.view('i8') == iNaT
    else:
        result = np.isnan(values)

    # box
    if isinstance(obj, ABCSeries):
        from pandas import Series
        result = Series(result, index=obj.index, name=obj.name, copy=False)

    return result


def _isna_ndarraylike_old(obj):
    values = getattr(obj, 'values', obj)
    dtype = values.dtype

    if is_string_dtype(dtype):
        # Working around NumPy ticket 1542
        shape = values.shape

        if is_string_like_dtype(dtype):
            result = np.zeros(values.shape, dtype=bool)
        else:
            result = np.empty(shape, dtype=bool)
            vec = libmissing.isnaobj_old(values.ravel())
            result[:] = vec.reshape(shape)

    elif is_datetime64_dtype(dtype):
        # this is the NaT pattern
        result = values.view('i8') == iNaT
    else:
        result = ~np.isfinite(values)

    # box
    if isinstance(obj, ABCSeries):
        from pandas import Series
        result = Series(result, index=obj.index, name=obj.name, copy=False)

    return result


def notna(obj):
    """
    Detect non-missing values for an array-like object.

    This function takes a scalar or array-like object and indictates
    whether values are valid (not missing, which is ``NaN`` in numeric
    arrays, ``None`` or ``NaN`` in object arrays, ``NaT`` in datetimelike).

    Parameters
    ----------
    obj : array-like or object value
        Object to check for *not* null or *non*-missing values.

    Returns
    -------
    bool or array-like of bool
        For scalar input, returns a scalar boolean.
        For array input, returns an array of boolean indicating whether each
        corresponding element is valid.

    See Also
    --------
    isna : boolean inverse of pandas.notna.
    Series.notna : Detetct valid values in a Series.
    DataFrame.notna : Detect valid values in a DataFrame.
    Index.notna : Detect valid values in an Index.

    Examples
    --------
    Scalar arguments (including strings) result in a scalar boolean.

    >>> pd.notna('dog')
    True

    >>> pd.notna(np.nan)
    False

    ndarrays result in an ndarray of booleans.

    >>> array = np.array([[1, np.nan, 3], [4, 5, np.nan]])
    >>> array
    array([[ 1., nan,  3.],
           [ 4.,  5., nan]])
    >>> pd.notna(array)
    array([[ True, False,  True],
           [ True,  True, False]])

    For indexes, an ndarray of booleans is returned.

    >>> index = pd.DatetimeIndex(["2017-07-05", "2017-07-06", None,
    ...                          "2017-07-08"])
    >>> index
    DatetimeIndex(['2017-07-05', '2017-07-06', 'NaT', '2017-07-08'],
                  dtype='datetime64[ns]', freq=None)
    >>> pd.notna(index)
    array([ True,  True, False,  True])

    For Series and DataFrame, the same type is returned, containing booleans.

    >>> df = pd.DataFrame([['ant', 'bee', 'cat'], ['dog', None, 'fly']])
    >>> df
         0     1    2
    0  ant   bee  cat
    1  dog  None  fly
    >>> pd.notna(df)
          0      1     2
    0  True   True  True
    1  True  False  True

    >>> pd.notna(df[1])
    0     True
    1    False
    Name: 1, dtype: bool
    """
    res = isna(obj)
    if is_scalar(res):
        return not res
    return ~res


notnull = notna


def is_null_datelike_scalar(other):
    """ test whether the object is a null datelike, e.g. Nat
    but guard against passing a non-scalar """
    if other is NaT or other is None:
        return True
    elif is_scalar(other):

        # a timedelta
        if hasattr(other, 'dtype'):
            return other.view('i8') == iNaT
        elif is_integer(other) and other == iNaT:
            return True
        return isna(other)
    return False


def _isna_compat(arr, fill_value=np.nan):
    """
    Parameters
    ----------
    arr: a numpy array
    fill_value: fill value, default to np.nan

    Returns
    -------
    True if we can fill using this fill_value
    """
    dtype = arr.dtype
    if isna(fill_value):
        return not (is_bool_dtype(dtype) or
                    is_integer_dtype(dtype))
    return True


def array_equivalent(left, right, strict_nan=False):
    """
    True if two arrays, left and right, have equal non-NaN elements, and NaNs
    in corresponding locations.  False otherwise. It is assumed that left and
    right are NumPy arrays of the same dtype. The behavior of this function
    (particularly with respect to NaNs) is not defined if the dtypes are
    different.

    Parameters
    ----------
    left, right : ndarrays
    strict_nan : bool, default False
        If True, consider NaN and None to be different.

    Returns
    -------
    b : bool
        Returns True if the arrays are equivalent.

    Examples
    --------
    >>> array_equivalent(
    ...     np.array([1, 2, np.nan]),
    ...     np.array([1, 2, np.nan]))
    True
    >>> array_equivalent(
    ...     np.array([1, np.nan, 2]),
    ...     np.array([1, 2, np.nan]))
    False
    """

    left, right = np.asarray(left), np.asarray(right)

    # shape compat
    if left.shape != right.shape:
        return False

    # Object arrays can contain None, NaN and NaT.
    # string dtypes must be come to this path for NumPy 1.7.1 compat
    if is_string_dtype(left) or is_string_dtype(right):

        if not strict_nan:
            # isna considers NaN and None to be equivalent.
            return lib.array_equivalent_object(
                _ensure_object(left.ravel()), _ensure_object(right.ravel()))

        for left_value, right_value in zip(left, right):
            if left_value is NaT and right_value is not NaT:
                return False

            elif isinstance(left_value, float) and np.isnan(left_value):
                if (not isinstance(right_value, float) or
                        not np.isnan(right_value)):
                    return False
            else:
                if left_value != right_value:
                    return False
        return True

    # NaNs can occur in float and complex arrays.
    if is_float_dtype(left) or is_complex_dtype(left):

        # empty
        if not (np.prod(left.shape) and np.prod(right.shape)):
            return True
        return ((left == right) | (isna(left) & isna(right))).all()

    # numpy will will not allow this type of datetimelike vs integer comparison
    elif is_datetimelike_v_numeric(left, right):
        return False

    # M8/m8
    elif needs_i8_conversion(left) and needs_i8_conversion(right):
        if not is_dtype_equal(left.dtype, right.dtype):
            return False

        left = left.view('i8')
        right = right.view('i8')

    # if we have structured dtypes, compare first
    if (left.dtype.type is np.void or
            right.dtype.type is np.void):
        if left.dtype != right.dtype:
            return False

    return np.array_equal(left, right)


def _infer_fill_value(val):
    """
    infer the fill value for the nan/NaT from the provided
    scalar/ndarray/list-like if we are a NaT, return the correct dtyped
    element to provide proper block construction
    """

    if not is_list_like(val):
        val = [val]
    val = np.array(val, copy=False)
    if is_datetimelike(val):
        return np.array('NaT', dtype=val.dtype)
    elif is_object_dtype(val.dtype):
        dtype = lib.infer_dtype(_ensure_object(val))
        if dtype in ['datetime', 'datetime64']:
            return np.array('NaT', dtype=_NS_DTYPE)
        elif dtype in ['timedelta', 'timedelta64']:
            return np.array('NaT', dtype=_TD_DTYPE)
    return np.nan


def _maybe_fill(arr, fill_value=np.nan):
    """
    if we have a compatible fill_value and arr dtype, then fill
    """
    if _isna_compat(arr, fill_value):
        arr.fill(fill_value)
    return arr


def na_value_for_dtype(dtype, compat=True):
    """
    Return a dtype compat na value

    Parameters
    ----------
    dtype : string / dtype
    compat : boolean, default True

    Returns
    -------
    np.dtype or a pandas dtype
    """
    dtype = pandas_dtype(dtype)

    if is_extension_array_dtype(dtype):
        return dtype.na_value
    if (is_datetime64_dtype(dtype) or is_datetime64tz_dtype(dtype) or
            is_timedelta64_dtype(dtype) or is_period_dtype(dtype)):
        return NaT
    elif is_float_dtype(dtype):
        return np.nan
    elif is_integer_dtype(dtype):
        if compat:
            return 0
        return np.nan
    elif is_bool_dtype(dtype):
        return False
    return np.nan


def remove_na_arraylike(arr):
    """
    Return array-like containing only true/non-NaN values, possibly empty.
    """
    if is_extension_array_dtype(arr):
        return arr[notna(arr)]
    else:
        return arr[notna(lib.values_from_object(arr))]
