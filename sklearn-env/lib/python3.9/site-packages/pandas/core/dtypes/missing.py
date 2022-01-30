"""
missing types & inference
"""
from __future__ import annotations

from decimal import Decimal
from functools import partial

import numpy as np

from pandas._config import get_option

from pandas._libs import lib
import pandas._libs.missing as libmissing
from pandas._libs.tslibs import (
    NaT,
    iNaT,
)
from pandas._typing import (
    ArrayLike,
    DtypeObj,
    npt,
)

from pandas.core.dtypes.common import (
    DT64NS_DTYPE,
    TD64NS_DTYPE,
    ensure_object,
    is_bool_dtype,
    is_categorical_dtype,
    is_complex_dtype,
    is_datetimelike_v_numeric,
    is_dtype_equal,
    is_extension_array_dtype,
    is_float_dtype,
    is_integer_dtype,
    is_object_dtype,
    is_scalar,
    is_string_or_object_np_dtype,
    needs_i8_conversion,
)
from pandas.core.dtypes.dtypes import (
    CategoricalDtype,
    ExtensionDtype,
    IntervalDtype,
    PeriodDtype,
)
from pandas.core.dtypes.generic import (
    ABCDataFrame,
    ABCExtensionArray,
    ABCIndex,
    ABCMultiIndex,
    ABCSeries,
)
from pandas.core.dtypes.inference import is_list_like

isposinf_scalar = libmissing.isposinf_scalar
isneginf_scalar = libmissing.isneginf_scalar

nan_checker = np.isnan
INF_AS_NA = False
_dtype_object = np.dtype("object")
_dtype_str = np.dtype(str)


def isna(obj):
    """
    Detect missing values for an array-like object.

    This function takes a scalar or array-like object and indicates
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
    notna : Boolean inverse of pandas.isna.
    Series.isna : Detect missing values in a Series.
    DataFrame.isna : Detect missing values in a DataFrame.
    Index.isna : Detect missing values in an Index.

    Examples
    --------
    Scalar arguments (including strings) result in a scalar boolean.

    >>> pd.isna('dog')
    False

    >>> pd.isna(pd.NA)
    True

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


def _isna(obj, inf_as_na: bool = False):
    """
    Detect missing values, treating None, NaN or NA as null. Infinite
    values will also be treated as null if inf_as_na is True.

    Parameters
    ----------
    obj: ndarray or object value
        Input array or scalar value.
    inf_as_na: bool
        Whether to treat infinity as null.

    Returns
    -------
    boolean ndarray or boolean
    """
    if is_scalar(obj):
        return libmissing.checknull(obj, inf_as_na=inf_as_na)
    elif isinstance(obj, ABCMultiIndex):
        raise NotImplementedError("isna is not defined for MultiIndex")
    elif isinstance(obj, type):
        return False
    elif isinstance(obj, (np.ndarray, ABCExtensionArray)):
        return _isna_array(obj, inf_as_na=inf_as_na)
    elif isinstance(obj, ABCIndex):
        # Try to use cached isna, which also short-circuits for integer dtypes
        #  and avoids materializing RangeIndex._values
        if not obj._can_hold_na:
            return obj.isna()
        return _isna_array(obj._values, inf_as_na=inf_as_na)

    elif isinstance(obj, ABCSeries):
        result = _isna_array(obj._values, inf_as_na=inf_as_na)
        # box
        result = obj._constructor(result, index=obj.index, name=obj.name, copy=False)
        return result
    elif isinstance(obj, ABCDataFrame):
        return obj.isna()
    elif isinstance(obj, list):
        return _isna_array(np.asarray(obj, dtype=object), inf_as_na=inf_as_na)
    elif hasattr(obj, "__array__"):
        return _isna_array(np.asarray(obj), inf_as_na=inf_as_na)
    else:
        return False


def _use_inf_as_na(key):
    """
    Option change callback for na/inf behaviour.

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

    * https://stackoverflow.com/questions/4859217/
      programmatically-creating-variables-in-python/4859312#4859312
    """
    inf_as_na = get_option(key)
    globals()["_isna"] = partial(_isna, inf_as_na=inf_as_na)
    if inf_as_na:
        globals()["nan_checker"] = lambda x: ~np.isfinite(x)
        globals()["INF_AS_NA"] = True
    else:
        globals()["nan_checker"] = np.isnan
        globals()["INF_AS_NA"] = False


def _isna_array(values: ArrayLike, inf_as_na: bool = False):
    """
    Return an array indicating which values of the input array are NaN / NA.

    Parameters
    ----------
    obj: ndarray or ExtensionArray
        The input array whose elements are to be checked.
    inf_as_na: bool
        Whether or not to treat infinite values as NA.

    Returns
    -------
    array-like
        Array of boolean values denoting the NA status of each element.
    """
    dtype = values.dtype

    if not isinstance(values, np.ndarray):
        # i.e. ExtensionArray
        if inf_as_na and is_categorical_dtype(dtype):
            result = libmissing.isnaobj(values.to_numpy(), inf_as_na=inf_as_na)
        else:
            # error: Incompatible types in assignment (expression has type
            # "Union[ndarray[Any, Any], ExtensionArraySupportsAnyAll]", variable has
            # type "ndarray[Any, dtype[bool_]]")
            result = values.isna()  # type: ignore[assignment]
    elif is_string_or_object_np_dtype(values.dtype):
        result = _isna_string_dtype(values, inf_as_na=inf_as_na)
    elif needs_i8_conversion(dtype):
        # this is the NaT pattern
        result = values.view("i8") == iNaT
    else:
        if inf_as_na:
            result = ~np.isfinite(values)
        else:
            result = np.isnan(values)

    return result


def _isna_string_dtype(values: np.ndarray, inf_as_na: bool) -> npt.NDArray[np.bool_]:
    # Working around NumPy ticket 1542
    dtype = values.dtype

    if dtype.kind in ("S", "U"):
        result = np.zeros(values.shape, dtype=bool)
    else:

        if values.ndim == 1:
            result = libmissing.isnaobj(values, inf_as_na=inf_as_na)
        elif values.ndim == 2:
            result = libmissing.isnaobj2d(values, inf_as_na=inf_as_na)
        else:
            # 0-D, reached via e.g. mask_missing
            result = libmissing.isnaobj(values.ravel(), inf_as_na=inf_as_na)
            result = result.reshape(values.shape)

    return result


def notna(obj):
    """
    Detect non-missing values for an array-like object.

    This function takes a scalar or array-like object and indicates
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
    isna : Boolean inverse of pandas.notna.
    Series.notna : Detect valid values in a Series.
    DataFrame.notna : Detect valid values in a DataFrame.
    Index.notna : Detect valid values in an Index.

    Examples
    --------
    Scalar arguments (including strings) result in a scalar boolean.

    >>> pd.notna('dog')
    True

    >>> pd.notna(pd.NA)
    False

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


def isna_compat(arr, fill_value=np.nan) -> bool:
    """
    Parameters
    ----------
    arr: a numpy array
    fill_value: fill value, default to np.nan

    Returns
    -------
    True if we can fill using this fill_value
    """
    if isna(fill_value):
        dtype = arr.dtype
        return not (is_bool_dtype(dtype) or is_integer_dtype(dtype))
    return True


def array_equivalent(
    left,
    right,
    strict_nan: bool = False,
    dtype_equal: bool = False,
) -> bool:
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
    dtype_equal : bool, default False
        Whether `left` and `right` are known to have the same dtype
        according to `is_dtype_equal`. Some methods like `BlockManager.equals`.
        require that the dtypes match. Setting this to ``True`` can improve
        performance, but will give different results for arrays that are
        equal but different dtypes.

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

    if dtype_equal:
        # fastpath when we require that the dtypes match (Block.equals)
        if left.dtype.kind in ["f", "c"]:
            return _array_equivalent_float(left, right)
        elif is_datetimelike_v_numeric(left.dtype, right.dtype):
            return False
        elif needs_i8_conversion(left.dtype):
            return _array_equivalent_datetimelike(left, right)
        elif is_string_or_object_np_dtype(left.dtype):
            # TODO: fastpath for pandas' StringDtype
            return _array_equivalent_object(left, right, strict_nan)
        else:
            return np.array_equal(left, right)

    # Slow path when we allow comparing different dtypes.
    # Object arrays can contain None, NaN and NaT.
    # string dtypes must be come to this path for NumPy 1.7.1 compat
    if left.dtype.kind in "OSU" or right.dtype.kind in "OSU":
        # Note: `in "OSU"` is non-trivially faster than `in ["O", "S", "U"]`
        #  or `in ("O", "S", "U")`
        return _array_equivalent_object(left, right, strict_nan)

    # NaNs can occur in float and complex arrays.
    if is_float_dtype(left.dtype) or is_complex_dtype(left.dtype):
        if not (left.size and right.size):
            return True
        return ((left == right) | (isna(left) & isna(right))).all()

    elif is_datetimelike_v_numeric(left, right):
        # GH#29553 avoid numpy deprecation warning
        return False

    elif needs_i8_conversion(left.dtype) or needs_i8_conversion(right.dtype):
        # datetime64, timedelta64, Period
        if not is_dtype_equal(left.dtype, right.dtype):
            return False

        left = left.view("i8")
        right = right.view("i8")

    # if we have structured dtypes, compare first
    if (
        left.dtype.type is np.void or right.dtype.type is np.void
    ) and left.dtype != right.dtype:
        return False

    return np.array_equal(left, right)


def _array_equivalent_float(left, right) -> bool:
    return bool(((left == right) | (np.isnan(left) & np.isnan(right))).all())


def _array_equivalent_datetimelike(left, right):
    return np.array_equal(left.view("i8"), right.view("i8"))


def _array_equivalent_object(left: np.ndarray, right: np.ndarray, strict_nan: bool):
    if not strict_nan:
        # isna considers NaN and None to be equivalent.

        if left.flags["F_CONTIGUOUS"] and right.flags["F_CONTIGUOUS"]:
            # we can improve performance by doing a copy-free ravel
            # e.g. in frame_methods.Equals.time_frame_nonunique_equal
            #  if we transposed the frames
            left = left.ravel("K")
            right = right.ravel("K")

        return lib.array_equivalent_object(
            ensure_object(left.ravel()), ensure_object(right.ravel())
        )

    for left_value, right_value in zip(left, right):
        if left_value is NaT and right_value is not NaT:
            return False

        elif left_value is libmissing.NA and right_value is not libmissing.NA:
            return False

        elif isinstance(left_value, float) and np.isnan(left_value):
            if not isinstance(right_value, float) or not np.isnan(right_value):
                return False
        else:
            try:
                if np.any(np.asarray(left_value != right_value)):
                    return False
            except TypeError as err:
                if "boolean value of NA is ambiguous" in str(err):
                    return False
                raise
    return True


def array_equals(left: ArrayLike, right: ArrayLike) -> bool:
    """
    ExtensionArray-compatible implementation of array_equivalent.
    """
    if not is_dtype_equal(left.dtype, right.dtype):
        return False
    elif isinstance(left, ABCExtensionArray):
        return left.equals(right)
    else:
        return array_equivalent(left, right, dtype_equal=True)


def infer_fill_value(val):
    """
    infer the fill value for the nan/NaT from the provided
    scalar/ndarray/list-like if we are a NaT, return the correct dtyped
    element to provide proper block construction
    """
    if not is_list_like(val):
        val = [val]
    val = np.array(val, copy=False)
    if needs_i8_conversion(val.dtype):
        return np.array("NaT", dtype=val.dtype)
    elif is_object_dtype(val.dtype):
        dtype = lib.infer_dtype(ensure_object(val), skipna=False)
        if dtype in ["datetime", "datetime64"]:
            return np.array("NaT", dtype=DT64NS_DTYPE)
        elif dtype in ["timedelta", "timedelta64"]:
            return np.array("NaT", dtype=TD64NS_DTYPE)
    return np.nan


def maybe_fill(arr: np.ndarray) -> np.ndarray:
    """
    Fill numpy.ndarray with NaN, unless we have a integer or boolean dtype.
    """
    if arr.dtype.kind not in ("u", "i", "b"):
        arr.fill(np.nan)
    return arr


def na_value_for_dtype(dtype: DtypeObj, compat: bool = True):
    """
    Return a dtype compat na value

    Parameters
    ----------
    dtype : string / dtype
    compat : bool, default True

    Returns
    -------
    np.dtype or a pandas dtype

    Examples
    --------
    >>> na_value_for_dtype(np.dtype('int64'))
    0
    >>> na_value_for_dtype(np.dtype('int64'), compat=False)
    nan
    >>> na_value_for_dtype(np.dtype('float64'))
    nan
    >>> na_value_for_dtype(np.dtype('bool'))
    False
    >>> na_value_for_dtype(np.dtype('datetime64[ns]'))
    numpy.datetime64('NaT')
    """

    if isinstance(dtype, ExtensionDtype):
        return dtype.na_value
    elif needs_i8_conversion(dtype):
        return dtype.type("NaT", "ns")
    elif is_float_dtype(dtype):
        return np.nan
    elif is_integer_dtype(dtype):
        if compat:
            return 0
        return np.nan
    elif is_bool_dtype(dtype):
        if compat:
            return False
        return np.nan
    return np.nan


def remove_na_arraylike(arr):
    """
    Return array-like containing only true/non-NaN values, possibly empty.
    """
    if is_extension_array_dtype(arr):
        return arr[notna(arr)]
    else:
        return arr[notna(np.asarray(arr))]


def is_valid_na_for_dtype(obj, dtype: DtypeObj) -> bool:
    """
    isna check that excludes incompatible dtypes

    Parameters
    ----------
    obj : object
    dtype : np.datetime64, np.timedelta64, DatetimeTZDtype, or PeriodDtype

    Returns
    -------
    bool
    """
    if not lib.is_scalar(obj) or not isna(obj):
        return False
    elif dtype.kind == "M":
        if isinstance(dtype, np.dtype):
            # i.e. not tzaware
            return not isinstance(obj, (np.timedelta64, Decimal))
        # we have to rule out tznaive dt64("NaT")
        return not isinstance(obj, (np.timedelta64, np.datetime64, Decimal))
    elif dtype.kind == "m":
        return not isinstance(obj, (np.datetime64, Decimal))
    elif dtype.kind in ["i", "u", "f", "c"]:
        # Numeric
        return obj is not NaT and not isinstance(obj, (np.datetime64, np.timedelta64))

    elif dtype == _dtype_str:
        # numpy string dtypes to avoid float np.nan
        return not isinstance(obj, (np.datetime64, np.timedelta64, Decimal, float))

    elif dtype == _dtype_object:
        # This is needed for Categorical, but is kind of weird
        return True

    elif isinstance(dtype, PeriodDtype):
        return not isinstance(obj, (np.datetime64, np.timedelta64, Decimal))

    elif isinstance(dtype, IntervalDtype):
        return lib.is_float(obj) or obj is None or obj is libmissing.NA

    elif isinstance(dtype, CategoricalDtype):
        return is_valid_na_for_dtype(obj, dtype.categories.dtype)

    # fallback, default to allowing NaN, None, NA, NaT
    return not isinstance(obj, (np.datetime64, np.timedelta64, Decimal))
