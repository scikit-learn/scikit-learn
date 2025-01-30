import numpy as np

# Define classes of supported dtypes and Python scalar types
# Variables ending in `_dtypes` only contain numpy.dtypes of the respective
# class; variables ending in `_types` additionally include Python scalar types.
signed_integer_dtypes = {np.int8, np.int16, np.int32, np.int64}
signed_integer_types = signed_integer_dtypes | {int}

unsigned_integer_dtypes = {np.uint8, np.uint16, np.uint32, np.uint64}

integer_dtypes = signed_integer_dtypes | unsigned_integer_dtypes
integer_types = signed_integer_types | unsigned_integer_dtypes

floating_dtypes = {np.float16, np.float32, np.float64}
floating_types = floating_dtypes | {float}

complex_dtypes = {np.complex64, np.complex128}
complex_types = complex_dtypes | {complex}

inexact_dtypes = floating_dtypes | complex_dtypes
inexact_types = floating_types | complex_types

bool_types = {np.dtype(bool), bool}

numeric_dtypes = integer_dtypes | inexact_dtypes | {np.bool_}
numeric_types = integer_types | inexact_types | bool_types


def numeric_dtype_min_max(dtype):
    """Return minimum and maximum representable value for a given dtype.

    A convenient wrapper around `numpy.finfo` and `numpy.iinfo` that
    additionally supports numpy.bool as well.

    Parameters
    ----------
    dtype : numpy.dtype
        The dtype. Tries to convert Python "types" such as int or float, to
        the corresponding NumPy dtype.

    Returns
    -------
    min, max : number
        Minimum and maximum of the given `dtype`. These scalars are themselves
        of the given `dtype`.

    Examples
    --------
    >>> import numpy as np
    >>> numeric_dtype_min_max(np.uint8)
    (0, 255)
    >>> numeric_dtype_min_max(bool)
    (False, True)
    >>> numeric_dtype_min_max(np.float64)
    (-1.7976931348623157e+308, 1.7976931348623157e+308)
    >>> numeric_dtype_min_max(int)
    (-9223372036854775808, 9223372036854775807)
    """
    dtype = np.dtype(dtype)
    if np.issubdtype(dtype, np.integer):
        info = np.iinfo(dtype)
        min_ = dtype.type(info.min)
        max_ = dtype.type(info.max)
    elif np.issubdtype(dtype, np.inexact):
        info = np.finfo(dtype)
        min_ = info.min
        max_ = info.max
    elif np.issubdtype(dtype, np.dtype(bool)):
        min_ = dtype.type(False)
        max_ = dtype.type(True)
    else:
        raise ValueError(f"unsupported dtype {dtype!r}")
    return min_, max_
