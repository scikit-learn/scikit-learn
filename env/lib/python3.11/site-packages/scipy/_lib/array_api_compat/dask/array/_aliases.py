from __future__ import annotations

from ...common import _aliases
from ...common._helpers import _check_device

from ..._internal import get_xp

import numpy as np
from numpy import (
    # Constants
    e,
    inf,
    nan,
    pi,
    newaxis,
    # Dtypes
    bool_ as bool,
    float32,
    float64,
    int8,
    int16,
    int32,
    int64,
    uint8,
    uint16,
    uint32,
    uint64,
    complex64,
    complex128,
    iinfo,
    finfo,
    can_cast,
    result_type,
)

from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from typing import Optional, Union

    from ...common._typing import Device, Dtype, Array

import dask.array as da

isdtype = get_xp(np)(_aliases.isdtype)
astype = _aliases.astype

# Common aliases

# This arange func is modified from the common one to
# not pass stop/step as keyword arguments, which will cause
# an error with dask

# TODO: delete the xp stuff, it shouldn't be necessary
def _dask_arange(
    start: Union[int, float],
    /,
    stop: Optional[Union[int, float]] = None,
    step: Union[int, float] = 1,
    *,
    xp,
    dtype: Optional[Dtype] = None,
    device: Optional[Device] = None,
    **kwargs,
) -> Array:
    _check_device(xp, device)
    args = [start]
    if stop is not None:
        args.append(stop)
    else:
        # stop is None, so start is actually stop
        # prepend the default value for start which is 0
        args.insert(0, 0)
    args.append(step)
    return xp.arange(*args, dtype=dtype, **kwargs)

arange = get_xp(da)(_dask_arange)
eye = get_xp(da)(_aliases.eye)

from functools import partial
asarray = partial(_aliases._asarray, namespace='dask.array')
asarray.__doc__ = _aliases._asarray.__doc__

linspace = get_xp(da)(_aliases.linspace)
eye = get_xp(da)(_aliases.eye)
UniqueAllResult = get_xp(da)(_aliases.UniqueAllResult)
UniqueCountsResult = get_xp(da)(_aliases.UniqueCountsResult)
UniqueInverseResult = get_xp(da)(_aliases.UniqueInverseResult)
unique_all = get_xp(da)(_aliases.unique_all)
unique_counts = get_xp(da)(_aliases.unique_counts)
unique_inverse = get_xp(da)(_aliases.unique_inverse)
unique_values = get_xp(da)(_aliases.unique_values)
permute_dims = get_xp(da)(_aliases.permute_dims)
std = get_xp(da)(_aliases.std)
var = get_xp(da)(_aliases.var)
empty = get_xp(da)(_aliases.empty)
empty_like = get_xp(da)(_aliases.empty_like)
full = get_xp(da)(_aliases.full)
full_like = get_xp(da)(_aliases.full_like)
ones = get_xp(da)(_aliases.ones)
ones_like = get_xp(da)(_aliases.ones_like)
zeros = get_xp(da)(_aliases.zeros)
zeros_like = get_xp(da)(_aliases.zeros_like)
reshape = get_xp(da)(_aliases.reshape)
matrix_transpose = get_xp(da)(_aliases.matrix_transpose)
vecdot = get_xp(da)(_aliases.vecdot)

nonzero = get_xp(da)(_aliases.nonzero)
sum = get_xp(np)(_aliases.sum)
prod = get_xp(np)(_aliases.prod)
ceil = get_xp(np)(_aliases.ceil)
floor = get_xp(np)(_aliases.floor)
trunc = get_xp(np)(_aliases.trunc)
matmul = get_xp(np)(_aliases.matmul)
tensordot = get_xp(np)(_aliases.tensordot)

from dask.array import (
    # Element wise aliases
    arccos as acos,
    arccosh as acosh,
    arcsin as asin,
    arcsinh as asinh,
    arctan as atan,
    arctan2 as atan2,
    arctanh as atanh,
    left_shift as bitwise_left_shift,
    right_shift as bitwise_right_shift,
    invert as bitwise_invert,
    power as pow,
    # Other
    concatenate as concat,
)

# exclude these from all since
_da_unsupported = ['sort', 'argsort']

common_aliases = [alias for alias in _aliases.__all__ if alias not in _da_unsupported]

__all__ = common_aliases + ['asarray', 'bool', 'acos',
                            'acosh', 'asin', 'asinh', 'atan', 'atan2',
                            'atanh', 'bitwise_left_shift', 'bitwise_invert',
                            'bitwise_right_shift', 'concat', 'pow',
                            'e', 'inf', 'nan', 'pi', 'newaxis', 'float32', 'float64', 'int8',
                            'int16', 'int32', 'int64', 'uint8', 'uint16', 'uint32', 'uint64',
                            'complex64', 'complex128', 'iinfo', 'finfo', 'can_cast', 'result_type']

_all_ignore = ['get_xp', 'da', 'partial', 'common_aliases', 'np']
