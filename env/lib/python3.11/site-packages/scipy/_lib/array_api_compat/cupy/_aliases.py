from __future__ import annotations

from functools import partial

import cupy as cp

from ..common import _aliases
from .._internal import get_xp

asarray = asarray_cupy = partial(_aliases._asarray, namespace='cupy')
asarray.__doc__ = _aliases._asarray.__doc__
del partial

bool = cp.bool_

# Basic renames
acos = cp.arccos
acosh = cp.arccosh
asin = cp.arcsin
asinh = cp.arcsinh
atan = cp.arctan
atan2 = cp.arctan2
atanh = cp.arctanh
bitwise_left_shift = cp.left_shift
bitwise_invert = cp.invert
bitwise_right_shift = cp.right_shift
concat = cp.concatenate
pow = cp.power

arange = get_xp(cp)(_aliases.arange)
empty = get_xp(cp)(_aliases.empty)
empty_like = get_xp(cp)(_aliases.empty_like)
eye = get_xp(cp)(_aliases.eye)
full = get_xp(cp)(_aliases.full)
full_like = get_xp(cp)(_aliases.full_like)
linspace = get_xp(cp)(_aliases.linspace)
ones = get_xp(cp)(_aliases.ones)
ones_like = get_xp(cp)(_aliases.ones_like)
zeros = get_xp(cp)(_aliases.zeros)
zeros_like = get_xp(cp)(_aliases.zeros_like)
UniqueAllResult = get_xp(cp)(_aliases.UniqueAllResult)
UniqueCountsResult = get_xp(cp)(_aliases.UniqueCountsResult)
UniqueInverseResult = get_xp(cp)(_aliases.UniqueInverseResult)
unique_all = get_xp(cp)(_aliases.unique_all)
unique_counts = get_xp(cp)(_aliases.unique_counts)
unique_inverse = get_xp(cp)(_aliases.unique_inverse)
unique_values = get_xp(cp)(_aliases.unique_values)
astype = _aliases.astype
std = get_xp(cp)(_aliases.std)
var = get_xp(cp)(_aliases.var)
permute_dims = get_xp(cp)(_aliases.permute_dims)
reshape = get_xp(cp)(_aliases.reshape)
argsort = get_xp(cp)(_aliases.argsort)
sort = get_xp(cp)(_aliases.sort)
nonzero = get_xp(cp)(_aliases.nonzero)
sum = get_xp(cp)(_aliases.sum)
prod = get_xp(cp)(_aliases.prod)
ceil = get_xp(cp)(_aliases.ceil)
floor = get_xp(cp)(_aliases.floor)
trunc = get_xp(cp)(_aliases.trunc)
matmul = get_xp(cp)(_aliases.matmul)
matrix_transpose = get_xp(cp)(_aliases.matrix_transpose)
tensordot = get_xp(cp)(_aliases.tensordot)

# These functions are completely new here. If the library already has them
# (i.e., numpy 2.0), use the library version instead of our wrapper.
if hasattr(cp, 'vecdot'):
    vecdot = cp.vecdot
else:
    vecdot = get_xp(cp)(_aliases.vecdot)
if hasattr(cp, 'isdtype'):
    isdtype = cp.isdtype
else:
    isdtype = get_xp(cp)(_aliases.isdtype)

__all__ = _aliases.__all__ + ['asarray', 'asarray_cupy', 'bool', 'acos',
                              'acosh', 'asin', 'asinh', 'atan', 'atan2',
                              'atanh', 'bitwise_left_shift', 'bitwise_invert',
                              'bitwise_right_shift', 'concat', 'pow']

_all_ignore = ['cp', 'get_xp']
