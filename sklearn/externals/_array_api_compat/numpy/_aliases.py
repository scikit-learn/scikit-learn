from __future__ import annotations

from functools import partial

from ..common import _aliases

from .._internal import get_xp

asarray = asarray_numpy = partial(_aliases._asarray, namespace='numpy')
asarray.__doc__ = _aliases._asarray.__doc__
del partial

import numpy as np
bool = np.bool_

# Basic renames
acos = np.arccos
acosh = np.arccosh
asin = np.arcsin
asinh = np.arcsinh
atan = np.arctan
atan2 = np.arctan2
atanh = np.arctanh
bitwise_left_shift = np.left_shift
bitwise_invert = np.invert
bitwise_right_shift = np.right_shift
concat = np.concatenate
pow = np.power

arange = get_xp(np)(_aliases.arange)
empty = get_xp(np)(_aliases.empty)
empty_like = get_xp(np)(_aliases.empty_like)
eye = get_xp(np)(_aliases.eye)
full = get_xp(np)(_aliases.full)
full_like = get_xp(np)(_aliases.full_like)
linspace = get_xp(np)(_aliases.linspace)
ones = get_xp(np)(_aliases.ones)
ones_like = get_xp(np)(_aliases.ones_like)
zeros = get_xp(np)(_aliases.zeros)
zeros_like = get_xp(np)(_aliases.zeros_like)
UniqueAllResult = get_xp(np)(_aliases.UniqueAllResult)
UniqueCountsResult = get_xp(np)(_aliases.UniqueCountsResult)
UniqueInverseResult = get_xp(np)(_aliases.UniqueInverseResult)
unique_all = get_xp(np)(_aliases.unique_all)
unique_counts = get_xp(np)(_aliases.unique_counts)
unique_inverse = get_xp(np)(_aliases.unique_inverse)
unique_values = get_xp(np)(_aliases.unique_values)
astype = _aliases.astype
std = get_xp(np)(_aliases.std)
var = get_xp(np)(_aliases.var)
permute_dims = get_xp(np)(_aliases.permute_dims)
reshape = get_xp(np)(_aliases.reshape)
argsort = get_xp(np)(_aliases.argsort)
sort = get_xp(np)(_aliases.sort)
sum = get_xp(np)(_aliases.sum)
prod = get_xp(np)(_aliases.prod)
ceil = get_xp(np)(_aliases.ceil)
floor = get_xp(np)(_aliases.floor)
trunc = get_xp(np)(_aliases.trunc)
matmul = get_xp(np)(_aliases.matmul)
matrix_transpose = get_xp(np)(_aliases.matrix_transpose)
tensordot = get_xp(np)(_aliases.tensordot)
vecdot = get_xp(np)(_aliases.vecdot)
isdtype = get_xp(np)(_aliases.isdtype)

__all__ = _aliases.__all__ + ['asarray', 'asarray_numpy', 'bool', 'acos',
                              'acosh', 'asin', 'asinh', 'atan', 'atan2',
                              'atanh', 'bitwise_left_shift', 'bitwise_invert',
                              'bitwise_right_shift', 'concat', 'pow']
