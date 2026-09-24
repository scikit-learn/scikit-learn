# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from libc.string cimport memcpy
from cython cimport floating

from sklearn.utils._typedefs cimport uint32_t, uint64_t


cdef inline bint inlinable_isnan(floating x) noexcept nogil:
    """Check whether x is NaN.

    Prefer this over libc.math.isnan in hot loops: unlike that libm call,
    which some compilers/libc fail to inline, this is guaranteed to be
    inlined. See https://github.com/scikit-learn/scikit-learn/issues/34869.

    TODO: remove this helper in favor of libc.math.isnan when conda-forge bumps
    its minimal glibc version to 2.28, see:
    https://github.com/conda-forge/conda-forge.github.io/issues/2383
    """
    cdef uint32_t bits32
    cdef uint64_t bits64
    if floating is float:
        memcpy(&bits32, &x, sizeof(bits32))
        return (bits32 & 0x7fffffff) > 0x7f800000
    else:
        memcpy(&bits64, &x, sizeof(bits64))
        return (bits64 & <uint64_t>0x7fffffffffffffff) > <uint64_t>0x7ff0000000000000
