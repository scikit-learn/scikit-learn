# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from ._typedefs cimport uint32_t


cdef inline uint32_t DEFAULT_SEED = 1

cdef enum:
    # Max value for our rand_r replacement (near the bottom).
    # We don't use RAND_MAX because it's different across platforms and
    # particularly tiny on Windows/MSVC.
    # It corresponds to the maximum representable value for
    # 32-bit signed integers (i.e. 2^31 - 1).
    RAND_R_MAX = 2147483647


# rand_r replacement using a 32bit XorShift generator
# See http://www.jstatsoft.org/v08/i14/paper for details
cdef inline uint32_t our_rand_r(uint32_t* seed) nogil:
    """Generate a pseudo-random np.uint32 from a np.uint32 seed"""
    # seed shouldn't ever be 0.
    if (seed[0] == 0):
        seed[0] = DEFAULT_SEED

    seed[0] ^= <uint32_t>(seed[0] << 13)
    seed[0] ^= <uint32_t>(seed[0] >> 17)
    seed[0] ^= <uint32_t>(seed[0] << 5)

    # Use the modulo to make sure that we don't return a values greater than the
    # maximum representable value for signed 32bit integers (i.e. 2^31 - 1).
    # Note that the parenthesis are needed to avoid overflow: here
    # RAND_R_MAX is cast to uint32_t before 1 is added.
    return seed[0] % ((<uint32_t>RAND_R_MAX) + 1)
