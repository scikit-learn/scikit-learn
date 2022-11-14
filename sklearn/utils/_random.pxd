# Authors: Arnaud Joly
#
# License: BSD 3 clause


cimport numpy as cnp
ctypedef cnp.npy_uint32 UINT32_t

cdef inline UINT32_t DEFAULT_SEED = 1

cpdef sample_without_replacement(cnp.int_t n_population,
                                 cnp.int_t n_samples,
                                 method=*,
                                 random_state=*)

# rand_r replacement using a 32bit XorShift generator
# See http://www.jstatsoft.org/v08/i14/paper for details
cdef inline UINT32_t our_rand_r(UINT32_t* seed) nogil:
    """Generate a pseudo-random np.uint32 from a np.uint32 seed"""
    # seed shouldn't ever be 0.
    if (seed[0] == 0): seed[0] = DEFAULT_SEED

    seed[0] ^= <UINT32_t>(seed[0] << 13)
    seed[0] ^= <UINT32_t>(seed[0] >> 17)
    seed[0] ^= <UINT32_t>(seed[0] << 5)

    return seed[0]
