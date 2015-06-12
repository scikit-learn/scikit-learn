"""Projection on the elastic-net ball (Cython version)
**References:**

J. Mairal, F. Bach, J. Ponce, G. Sapiro, 2009: Online dictionary learning
for sparse coding (http://www.di.ens.fr/sierra/pdfs/icml09.pdf)
"""

# Author: Arthur Mensch
# License: BSD

cimport cython
from cython cimport view

cimport numpy as np
from libc.math cimport sqrt

from libc.stdlib cimport malloc, free

ctypedef np.float64_t DOUBLE
ctypedef np.uint32_t UINT32_t
ctypedef np.uint8_t UINT8_t

cdef enum:
    MASKED = 0
    UNMASKED = 1
    GREATER = 2
    LOWER = 3

cdef enum:
    # Max value for our rand_r replacement (near the bottom).
    # We don't use RAND_MAX because it's different across platforms and
    # particularly tiny on Windows/MSVC.
    RAND_R_MAX = 0x7FFFFFFF


# The following two functions are shamelessly copied from the tree code.
@cython.cdivision(True)
cdef inline UINT32_t our_rand_r(UINT32_t* seed) nogil:
    seed[0] ^= <UINT32_t>(seed[0] << 13)
    seed[0] ^= <UINT32_t>(seed[0] >> 17)
    seed[0] ^= <UINT32_t>(seed[0] << 5)

    return seed[0] % (<UINT32_t>RAND_R_MAX + 1)


@cython.cdivision(True)
cdef inline UINT32_t _rand_int(UINT32_t end, UINT32_t* random_state) nogil:
    """Generate a random integer in [0; end)."""
    return our_rand_r(random_state) % end


cdef inline double _positive(double a) nogil:
    if a > 0:
        return a
    else:
        return 0


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.initializedcheck(False)
cdef double _masked_enet_norm(DOUBLE[:] v,
                                      UINT8_t * mask, int size,
                                      int masker, double gamma) nogil:
    cdef double res = 0
    for i in range(size):
        if mask[i] == masker:
            res += v[i] * (1. + gamma * v[i] / 2.)
    return res


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.initializedcheck(False)
cdef int _choose_index_in_mask(UINT8_t * mask, int size,
                               UINT32_t* random_state) nogil:
    cdef int rand_i = _rand_int(size, random_state)
    cdef int i
    for i in range(size):
        if mask[i]:
            if not rand_i:
                return i
            rand_i -= 1


# XXX: Implementation should be faster
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.initializedcheck(False)
@cython.nonecheck(False)
cpdef enet_projection(DOUBLE[:] v, double radius, double l1_gamma):
    """
    Project a vector on the elastic-net ball

    Parameters
    ---------------------------------------------
    res: double memory-view,
        Vector in which to store the projection
    v: double memory-view,
        Vector to project
    radius: float,
        Radius of the elastic-net ball
    l1_gamma: float,
        Ratio of L1 norm (between 0 and 1)
    """
    cdef int n = v.shape[0]
    cdef double gamma
    cdef double norm
    cdef double s
    cdef double rho
    cdef int true_size
    cdef int G_true_size
    cdef int L_true_size
    cdef double l
    cdef double a
    cdef double b_
    cdef double c
    cdef double d_s
    cdef int i
    cdef double norm2
    cdef double sign
    cdef UINT32_t random_state = 0
    cdef UINT8_t * mask = <UINT8_t *>malloc(n * sizeof(UINT8_t))
    cdef DOUBLE[:] res = view.array(shape=(n,),
                                    itemsize=sizeof(DOUBLE), format="d")

    # projection onto the l2 ball
    if l1_gamma == 0.:
        norm = 0.
        for i in range(n):
            norm += v[i] ** 2
        norm = sqrt(norm) / radius
        if norm <= 1:
            res[:] = v[:]
        else:
            for i in range(n):
                res[i] = v[i] / norm
        return res

    gamma = 2 / l1_gamma

    # res = np.sign(v)
    # v = np.abs(v)
    for i in range(n):
        if v[i] >= 0:
            res[i] = 1
        else:
            res[i] = -1
            v[i] *= -1
        mask[i] = UNMASKED
    norm = _masked_enet_norm(v, mask, n, UNMASKED, gamma)
    radius /= l1_gamma

    # Return res if it is already in the elastic-net ball
    if norm <= radius:
        res[:] = v[:]
        return res

    # Main loop for finding norm
    s = 0
    rho = 0
    true_size = n
    while true_size:
        idx = _choose_index_in_mask(mask, n, &random_state)
        G_true_size = 0
        L_true_size = 0
        for i in range(n):
            if mask[i]:
                if v[i] >= v[idx]:
                    mask[i] = GREATER
                    G_true_size += 1
                else:
                    mask[i] = LOWER
                    L_true_size += 1
        d_s = _masked_enet_norm(v, mask, n, GREATER, gamma)
        if s + d_s - (rho + G_true_size) * (1. + gamma / 2. * v[idx]) * v[idx]\
                < radius * (1. + gamma * v[idx]) ** 2:
            s += d_s
            rho += G_true_size
            for i in range(n):
                if mask[i] == GREATER:
                    mask[i] = MASKED
            true_size = L_true_size
        else:
            for i in range(n):
                if mask[i] == LOWER:
                    mask[i] = MASKED
            mask[idx] = MASKED
            true_size = G_true_size - 1
    # if l1_gamma = np.inf
    if gamma != 0.:
        a = gamma ** 2 * radius + gamma * rho * 0.5
        b_ = 2 * radius * gamma + rho
        c = radius - s
        l = (-b_ + sqrt(b_ ** 2 - 4 * a * c)) / (2*a)
    else:
        l = (s - radius) / rho
    for i in range(n):
        sign = res[i]
        res[i] *= _positive(v[i] - l) / ( 1 + l * gamma)
        v[i] *= sign
    free(mask)

    return res


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.initializedcheck(False)
cpdef double enet_norm(DOUBLE[:] v, double l1_gamma=0.1) nogil:
    """Returns the elastic net norm of a vector

    Parameters
    -----------------------------------------
    v: double memory-view,
        Vector

    l1_gamma: float,
        Ratio of l1 norm (between 0 and 1)

    Returns
    ------------------------------------------
    norm: float,
        Elastic-net norm
    """
    cdef int n = v.shape[0]
    cdef double res = 0
    cdef double b_ = 0
    cdef int i
    for i in range(n):
        b_ = v[i] if v[i] > 0. else - v[i]
        res += b_ * (l1_gamma + b_)
    return res