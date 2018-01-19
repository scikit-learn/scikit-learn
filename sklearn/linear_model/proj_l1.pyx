# encoding: utf-8
# cython: cdivision=True
# cython: boundscheck=False
# cython: wraparound=False
#
# Author: Arthur Mensch <arthur.mensch@inria.fr>,
#         Elvis Dohmatob <gmdopp@gmail.com>
# License: BSD

import warnings
from libc.math cimport fabs, sqrt
cimport numpy as np
import numpy as np
from cython cimport floating
from utils cimport fsign
from blas_api cimport fused_scal, fused_copy


cdef inline void swap(floating *b, unsigned int i, unsigned int j,
                      floating * buf) nogil:
    buf[0] = b[i]
    b[i] = b[j]
    b[j] = buf[0]
    return


cdef void enet_projection(unsigned int m, floating *v, floating *out, floating radius,
                          floating l1_ratio) nogil:
    """I borrowed and adapted this piece of projection code from Arthur's modl.utils.
    Hash-tag don't reinvent the wheel!"""
    cdef unsigned int i
    cdef unsigned int j
    cdef unsigned int size_U
    cdef unsigned int start_U
    cdef floating buf = 0
    cdef floating gamma
    cdef unsigned int pivot
    cdef unsigned int rho
    cdef unsigned int drho
    cdef floating s
    cdef floating ds
    cdef floating a
    cdef floating d
    cdef floating c
    cdef floating l
    cdef floating norm = 0
    cdef floating tmp

    # XXX this is nasty, we should be doing such corner-case handling explicitly :/
    if radius == 0:
        for i in range(m):
            out[i] = 0

    # L2 projection
    if l1_ratio == 0:
        for i in range(m):
            norm += v[i] ** 2
        if norm <= radius:
            norm = 1
        else:
            norm = sqrt(norm / radius)
        for i in range(m):
            out[i] = v[i] / norm
    else:
        # Scaling by 1 / l1_ratio
        gamma = 2 / l1_ratio - 2
        radius /= l1_ratio
        # Preparing data
        for j in range(m):
            out[j] = fabs(v[j])
            norm += out[j] * (1 + gamma / 2 * out[j])
        if norm <= radius:
            for i in range(m):
                out[i] = v[i]
        else:
            # s and rho computation
            s = 0
            rho = 0
            start_U = 0
            size_U = m
            while size_U > 0:
                pivot = start_U + size_U / 2
                # Putting pivot at the beginning
                swap(out, pivot, start_U, &buf)
                pivot = start_U
                drho = 1
                ds = out[pivot] * (1 + gamma / 2 * out[pivot])
                # Ordering : [pivot, >=, <], using Lobato quicksort
                for i in range(start_U + 1, start_U + size_U):
                    if out[i] >= out[pivot]:
                        ds += out[i] * (1 + gamma / 2 * out[i])
                        swap(out, i, start_U + drho, &buf)
                        drho += 1
                if s + ds - (rho + drho) * (1 + gamma / 2 * out[pivot])\
                        * out[pivot] < radius * (1 + gamma * out[pivot]) ** 2:
                    # U <- L : [<]
                    start_U += drho
                    size_U -= drho
                    rho += drho
                    s += ds
                else:
                    # U <- G \ k : [>=]
                    start_U += 1
                    size_U = drho - 1

            # Projection
            if gamma != 0:
                a = gamma ** 2 * radius + gamma * rho * .5
                d = 2 * radius * gamma + rho
                c = radius - s
                l = (-d + sqrt(d ** 2 - 4 * a * c)) / (2 * a)
            else:
                l = (s - radius) / rho
            for i in range(m):
                tmp = fabs(v[i]) - l
                if tmp > 0.:
                    out[i] = fsign(v[i]) * tmp / (1. + l * gamma)
                else:
                    out[i] = 0.


cdef inline void proj_l1(int n, floating *w, floating reg,
                         floating ajj) nogil:
    """Computes (in-place)

        argmin .5 * ||z - w / ajj||_2^2 subject to ||z||_1 <= reg
    """
    cdef int k
    cdef floating scaling
    cdef floating[:] out
    with gil:
        if floating is float:
            dtype = np.float32
        else:
            dtype = np.float64

    # some scaling tricks
    if ajj == 0. or reg == 0.:
        for k in range(n):
            w[k] = 0.
    else:
        if ajj != 1.:
            scaling = 1. / ajj
            fused_scal(n, scaling, w, 1)
        with gil:
            out = np.zeros(n, dtype=dtype)
        enet_projection(n, w, &out[0], reg, 1.)
        fused_copy(n, &out[0], 1, w, 1)

