# cython: boundscheck=False
# cython: wraparound=False
# cython: cdivision=True
#
# Author: Christopher Moody <chrisemoody@gmail.com>
# Author: Nick Travers <nickt@squareup.com>
# Implementation by Chris Moody & Nick Travers
# See http://homepage.tudelft.nl/19j49/t-SNE.html for reference
# implementations and papers describing the technique


import numpy as np
cimport numpy as np
from libc.stdio cimport printf
from libc.math cimport sqrt, log
from libc.stdlib cimport malloc, free
from cython.parallel cimport prange, parallel

from ..neighbors._quad_tree cimport _QuadTree

np.import_array()


cdef char* EMPTY_STRING = ""

cdef extern from "math.h":
    float fabsf(float x) nogil

# Smallest strictly positive value that can be represented by floating
# point numbers for different precision levels. This is useful to avoid
# taking the log of zero when computing the KL divergence.
cdef float FLOAT32_TINY = np.finfo(np.float32).tiny

# Useful to void division by zero or divergence to +inf.
cdef float FLOAT64_EPS = np.finfo(np.float64).eps

# This is effectively an ifdef statement in Cython
# It allows us to write printf debugging lines
# and remove them at compile time
cdef enum:
    DEBUGFLAG = 0

cdef extern from "time.h":
    # Declare only what is necessary from `tm` structure.
    ctypedef long clock_t
    clock_t clock() nogil
    double CLOCKS_PER_SEC


cdef float compute_gradient(float[:] val_P,
                            float[:, :] pos_reference,
                            np.int64_t[:] neighbors,
                            np.int64_t[:] indptr,
                            float[:, :] tot_force,
                            _QuadTree qt,
                            float theta,
                            int dof,
                            long start,
                            long stop,
                            bint compute_error,
                            int num_threads) nogil:
    # Having created the tree, calculate the gradient
    # in two components, the positive and negative forces
    cdef:
        long i, coord
        int ax
        long n_samples = pos_reference.shape[0]
        int n_dimensions = qt.n_dimensions
        clock_t t1 = 0, t2 = 0
        double sQ
        float error
        int take_timing = 1 if qt.verbose > 15 else 0

    if qt.verbose > 11:
        printf("[t-SNE] Allocating %li elements in force arrays\n",
                n_samples * n_dimensions * 2)
    cdef float* neg_f = <float*> malloc(sizeof(float) * n_samples * n_dimensions)
    cdef float* pos_f = <float*> malloc(sizeof(float) * n_samples * n_dimensions)

    if take_timing:
        t1 = clock()
    sQ = compute_gradient_negative(pos_reference, neg_f, qt, dof, theta, start,
                                   stop, num_threads)
    if take_timing:
        t2 = clock()
        printf("[t-SNE] Computing negative gradient: %e ticks\n", ((float) (t2 - t1)))

    if take_timing:
        t1 = clock()
    error = compute_gradient_positive(val_P, pos_reference, neighbors, indptr,
                                      pos_f, n_dimensions, dof, sQ, start,
                                      qt.verbose, compute_error, num_threads)
    if take_timing:
        t2 = clock()
        printf("[t-SNE] Computing positive gradient: %e ticks\n",
               ((float) (t2 - t1)))
    for i in prange(start, n_samples, nogil=True, num_threads=num_threads,
                    schedule='static'):
        for ax in range(n_dimensions):
            coord = i * n_dimensions + ax
            tot_force[i, ax] = pos_f[coord] - (neg_f[coord] / sQ)

    free(neg_f)
    free(pos_f)
    return error


cdef float compute_gradient_positive(float[:] val_P,
                                     float[:, :] pos_reference,
                                     np.int64_t[:] neighbors,
                                     np.int64_t[:] indptr,
                                     float* pos_f,
                                     int n_dimensions,
                                     int dof,
                                     double sum_Q,
                                     np.int64_t start,
                                     int verbose,
                                     bint compute_error,
                                     int num_threads) nogil:
    # Sum over the following expression for i not equal to j
    # grad_i = p_ij (1 + ||y_i - y_j||^2)^-1 (y_i - y_j)
    # This is equivalent to compute_edge_forces in the authors' code
    # It just goes over the nearest neighbors instead of all the data points
    # (unlike the non-nearest neighbors version of `compute_gradient_positive')
    cdef:
        int ax
        long i, j, k
        long n_samples = indptr.shape[0] - 1
        float C = 0.0
        float dij, qij, pij
        float exponent = (dof + 1.0) / 2.0
        float float_dof = (float) (dof)
        float* buff
        clock_t t1 = 0, t2 = 0
        float dt

    if verbose > 10:
        t1 = clock()

    with nogil, parallel(num_threads=num_threads):
        # Define private buffer variables
        buff = <float *> malloc(sizeof(float) * n_dimensions)

        for i in prange(start, n_samples, schedule='static'):
            # Init the gradient vector
            for ax in range(n_dimensions):
                pos_f[i * n_dimensions + ax] = 0.0
            # Compute the positive interaction for the nearest neighbors
            for k in range(indptr[i], indptr[i+1]):
                j = neighbors[k]
                dij = 0.0
                pij = val_P[k]
                for ax in range(n_dimensions):
                    buff[ax] = pos_reference[i, ax] - pos_reference[j, ax]
                    dij += buff[ax] * buff[ax]
                qij = float_dof / (float_dof + dij)
                if dof != 1:  # i.e. exponent != 1
                    qij = qij ** exponent
                dij = pij * qij

                # only compute the error when needed
                if compute_error:
                    qij = qij / sum_Q
                    C += pij * log(max(pij, FLOAT32_TINY) \
                        / max(qij, FLOAT32_TINY))
                for ax in range(n_dimensions):
                    pos_f[i * n_dimensions + ax] += dij * buff[ax]

        free(buff)
    if verbose > 10:
        t2 = clock()
        dt = ((float) (t2 - t1))
        printf("[t-SNE] Computed error=%1.4f in %1.1e ticks\n", C, dt)
    return C


cdef double compute_gradient_negative(float[:, :] pos_reference,
                                      float* neg_f,
                                      _QuadTree qt,
                                      int dof,
                                      float theta,
                                      long start,
                                      long stop,
                                      int num_threads) nogil:
    if stop == -1:
        stop = pos_reference.shape[0]
    cdef:
        int ax
        int n_dimensions = qt.n_dimensions
        int offset = n_dimensions + 2
        long i, j, idx
        long n = stop - start
        long dta = 0
        long dtb = 0
        float size, dist2s, mult
        float exponent = (dof + 1.0) / 2.0
        float float_dof = (float) (dof)
        double qijZ, sum_Q = 0.0
        float* force
        float* neg_force
        float* pos
        clock_t t1 = 0, t2 = 0, t3 = 0
        int take_timing = 1 if qt.verbose > 20 else 0


    with nogil, parallel(num_threads=num_threads):
        # Define thread-local buffers
        summary = <float*> malloc(sizeof(float) * n * offset)
        pos = <float *> malloc(sizeof(float) * n_dimensions)
        force = <float *> malloc(sizeof(float) * n_dimensions)
        neg_force = <float *> malloc(sizeof(float) * n_dimensions)

        for i in prange(start, stop, schedule='static'):
            # Clear the arrays
            for ax in range(n_dimensions):
                force[ax] = 0.0
                neg_force[ax] = 0.0
                pos[ax] = pos_reference[i, ax]

            # Find which nodes are summarizing and collect their centers of mass
            # deltas, and sizes, into vectorized arrays
            if take_timing:
                t1 = clock()
            idx = qt.summarize(pos, summary, theta*theta)
            if take_timing:
                t2 = clock()
            # Compute the t-SNE negative force
            # for the digits dataset, walking the tree
            # is about 10-15x more expensive than the
            # following for loop
            for j in range(idx // offset):

                dist2s = summary[j * offset + n_dimensions]
                size = summary[j * offset + n_dimensions + 1]
                qijZ = float_dof / (float_dof + dist2s)  # 1/(1+dist)
                if dof != 1:  # i.e. exponent != 1
                    qijZ = qijZ ** exponent

                sum_Q += size * qijZ   # size of the node * q
                mult = size * qijZ * qijZ
                for ax in range(n_dimensions):
                    neg_force[ax] += mult * summary[j * offset + ax]
            if take_timing:
                t3 = clock()
            for ax in range(n_dimensions):
                neg_f[i * n_dimensions + ax] = neg_force[ax]
            if take_timing:
                dta += t2 - t1
                dtb += t3 - t2
        free(pos)
        free(force)
        free(neg_force)
        free(summary)
    if take_timing:
        printf("[t-SNE] Tree: %li clock ticks | ", dta)
        printf("Force computation: %li clock ticks\n", dtb)

    # Put sum_Q to machine EPSILON to avoid divisions by 0
    sum_Q = max(sum_Q, FLOAT64_EPS)
    return sum_Q


def gradient(float[:] val_P,
             float[:, :] pos_output,
             np.int64_t[:] neighbors,
             np.int64_t[:] indptr,
             float[:, :] forces,
             float theta,
             int n_dimensions,
             int verbose,
             int dof=1,
             long skip_num_points=0,
             bint compute_error=1,
             int num_threads=1):
    # This function is designed to be called from external Python
    # it passes the 'forces' array by reference and fills thats array
    # up in-place
    cdef float C
    cdef int n
    n = pos_output.shape[0]
    assert val_P.itemsize == 4
    assert pos_output.itemsize == 4
    assert forces.itemsize == 4
    m = "Forces array and pos_output shapes are incompatible"
    assert n == forces.shape[0], m
    m = "Pij and pos_output shapes are incompatible"
    assert n == indptr.shape[0] - 1, m
    if verbose > 10:
        printf("[t-SNE] Initializing tree of n_dimensions %i\n", n_dimensions)
    cdef _QuadTree qt = _QuadTree(pos_output.shape[1], verbose)
    if verbose > 10:
        printf("[t-SNE] Inserting %li points\n", pos_output.shape[0])
    qt.build_tree(pos_output)
    if verbose > 10:
        # XXX: format hack to workaround lack of `const char *` type
        # in the generated C code that triggers error with gcc 4.9
        # and -Werror=format-security
        printf("[t-SNE] Computing gradient\n%s", EMPTY_STRING)

    C = compute_gradient(val_P, pos_output, neighbors, indptr, forces,
                         qt, theta, dof, skip_num_points, -1, compute_error,
                         num_threads)

    if verbose > 10:
        # XXX: format hack to workaround lack of `const char *` type
        # in the generated C code
        # and -Werror=format-security
        printf("[t-SNE] Checking tree consistency\n%s", EMPTY_STRING)
    m = "Tree consistency failed: unexpected number of points on the tree"
    assert qt.cells[0].cumulative_size == qt.n_points, m
    if not compute_error:
        C = np.nan
    return C
