# encoding: utf-8
# cython: boundscheck=False
# cython: wraparound=False
# cython: cdivision=True
#
# Author: Elvis Dohmatob <gmdopp@gmail.com>
# License: BSD
#
# Synopsis: Implementation of block coordinate-descent (BCD) algorithm
#           for general penalized / constrained multi-task least-squares
#           regression problems. Out-of-the-box, it supports
#           - real and complex-valued data
#           - Gram mode (crucial when p << n) and non-Gram mode
#             (crucial when p >> n)
#           - a variety of penalties and constraints, including:
#             * L11 penalty (sum of L1 norms), as in multi-task Lasso;
#             * L21 penalty (sum of L2 norms), as in G-Lasso;
#             * L2INF constraint (L2 constraint on each block), as in vanilla
#               dictionary-learning;
#             * L1INF constraint (L1 constraint on each block), e.g as in
#               sparse PCA
#
#           It also supports user-defined penalties (via their prox operators).
#           For example, with this technology, it is straightforward to
#           implement otherwise non-trivial models like:
#           - BCD dictionary update in online structured dictionary-learning
#             with general (e.g see Dohmatob el al. "Learning brain regions
#             via large-scale online structured sparse dictionary-learning",
#             In NIPS 2016")
#           - Certain problems arising in data-driven compressive sensing
#
# XXX TODO:
#    - Add single-precision (32-bit) support. Not easy due to issues with
#      specialization of fused types. Long story short...
#    - Add support for sparse matrices and vectors
#    - Screening rules for LASSO-type problems (inexact Stanford, GapSafe,
#      etc.)

import warnings
from numbers import Number
cimport cython
cimport numpy as np
import numpy as np
ctypedef np.uint32_t UINT32_t
from types cimport floating, complexing
from blas_api cimport (CblasColMajor, fused_scal, fused_copy, fused_nrm2,
                       fused_dotu, fused_dotc, fused_geru, fused_axpy)
from dual_gap cimport _compute_dual_gap
from utils cimport fmax, abs_max, diff_abs_max, relu
from proj_l2 cimport proj_l2
from proj_l1 cimport proj_l1
from prox_l1 cimport prox_l1
from prox_l2 cimport prox_l2

np.import_array()

cdef extern from "complex.h" nogil:
    double creal(double complex)
    float crealf(float complex)


cdef inline UINT32_t our_rand_r(UINT32_t* seed) nogil:
    seed[0] ^= <UINT32_t>(seed[0] << 13)
    seed[0] ^= <UINT32_t>(seed[0] >> 17)
    seed[0] ^= <UINT32_t>(seed[0] << 5)
    return seed[0] % (<UINT32_t>RAND_R_MAX + 1)


cdef inline UINT32_t rand_int(UINT32_t end, UINT32_t* random_state) nogil:
    """Generate a random integer in [0; end)."""
    return our_rand_r(random_state) % end


cdef enum:
    # Max value for our rand_r replacement (near the bottom).
    # We don't use RAND_MAX because it's different across platforms and
    # particularly tiny on Windows/MSVC.
    RAND_R_MAX = 0x7FFFFFFF


# What's a proximal operator ?
# := argmin_z .5 * ||z - w / ajj||_2^2 + (reg / ajj) * pen(z)
# where ajj > 0 is a constant. ajj = 0 corresponds armin_z pen(z),  while
# ajj = 1 corresponds to classical definition of the prox.
ctypedef void (*PROX)(int n,
                      complexing *w,
                      floating reg,
                      floating ajj) nogil except *


def coordescendant(np.ndarray[complexing, ndim=2, mode="c"] W,
                   floating reg,
                   floating l2_reg,
                   np.ndarray[complexing, ndim=2, mode="fortran"] X_or_Gram,
                   np.ndarray[complexing, ndim=2, mode="fortran"] Y_or_Cov,
                   bint precomputed=1,
                   floating Y_norm2=np.nan,
                   object penalty_model=L11_PENALTY,
                   unsigned int max_iter=100,
                   floating tol=0,
                   bint random=0,
                   object rng=None,
                   bint pos=0,
                   bint emulate_sklearn_dl=False):
    """Multi-output penalized least-squares solver for both real- and complex-valued
    data.

    The optimization problem is of the form

        minimize .5 ||Y - XW||_F^2 + l2_reg ||W||_F^2 + reg \sum_j g(W_j),
           W

    where
    - g is l.s.c convex function with known prox (see .pxd header file),
    - If `precomputed` is True, then
          X_or_Gram is X^H X (shape n_features x n_features) and Y_or_Cov is X^HY
          (^H denotes conjugate transposition)
      else
          X_or_Gram is X (shape n_samples x n_features) and Y_or_Cov is Y
          (shape n_samples x n_targets)

    Parameters
    ----------
    emulate_sklearn_dl: optional, default False
        If True, then a weird update is used to update W[j] in `precomputed` mode,
        so as to emulate sklearn dictionary-learning (buggy?) algebraic manips
 
    Notes
    -----
    A variety of problems belong to this category in both their penalized and
    constrained versions. Just to name a few:
    - multi-task elastic-net
    - BCD dictionary update in online structured dictionary-learning
    - certain problems arising in data-driven compressive sensing

    """
    # some sanity checks
    if pos and not (complexing is double or complexing is float):
            raise TypeError("pos=True for complex data makes no sense")

    # specialization of fuzed types / functions
    if complexing is float:
        dtype = np.float32
        real_part = crealf
    elif complexing is double:
        dtype = np.float64
        real_part = creal
    elif complexing is complex:
        dtype = np.complex128
        real_part = creal
    else:
        dtype = np.complex64
        real_part = crealf

    # select an appropriate prox handle by model
    cdef PROX prox
    if hasattr(penalty_model, "__call__"):
        # user supplied their a callable for the prox operator
        user_prox = penalty_model
        prox = <PROX>0
    else:
        user_prox = None
        if penalty_model == NOP:
            prox = <PROX>0
        elif penalty_model == L11_PENALTY:
            prox = prox_l1
        elif penalty_model == L21_PENALTY:
            prox = prox_l2
        elif penalty_model == L2INF_CONSTRAINT:
            prox = proj_l2
        elif penalty_model == L1INF_CONSTRAINT:
            prox = proj_l1
        else:
            raise NotImplementedError("penalty_model=%s" % penalty_model)

        # sklearn only supports contraints (not penalties)
        if penalty_model > 0:
            emulate_sklearn_dl = False

    # misc
    cdef int penalty_model_int
    cdef bint dual_gap_available = False
    if penalty_model in [L11_PENALTY, L21_PENALTY]:
        dual_gap_available = True
        penalty_model_int = penalty_model
    cdef int n_samples = X_or_Gram.shape[0]
    cdef int n_features = X_or_Gram.shape[1]
    cdef int n_targets = Y_or_Cov.shape[1]
    cdef complexing[:] Wj = np.zeros(n_targets, dtype=dtype)
    cdef floating ajj
    cdef complexing alpha, beta
    cdef int X_size = n_samples * n_features
    cdef int W_size = n_features * n_targets
    cdef int R_size = n_samples * n_targets
    cdef int inc = 1
    cdef floating Wj_abs_max, d_Wj_abs_max, W_abs_max, d_W_abs_max
    cdef floating d_W_abs_tol = tol
    cdef int j, k
    cdef unsigned int n_iter = 0
    cdef UINT32_t rand_r_state_seed
    cdef UINT32_t *rand_r_state
    if random and rng is not None:
        rand_r_state_seed = rng.randint(0, RAND_R_MAX)
        rand_r_state = &rand_r_state_seed
    else:
        random = 0

    # pre-compute conjugate of X_or_Gram for latter use in *gemm
    cdef np.ndarray[complexing, ndim=2, mode="fortran"] X_or_Gram_conj
    if complexing is float or complexing is double:
        X_or_Gram_conj = X_or_Gram
    else:
        X_or_Gram_conj = X_or_Gram.conjugate()

    # pointers are faster than numpy array addresses
    cdef complexing *W_ptr = &W[0, 0]
    cdef complexing *Wj_ptr = &Wj[0]
    cdef complexing *X_or_Gram_ptr = &X_or_Gram[0, 0]
    cdef complexing *Y_or_Cov_ptr = &Y_or_Cov[0, 0]
    cdef complexing *X_or_Gram_conj_ptr = &X_or_Gram_conj[0, 0]

    # some sanity checking
    if reg == l2_reg == 0:
        warnings.warn("Coordinate descent with reg = l2_reg = 0 may lead"
                      " to unexpected results and is discouraged.")
    if len(W) != n_features or len(W[0]) != n_targets or len(Y_or_Cov) != n_samples or \
       len(Y_or_Cov[0]) != n_targets:
        raise ValueError("Inconsistent matrix dimensions")

    # initialize residuals
    if precomputed and n_samples != n_features:
        raise ValueError("You specified precomputed=True. But X_or_Gram "
                         "has shape (%i, %i) which is non-square" % (n_samples,
                                                                     n_features))
    cdef np.ndarray[complexing, ndim=2, mode="fortran"] R
    R = np.asfortranarray(Y_or_Cov - np.dot(X_or_Gram, W))
    cdef complexing *R_ptr = &R[0, 0]

    # stuff for computing dual gap
    cdef floating gap = tol + 1.
    cdef np.ndarray[complexing, ndim=2, mode="fortran"] Grad = np.zeros_like(W, order="F")
    cdef complexing *Grad_ptr = &Grad[0, 0]

    cdef np.ndarray[double, ndim=1] X_col_norms_squared = np.zeros(n_features)
    for j in range(n_features):
        if precomputed:
            X_col_norms_squared[j] = X_or_Gram[j, j].real
        else:
            X_col_norms_squared[j] = real_part(fused_dotc(n_samples,
                                                         X_or_Gram_ptr + j * n_samples,
                                                         inc,
                                                         X_or_Gram_ptr + j * n_samples,
                                                         inc))

    # main loop: the rest of code doesn't need the GIL
    with nogil:
        for n_iter in range(max_iter):
            W_abs_max = 0.
            d_W_abs_max = 0.
            for j in range(n_features):
                if random:
                    j = rand_int(n_features, rand_r_state)
                ajj = X_col_norms_squared[j]
                ajj += l2_reg

                # check for dead features
                if ajj == 0. and not emulate_sklearn_dl:
                    continue

                # store previous value of coefficients for this feature
                if n_targets == 1:
                    Wj_ptr[0] = W_ptr[j]
                else:
                    fused_copy(n_targets,
                               W_ptr + j * n_targets,
                               inc,
                               Wj_ptr,
                               inc)

                # rank-1 update: R += np.outer(X_or_Gram[:, j], W[j])
                if n_targets == 1:
                    alpha = W_ptr[j]
                    if alpha != 0.:
                        fused_axpy(n_samples,
                                   alpha,
                                   X_or_Gram_ptr + j * n_samples,
                                   inc,
                                   R_ptr,
                                   inc)
                else:
                    alpha = 1
                    if fused_nrm2(n_targets,
                                  W_ptr + j * n_targets,
                                  inc):
                        fused_geru(CblasColMajor,
                                   n_samples,
                                   n_targets,
                                   alpha,
                                   X_or_Gram_ptr + j * n_samples,
                                   inc,
                                   W_ptr + j * n_targets,
                                   inc,
                                   R_ptr,
                                   n_samples)


                # XXX The following step only makes sense if we're not running in
                # "precomputed" mode. However, for some reason, sklearn's
                # implementation of dict-learning BCD update does it regardless.
                if (not precomputed) or emulate_sklearn_dl:
                    # W[j] <- ajj * np.dot(X_or_Gram_conj[:, j], R)
                    if emulate_sklearn_dl:
                        alpha = ajj
                    else:
                        alpha = 1.
                    for k in range(n_targets):
                        W_ptr[j * n_targets + k] = alpha * fused_dotc(n_samples,
                                                                     X_or_Gram_ptr + j * n_samples,
                                                                     inc,
                                                                     R_ptr + k * n_samples,
                                                                     inc)
                else:
                    # copy: W[j] <- R[:, j]
                    if n_targets == 1:
                        W_ptr[j] = R_ptr[j]
                    else:
                        fused_copy(n_targets,
                                   R_ptr + j,
                                   n_samples,
                                   W_ptr + j * n_targets,
                                   inc)

                # proximal update
                if user_prox is not None:
                    # invoke user-supplied prox operator
                    # XXX it would be nice to do this without requiring the GIL
                    if pos:
                        relu(n_targets,
                             <floating *>W_ptr + j * n_targets)
                    with gil:
                        user_prox(W[j],
                                  reg,
                                  ajj)
                else:
                    # invoke a standard prox operator
                    if pos:
                        relu(n_targets, <floating *>W_ptr + j * n_targets)
                    if prox:
                        prox(n_targets,
                             W_ptr + j * n_targets,
                             reg,
                             ajj)
                    else:
                        # N.B.: W[j] /= ajj
                        if ajj != 0.:
                            alpha = 1. / ajj
                        else:
                            alpha = 0.
                        fused_scal(n_targets,
                                   alpha,
                                   W_ptr + j * n_targets,
                                   inc)

                # rank-1 update: R -= np.outer(X_or_Gram[:, j], W[j])
                if n_targets == 1:
                    alpha = -W_ptr[j]
                    if alpha != 0.:
                        fused_axpy(n_samples,
                                   alpha,
                                   X_or_Gram_ptr + j * n_samples,
                                   inc,
                                   R_ptr,
                                   inc)
                else:
                    alpha = -1
                    if fused_nrm2(n_targets,
                                  W_ptr + j * n_targets,
                                  inc):
                        fused_geru(CblasColMajor,
                                   n_samples,
                                   n_targets,
                                   alpha,
                                   X_or_Gram_ptr + j * n_samples,
                                   inc,
                                   W_ptr + j * n_targets,
                                   inc,
                                   R_ptr,
                                   n_samples)

                # update the maximum absolute coefficient
                d_Wj_abs_max = diff_abs_max(n_targets, W_ptr + j * n_targets, Wj_ptr)
                d_W_abs_max = fmax(d_W_abs_max, d_Wj_abs_max)
                Wj_abs_max = abs_max(n_targets, W_ptr + j * n_targets)
                W_abs_max = fmax(W_abs_max, Wj_abs_max)

            # check convergence
            if (W_abs_max == 0. or
                d_W_abs_max / W_abs_max < d_W_abs_tol or
                n_iter == max_iter - 1):
                # the biggest coordinate update of this iteration was smaller than
                # check the duality gap as ultimate stopping criterion
                if dual_gap_available:
                    gap = _compute_dual_gap(n_samples,
                                            n_features,
                                            n_targets,
                                            W_ptr,
                                            reg,
                                            l2_reg,
                                            X_or_Gram_conj_ptr,
                                            Y_or_Cov_ptr,
                                            R_ptr,
                                            Grad_ptr,
                                            Y_norm2,
                                            precomputed,
                                            penalty_model_int)
                else:
                    gap = tol

                # XXX TODO: check for very bad things like negative gap, etc.
                pass

                # exit if we reached desired tolerance
                if gap < tol:
                    break

    return W, gap, tol, n_iter + 1
