# Synopsis: Python wrappers for some cython functions for testing
# Author: Elvis Dohmatob <gmdopp@gmail.com>

cimport numpy as np
import numpy as np
from cython cimport floating
from blas_api cimport (fused_asum, fused_dot, fused_nrm2, fused_copy)
from utils cimport abs_max, arr_max, diff_abs_max, relu, fused_nrm2_squared
from dual_gap cimport compute_dual_gap
from coordescendant cimport L21_PENALTY
from prox_l1 cimport prox_l1
from proj_l1 cimport proj_l1, enet_projection
from prox_l2 cimport prox_l2
from proj_l2 cimport proj_l2


def _py_fused_dot(np.ndarray[floating, ndim=1] X,
                  np.ndarray[floating, ndim=1] Y):
    return fused_dot(len(X), &X[0], 1, &Y[0], 1)


def _py_fused_asum(np.ndarray[floating, ndim=1] w):
    return fused_asum(len(w), &w[0], 1)


def _py_fused_nrm2(np.ndarray[floating, ndim=1] w):
    return fused_nrm2(len(w), &w[0], 1)


def _py_fused_nrm2_squared(np.ndarray[floating, ndim=1] w):
    return fused_nrm2_squared(len(w), &w[0], 1)


def _py_arr_max(np.ndarray[floating, ndim=1] a):
    """For testing"""
    return arr_max(len(a), &a[0], 1)


def _py_diff_abs_max(np.ndarray[floating, ndim=1] a,
                     np.ndarray[floating, ndim=1] b):
    return diff_abs_max(len(a), &a[0], 1, &b[0], 1)


def _py_abs_max(np.ndarray[floating, ndim=1] a):
    return abs_max(len(a), &a[0], 1)
 

def _py_relu(np.ndarray[floating, ndim=1] a):
    relu(len(a), &a[0], 1)


def _py_fused_copy(np.ndarray[floating, ndim=1] a,
                   np.ndarray[floating, ndim=1] b):
    fused_copy(len(a), &a[0], 1, &b[0], 1)


def _py_compute_dual_gap(np.ndarray[floating, ndim=2, mode="c"] W,
                         floating reg, floating l2_reg,
                         np.ndarray[floating, ndim=2, mode="fortran"] X_or_Gram,
                         np.ndarray[floating, ndim=2, mode="fortran"] Y_or_Cov,
                         np.ndarray[floating, ndim=2, mode="fortran"] R,
                         np.ndarray[floating, ndim=2, mode="fortran"] Grad,
                         floating Y_norm2=np.nan, int precomputed=True,
                         int penalty_model=L21_PENALTY, bint positive=False):
    return compute_dual_gap(
        len(X_or_Gram), len(W), len(W[0]), &W[0, 0], reg, l2_reg,
        &X_or_Gram[0, 0], &Y_or_Cov[0, 0], &R[0, 0], &Grad[0, 0], Y_norm2,
        precomputed, penalty_model, positive)


def _py_prox_l1(np.ndarray[floating, ndim=1] w, floating reg,
                floating ajj=1.):
    prox_l1(len(w), &w[0], reg, ajj)


def _py_proj_l1(np.ndarray[floating, ndim=1] w, floating reg,
                floating ajj=1.):
    proj_l1(len(w), &w[0], reg, ajj)


def _py_prox_l2(np.ndarray[floating, ndim=1] w, floating reg,
                floating ajj=1.):
    prox_l2(len(w), &w[0], reg, ajj)


def _py_proj_l2(np.ndarray[floating, ndim=1] w, floating reg,
                floating ajj=1.):
    proj_l2(len(w), &w[0], reg, ajj)


def _py_enet_projection(floating[:] v, floating[:] out, floating radius,
                        floating l1_ratio):
    cdef unsigned int m = len(v)
    enet_projection(m, &v[0], &out[0], radius, l1_ratio)
