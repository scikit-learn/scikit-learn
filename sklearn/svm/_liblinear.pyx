"""
Wrapper for liblinear

Author: fabian.pedregosa@inria.fr
"""

import  numpy as np
cimport numpy as np

from ..utils._cython_blas cimport _dot, _axpy, _scal, _nrm2

include "_liblinear.pxi"

np.import_array()


def train_wrap(X, np.ndarray[np.float64_t, ndim=1, mode='c'] Y,
               bint is_sparse, int solver_type, double eps, double bias,
               double C, np.ndarray[np.float64_t, ndim=1] class_weight,
               int max_iter, unsigned random_seed, double epsilon,
               np.ndarray[np.float64_t, ndim=1, mode='c'] sample_weight):
    cdef parameter *param
    cdef problem *problem
    cdef model *model
    cdef char_const_ptr error_msg
    cdef int len_w

    if is_sparse:
        problem = csr_set_problem(
                (<np.ndarray>X.data).data, X.dtype == np.float64,
                (<np.ndarray[np.int32_t,   ndim=1, mode='c']>X.indices).data,
                (<np.ndarray[np.int32_t,   ndim=1, mode='c']>X.indptr).data,
                (<np.int32_t>X.shape[0]), (<np.int32_t>X.shape[1]),
                (<np.int32_t>X.nnz), bias, sample_weight.data, Y.data)
    else:
        problem = set_problem(
                (<np.ndarray>X).data, X.dtype == np.float64,
                (<np.int32_t>X.shape[0]), (<np.int32_t>X.shape[1]),
                (<np.int32_t>np.count_nonzero(X)), bias, sample_weight.data,
                Y.data)

    cdef np.ndarray[np.int32_t, ndim=1, mode='c'] \
        class_weight_label = np.arange(class_weight.shape[0], dtype=np.intc)
    param = set_parameter(solver_type, eps, C, class_weight.shape[0],
                          class_weight_label.data, class_weight.data,
                          max_iter, random_seed, epsilon)

    error_msg = check_parameter(problem, param)
    if error_msg:
        free_problem(problem)
        free_parameter(param)
        raise ValueError(error_msg)
    
    cdef BlasFunctions blas_functions
    blas_functions.dot = _dot[double]
    blas_functions.axpy = _axpy[double]
    blas_functions.scal = _scal[double]
    blas_functions.nrm2 = _nrm2[double]

    # early return
    with nogil:
        model = train(problem, param, &blas_functions)

    ### FREE
    free_problem(problem)
    free_parameter(param)
    # destroy_param(param)  don't call this or it will destroy class_weight_label and class_weight

    # coef matrix holder created as fortran since that's what's used in liblinear
    cdef np.ndarray[np.float64_t, ndim=2, mode='fortran'] w
    cdef int nr_class = get_nr_class(model)

    cdef int labels_ = nr_class
    if nr_class == 2:
        labels_ = 1
    cdef np.ndarray[np.int32_t, ndim=1, mode='c'] n_iter = np.zeros(labels_, dtype=np.intc)
    get_n_iter(model, <int *>n_iter.data)

    cdef int nr_feature = get_nr_feature(model)
    if bias > 0: nr_feature = nr_feature + 1
    if nr_class == 2 and solver_type != 4:  # solver is not Crammer-Singer
        w = np.empty((1, nr_feature),order='F')
        copy_w(w.data, model, nr_feature)
    else:
        len_w = (nr_class) * nr_feature
        w = np.empty((nr_class, nr_feature),order='F')
        copy_w(w.data, model, len_w)

    free_and_destroy_model(&model)

    return w, n_iter


def set_verbosity_wrap(int verbosity):
    """
    Control verbosity of libsvm library
    """
    set_verbosity(verbosity)
