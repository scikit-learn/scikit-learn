"""
Wrapper for liblinear

Author: fabian.pedregosa@inria.fr
"""

import  numpy as np

from sklearn.utils._cython_blas cimport _dot, _axpy, _scal, _nrm2
from sklearn.utils._typedefs cimport float32_t, float64_t, int32_t

include "_liblinear.pxi"


def train_wrap(
    object X,
    const float64_t[::1] Y,
    bint is_sparse,
    int solver_type,
    double eps,
    double bias,
    double C,
    const float64_t[:] class_weight,
    int max_iter,
    unsigned random_seed,
    double epsilon,
    const float64_t[::1] sample_weight
):
    cdef parameter *param
    cdef problem *problem
    cdef model *model
    cdef char_const_ptr error_msg
    cdef int len_w
    cdef bint X_has_type_float64 = X.dtype == np.float64
    cdef char * X_data_bytes_ptr
    cdef const float64_t[::1] X_data_64
    cdef const float32_t[::1] X_data_32
    cdef const int32_t[::1] X_indices
    cdef const int32_t[::1] X_indptr

    if is_sparse:
        X_indices = X.indices
        X_indptr = X.indptr
        if X_has_type_float64:
            X_data_64 = X.data
            X_data_bytes_ptr = <char *> &X_data_64[0]
        else:
            X_data_32 = X.data
            X_data_bytes_ptr = <char *> &X_data_32[0]

        problem = csr_set_problem(
            X_data_bytes_ptr,
            X_has_type_float64,
            <char *> &X_indices[0],
            <char *> &X_indptr[0],
            (<int32_t>X.shape[0]),
            (<int32_t>X.shape[1]),
            (<int32_t>X.nnz),
            bias,
            <char *> &sample_weight[0],
            <char *> &Y[0]
        )
    else:
        X_as_1d_array = X.reshape(-1)
        if X_has_type_float64:
            X_data_64 = X_as_1d_array
            X_data_bytes_ptr = <char *> &X_data_64[0]
        else:
            X_data_32 = X_as_1d_array
            X_data_bytes_ptr = <char *> &X_data_32[0]

        problem = set_problem(
            X_data_bytes_ptr,
            X_has_type_float64,
            (<int32_t>X.shape[0]),
            (<int32_t>X.shape[1]),
            (<int32_t>np.count_nonzero(X)),
            bias,
            <char *> &sample_weight[0],
            <char *> &Y[0]
        )

    cdef int32_t[::1] class_weight_label = np.arange(class_weight.shape[0], dtype=np.intc)
    param = set_parameter(
        solver_type,
        eps,
        C,
        class_weight.shape[0],
        <char *> &class_weight_label[0] if class_weight_label.size > 0 else NULL,
        <char *> &class_weight[0] if class_weight.size > 0 else NULL,
        max_iter,
        random_seed,
        epsilon
    )

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

    # FREE
    free_problem(problem)
    free_parameter(param)
    # destroy_param(param)  don't call this or it will destroy class_weight_label and class_weight

    # coef matrix holder created as fortran since that's what's used in liblinear
    cdef float64_t[::1, :] w
    cdef int nr_class = get_nr_class(model)

    cdef int labels_ = nr_class
    if nr_class == 2:
        labels_ = 1
    cdef int32_t[::1] n_iter = np.zeros(labels_, dtype=np.intc)
    get_n_iter(model, <int *> &n_iter[0])

    cdef int nr_feature = get_nr_feature(model)
    if bias > 0:
        nr_feature = nr_feature + 1
    if nr_class == 2 and solver_type != 4:  # solver is not Crammer-Singer
        w = np.empty((1, nr_feature), order='F')
        copy_w(&w[0, 0], model, nr_feature)
    else:
        len_w = (nr_class) * nr_feature
        w = np.empty((nr_class, nr_feature), order='F')
        copy_w(&w[0, 0], model, len_w)

    free_and_destroy_model(&model)

    return w.base, n_iter.base


def set_verbosity_wrap(int verbosity):
    """
    Control verbosity of libsvm library
    """
    set_verbosity(verbosity)
