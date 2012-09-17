"""
Wrapper for liblinear

Author: fabian.pedregosa@inria.fr
"""

import  numpy as np
cimport numpy as np
cimport liblinear


def train_wrap ( np.ndarray[np.float64_t, ndim=2, mode='c'] X,
                 np.ndarray[np.float64_t, ndim=1, mode='c'] Y,
                 int solver_type, double eps, double bias, double C, 
                 np.ndarray[np.int32_t, ndim=1] weight_label,
                 np.ndarray[np.float64_t, ndim=1] weight):
    """
    Wrapper for train methd.
    """
    cdef parameter *param
    cdef problem *problem
    cdef model *model
    cdef char_const_ptr error_msg
    cdef int len_w

    problem = set_problem(X.data, Y.data, X.shape, bias)

    param = set_parameter(solver_type, eps, C, weight.shape[0], weight_label.data, weight.data)

    error_msg = check_parameter(problem, param)
    if error_msg:
        free_problem(problem)
        free_parameter(param)
        raise ValueError(error_msg)
 
    # early return
    model = train(problem, param)

    # coef matrix holder created as fortran since that's what's used in liblinear
    cdef np.ndarray[np.float64_t, ndim=2, mode='fortran'] w  
    cdef int nr_class = get_nr_class(model)
    cdef int nr_feature = get_nr_feature(model)
    if bias > 0: nr_feature = nr_feature + 1
    if nr_class == 2:
        w = np.empty((1, nr_feature),order='F')
        copy_w(w.data, model, nr_feature)
    else:
        len_w = (nr_class) * nr_feature
        w = np.empty((nr_class, nr_feature),order='F') 
        copy_w(w.data, model, len_w)

    ### FREE
    free_and_destroy_model(&model)
    free_problem(problem)
    free_parameter(param)
    # destroy_param(param)  don't call this or it will destroy weight_label and weight

    return w


cdef _csr_train_wrap(np.int32_t n_features,
                 np.ndarray[np.float64_t, ndim=1, mode='c'] X_values,
                 np.ndarray[np.int32_t,   ndim=1, mode='c'] X_indices,
                 np.ndarray[np.int32_t,   ndim=1, mode='c'] X_indptr,
                 np.ndarray[np.float64_t, ndim=1, mode='c'] Y,
                 int solver_type, double eps, double bias, double C,
                 np.ndarray[np.int32_t, ndim=1] weight_label,
                 np.ndarray[np.float64_t, ndim=1] weight):
    cdef parameter *param
    cdef problem *problem
    cdef model *model
    cdef char_const_ptr error_msg
    cdef int len_w

    problem = csr_set_problem(X_values.data, X_indices.shape,
                              X_indices.data, X_indptr.shape,
                              X_indptr.data, Y.data, n_features, bias)

    param = set_parameter(solver_type, eps, C, weight.shape[0],
                          weight_label.data, weight.data)

    error_msg = check_parameter(problem, param)
    if error_msg:
        free_problem(problem)
        free_parameter(param)
        raise ValueError(error_msg)

    # early return
    model = train(problem, param)

    # fortran order since that's what liblinear does
    cdef np.ndarray[np.float64_t, ndim=2, mode='fortran'] w
    cdef int nr_class = get_nr_class(model)
    cdef int nr_feature = n_features
    if bias > 0: nr_feature = nr_feature + 1
    if nr_class == 2:
        w = np.empty((1, nr_feature),order='F')
        copy_w(w.data, model, nr_feature)
    else:
        len_w = (nr_class * nr_feature)
        w = np.empty((nr_class, nr_feature),order='F')
        copy_w(w.data, model, len_w)

    ### FREE
    free_and_destroy_model(&model)
    free_problem(problem)
    free_parameter(param)
    # destroy_param(param)  don't call this or it will destroy weight_label and weight

    return w


def csr_train_wrap(X, Y, solver_type, eps, bias, C, weight_label, weight):
    """
    Wrapper for train.

    X matrix is given in CSR sparse format.
    """
    return _csr_train_wrap(X.shape[1], X.data, X.indices, X.indptr, Y,
                           solver_type, eps, bias, C, weight_label, weight)


def set_verbosity_wrap(int verbosity):
    """
    Control verbosity of libsvm library
    """
    set_verbosity(verbosity)
