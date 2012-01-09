"""
Wrapper for liblinear

Author: fabian.pedregosa@inria.fr
"""

import  numpy as np
cimport numpy as np

cdef extern from "src/liblinear/linear.h":
    cdef struct feature_node
    cdef struct problem
    cdef struct model
    cdef struct parameter
    char *check_parameter (problem *prob, parameter *param)
    model *train (problem *prob, parameter *param)
    int get_nr_feature (model *model)
    int get_nr_class (model *model)
    void free_and_destroy_model (model **)
    void destroy_param (parameter *)

cdef extern from "src/liblinear/liblinear_helper.c":
    void copy_w(char *, model *, int)
    parameter *set_parameter (int, double, double, int,
                             char *, char *)
    problem *set_problem (char *, char *, np.npy_intp *, double)
    problem *csr_set_problem (char *values, np.npy_intp *n_indices,
        char *indices, np.npy_intp *n_indptr, char *indptr, char *Y,
        np.npy_intp n_features, double bias)
    parameter *set_parameter(int, double, double, int, char *, char *)

    model *set_model(parameter *, char *, np.npy_intp *, char *, double)
    int copy_predict(char *, model *, np.npy_intp *, char *)

    int csr_copy_predict(
        np.npy_intp n_features, np.npy_intp *data_size, char *data,
        np.npy_intp *index_size, char *index, np.npy_intp
        *intptr_size, char *intptr, model *model, char *dec_values)

    int csr_copy_predict_values(
        np.npy_intp n_features, np.npy_intp *data_size, char *data, np.npy_intp
        *index_size, char *index, np.npy_intp *indptr_shape, char
        *intptr, model *model_, char *dec_values, int nr_class)


    int csr_copy_predict_proba(
        np.npy_intp n_features, np.npy_intp *data_size, char *data,
        np.npy_intp *index_size, char *index, np.npy_intp
        *indptr_shape, char *indptr, model *model_, char *dec_values)

    int copy_prob_predict(char *, model *, np.npy_intp *, char *)
    int copy_predict_values(char *, model *, np.npy_intp *, char *, int)
    int copy_label(char *, model *, int)
    double get_bias(model *)
    void free_problem (problem *)
    void free_parameter (parameter *)


def train_wrap ( np.ndarray[np.float64_t, ndim=2, mode='c'] X,
                 np.ndarray[np.int32_t, ndim=1, mode='c'] Y, int
                 solver_type, double eps, double bias, double C, 
                 np.ndarray[np.int32_t, ndim=1] weight_label,
                 np.ndarray[np.float64_t, ndim=1] weight):
    """
    Wrapper for train methd.
    """
    cdef parameter *param
    cdef problem *problem
    cdef model *model
    cdef char *error_msg
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

    cdef np.ndarray[np.int32_t, ndim=1, mode='c'] label
    label = np.empty(nr_class, dtype=np.int32)
    copy_label(label.data, model, nr_class)

    ### FREE
    free_and_destroy_model(&model)
    free_problem(problem)
    free_parameter(param)
    # destroy_param(param)  don't call this or it will destroy weight_label and weight

    return w, label

def csr_train_wrap ( int n_features,
                 np.ndarray[np.float64_t, ndim=1, mode='c'] X_values,
                 np.ndarray[np.int32_t,   ndim=1, mode='c'] X_indices,
                 np.ndarray[np.int32_t,   ndim=1, mode='c'] X_indptr,
                 np.ndarray[np.int32_t, ndim=1, mode='c'] Y, int
                 solver_type, double eps, double bias, double C,
                 np.ndarray[np.int32_t, ndim=1] weight_label,
                 np.ndarray[np.float64_t, ndim=1] weight):
    """
    Wrapper for train.

    X matrix is given in CSR sparse format.
    """
    cdef parameter *param
    cdef problem *problem
    cdef model *model
    cdef char *error_msg
    cdef int len_w

    problem = csr_set_problem(X_values.data, X_indices.shape,
              X_indices.data, X_indptr.shape, X_indptr.data, Y.data, n_features, bias)

    param = set_parameter(solver_type, eps, C, weight.shape[0], weight_label.data, weight.data)

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

    cdef np.ndarray[np.int32_t, ndim=1, mode='c'] label
    label = np.empty((nr_class), dtype=np.int32)
    copy_label(label.data, model, nr_class)

    ### FREE
    free_and_destroy_model(&model)
    free_problem(problem)
    free_parameter(param)
    # destroy_param(param)  don't call this or it will destroy weight_label and weight

    return w, label


def decision_function_wrap(
    np.ndarray[np.float64_t, ndim=2, mode='c'] T,
    np.ndarray[np.float64_t, ndim=2, mode='fortran'] coef_,   
    int solver_type, double eps, double C,
    np.ndarray[np.int32_t, ndim=1, mode='c'] weight_label,
    np.ndarray[np.float64_t, ndim=1, mode='c'] weight,
    np.ndarray[np.int32_t, ndim=1, mode='c'] label,
    double bias):

    cdef np.ndarray[np.float64_t, ndim=2, mode='c'] dec_values
    cdef parameter *param
    cdef model *model

    param = set_parameter(
        solver_type, eps, C, weight.shape[0], weight_label.data, weight.data)

    model = set_model(param, coef_.data, coef_.shape, label.data, bias)

    n_class = label.shape[0]
    if n_class <= 2: n_class = 1
    dec_values = np.empty((T.shape[0], n_class), dtype=np.float64)

    if copy_predict_values(T.data, model, T.shape, dec_values.data, n_class) < 0:
        raise MemoryError("We've run out of of memory")

    ### FREE
    free_parameter(param)
    free_and_destroy_model(&model)
    return dec_values



def csr_decision_function_wrap(
    int n_features,
    np.ndarray[np.float64_t, ndim=1, mode='c'] T_values,
    np.ndarray[np.int32_t,   ndim=1, mode='c'] T_indices,
    np.ndarray[np.int32_t,   ndim=1, mode='c'] T_indptr,
    np.ndarray[np.float64_t, ndim=2, mode='fortran'] coef_,
    int solver_type, double eps, double C,
    np.ndarray[np.int32_t, ndim=1, mode='c'] weight_label,
    np.ndarray[np.float64_t, ndim=1, mode='c'] weight,
    np.ndarray[np.int32_t, ndim=1, mode='c'] label,
    double bias):
    """
    Predict from model

    Test data given in CSR format
    """

    cdef np.ndarray[np.float64_t, ndim=2, mode='c'] dec_values
    cdef parameter *param
    cdef model *model

    param = set_parameter(
        solver_type, eps, C, weight.shape[0], weight_label.data, weight.data)

    model = set_model(param, coef_.data, coef_.shape, label.data, bias)

    n_class = label.shape[0]
    if n_class <= 2: n_class = 1

    dec_values = np.empty((T_indptr.shape[0] - 1, n_class), dtype=np.float64)

    if csr_copy_predict_values(
        n_features, T_values.shape, T_values.data, T_indices.shape,
        T_indices.data, T_indptr.shape, T_indptr.data, model,
        dec_values.data, n_class) < 0:
        raise MemoryError("We've run out of of memory")

    ### FREE
    free_parameter(param)
    free_and_destroy_model(&model)
    return dec_values
                          

def predict_wrap(
    np.ndarray[np.float64_t, ndim=2, mode='c'] T,
    np.ndarray[np.float64_t, ndim=2, mode='fortran'] coef_,
    int solver_type, double eps, double C,
    np.ndarray[np.int32_t, ndim=1, mode='c'] weight_label,
    np.ndarray[np.float64_t, ndim=1, mode='c'] weight,
    np.ndarray[np.int32_t, ndim=1, mode='c'] label,
    double bias):

    cdef np.ndarray[np.int32_t, ndim=1, mode='c'] dec_values
    cdef parameter *param
    cdef model *model

    param = set_parameter(solver_type, eps, C, weight.shape[0], weight_label.data, weight.data)

    model = set_model(param, coef_.data, coef_.shape, label.data, bias)

    dec_values = np.empty(T.shape[0], dtype=np.int32)
    if copy_predict(T.data, model, T.shape, dec_values.data) < 0:
        raise MemoryError("We've run out of of memory")

    ### FREE
    free_parameter(param)
    free_and_destroy_model(&model)
    return dec_values


def csr_predict_wrap(
        int n_features,
        np.ndarray[np.float64_t, ndim=1, mode='c'] T_values,
        np.ndarray[np.int32_t,   ndim=1, mode='c'] T_indices,
        np.ndarray[np.int32_t,   ndim=1, mode='c'] T_indptr,
        np.ndarray[np.float64_t, ndim=2, mode='fortran'] coef_,
        int solver_type, double eps, double C,
        np.ndarray[np.int32_t, ndim=1, mode='c'] weight_label,
        np.ndarray[np.float64_t, ndim=1, mode='c'] weight,
        np.ndarray[np.int32_t, ndim=1, mode='c'] label,
        double bias):
    """
    Predict from model

    Test data given in CSR format
    """

    cdef np.ndarray[np.int32_t, ndim=1, mode='c'] dec_values
    cdef parameter *param
    cdef model *model

    param = set_parameter(solver_type, eps, C, weight.shape[0], weight_label.data, weight.data)

    model = set_model(param, coef_.data, coef_.shape, label.data, bias)

    dec_values = np.empty(T_indptr.shape[0] - 1, dtype=np.int32)
    if csr_copy_predict(n_features, T_values.shape, T_values.data,
                        T_indices.shape, T_indices.data,
                        T_indptr.shape, T_indptr.data,
                        model, dec_values.data) < 0:
        raise MemoryError("We've run out of of memory")

    ### FREE
    free_parameter(param)
    free_and_destroy_model(&model)
    return dec_values

    


def predict_prob_wrap(np.ndarray[np.float64_t, ndim=2, mode='c'] T,
                 np.ndarray[np.float64_t, ndim=2, mode='fortran'] coef_,
                 int solver_type, double eps, double C,
                 np.ndarray[np.int32_t, ndim=1, mode='c'] weight_label,
                 np.ndarray[np.float64_t, ndim=1, mode='c'] weight,
                 np.ndarray[np.int32_t, ndim=1, mode='c'] label,
                 double bias):
    """
    Predict probabilities

    svm_model stores all parameters needed to predict a given value.

    For speed, all real work is done at the C level in function
    copy_predict (libsvm_helper.c).

    We have to reconstruct model and parameters to make sure we stay
    in sync with the python object. predict_wrap skips this step.

    See scikits.learn.svm.predict for a complete list of parameters.

    Parameters
    ----------
    X: array-like, dtype=float
    Y: array
        target vector

    Returns
    -------
    dec_values : array
        predicted values.
    """
    cdef np.ndarray[np.float64_t, ndim=2, mode='c'] dec_values
    cdef parameter *param
    cdef model *model

    param = set_parameter(solver_type, eps, C, weight.shape[0], weight_label.data, weight.data)

    model = set_model(param, coef_.data, coef_.shape, label.data, bias)

    cdef int nr_class = get_nr_class(model)
    dec_values = np.empty((T.shape[0], nr_class), dtype=np.float64)
    if copy_prob_predict(T.data, model, T.shape, dec_values.data) < 0:
        raise MemoryError("We've run out of of memory")

    ### FREE
    free_parameter(param)
    free_and_destroy_model(&model)

    return dec_values



def csr_predict_prob(
        int n_features,
        np.ndarray[np.float64_t, ndim=1, mode='c'] T_values,
        np.ndarray[np.int32_t,   ndim=1, mode='c'] T_indices,
        np.ndarray[np.int32_t,   ndim=1, mode='c'] T_indptr,
        np.ndarray[np.float64_t, ndim=2, mode='fortran'] coef_,
        int solver_type, double eps, double C,
        np.ndarray[np.int32_t, ndim=1, mode='c'] weight_label,
        np.ndarray[np.float64_t, ndim=1, mode='c'] weight,
        np.ndarray[np.int32_t, ndim=1, mode='c'] label,
        double bias):
    """
    Predict probability from model

    Test data given in CSR format
    """

    cdef np.ndarray[np.float64_t, ndim=2, mode='c'] dec_values
    cdef parameter *param
    cdef model *model

    param = set_parameter(solver_type, eps, C, weight.shape[0],
                          weight_label.data, weight.data)

    model = set_model(param, coef_.data, coef_.shape, label.data, bias)
    cdef int nr_class = get_nr_class(model)
    dec_values = np.empty((T_indptr.shape[0]-1, nr_class), dtype=np.float64)

    if csr_copy_predict_proba(n_features, T_values.shape, T_values.data,
                        T_indices.shape, T_indices.data,
                        T_indptr.shape, T_indptr.data,
                        model, dec_values.data) < 0:
        raise MemoryError("We've run out of of memory")

    ### FREE
    free_parameter(param)
    free_and_destroy_model(&model)
    return dec_values
