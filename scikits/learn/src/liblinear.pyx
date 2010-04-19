"""
Wrapper for liblinear

Author: fabian.pedregosa@inria.fr
"""

import  numpy as np
cimport numpy as np

cdef extern from "linear.h":
    cdef struct feature_node
    cdef struct problem
    cdef struct model
    cdef struct parameter
    char *check_parameter (problem *prob, parameter *param)
    model *train (problem *prob, parameter *param)
    int get_nr_feature (model *model)
    int get_nr_class (model *model)
    void destroy_model (model *)
    void destroy_param (parameter *)

cdef extern from "liblinear_helper.c":
    void copy_w(char *, model *, int)
    parameter *set_parameter (int, double, double, int,
                             char *, char *)
    problem *set_problem (char *, char *, np.npy_intp *, double)
    parameter *set_parameter(int, double, double, int, char *, char *)
                          
    model *set_model(parameter *, char *, np.npy_intp *, char *, double)
    int copy_predict(char *, model *, np.npy_intp *, char *)
    int copy_label(char *, model *, int)
    double get_bias(model *)
    void free_problem (problem *)
    void free_parameter (parameter *)


def train_wrap ( np.ndarray[np.float64_t, ndim=2, mode='c'] X,
                 np.ndarray[np.int_t, ndim=1, mode='c'] Y, int
                 solver_type, double eps, double bias, double C, int nr_weight,
                 np.ndarray[np.int_t, ndim=1] weight_label,
                 np.ndarray[np.float64_t, ndim=1] weight):
    """
    Wrapper for train
    """
    cdef parameter *param
    cdef problem *problem
    cdef model *model
    cdef char *error_msg

    problem = set_problem(X.data, Y.data, X.shape, bias)

    param = set_parameter(solver_type, eps, C, nr_weight, weight_label.data, weight.data)

    error_msg = check_parameter(problem, param)
    if error_msg:
        free_problem(problem)
        free_parameter(param)
        raise ValueError(error_msg)

    # early return
    model = train(problem, param)

    cdef np.ndarray[np.float64_t, ndim=2, mode='c'] w
    cdef int nr_class = get_nr_class(model)
    cdef int nr_feature = get_nr_feature(model)
    w = np.empty((nr_class, nr_feature))
    copy_w(w.data, model, nr_class * nr_feature)

    bias = get_bias(model)

    cdef np.ndarray[np.int_t, ndim=1, mode='c'] label
    label = np.empty((nr_class), dtype=np.int)
    copy_label(label.data, model, nr_class)

    ### FREE
    destroy_model(model)
    free_problem(problem)
    free_parameter(param)
    # destroy_param(param)  don't call this or it will destroy weight_label and weight

    return w, label, bias

def predict_wrap(np.ndarray[np.float64_t, ndim=2, mode='c'] T,
                 np.ndarray[np.float64_t, ndim=2, mode='c'] coef_,
                 int solver_type, double eps, double C,
                 np.ndarray[np.int_t, ndim=1, mode='c'] weight_label,
                 np.ndarray[np.float64_t, ndim=1, mode='c'] weight,
                 np.ndarray[np.int_t, ndim=1, mode='c'] label,
                 double bias):

    cdef np.ndarray[np.int32_t, ndim=1, mode='c'] dec_values
    cdef parameter *param
    cdef model *model

    param = set_parameter(solver_type, eps, C, weight.shape[0], weight_label.data, weight.data)

    model = set_model(param, coef_.data, coef_.shape, label.data, bias)

    dec_values = np.empty(T.shape[0], dtype=np.int)
    if copy_predict(T.data, model, T.shape, dec_values.data) < 0:
        raise MemoryError("We've run out of of memory")

    ### FREE
    return dec_values
                          
    
