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
    void set_verbosity(int)
