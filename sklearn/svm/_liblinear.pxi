from sklearn.utils._typedefs cimport intp_t

cdef extern from "_cython_blas_helpers.h":
    ctypedef double (*dot_func)(int, const double*, int, const double*, int)
    ctypedef void (*axpy_func)(int, double, const double*, int, double*, int)
    ctypedef void (*scal_func)(int, double, const double*, int)
    ctypedef double (*nrm2_func)(int, const double*, int)
    cdef struct BlasFunctions:
        dot_func dot
        axpy_func axpy
        scal_func scal
        nrm2_func nrm2


cdef extern from "linear.h":
    cdef struct feature_node
    cdef struct problem
    cdef struct model
    cdef struct parameter
    ctypedef problem* problem_const_ptr "problem const *"
    ctypedef parameter* parameter_const_ptr "parameter const *"
    ctypedef char* char_const_ptr "char const *"
    char_const_ptr check_parameter(problem_const_ptr prob, parameter_const_ptr param)
    model *train(problem_const_ptr prob, parameter_const_ptr param, BlasFunctions *blas_functions) nogil
    int get_nr_feature (model *model)
    int get_nr_class (model *model)
    void get_n_iter (model *model, int *n_iter)
    void free_and_destroy_model (model **)
    void destroy_param (parameter *)


cdef extern from "liblinear_helper.c":
    void copy_w(void *, model *, int)
    parameter *set_parameter(int, double, double, int, char *, char *, int, int, double)
    problem *set_problem (char *, int, int, int, int, double, char *, char *)
    problem *csr_set_problem (char *, int, char *, char *, int, int, int, double, char *, char *)

    model *set_model(parameter *, char *, intp_t *, char *, double)

    double get_bias(model *)
    void free_problem (problem *)
    void free_parameter (parameter *)
    void set_verbosity(int)
