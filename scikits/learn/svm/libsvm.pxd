cimport numpy as np

################################################################################
# Includes

cdef extern from "svm.h":
    cdef struct svm_node
    cdef struct svm_model
    cdef struct svm_parameter
    cdef struct svm_problem
    char *svm_check_parameter(svm_problem *, svm_parameter *)
    svm_model *svm_train(svm_problem *, svm_parameter *)
    void svm_free_and_destroy_model(svm_model** model_ptr_ptr)
    void svm_cross_validation(svm_problem *, svm_parameter *, int nr_fold, double *target)


cdef extern from "libsvm_helper.c":
    # this file contains methods for accessing libsvm 'hidden' fields
    svm_node **dense_to_sparse (char *, np.npy_intp *)
    svm_parameter *set_parameter (int , int , int , double, double ,
                                  double , double , double , double,
                                  double, int, int, int, char *, char *)
    svm_problem * set_problem (char *, char *, char *, np.npy_intp *, int)

    svm_model *set_model (svm_parameter *, int, char *, np.npy_intp *,
                         char *, np.npy_intp *, np.npy_intp *, char *,
                         char *, char *, char *, char *, char *)

    void copy_sv_coef   (char *, svm_model *)
    void copy_intercept (char *, svm_model *, np.npy_intp *)
    void copy_SV        (char *, svm_model *, np.npy_intp *)
    int copy_support (char *data, svm_model *model)
    int copy_predict (char *, svm_model *, np.npy_intp *, char *)
    int copy_predict_proba (char *, svm_model *, np.npy_intp *, char *)
    int copy_predict_values(char *, svm_model *, np.npy_intp *, char *, int)
    void copy_nSV     (char *, svm_model *)
    void copy_label   (char *, svm_model *)
    void copy_probA   (char *, svm_model *, np.npy_intp *)
    void copy_probB   (char *, svm_model *, np.npy_intp *)
    np.npy_intp  get_l  (svm_model *)
    np.npy_intp  get_nr (svm_model *)
    int  free_problem   (svm_problem *)
    int  free_model     (svm_model *)
    int  free_param     (svm_parameter *)
    void set_verbosity(int)


