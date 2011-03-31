"""
Binding for libdecisiontree
---------------------------

Notes
-----

Authors
-------
2011: Noel Dawe <end1@sfu.ca>
"""

import  numpy as np
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

################################################################################
# Wrapper functions

def train(np.ndarray[np.float64_t, ndim=2, mode='c'] X,
          np.ndarray[np.float64_t, ndim=1, mode='c'] Y, int
          svm_type, int kernel_type, int degree, double gamma,
          double coef0, double eps, double C, double nu,
          double cache_size, double p,
          np.ndarray[np.int32_t, ndim=1, mode='c']
              class_weight_label=np.empty(0, dtype=np.int32),
          np.ndarray[np.float64_t, ndim=1, mode='c']
              class_weight=np.empty(0),
          np.ndarray[np.float64_t, ndim=1, mode='c']
              sample_weight=np.empty(0),
          int shrinking=0, int probability=0):

    """
    Train the model using libsvm (low-level method)

    Parameters
    ----------
    X: array-like, dtype=float, size=[n_samples, n_features]

    Y: array, dtype=float, size=[n_samples]
        target vector

    svm_type : {0, 1, 2, 3, 4}
        Type of SVM: C SVC, nu SVC, one class, epsilon SVR, nu SVR

    kernel_type : {0, 1, 2, 3, 4}
        Kernel to use in the model: linear, polynomial, RBF, sigmoid
        or precomputed.

    degree : int
        Degree of the polynomial kernel (only relevant if kernel is
        set to polynomial)

    gamma : float
        Gamma parameter in RBF kernel (only relevant if kernel is set
        to RBF)

    coef0 : float
        Independent parameter in poly/sigmoid kernel.

    eps : float
        Stopping criteria.

    C : float
        C parameter in C-Support Vector Classification

    nu : float

    cache_size : float

    Return
    ------
    support : array, shape=[n_support]
        index of support vectors

    support_vectors : array, shape=[n_support, n_features]
        support vectors (equivalent to X[support]). Will return an
        empty array in the case of precomputed kernel.

    n_class_SV : array
        number of support vectors in each class.

    sv_coef : array
        coefficients of support vectors in decision function.

    intercept : array
        intercept in decision function

    label : labels for different classes (only relevant in classification).

    probA, probB : array
        probability estimates, empty array for probability=False
    """

    cdef svm_parameter *param
    cdef svm_problem *problem
    cdef svm_model *model
    cdef char *error_msg
    cdef np.npy_intp SV_len    
    cdef np.npy_intp nr


    if len(sample_weight) == 0:
        sample_weight = np.ones(X.shape[0], dtype=np.float64)
    else:
        assert sample_weight.shape[0] == X.shape[0], \
               "sample_weight and X have incompatible shapes: " + \
               "sample_weight has %s samples while X has %s" % \
               (sample_weight.shape[0], X.shape[0])

    # set libsvm problem
    problem = set_problem(X.data, Y.data, sample_weight.data,
                          X.shape, kernel_type)

    param = set_parameter(svm_type, kernel_type, degree, gamma, coef0,
                          nu, cache_size, C, eps, p, shrinking,
                          probability, <int> class_weight.shape[0],
                          class_weight_label.data, class_weight.data)

    # check parameters
    if (param == NULL or problem == NULL):
        raise MemoryError("Seems we've run out of of memory")
    error_msg = svm_check_parameter(problem, param);
    if error_msg:
        free_problem(problem)
        free_param(param)
        raise ValueError(error_msg)

    # this does the real work
    model = svm_train(problem, param)

    # from here until the end, we just copy the data returned by
    # svm_train
    SV_len  = get_l(model)
    n_class = get_nr(model)

    # copy model.sv_coef
    cdef np.ndarray[np.float64_t, ndim=2, mode='c'] sv_coef
    sv_coef = np.empty((n_class-1, SV_len), dtype=np.float64)
    copy_sv_coef (sv_coef.data, model)

    # copy model.rho into the intercept
    # the intercept is just model.rho but with sign changed
    cdef np.ndarray[np.float64_t, ndim=1, mode='c'] intercept
    intercept = np.empty(n_class*(n_class-1)/2, dtype=np.float64)
    copy_intercept (intercept.data, model, intercept.shape)

    cdef np.ndarray[np.int32_t, ndim=1, mode='c'] support
    support = np.empty (SV_len, dtype=np.int32)
    copy_support (support.data, model)

    # copy model.SV
    cdef np.ndarray[np.float64_t, ndim=2, mode='c'] support_vectors
    if kernel_type == 4:
        support_vectors = np.empty((0, 0), dtype=np.float64)
    else:
        support_vectors = np.empty((SV_len, X.shape[1]), dtype=np.float64)
        copy_SV(support_vectors.data, model, support_vectors.shape)

    # copy model.nSV
    # TODO: do only in classification
    cdef np.ndarray[np.int32_t, ndim=1, mode='c'] n_class_SV
    n_class_SV = np.empty(n_class, dtype=np.int32)
    copy_nSV(n_class_SV.data, model)

    # copy label
    cdef np.ndarray[np.int32_t, ndim=1, mode='c'] label
    label = np.empty((n_class), dtype=np.int32)
    copy_label(label.data, model)

    # copy probabilities
    cdef np.ndarray[np.float64_t, ndim=1, mode='c'] probA
    cdef np.ndarray[np.float64_t, ndim=1, mode='c'] probB
    if probability != 0:
        if svm_type < 2: # SVC and NuSVC
            probA = np.empty(n_class*(n_class-1)/2, dtype=np.float64)
            probB = np.empty(n_class*(n_class-1)/2, dtype=np.float64)
            copy_probB(probB.data, model, probB.shape)
        else:
            probA = np.empty(1, dtype=np.float64)
            probB = np.empty(0, dtype=np.float64)
        copy_probA(probA.data, model, probA.shape)

    # memory deallocation
    svm_free_and_destroy_model(&model)
    free_problem(problem)
    free_param(param)

    return support, support_vectors, n_class_SV, sv_coef, intercept, label, \
           probA, probB


def predict(np.ndarray[np.float64_t, ndim=2, mode='c'] X,
            np.ndarray[np.float64_t, ndim=2, mode='c'] SV,
            np.ndarray[np.float64_t, ndim=2, mode='c'] sv_coef,
            np.ndarray[np.float64_t, ndim=1, mode='c'] intercept,
            int svm_type, int kernel_type, int degree,
            double gamma, double coef0, double eps, double C, 
            double nu, double cache_size, double p,
            np.ndarray[np.int32_t, ndim=1, mode='c'] nSV,
            np.ndarray[np.int32_t, ndim=1, mode='c'] support,
            np.ndarray[np.int32_t, ndim=1, mode='c'] label,
            np.ndarray[np.int32_t, ndim=1]
                class_weight_label=np.empty(0, dtype=np.int32),
          np.ndarray[np.float64_t, ndim=1, mode='c']
              class_weight=np.empty(0),
            np.ndarray[np.float64_t, ndim=1, mode='c'] probA=np.empty(0),
            np.ndarray[np.float64_t, ndim=1, mode='c'] probB=np.empty(0),
            int shrinking=0, int probability=0):
    """
    Predict target values of X given a model (low-level method)

    Parameters
    ----------
    X: array-like, dtype=float, size=[n_samples, n_features]

    svm_type : {0, 1, 2, 3, 4}
        Type of SVM: C SVC, nu SVC, one class, epsilon SVR, nu SVR

    kernel_type : {0, 1, 2, 3, 4}
        Kernel to use in the model: linear, polynomial, RBF, sigmoid
        or precomputed.

    degree : int
        Degree of the polynomial kernel (only relevant if kernel is
        set to polynomial)

    gamma : float
        Gamma parameter in RBF kernel (only relevant if kernel is set
        to RBF)

    coef0 : float
        Independent parameter in poly/sigmoid kernel.

    eps : float
        Stopping criteria.

    C : float
        C parameter in C-Support Vector Classification


    Return
    ------
    dec_values : array
        predicted values.
    """
    cdef np.ndarray[np.float64_t, ndim=1, mode='c'] dec_values
    cdef svm_parameter *param
    cdef svm_model *model

    param = set_parameter(svm_type, kernel_type, degree, gamma, coef0,
                          nu, cache_size, C, eps, p, shrinking,
                          probability, <int> class_weight.shape[0],
                          class_weight_label.data, class_weight.data)

    model = set_model(param, <int> nSV.shape[0], SV.data, SV.shape,
                      support.data, support.shape, sv_coef.strides,
                      sv_coef.data, intercept.data, nSV.data,
                      label.data, probA.data, probB.data)
    
    #TODO: use check_model
    dec_values = np.empty(X.shape[0])
    if copy_predict(X.data, model, X.shape, dec_values.data) < 0:
        raise MemoryError("We've run out of of memory")
    free_model(model)
    free_param(param)
    return dec_values
