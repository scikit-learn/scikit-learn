"""
Binding for libsvm[1]
---------------------
We do not use the binding that ships with libsvm because we need to
access svm_model.sv_coeff (and other fields), but libsvm does not
provide an accessor. Our solution is to export svm_model and access it
manually, this is done un function see svm_train_wrap.

libsvm uses an sparse representation for the training vectors. In
method dense_to_sparse we translate a dense matrix representation as
those produced by NumPy to a sparse representation that libsvm can
understand.

Low-level memory management is done in libsvm_helper.c. If we happen
to run out of memory a MemoryError will be raised. In practice this is
not very helpful since hight changes are malloc fails inside svm.cpp,
where no sort of memory checks are done.

These are low-level routines, not meant to be used directly. See
scikits.learn.svm for a higher-level API.

[1] http://www.csie.ntu.edu.tw/~cjlin/libsvm/

Notes
-----
Maybe we could speed it a bit further by decorating functions with
@cython.boundscheck(False), but probably it is not worth since all
work is done in lisvm_helper.c

TODO
----
Coding style, eliminate SvmModelPtr

Authors
-------
2010: Fabian Pedregosa <fabian.pedregosa@inria.fr>
      Gael Varoquaux <gael.varoquaux@normalesup.org>
"""

import  numpy as np
cimport numpy as np

################################################################################
# Default Parameters

DEF NR_WEIGHT   = 0
DEF NU          = 0.5
DEF P           = 0.1
DEF SHRINKING   = 1
DEF PROBABILITY = 0

################################################################################
# Includes

cdef extern from "svm.h":
    cdef struct svm_node
    cdef struct svm_model
    cdef struct svm_parameter
    cdef struct svm_problem
    char *svm_check_parameter(svm_problem *, svm_parameter *)
    svm_model *svm_train(svm_problem *, svm_parameter *)
    double svm_predict(svm_model *, svm_node *)

cdef extern from "libsvm_helper.c":
    # this file contains methods for accessing libsvm 'hidden' fields
    svm_node **dense_to_sparse (char *, int *)
    svm_parameter *set_parameter(int , int , int , double, double ,
                                  double , double , double , double,
                                  double, int, int, int, char *, char *)
    svm_problem *set_problem(char *, char *, int *)
    void copy_sv_coef (char *, svm_model *, int *)
    void copy_rho     (char *, svm_model *, int *)
    void copy_SV      (char *, svm_model *, int *)
    int  copy_predict (char *, svm_model *, int *, char *)
    int  get_l  (svm_model *)
    int  get_nr (svm_model *)
    int  free_problem(svm_problem *)
    int  free_model  (svm_model *)
    int  free_param  (svm_parameter *)


################################################################################
# Object definition

cdef class SvmModelPtr:
    """
    Container for two pointers, one to svm_problem and another one to
    svm_model.

    This lets us efficiently pass a svm_model from train_wrap to
    predict_from_model_wrap. See the internals of function
    scikits.learn.svm.SVM.predict for an example on how this is
    handled.
    """
    cdef svm_model *model
    cdef svm_problem *problem

    def __dealloc__(self):
        free_problem(self.problem)
        free_model(self.model)


################################################################################
# Wrapper functions


def predict_wrap( np.ndarray[np.float_t, ndim=2, mode='c'] X,
                  np.ndarray[np.float_t, ndim=1, mode='c'] Y,
                  np.ndarray[np.float_t, ndim=2, mode='c'] T,
                  int svm_type, int kernel_type,
                  int degree, double gamma,
                  double coef0, double eps, double C,
                  int nr_weight,
                  np.ndarray[np.int_t, ndim=1] weight_label,
                  np.ndarray[np.float_t, ndim=1] weight,
                  double nu, double cache_size,
                  double p, int shrinking,
                  int probability):
    """
    Wrapper for svm_train_predict in libsvm.
    Predict T learning from X, Y, where X are data points and Y are
    labels.
    
    Parameters
    ----------
    X: array-like, dtype=float, size=[N, D]

    Y: array, dtype=float, size=[N]
        target vector

    T: array-like, dtype=float, size=[M, D]
        test vector where M = number of test samples, D = dimension of
        sample space.

    Optional Parameters
    -------------------
    See scikits.learn.svm.predict for a complete list of parameters.

    Return
    ------
    dec_values : array-like, dtype=float, size=[N]
        predicted values

    Notes
    -----
    In case of dimension mismatch, no error will be reported but
    meaningless output will be given.
    """

    cdef svm_problem *problem
    cdef svm_parameter *param
    cdef svm_model *model
    cdef char *error_msg

    # set libsvm problem
    problem = set_problem(X.data, Y.data, X.shape)

    # set parameters
    if (gamma == 0): gamma = 1.0/X.shape[0]
    param = set_parameter(svm_type, kernel_type, degree, gamma, coef0, nu,
                          cache_size, C, eps, p, shrinking,
                          probability, nr_weight, weight_label.data,
                          weight.data)


    # check parameters
    if (problem == NULL or param == NULL):
        raise MemoryError("We've run out of of memory in predict_wrap")
    error_msg = svm_check_parameter(problem, param)
    if error_msg:
        free_problem(problem)
        free_param(param)
        raise ValueError(error_msg)

    # call svm_train, this does the real work
    # maybe we should check that the return model is not null
    model = svm_train(problem, param)
    if model == NULL: raise MemoryError

    # predict values
    cdef np.ndarray[np.float_t, ndim=1, mode='c'] dec_values
    dec_values = np.empty(T.shape[0])
    if copy_predict(T.data, model, T.shape, dec_values.data) < 0:
        raise MemoryError("We've run out of of memory in svm_predict")

    # free memory
    free_model(model)
    free_problem(problem)
    free_param(param)

    return dec_values



def train_wrap (  np.ndarray[np.float_t, ndim=2, mode='c'] X, 
                  np.ndarray[np.float_t, ndim=1, mode='c'] Y,
                  int svm_type, int kernel_type,
                  int degree, double gamma,
                  double coef0, double eps, double C,
                  int nr_weight,
                  np.ndarray[np.int_t, ndim=1] weight_label,
                  np.ndarray[np.float_t, ndim=1] weight,
                  double nu, double cache_size,
                  double p, int shrinking,
                  int probability):
    """
    Wrapper for svm_train in libsvm

    Parameters
    ----------
    X: array-like, dtype=float, size=[N, D]

    Y: array, dtype=float, size=[N]
        target vector

    Optional Parameters
    -------------------
    See scikits.learn.svm.predict for a complete list of parameters.

    Return
    ------
    sv_coef: array of coeficients for support vector in decision
            function (aka alphas)
    rho : array
        constants in decision functions
    SV : array-like
        support vectors
    ptr : SvmModelPtr
        container for efficient C/Python communication
    """

    cdef svm_parameter *param
    cdef svm_problem *problem
    cdef svm_model *model
    cdef char *error_msg

    # set libsvm problem
    problem = set_problem(X.data, Y.data, X.shape)

    # set parameters
    if (gamma == 0): gamma = 1.0/X.shape[0]
    param = set_parameter(svm_type, kernel_type, degree, gamma,
                          coef0, nu, cache_size,
                          C, eps, p, shrinking, probability,
                          nr_weight, weight_label.data, weight.data)

    # check parameters
    if (param == NULL or problem == NULL):
        raise MemoryError("Seems we've run out of of memory")
    error_msg = svm_check_parameter(problem, param);
    if error_msg:
        free_problem(problem)
        free_param(param)
        raise ValueError(error_msg)

    # call svm_train, this does the real work
    model = svm_train(problem, param)

    cdef int nSV = get_l(model)
    cdef int nclass = get_nr(model)

    # copy model.sv_coef 
    cdef np.ndarray[np.float_t, ndim=2, mode='c'] sv_coef
    sv_coef = np.empty((nclass-1, nSV))
    copy_sv_coef(sv_coef.data, model, sv_coef.strides)

    # copy model.rho
    cdef np.ndarray[np.float_t, ndim=1, mode='c'] rho
    rho = np.empty((nclass*nclass-1)/2)
    copy_rho(rho.data, model, rho.strides)

    # copy model.SV
    cdef np.ndarray[np.float_t, ndim=2, mode='c'] SV
    SV = np.zeros((nSV, X.shape[1]))
    copy_SV(SV.data, model, SV.strides)

    cdef SvmModelPtr ptr
    ptr = SvmModelPtr()
    ptr.model = model
    ptr.problem = problem
    return sv_coef, rho, SV, ptr


def predict_from_model_wrap(np.ndarray[np.float_t, ndim=2] T,
                                SvmModelPtr ptr):
    """
    Predict values T given a pointer to svm_model.

    svm_model stores all parameters needed to predict a given value.

    For speed, all real work is done at the C level in function
    copy_predict (libsvm_helper.c).

    Parameters
    ----------
    X: array-like, dtype=float
    Y: array
        target vector

    Optional Parameters
    -------------------
    See scikits.learn.svm.predict for a complete list of parameters.

    Return
    ------
    dec_values : array
        predicted values.
    """
    cdef np.ndarray[np.float_t, ndim=1, mode='c'] dec_values
    dec_values = np.empty(T.shape[0])
    if copy_predict(T.data, ptr.model, T.shape, dec_values.data) < 0:
        raise MemoryError("We've run out of of memory")
    return dec_values
