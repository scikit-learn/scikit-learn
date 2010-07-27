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

We define arrays to be the same type as those in libsvm, usually of 
type C double and C int.

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
Also, the signature mode='c' is somewhat superficial, since we already
check that arrays are C-contiguous in svm.py

Authors
-------
2010: Fabian Pedregosa <fabian.pedregosa@inria.fr>
      Gael Varoquaux <gael.varoquaux@normalesup.org>
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


cdef extern from "libsvm_helper.c":
    # this file contains methods for accessing libsvm 'hidden' fields
    svm_node **dense_to_sparse (char *, np.npy_intp *)
    svm_parameter *set_parameter (int , int , int , double, double ,
                                  double , double , double , double,
                                  double, int, int, int, char *, char *)
    svm_problem * set_problem (char *, char *, np.npy_intp *, int)
    svm_problem * csr_set_problem (char *, np.npy_intp *,
        char *, np.npy_intp *, char *, char *, int )
    svm_model *set_model (svm_parameter *, int, char *, np.npy_intp *, np.npy_intp *,
                         char *, char *, char *, char *, char *, char *)
    svm_model *csr_set_model(svm_parameter *param, int nr_class,
                            char *SV_data, np.npy_intp *SV_indices_dims,
                            char *SV_indices, np.npy_intp *SV_intptr_dims,
                            char *SV_intptr,
                            char *sv_coef, char *rho, char *nSV, char *label,
                            char *probA, char *probB)
    void copy_sv_coef   (char *, svm_model *)
    void copy_intercept (char *, svm_model *, np.npy_intp *)
    void copy_SV        (char *, svm_model *, np.npy_intp *)
    int copy_predict (char *, svm_model *, np.npy_intp *, char *)
    int csr_copy_predict (np.npy_intp *data_size, char *data, np.npy_intp *index_size,
		char *index, np.npy_intp *intptr_size, char *size,
                svm_model *model, char *dec_values)
    int  copy_predict_proba (char *, svm_model *, np.npy_intp *, char *)
    int  copy_predict_values(char *, svm_model *, np.npy_intp *, char *, int)
    int  csr_copy_SV (char *values, np.npy_intp *n_indices,
		char *indices, np.npy_intp *n_indptr, char *indptr,
                svm_model *model, int n_features)
    np.npy_intp get_nonzero_SV ( svm_model *)
    void copy_nSV     (char *, svm_model *)
    void copy_label   (char *, svm_model *)
    void copy_probA   (char *, svm_model *, np.npy_intp *)
    void copy_probB   (char *, svm_model *, np.npy_intp *)
    np.npy_intp  get_l  (svm_model *)
    np.npy_intp  get_nr (svm_model *)
    int  free_problem   (svm_problem *)
    int  free_model     (svm_model *)
    int  free_model_SV  (svm_model *)
    int  free_param     (svm_parameter *)


################################################################################
# Wrapper functions

def train_wrap (  np.ndarray[np.float64_t, ndim=2, mode='c'] X, 
                  np.ndarray[np.float64_t, ndim=1, mode='c'] Y, int
                  svm_type, int kernel_type, int degree, double gamma,
                  double coef0, double eps, double C, 
                  np.ndarray[np.float64_t, ndim=2, mode='c'] SV,
                  np.ndarray[np.float64_t, ndim=2, mode='c'] sv_coef,
                  np.ndarray[np.float64_t, ndim=1, mode='c'] intercept,
                  np.ndarray[np.int32_t,   ndim=1, mode='c'] weight_label,
                  np.ndarray[np.float64_t, ndim=1, mode='c'] weight,
                  np.ndarray[np.int32_t,   ndim=1, mode='c'] nclass_SV,
                  double nu, double cache_size, double p, int
                  shrinking, int probability):
    """
    Wrap svm_train from libsvm

    Parameters
    ----------
    X: array-like, dtype=float, size=[N, D]

    Y: array, dtype=float, size=[N]
        target vector

    ...

    Notes
    -------------------
    See scikits.learn.svm.predict for a complete list of parameters.

    Return
    ------
    sv_coef: array of coeficients for support vector in decision
            function (aka alphas)
    intercept : array
        constants in decision functions
    SV : array-like
        support vectors
    TODO
    """

    cdef svm_parameter *param
    cdef svm_problem *problem
    cdef svm_model *model
    cdef char *error_msg

    # set libsvm problem
    problem = set_problem(X.data, Y.data, X.shape, kernel_type)

    # set parameters
    param = set_parameter(svm_type, kernel_type, degree, gamma,
                          coef0, nu, cache_size,
                          C, eps, p, shrinking, probability,
                          <int> weight.shape[0], weight_label.data, weight.data)

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

    cdef np.npy_intp SV_len = get_l(model)
    cdef np.npy_intp nr     = get_nr(model)

    # copy model.sv_coef
    # we create a new array instead of resizing, otherwise
    # it would not erase previous information
    sv_coef.resize ((nr-1, SV_len), refcheck=False)
    copy_sv_coef (sv_coef.data, model)

    # copy model.rho into the intercept
    # the intercept is just model.rho but with sign changed
    intercept.resize (nr*(nr-1)/2, refcheck=False)
    copy_intercept (intercept.data, model, intercept.shape)

    # copy model.SV
    # we erase any previous information in SV
    SV.resize((0,0), refcheck=False)
    if kernel_type == 4:
        SV.resize((SV_len, 1), refcheck=False)
    else:
        SV.resize((SV_len, X.shape[1]), refcheck=False)
    copy_SV(SV.data, model, SV.shape)

    # copy model.nSV
    # TODO: do only in classification
    nclass_SV.resize(nr, refcheck=False)
    copy_nSV(nclass_SV.data, model)

    # copy label
    cdef np.ndarray[np.int32_t, ndim=1, mode='c'] label
    label = np.empty((nr), dtype=np.int32)
    copy_label(label.data, model)

    # copy probabilities
    cdef np.ndarray[np.float64_t, ndim=1, mode='c'] probA
    cdef np.ndarray[np.float64_t, ndim=1, mode='c'] probB
    if probability != 0:
        # this is only valid for SVC
        probA = np.empty(nr*(nr-1)/2, dtype=np.float64)
        probB = np.empty(nr*(nr-1)/2, dtype=np.float64)
        copy_probA(probA.data, model, probA.shape)
        copy_probB(probB.data, model, probB.shape)

    free_model(model)
    free_problem(problem)
    free_param(param)

    return label, probA, probB

def csr_train_wrap ( int n_features,
                     np.ndarray[np.float64_t, ndim=1, mode='c'] values,
                     np.ndarray[np.int32_t,   ndim=1, mode='c'] indices,
                     np.ndarray[np.int32_t,   ndim=1, mode='c'] indptr,
                     np.ndarray[np.float64_t, ndim=1, mode='c'] Y, 
                     int svm_type, int kernel_type, int degree, double gamma,
                     double coef0, double eps, double C, 
                     np.ndarray[np.float64_t, ndim=1, mode='c'] SV_data,
                     np.ndarray[np.int32_t,   ndim=1, mode='c'] SV_indices,
                     np.ndarray[np.int32_t,   ndim=1, mode='c'] SV_indptr,            
                     np.ndarray[np.float64_t, ndim=1, mode='c'] sv_coef_data,
                     np.ndarray[np.float64_t, ndim=1, mode='c'] intercept,
                     np.ndarray[np.int32_t,   ndim=1, mode='c'] weight_label,
                     np.ndarray[np.float64_t, ndim=1, mode='c'] weight,
                     np.ndarray[np.int32_t,   ndim=1, mode='c'] nclass_SV,
                     double nu, double cache_size, double p, int
                     shrinking, int probability):
    """
    Wrap svm_train from libsvm using a scipy.sparse.csr matrix 

    Work in progress.

    Parameters
    ----------
    n_features : number of features.
        XXX: can we retrieve this from any other parameter ?

    X: array-like, dtype=float, size=[N, D]

    Y: array, dtype=float, size=[N]
        target vector

    ...

    Notes
    -------------------
    See scikits.learn.svm.predict for a complete list of parameters.

    """

    cdef svm_parameter *param
    cdef svm_problem *problem
    cdef svm_model *model
    cdef char *error_msg

    # set libsvm problem
    problem = csr_set_problem(values.data, indices.shape, indices.data, 
                              indptr.shape, indptr.data, Y.data, kernel_type)

    # set parameters
    param = set_parameter(svm_type, kernel_type, degree, gamma,
                          coef0, nu, cache_size,
                          C, eps, p, shrinking, probability,
                          <int> weight.shape[0], weight_label.data, weight.data)

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

    cdef np.npy_intp SV_len = get_l(model)
    cdef np.npy_intp nr     = get_nr(model)

    # copy model.sv_coef
    # we create a new array instead of resizing, otherwise
    # it would not erase previous information
    sv_coef_data.resize ((nr-1)*SV_len, refcheck=False)
    copy_sv_coef (sv_coef_data.data, model)

    # copy model.rho into the intercept
    # the intercept is just model.rho but with sign changed
    intercept.resize (nr*(nr-1)/2, refcheck=False)
    copy_intercept (intercept.data, model, intercept.shape)

    # copy model.SV
    # we erase any previous information in SV
    # TODO: custom kernel
    cdef np.npy_intp nonzero_SV
    nonzero_SV = get_nonzero_SV (model)

    # SV_data.resize((0,0), refcheck=False) # why is this needed ?
    SV_data.resize (nonzero_SV, refcheck=False)
    SV_indices.resize (nonzero_SV, refcheck=False)
    SV_indptr.resize (<np.npy_intp> SV_len + 1, refcheck=False)
    csr_copy_SV(SV_data.data, SV_indices.shape, SV_indices.data,
                SV_indptr.shape, SV_indptr.data, model, n_features)

    # copy model.nSV
    # TODO: do only in classification
    nclass_SV.resize(nr, refcheck=False)
    copy_nSV(nclass_SV.data, model)

    # # copy label
    cdef np.ndarray[np.int32_t, ndim=1, mode='c'] label
    label = np.empty((nr), dtype=np.int32)
    copy_label(label.data, model)

    # # copy probabilities
    cdef np.ndarray[np.float64_t, ndim=1, mode='c'] probA
    cdef np.ndarray[np.float64_t, ndim=1, mode='c'] probB
    if probability != 0:
        # this is only valid for SVC
        probA = np.empty(nr*(nr-1)/2, dtype=np.float64)
        probB = np.empty(nr*(nr-1)/2, dtype=np.float64)
        copy_probA(probA.data, model, probA.shape)
        copy_probB(probB.data, model, probB.shape)

    free_model(model)
    free_problem(problem)
    free_param(param)

    return label, probA, probB


def predict_from_model_wrap(np.ndarray[np.float64_t, ndim=2, mode='c'] T,
                            np.ndarray[np.float64_t, ndim=2, mode='c'] SV,
                            np.ndarray[np.float64_t, ndim=2, mode='c'] sv_coef,
                            np.ndarray[np.float64_t, ndim=1, mode='c']
                            intercept, int svm_type, int kernel_type, int
                            degree, double gamma, double coef0, double
                            eps, double C, 
                            np.ndarray[np.int32_t, ndim=1] weight_label,
                            np.ndarray[np.float64_t, ndim=1] weight,
                            double nu, double cache_size, double p, int
                            shrinking, int probability,
                            np.ndarray[np.int32_t, ndim=1, mode='c'] nSV,
                            np.ndarray[np.int32_t, ndim=1, mode='c'] label,
                            np.ndarray[np.float64_t, ndim=1, mode='c'] probA,
                            np.ndarray[np.float64_t, ndim=1, mode='c'] probB):
    """
    Predict values T given a model.

    For speed, all real work is done at the C level in function
    copy_predict (libsvm_helper.c).

    We have to reconstruct model and parameters to make sure we stay
    in sync with the python object.

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
    cdef np.ndarray[np.float64_t, ndim=1, mode='c'] dec_values
    cdef svm_parameter *param
    cdef svm_model *model
    param = set_parameter(svm_type, kernel_type, degree, gamma,
                          coef0, nu, cache_size, C, eps, p, shrinking,
                          probability, <int> weight.shape[0], weight_label.data,
                          weight.data)

    model = set_model(param, <int> nSV.shape[0], SV.data, SV.shape,
                      sv_coef.strides, sv_coef.data, intercept.data,
                      nSV.data, label.data, probA.data, probB.data)
    #TODO: use check_model
    dec_values = np.empty(T.shape[0])
    if copy_predict(T.data, model, T.shape, dec_values.data) < 0:
        raise MemoryError("We've run out of of memory")
    # free model and param
    free_model_SV(model)
    free_model(model)
    free_param(param)
    return dec_values




def csr_predict_from_model_wrap(np.ndarray[np.float64_t, ndim=1, mode='c'] T_data,
                            np.ndarray[np.int32_t,   ndim=1, mode='c'] T_indices,
                            np.ndarray[np.int32_t,   ndim=1, mode='c'] T_indptr,
                            np.ndarray[np.float64_t, ndim=1, mode='c'] SV_data,
                            np.ndarray[np.int32_t,   ndim=1, mode='c'] SV_indices,
                            np.ndarray[np.int32_t,   ndim=1, mode='c'] SV_indptr,
                            np.ndarray[np.float64_t, ndim=1, mode='c'] sv_coef,
                            np.ndarray[np.float64_t, ndim=1, mode='c']
                            intercept, int svm_type, int kernel_type, int
                            degree, double gamma, double coef0, double
                            eps, double C, 
                            np.ndarray[np.int32_t, ndim=1] weight_label,
                            np.ndarray[np.float64_t, ndim=1] weight,
                            double nu, double cache_size, double p, int
                            shrinking, int probability,
                            np.ndarray[np.int32_t, ndim=1, mode='c'] nSV,
                            np.ndarray[np.int32_t, ndim=1, mode='c'] label,
                            np.ndarray[np.float64_t, ndim=1, mode='c'] probA,
                            np.ndarray[np.float64_t, ndim=1, mode='c'] probB):
    """
    Predict values T given a model.

    For speed, all real work is done at the C level in function
    copy_predict (libsvm_helper.c).

    We have to reconstruct model and parameters to make sure we stay
    in sync with the python object.

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
    cdef np.ndarray[np.float64_t, ndim=1, mode='c'] dec_values
    cdef svm_parameter *param
    cdef svm_model *model
    param = set_parameter(svm_type, kernel_type, degree, gamma,
                          coef0, nu, cache_size, C, eps, p, shrinking,
                          probability, <int> weight.shape[0], weight_label.data,
                          weight.data)

    model = csr_set_model(param, <int> nSV.shape[0], SV_data.data,
                          SV_indices.shape, SV_indices.data,
                          SV_indptr.shape, SV_indptr.data,
                          sv_coef.data, intercept.data,
                          nSV.data, label.data, probA.data, probB.data)
    #TODO: use check_model
    dec_values = np.empty(T_indptr.shape[0]-1)
    if csr_copy_predict(T_data.shape, T_data.data,
                        T_indices.shape, T_indices.data,
                        T_indptr.shape, T_indptr.data,
                        model, dec_values.data) < 0:
        raise MemoryError("We've run out of of memory")
    # free model and param
    free_model_SV(model)
    free_model(model)
    free_param(param)
    return dec_values



def predict_prob_from_model_wrap(np.ndarray[np.float64_t, ndim=2, mode='c'] T,
                            np.ndarray[np.float64_t, ndim=2, mode='c'] SV,
                            np.ndarray[np.float64_t, ndim=2, mode='c'] sv_coef,
                            np.ndarray[np.float64_t, ndim=1, mode='c']
                            intercept, int svm_type, int kernel_type, int
                            degree, double gamma, double coef0, double
                            eps, double C, 
                            np.ndarray[np.int32_t, ndim=1] weight_label,
                            np.ndarray[np.float_t, ndim=1] weight,
                            double nu, double cache_size, double p, int
                            shrinking, int probability,
                            np.ndarray[np.int32_t, ndim=1, mode='c'] nSV,
                            np.ndarray[np.int32_t, ndim=1, mode='c'] label,
                            np.ndarray[np.float64_t, ndim=1, mode='c'] probA,
                            np.ndarray[np.float64_t, ndim=1, mode='c'] probB):
    """
    Predict probabilities

    svm_model stores all parameters needed to predict a given value.

    For speed, all real work is done at the C level in function
    copy_predict (libsvm_helper.c).

    We have to reconstruct model and parameters to make sure we stay
    in sync with the python object.

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
    cdef np.ndarray[np.float64_t, ndim=2, mode='c'] dec_values
    cdef svm_parameter *param
    cdef svm_model *model
    param = set_parameter(svm_type, kernel_type, degree, gamma,
                          coef0, nu, cache_size, C, eps, p, shrinking,
                          probability, <int> weight.shape[0], weight_label.data,
                          weight.data)
    model = set_model(param, <int> nSV.shape[0], SV.data, SV.shape, sv_coef.strides,
                      sv_coef.data, intercept.data, nSV.data, label.data,
                      probA.data, probB.data)

    cdef np.npy_intp nr = get_nr(model)    
    dec_values = np.empty((T.shape[0], nr), dtype=np.float64)
    if copy_predict_proba(T.data, model, T.shape, dec_values.data) < 0:
        raise MemoryError("We've run out of of memory")
    # free model and param
    free_model_SV(model)
    free_model(model)
    free_param(param)
    return dec_values


def predict_margin_from_model_wrap(np.ndarray[np.float64_t, ndim=2, mode='c'] T,
                            np.ndarray[np.float64_t, ndim=2, mode='c'] SV,
                            np.ndarray[np.float64_t, ndim=2, mode='c'] sv_coef,
                            np.ndarray[np.float64_t, ndim=1, mode='c']
                            intercept, int svm_type, int kernel_type, int
                            degree, double gamma, double coef0, double
                            eps, double C, 
                            np.ndarray[np.int32_t, ndim=1] weight_label,
                            np.ndarray[np.float_t, ndim=1] weight,
                            double nu, double cache_size, double p, int
                            shrinking, int probability,
                            np.ndarray[np.int32_t, ndim=1, mode='c'] nSV,
                            np.ndarray[np.int32_t, ndim=1, mode='c'] label,
                            np.ndarray[np.float64_t, ndim=1, mode='c'] probA,
                            np.ndarray[np.float64_t, ndim=1, mode='c'] probB):
    """
    Predict margin (libsvm name for this is predict_values)

    We have to reconstruct model and parameters to make sure we stay
    in sync with the python object.
    """
    cdef np.ndarray[np.float64_t, ndim=2, mode='c'] dec_values
    cdef svm_parameter *param
    cdef svm_model *model
    cdef np.npy_intp nr

    param = set_parameter(svm_type, kernel_type, degree, gamma,
                          coef0, nu, cache_size, C, eps, p, shrinking,
                          probability, <int> weight.shape[0], weight_label.data,
                          weight.data)
    model = set_model(param, <int> nSV.shape[0], SV.data, SV.shape, sv_coef.strides,
                      sv_coef.data, intercept.data, nSV.data, label.data,
                      probA.data, probB.data)

    if svm_type > 1:
        nr = 1
    else:
        nr = get_nr(model)
        nr = nr * (nr - 1) / 2
    
    dec_values = np.empty((T.shape[0], nr), dtype=np.float64)
    if copy_predict_values(T.data, model, T.shape, dec_values.data, nr) < 0:
        raise MemoryError("We've run out of of memory")
    # free model and param
    free_model_SV(model)
    free_model(model)
    free_param(param)
    return dec_values

