import warnings
import  numpy as np
cimport numpy as np
from scipy import sparse
from ..exceptions import ConvergenceWarning
from ..utils._cython_blas cimport _dot
np.import_array()

cdef extern from *:
    ctypedef char* const_char_p "const char*"

################################################################################
# Includes

cdef extern from "_svm_cython_blas_helpers.h":
    ctypedef double (*dot_func)(int, double*, int, double*, int)
    cdef struct BlasFunctions:
        dot_func dot

cdef extern from "svm.h":
    cdef struct svm_csr_node
    cdef struct svm_csr_model
    cdef struct svm_parameter
    cdef struct svm_csr_problem
    char *svm_csr_check_parameter(svm_csr_problem *, svm_parameter *)
    svm_csr_model *svm_csr_train(svm_csr_problem *, svm_parameter *, int *, BlasFunctions *) nogil
    void svm_csr_free_and_destroy_model(svm_csr_model** model_ptr_ptr)

cdef extern from "libsvm_sparse_helper.c":
    # this file contains methods for accessing libsvm 'hidden' fields
    svm_csr_problem * csr_set_problem (char *, np.npy_intp *,
         char *, np.npy_intp *, char *, char *, char *, int )
    svm_csr_model *csr_set_model(svm_parameter *param, int nr_class,
                            char *SV_data, np.npy_intp *SV_indices_dims,
                            char *SV_indices, np.npy_intp *SV_intptr_dims,
                            char *SV_intptr,
                            char *sv_coef, char *rho, char *nSV,
                            char *probA, char *probB)
    svm_parameter *set_parameter (int , int , int , double, double ,
                                  double , double , double , double,
                                  double, int, int, int, char *, char *, int,
                                  int)
    void copy_sv_coef   (char *, svm_csr_model *)
    void copy_support   (char *, svm_csr_model *)
    void copy_intercept (char *, svm_csr_model *, np.npy_intp *)
    int copy_predict (char *, svm_csr_model *, np.npy_intp *, char *, BlasFunctions *)
    int csr_copy_predict_values (np.npy_intp *data_size, char *data, np.npy_intp *index_size,
        	char *index, np.npy_intp *intptr_size, char *size,
                svm_csr_model *model, char *dec_values, int nr_class, BlasFunctions *)
    int csr_copy_predict (np.npy_intp *data_size, char *data, np.npy_intp *index_size,
        	char *index, np.npy_intp *intptr_size, char *size,
                svm_csr_model *model, char *dec_values, BlasFunctions *) nogil
    int csr_copy_predict_proba (np.npy_intp *data_size, char *data, np.npy_intp *index_size,
        	char *index, np.npy_intp *intptr_size, char *size,
                svm_csr_model *model, char *dec_values, BlasFunctions *) nogil

    int  copy_predict_values(char *, svm_csr_model *, np.npy_intp *, char *, int, BlasFunctions *)
    int  csr_copy_SV (char *values, np.npy_intp *n_indices,
        	char *indices, np.npy_intp *n_indptr, char *indptr,
                svm_csr_model *model, int n_features)
    np.npy_intp get_nonzero_SV ( svm_csr_model *)
    void copy_nSV     (char *, svm_csr_model *)
    void copy_probA   (char *, svm_csr_model *, np.npy_intp *)
    void copy_probB   (char *, svm_csr_model *, np.npy_intp *)
    np.npy_intp  get_l  (svm_csr_model *)
    np.npy_intp  get_nr (svm_csr_model *)
    int  free_problem   (svm_csr_problem *)
    int  free_model     (svm_csr_model *)
    int  free_param     (svm_parameter *)
    int free_model_SV(svm_csr_model *model)
    void set_verbosity(int)


np.import_array()


def libsvm_sparse_train ( int n_features,
                     np.ndarray[np.float64_t, ndim=1, mode='c'] values,
                     np.ndarray[np.int32_t,   ndim=1, mode='c'] indices,
                     np.ndarray[np.int32_t,   ndim=1, mode='c'] indptr,
                     np.ndarray[np.float64_t, ndim=1, mode='c'] Y,
                     int svm_type, int kernel_type, int degree, double gamma,
                     double coef0, double eps, double C,
                     np.ndarray[np.float64_t, ndim=1, mode='c'] class_weight,
                     np.ndarray[np.float64_t, ndim=1, mode='c'] sample_weight,
                     double nu, double cache_size, double p, int
                     shrinking, int probability, int max_iter,
                     int random_seed):
    """
    Wrap svm_train from libsvm using a scipy.sparse.csr matrix

    Work in progress.

    Parameters
    ----------
    n_features : number of features.
        XXX: can we retrieve this from any other parameter ?

    X : array-like, dtype=float, size=[N, D]

    Y : array, dtype=float, size=[N]
        target vector

    ...

    Notes
    -------------------
    See sklearn.svm.predict for a complete list of parameters.

    """

    cdef svm_parameter *param
    cdef svm_csr_problem *problem
    cdef svm_csr_model *model
    cdef const_char_p error_msg

    if len(sample_weight) == 0:
        sample_weight = np.ones(Y.shape[0], dtype=np.float64)
    else:
        assert sample_weight.shape[0] == indptr.shape[0] - 1, \
               "sample_weight and X have incompatible shapes: " + \
               "sample_weight has %s samples while X has %s" % \
               (sample_weight.shape[0], indptr.shape[0] - 1)

    # we should never end up here with a precomputed kernel matrix,
    # as this is always dense.
    assert(kernel_type != 4)

    # set libsvm problem
    problem = csr_set_problem(values.data, indices.shape, indices.data,
                              indptr.shape, indptr.data, Y.data,
                              sample_weight.data, kernel_type)

    cdef np.ndarray[np.int32_t, ndim=1, mode='c'] \
        class_weight_label = np.arange(class_weight.shape[0], dtype=np.int32)

    # set parameters
    param = set_parameter(svm_type, kernel_type, degree, gamma, coef0,
                          nu, cache_size, C, eps, p, shrinking,
                          probability, <int> class_weight.shape[0],
                          class_weight_label.data, class_weight.data, max_iter,
                          random_seed)

    # check parameters
    if (param == NULL or problem == NULL):
        raise MemoryError("Seems we've run out of memory")
    error_msg = svm_csr_check_parameter(problem, param);
    if error_msg:
        free_problem(problem)
        free_param(param)
        raise ValueError(error_msg)
    cdef BlasFunctions blas_functions
    blas_functions.dot = _dot[double]
    # call svm_train, this does the real work
    cdef int fit_status = 0
    with nogil:
        model = svm_csr_train(problem, param, &fit_status, &blas_functions)

    cdef np.npy_intp SV_len = get_l(model)
    cdef np.npy_intp n_class = get_nr(model)

    # copy model.sv_coef
    # we create a new array instead of resizing, otherwise
    # it would not erase previous information
    cdef np.ndarray sv_coef_data
    sv_coef_data = np.empty((n_class-1)*SV_len, dtype=np.float64)
    copy_sv_coef (sv_coef_data.data, model)

    cdef np.ndarray[np.int32_t, ndim=1, mode='c'] support
    support = np.empty(SV_len, dtype=np.int32)
    copy_support(support.data, model)

    # copy model.rho into the intercept
    # the intercept is just model.rho but with sign changed
    cdef np.ndarray intercept
    intercept = np.empty(n_class*(n_class-1)//2, dtype=np.float64)
    copy_intercept (intercept.data, model, intercept.shape)

    # copy model.SV
    # we erase any previous information in SV
    # TODO: custom kernel
    cdef np.npy_intp nonzero_SV
    nonzero_SV = get_nonzero_SV (model)

    cdef np.ndarray SV_data, SV_indices, SV_indptr
    SV_data = np.empty(nonzero_SV, dtype=np.float64)
    SV_indices = np.empty(nonzero_SV, dtype=np.int32)
    SV_indptr = np.empty(<np.npy_intp>SV_len + 1, dtype=np.int32)
    csr_copy_SV(SV_data.data, SV_indices.shape, SV_indices.data,
                SV_indptr.shape, SV_indptr.data, model, n_features)
    support_vectors_ = sparse.csr_matrix(
	(SV_data, SV_indices, SV_indptr), (SV_len, n_features))

    # copy model.nSV
    # TODO: do only in classification
    cdef np.ndarray n_class_SV
    n_class_SV = np.empty(n_class, dtype=np.int32)
    copy_nSV(n_class_SV.data, model)

    # # copy probabilities
    cdef np.ndarray probA, probB
    if probability != 0:
        if svm_type < 2: # SVC and NuSVC
            probA = np.empty(n_class*(n_class-1)//2, dtype=np.float64)
            probB = np.empty(n_class*(n_class-1)//2, dtype=np.float64)
            copy_probB(probB.data, model, probB.shape)
        else:
            probA = np.empty(1, dtype=np.float64)
            probB = np.empty(0, dtype=np.float64)
        copy_probA(probA.data, model, probA.shape)
    else:
        probA = np.empty(0, dtype=np.float64)
        probB = np.empty(0, dtype=np.float64)

    svm_csr_free_and_destroy_model (&model)
    free_problem(problem)
    free_param(param)

    return (support, support_vectors_, sv_coef_data, intercept, n_class_SV,
            probA, probB, fit_status)


def libsvm_sparse_predict (np.ndarray[np.float64_t, ndim=1, mode='c'] T_data,
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
                            np.ndarray[np.float64_t, ndim=1] class_weight,
                            double nu, double p, int
                            shrinking, int probability,
                            np.ndarray[np.int32_t, ndim=1, mode='c'] nSV,
                            np.ndarray[np.float64_t, ndim=1, mode='c'] probA,
                            np.ndarray[np.float64_t, ndim=1, mode='c'] probB):
    """
    Predict values T given a model.

    For speed, all real work is done at the C level in function
    copy_predict (libsvm_helper.c).

    We have to reconstruct model and parameters to make sure we stay
    in sync with the python object.

    See sklearn.svm.predict for a complete list of parameters.

    Parameters
    ----------
    X : array-like, dtype=float
    Y : array
        target vector

    Returns
    -------
    dec_values : array
        predicted values.
    """
    cdef np.ndarray[np.float64_t, ndim=1, mode='c'] dec_values
    cdef svm_parameter *param
    cdef svm_csr_model *model
    cdef np.ndarray[np.int32_t, ndim=1, mode='c'] \
        class_weight_label = np.arange(class_weight.shape[0], dtype=np.int32)
    cdef int rv
    param = set_parameter(svm_type, kernel_type, degree, gamma,
                          coef0, nu,
                          100., # cache size has no effect on predict
                          C, eps, p, shrinking,
                          probability, <int> class_weight.shape[0], class_weight_label.data,
                          class_weight.data, -1,
                          -1) # random seed has no effect on predict either

    model = csr_set_model(param, <int> nSV.shape[0], SV_data.data,
                          SV_indices.shape, SV_indices.data,
                          SV_indptr.shape, SV_indptr.data,
                          sv_coef.data, intercept.data,
                          nSV.data, probA.data, probB.data)
    #TODO: use check_model
    dec_values = np.empty(T_indptr.shape[0]-1)
    cdef BlasFunctions blas_functions
    blas_functions.dot = _dot[double]
    with nogil:
        rv = csr_copy_predict(T_data.shape, T_data.data,
                              T_indices.shape, T_indices.data,
                              T_indptr.shape, T_indptr.data,
                              model, dec_values.data,
                              &blas_functions)
    if rv < 0:
        raise MemoryError("We've run out of memory")
    # free model and param
    free_model_SV(model)
    free_model(model)
    free_param(param)
    return dec_values


def libsvm_sparse_predict_proba(
    np.ndarray[np.float64_t, ndim=1, mode='c'] T_data,
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
    np.ndarray[np.float64_t, ndim=1] class_weight,
    double nu, double p, int shrinking, int probability,
    np.ndarray[np.int32_t, ndim=1, mode='c'] nSV,
    np.ndarray[np.float64_t, ndim=1, mode='c'] probA,
    np.ndarray[np.float64_t, ndim=1, mode='c'] probB):
    """
    Predict values T given a model.
    """
    cdef np.ndarray[np.float64_t, ndim=2, mode='c'] dec_values
    cdef svm_parameter *param
    cdef svm_csr_model *model
    cdef np.ndarray[np.int32_t, ndim=1, mode='c'] \
        class_weight_label = np.arange(class_weight.shape[0], dtype=np.int32)
    param = set_parameter(svm_type, kernel_type, degree, gamma,
                          coef0, nu,
                          100., # cache size has no effect on predict
                          C, eps, p, shrinking,
                          probability, <int> class_weight.shape[0], class_weight_label.data,
                          class_weight.data, -1,
                          -1) # random seed has no effect on predict either

    model = csr_set_model(param, <int> nSV.shape[0], SV_data.data,
                          SV_indices.shape, SV_indices.data,
                          SV_indptr.shape, SV_indptr.data,
                          sv_coef.data, intercept.data,
                          nSV.data, probA.data, probB.data)
    #TODO: use check_model
    cdef np.npy_intp n_class = get_nr(model)
    cdef int rv
    dec_values = np.empty((T_indptr.shape[0]-1, n_class), dtype=np.float64)
    cdef BlasFunctions blas_functions
    blas_functions.dot = _dot[double]
    with nogil:
        rv = csr_copy_predict_proba(T_data.shape, T_data.data,
                                    T_indices.shape, T_indices.data,
                                    T_indptr.shape, T_indptr.data,
                                    model, dec_values.data,
                                    &blas_functions)
    if rv < 0:
        raise MemoryError("We've run out of memory")
    # free model and param
    free_model_SV(model)
    free_model(model)
    free_param(param)
    return dec_values




def libsvm_sparse_decision_function(
    np.ndarray[np.float64_t, ndim=1, mode='c'] T_data,
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
    np.ndarray[np.float64_t, ndim=1] class_weight,
    double nu, double p, int shrinking, int probability,
    np.ndarray[np.int32_t, ndim=1, mode='c'] nSV,
    np.ndarray[np.float64_t, ndim=1, mode='c'] probA,
    np.ndarray[np.float64_t, ndim=1, mode='c'] probB):
    """
    Predict margin (libsvm name for this is predict_values)

    We have to reconstruct model and parameters to make sure we stay
    in sync with the python object.
    """
    cdef np.ndarray[np.float64_t, ndim=2, mode='c'] dec_values
    cdef svm_parameter *param
    cdef np.npy_intp n_class

    cdef svm_csr_model *model
    cdef np.ndarray[np.int32_t, ndim=1, mode='c'] \
        class_weight_label = np.arange(class_weight.shape[0], dtype=np.int32)
    param = set_parameter(svm_type, kernel_type, degree, gamma,
                          coef0, nu,
                          100., # cache size has no effect on predict
                          C, eps, p, shrinking,
                          probability, <int> class_weight.shape[0],
                          class_weight_label.data, class_weight.data, -1, -1)

    model = csr_set_model(param, <int> nSV.shape[0], SV_data.data,
                          SV_indices.shape, SV_indices.data,
                          SV_indptr.shape, SV_indptr.data,
                          sv_coef.data, intercept.data,
                          nSV.data, probA.data, probB.data)

    if svm_type > 1:
        n_class = 1
    else:
        n_class = get_nr(model)
        n_class = n_class * (n_class - 1) // 2

    dec_values = np.empty((T_indptr.shape[0] - 1, n_class), dtype=np.float64)
    cdef BlasFunctions blas_functions
    blas_functions.dot = _dot[double]
    if csr_copy_predict_values(T_data.shape, T_data.data,
                        T_indices.shape, T_indices.data,
                        T_indptr.shape, T_indptr.data,
                        model, dec_values.data, n_class,
                        &blas_functions) < 0:
        raise MemoryError("We've run out of memory")
    # free model and param
    free_model_SV(model)
    free_model(model)
    free_param(param)

    return dec_values


def set_verbosity_wrap(int verbosity):
    """
    Control verbosity of libsvm library
    """
    set_verbosity(verbosity)
