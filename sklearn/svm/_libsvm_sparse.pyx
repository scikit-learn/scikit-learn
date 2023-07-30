import  numpy as np
cimport numpy as cnp
from scipy import sparse
from ..utils._cython_blas cimport _dot
cnp.import_array()

cdef extern from *:
    ctypedef char* const_char_p "const char*"

################################################################################
# Includes

cdef extern from "_svm_cython_blas_helpers.h":
    ctypedef double (*dot_func)(int, const double*, int, const double*, int)
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
    svm_csr_problem * csr_set_problem (
        char *, cnp.npy_intp *, char *, cnp.npy_intp *, char *, char *, char *, int)
    svm_csr_model *csr_set_model(svm_parameter *param, int nr_class,
                                 char *SV_data, cnp.npy_intp *SV_indices_dims,
                                 char *SV_indices, cnp.npy_intp *SV_intptr_dims,
                                 char *SV_intptr,
                                 char *sv_coef, char *rho, char *nSV,
                                 char *probA, char *probB)
    svm_parameter *set_parameter (int , int , int , double, double ,
                                  double , double , double , double,
                                  double, int, int, int, char *, char *, int,
                                  int)
    void copy_sv_coef   (char *, svm_csr_model *)
    void copy_n_iter  (char *, svm_csr_model *)
    void copy_support   (char *, svm_csr_model *)
    void copy_intercept (char *, svm_csr_model *, cnp.npy_intp *)
    int copy_predict (char *, svm_csr_model *, cnp.npy_intp *, char *, BlasFunctions *)
    int csr_copy_predict_values (cnp.npy_intp *data_size, char *data, cnp.npy_intp *index_size,
                                 char *index, cnp.npy_intp *intptr_size, char *size,
                                 svm_csr_model *model, char *dec_values, int nr_class, BlasFunctions *)
    int csr_copy_predict (cnp.npy_intp *data_size, char *data, cnp.npy_intp *index_size,
                          char *index, cnp.npy_intp *intptr_size, char *size,
                          svm_csr_model *model, char *dec_values, BlasFunctions *) nogil
    int csr_copy_predict_proba (cnp.npy_intp *data_size, char *data, cnp.npy_intp *index_size,
                                char *index, cnp.npy_intp *intptr_size, char *size,
                                svm_csr_model *model, char *dec_values, BlasFunctions *) nogil

    int  copy_predict_values(char *, svm_csr_model *, cnp.npy_intp *, char *, int, BlasFunctions *)
    int  csr_copy_SV (char *values, cnp.npy_intp *n_indices,
                      char *indices, cnp.npy_intp *n_indptr, char *indptr,
                      svm_csr_model *model, int n_features)
    cnp.npy_intp get_nonzero_SV (svm_csr_model *)
    void copy_nSV     (char *, svm_csr_model *)
    void copy_probA   (char *, svm_csr_model *, cnp.npy_intp *)
    void copy_probB   (char *, svm_csr_model *, cnp.npy_intp *)
    cnp.npy_intp  get_l  (svm_csr_model *)
    cnp.npy_intp  get_nr (svm_csr_model *)
    int  free_problem   (svm_csr_problem *)
    int  free_model     (svm_csr_model *)
    int  free_param     (svm_parameter *)
    int free_model_SV(svm_csr_model *model)
    void set_verbosity(int)


def libsvm_sparse_train (int n_features,
                         const cnp.float64_t[::1] values,
                         const cnp.int32_t[::1] indices,
                         const cnp.int32_t[::1] indptr,
                         const cnp.float64_t[::1] Y,
                         int svm_type, int kernel_type, int degree, double gamma,
                         double coef0, double eps, double C,
                         const cnp.float64_t[::1] class_weight,
                         const cnp.float64_t[::1] sample_weight,
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
    problem = csr_set_problem(
        <char *> &values[0],
        <cnp.npy_intp *> indices.shape,
        <char *> &indices[0],
        <cnp.npy_intp *> indptr.shape,
        <char *> &indptr[0],
        <char *> &Y[0],
        <char *> &sample_weight[0],
        kernel_type,
    )

    cdef cnp.int32_t[::1] \
        class_weight_label = np.arange(class_weight.shape[0], dtype=np.int32)

    # set parameters
    param = set_parameter(
        svm_type,
        kernel_type,
        degree,
        gamma,
        coef0,
        nu,
        cache_size,
        C,
        eps,
        p,
        shrinking,
        probability,
        <int> class_weight.shape[0],
        <char *> &class_weight_label[0] if class_weight_label.size > 0 else NULL,
        <char *> &class_weight[0] if class_weight.size > 0 else NULL, max_iter,
        random_seed,
    )

    # check parameters
    if (param == NULL or problem == NULL):
        raise MemoryError("Seems we've run out of memory")
    error_msg = svm_csr_check_parameter(problem, param)
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

    cdef cnp.npy_intp SV_len = get_l(model)
    cdef cnp.npy_intp n_class = get_nr(model)

    cdef int[::1] n_iter
    n_iter = np.empty(max(1, n_class * (n_class - 1) // 2), dtype=np.intc)
    copy_n_iter(<char *> &n_iter[0], model)

    # copy model.sv_coef
    # we create a new array instead of resizing, otherwise
    # it would not erase previous information
    cdef cnp.float64_t[::1] sv_coef_data
    sv_coef_data = np.empty((n_class-1)*SV_len, dtype=np.float64)
    copy_sv_coef (<char *> &sv_coef_data[0] if sv_coef_data.size > 0 else NULL, model)

    cdef cnp.int32_t[::1] support
    support = np.empty(SV_len, dtype=np.int32)
    copy_support(<char *> &support[0] if support.size > 0 else NULL, model)

    # copy model.rho into the intercept
    # the intercept is just model.rho but with sign changed
    cdef cnp.float64_t[::1]intercept
    intercept = np.empty(n_class*(n_class-1)//2, dtype=np.float64)
    copy_intercept (<char *> &intercept[0], model, <cnp.npy_intp *> intercept.shape)

    # copy model.SV
    # we erase any previous information in SV
    # TODO: custom kernel
    cdef cnp.npy_intp nonzero_SV
    nonzero_SV = get_nonzero_SV (model)

    cdef cnp.float64_t[::1] SV_data
    cdef cnp.int32_t[::1] SV_indices, SV_indptr
    SV_data = np.empty(nonzero_SV, dtype=np.float64)
    SV_indices = np.empty(nonzero_SV, dtype=np.int32)
    SV_indptr = np.empty(<cnp.npy_intp>SV_len + 1, dtype=np.int32)
    csr_copy_SV(
        <char *> &SV_data[0] if SV_data.size > 0 else NULL,
        <cnp.npy_intp *> SV_indices.shape,
        <char *> &SV_indices[0] if SV_indices.size > 0 else NULL,
        <cnp.npy_intp *> SV_indptr.shape,
        <char *> &SV_indptr[0] if SV_indptr.size > 0 else NULL,
        model,
        n_features,
    )
    support_vectors_ = sparse.csr_matrix(
        (SV_data, SV_indices, SV_indptr), (SV_len, n_features)
    )

    # copy model.nSV
    # TODO: do only in classification
    cdef cnp.int32_t[::1]n_class_SV
    n_class_SV = np.empty(n_class, dtype=np.int32)
    copy_nSV(<char *> &n_class_SV[0], model)

    # # copy probabilities
    cdef cnp.float64_t[::1] probA, probB
    if probability != 0:
        if svm_type < 2:  # SVC and NuSVC
            probA = np.empty(n_class*(n_class-1)//2, dtype=np.float64)
            probB = np.empty(n_class*(n_class-1)//2, dtype=np.float64)
            copy_probB(<char *> &probB[0], model, <cnp.npy_intp *> probB.shape)
        else:
            probA = np.empty(1, dtype=np.float64)
            probB = np.empty(0, dtype=np.float64)
        copy_probA(<char *> &probA[0], model, <cnp.npy_intp *> probA.shape)
    else:
        probA = np.empty(0, dtype=np.float64)
        probB = np.empty(0, dtype=np.float64)

    svm_csr_free_and_destroy_model (&model)
    free_problem(problem)
    free_param(param)

    return (
        support.base,
        support_vectors_,
        sv_coef_data.base,
        intercept.base,
        n_class_SV.base,
        probA.base,
        probB.base,
        fit_status,
        n_iter.base,
    )


def libsvm_sparse_predict (const cnp.float64_t[::1] T_data,
                           const cnp.int32_t[::1] T_indices,
                           const cnp.int32_t[::1] T_indptr,
                           const cnp.float64_t[::1] SV_data,
                           const cnp.int32_t[::1] SV_indices,
                           const cnp.int32_t[::1] SV_indptr,
                           const cnp.float64_t[::1] sv_coef,
                           const cnp.float64_t[::1]
                           intercept, int svm_type, int kernel_type, int
                           degree, double gamma, double coef0, double
                           eps, double C,
                           const cnp.float64_t[:] class_weight,
                           double nu, double p, int
                           shrinking, int probability,
                           const cnp.int32_t[::1] nSV,
                           const cnp.float64_t[::1] probA,
                           const cnp.float64_t[::1] probB):
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
    cdef cnp.float64_t[::1] dec_values
    cdef svm_parameter *param
    cdef svm_csr_model *model
    cdef cnp.int32_t[::1] \
        class_weight_label = np.arange(class_weight.shape[0], dtype=np.int32)
    cdef int rv
    param = set_parameter(
        svm_type,
        kernel_type,
        degree,
        gamma,
        coef0,
        nu,
        100.0,  # cache size has no effect on predict
        C,
        eps,
        p,
        shrinking,
        probability,
        <int> class_weight.shape[0],
        <char *> &class_weight_label[0] if class_weight_label.size > 0 else NULL,
        <char *> &class_weight[0] if class_weight.size > 0 else NULL,
        -1,
        -1,  # random seed has no effect on predict either
    )

    model = csr_set_model(
        param, <int> nSV.shape[0],
        <char *> &SV_data[0] if SV_data.size > 0 else NULL,
        <cnp.npy_intp *>SV_indices.shape,
        <char *> &SV_indices[0] if SV_indices.size > 0 else NULL,
        <cnp.npy_intp *> SV_indptr.shape,
        <char *> &SV_indptr[0] if SV_indptr.size > 0 else NULL,
        <char *> &sv_coef[0] if sv_coef.size > 0 else NULL,
        <char *> &intercept[0],
        <char *> &nSV[0],
        <char *> &probA[0] if probA.size > 0 else NULL,
        <char *> &probB[0] if probB.size > 0 else NULL,
    )
    # TODO: use check_model
    dec_values = np.empty(T_indptr.shape[0]-1)
    cdef BlasFunctions blas_functions
    blas_functions.dot = _dot[double]
    with nogil:
        rv = csr_copy_predict(
            <cnp.npy_intp *> T_data.shape,
            <char *> &T_data[0],
            <cnp.npy_intp *> T_indices.shape,
            <char *> &T_indices[0],
            <cnp.npy_intp *> T_indptr.shape,
            <char *> &T_indptr[0],
            model,
            <char *> &dec_values[0],
            &blas_functions,
        )
    if rv < 0:
        raise MemoryError("We've run out of memory")
    # free model and param
    free_model_SV(model)
    free_model(model)
    free_param(param)
    return dec_values.base


def libsvm_sparse_predict_proba(
    const cnp.float64_t[::1] T_data,
    const cnp.int32_t[::1] T_indices,
    const cnp.int32_t[::1] T_indptr,
    const cnp.float64_t[::1] SV_data,
    const cnp.int32_t[::1] SV_indices,
    const cnp.int32_t[::1] SV_indptr,
    const cnp.float64_t[::1] sv_coef,
    const cnp.float64_t[::1]
    intercept, int svm_type, int kernel_type, int
    degree, double gamma, double coef0, double
    eps, double C,
    const cnp.float64_t[:] class_weight,
    double nu, double p, int shrinking, int probability,
    const cnp.int32_t[::1] nSV,
    const cnp.float64_t[::1] probA,
    const cnp.float64_t[::1] probB,
):
    """
    Predict values T given a model.
    """
    cdef cnp.float64_t[:, ::1] dec_values
    cdef svm_parameter *param
    cdef svm_csr_model *model
    cdef cnp.int32_t[::1] \
        class_weight_label = np.arange(class_weight.shape[0], dtype=np.int32)
    param = set_parameter(
        svm_type,
        kernel_type,
        degree,
        gamma,
        coef0,
        nu,
        100.0,  # cache size has no effect on predict
        C,
        eps,
        p,
        shrinking,
        probability,
        <int> class_weight.shape[0],
        <char *> &class_weight_label[0] if class_weight_label.size > 0 else NULL,
        <char *> &class_weight[0] if class_weight.size > 0 else NULL,
        -1,
        -1,  # random seed has no effect on predict either
    )

    model = csr_set_model(
        param,
        <int> nSV.shape[0],
        <char *> &SV_data[0] if SV_data.size > 0 else NULL,
        <cnp.npy_intp *> SV_indices.shape,
        <char *> &SV_indices[0] if SV_indices.size > 0 else NULL,
        <cnp.npy_intp *> SV_indptr.shape,
        <char *> &SV_indptr[0] if SV_indptr.size > 0 else NULL,
        <char *> &sv_coef[0] if sv_coef.size > 0 else NULL,
        <char *> &intercept[0],
        <char *> &nSV[0],
        <char *> &probA[0] if probA.size > 0 else NULL,
        <char *> &probB[0] if probB.size > 0 else NULL,
    )
    # TODO: use check_model
    cdef cnp.npy_intp n_class = get_nr(model)
    cdef int rv
    dec_values = np.empty((T_indptr.shape[0]-1, n_class), dtype=np.float64)
    cdef BlasFunctions blas_functions
    blas_functions.dot = _dot[double]
    with nogil:
        rv = csr_copy_predict_proba(
            <cnp.npy_intp *> T_data.shape,
            <char *> &T_data[0],
            <cnp.npy_intp *> T_indices.shape,
            <char *> &T_indices[0],
            <cnp.npy_intp *> T_indptr.shape,
            <char *> &T_indptr[0],
            model,
            <char *> &dec_values[0, 0],
            &blas_functions,
        )
    if rv < 0:
        raise MemoryError("We've run out of memory")
    # free model and param
    free_model_SV(model)
    free_model(model)
    free_param(param)
    return dec_values.base


def libsvm_sparse_decision_function(
    const cnp.float64_t[::1] T_data,
    const cnp.int32_t[::1] T_indices,
    const cnp.int32_t[::1] T_indptr,
    const cnp.float64_t[::1] SV_data,
    const cnp.int32_t[::1] SV_indices,
    const cnp.int32_t[::1] SV_indptr,
    const cnp.float64_t[::1] sv_coef,
    const cnp.float64_t[::1]
    intercept, int svm_type, int kernel_type, int
    degree, double gamma, double coef0, double
    eps, double C,
    const cnp.float64_t[:] class_weight,
    double nu, double p, int shrinking, int probability,
    const cnp.int32_t[::1] nSV,
    const cnp.float64_t[::1] probA,
    const cnp.float64_t[::1] probB,
):
    """
    Predict margin (libsvm name for this is predict_values)

    We have to reconstruct model and parameters to make sure we stay
    in sync with the python object.
    """
    cdef cnp.float64_t[:, ::1] dec_values
    cdef svm_parameter *param
    cdef cnp.npy_intp n_class

    cdef svm_csr_model *model
    cdef cnp.int32_t[::1] \
        class_weight_label = np.arange(class_weight.shape[0], dtype=np.int32)
    param = set_parameter(
        svm_type,
        kernel_type,
        degree,
        gamma,
        coef0,
        nu,
        100.0,  # cache size has no effect on predict
        C,
        eps,
        p,
        shrinking,
        probability,
        <int> class_weight.shape[0],
        <char *> &class_weight_label[0] if class_weight_label.size > 0 else NULL,
        <char *> &class_weight[0] if class_weight.size > 0 else NULL,
        -1,
        -1,
    )

    model = csr_set_model(
        param,
        <int> nSV.shape[0],
        <char *> &SV_data[0] if SV_data.size > 0 else NULL,
        <cnp.npy_intp *> SV_indices.shape,
        <char *> &SV_indices[0] if SV_indices.size > 0 else NULL,
        <cnp.npy_intp *> SV_indptr.shape,
        <char *> &SV_indptr[0] if SV_indptr.size > 0 else NULL,
        <char *> &sv_coef[0] if sv_coef.size > 0 else NULL,
        <char *> &intercept[0],
        <char *> &nSV[0],
        <char *> &probA[0] if probA.size > 0 else NULL,
        <char *> &probB[0] if probB.size > 0 else NULL,
    )

    if svm_type > 1:
        n_class = 1
    else:
        n_class = get_nr(model)
        n_class = n_class * (n_class - 1) // 2

    dec_values = np.empty((T_indptr.shape[0] - 1, n_class), dtype=np.float64)
    cdef BlasFunctions blas_functions
    blas_functions.dot = _dot[double]
    if csr_copy_predict_values(
            <cnp.npy_intp *> T_data.shape,
            <char *> &T_data[0],
            <cnp.npy_intp *> T_indices.shape,
            <char *> &T_indices[0],
            <cnp.npy_intp *> T_indptr.shape,
            <char *> &T_indptr[0],
            model,
            <char *> &dec_values[0, 0],
            n_class,
            &blas_functions,
    ) < 0:
        raise MemoryError("We've run out of memory")
    # free model and param
    free_model_SV(model)
    free_model(model)
    free_param(param)

    return dec_values.base


def set_verbosity_wrap(int verbosity):
    """
    Control verbosity of libsvm library
    """
    set_verbosity(verbosity)
